// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <vector>
#include <thread> 

using namespace Rcpp;
using namespace RcppParallel;

// Structure to store prediction scores and actual labels
struct Prediction {
  double score;
  int label;
};

// Comparison function for sorting
bool comparePredictions(const Prediction &a, const Prediction &b) {
  return a.score > b.score;
}

// Worker for calculating AUC
struct AUCWorker : public Worker {
  const RMatrix<double> matrix;
  const RVector<int> group;
  RMatrix<double> aucResults;
  int numGroups;

  AUCWorker(const NumericMatrix matrix, const IntegerVector group, NumericMatrix aucResults, int numGroups)
    : matrix(matrix), group(group), aucResults(aucResults), numGroups(numGroups) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      std::vector<Prediction> predictions(matrix.nrow());
      for (std::size_t j = 0; j < matrix.nrow(); ++j) {
        predictions[j] = {matrix(j, i), group[j]};
      }

      std::sort(predictions.begin(), predictions.end(), comparePredictions);

      // Calculate AUC for each group
      for (int g = 0; g < numGroups; ++g) {
        double auc = 0.0;
        int numPositives = 0;
        int numNegatives = 0;
        double cumulativePositives = 0.0;

        for (const auto &pred : predictions) {
          if (pred.label == g) {
            numPositives++;
          } else {
            numNegatives++;
          }
        }

        if (numPositives > 0 && numNegatives > 0) {
          for (const auto &pred : predictions) {
            if (pred.label == g) {
              cumulativePositives+=1.0;
            } else {
              auc += cumulativePositives/numPositives;
            }
          }

          auc /= numNegatives;
        }  else{
          auc = 0.0;
        }
        aucResults(i, g) = auc;
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix calculateAUCParallel(NumericMatrix matrix, IntegerVector group) {
  int minLabel = Rcpp::min(group);
  IntegerVector adjustedGroup = group - minLabel;
  int numGroups = Rcpp::max(adjustedGroup) + 1; // Assuming adjustedGroup labels are 0-indexed
  NumericMatrix aucResults(matrix.ncol(), numGroups);
  AUCWorker worker(matrix, adjustedGroup, aucResults, numGroups);
  parallelFor(0, matrix.ncol(), worker);
  return aucResults;
}

struct RatioWorker : public Worker {
  const NumericMatrix& matrix;
  const IntegerMatrix& combinations;
  NumericMatrix& ratio_matrix;

  RatioWorker(const NumericMatrix& matrix,
              const IntegerMatrix& combinations,
              NumericMatrix& ratio_matrix)
    : matrix(matrix), combinations(combinations), ratio_matrix(ratio_matrix) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (int j = 0; j < matrix.nrow(); j++) {
        ratio_matrix(j, i) = matrix(j, combinations(i, 0)) /  matrix(j, combinations(i, 1));
      }
    }
  }
};

struct ResultProcessor : public Worker {
  const IntegerMatrix& combinations;
  const NumericMatrix& auc_results;
  const double threshold;
  std::vector<std::vector<std::tuple<int,int,int,double>>>& results_vector;

  ResultProcessor(const IntegerMatrix& combinations,
                  const NumericMatrix& auc_results,
                  const double threshold,
                  std::vector<std::vector<std::tuple<int,int,int,double>>>& results_vector)
    : combinations(combinations), auc_results(auc_results), threshold(threshold),
      results_vector(results_vector) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<std::tuple<int,int,int,double>> local_results;
    local_results.reserve(end - begin);
    
    for (std::size_t i = begin; i < end; i++) {
      for (int g = 0; g < auc_results.ncol(); g++) {
        double auc = auc_results(i, g);
        if (auc > threshold) {
          local_results.emplace_back(
              combinations(i, 0) + 1,
              combinations(i, 1) + 1,
              g + 1,
              auc
          );
        }
      }
    }
    static thread_local int thread_id = -1;
    if (thread_id == -1)
      thread_id = std::hash<std::thread::id>{}(std::this_thread::get_id()) % results_vector.size();

    auto& bucket = results_vector[thread_id];
    bucket.insert(bucket.end(), local_results.begin(), local_results.end());
  }
};

// [[Rcpp::export]]
DataFrame calculatePairedMarkersAUC(NumericMatrix matrix,
                                    IntegerVector group,
                                    double threshold,
                                    int batch_size = 1000,
                                    int n_threads = 0) {
  int n = matrix.ncol();
  int total_pairs = (n * (n - 1)) / 2;
  int num_chunks = ceil(total_pairs / (double)batch_size);

  if (n_threads <= 0) {
    n_threads = std::max(1u, std::thread::hardware_concurrency() - 4);
  }
  std::vector<std::vector<std::tuple<int,int,int,double>>> results_vector(n_threads);

  std::vector<int> final_col1, final_col2, final_groups;
  std::vector<double> final_aucs;

  int current_pair = 0;

  for (int chunk = 0; chunk < num_chunks; chunk++) {
    int current_batch_size = std::min(batch_size, total_pairs - current_pair);
    IntegerMatrix combinations(current_batch_size, 2);
    int pair_idx = 0;

    // Generate combinations
    for (int i = 0; i < n && pair_idx < current_batch_size; i++) {
      for (int j = i + 1; j < n && pair_idx < current_batch_size; j++) {
        if (current_pair >= (chunk * batch_size)) {
          combinations(pair_idx, 0) = i;
          combinations(pair_idx, 1) = j;
          pair_idx++;
        }
        current_pair++;
      }
    }

    // Parallel calculation
    NumericMatrix ratio_matrix(matrix.nrow(), current_batch_size);
    RatioWorker ratio_worker(matrix, combinations, ratio_matrix);
    parallelFor(0, current_batch_size, ratio_worker);

    NumericMatrix auc_results = calculateAUCParallel(ratio_matrix, group);

    ResultProcessor result_processor(combinations, auc_results, threshold, results_vector);
    parallelFor(0, current_batch_size, result_processor);
  }

  for (const auto& local_results : results_vector){
    for (const auto& result : local_results) {
      final_col1.push_back(std::get<0>(result));
      final_col2.push_back(std::get<1>(result));
      final_groups.push_back(std::get<2>(result));
      final_aucs.push_back(std::get<3>(result));
    }
  }


  return DataFrame::create(
    Named("col1") = wrap(final_col1),
    Named("col2") = wrap(final_col2),
    Named("group") = wrap(final_groups),
    Named("AUC") = wrap(final_aucs)
  );
}

