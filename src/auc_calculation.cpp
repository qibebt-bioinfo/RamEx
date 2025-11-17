// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <vector>
#include <thread> 

using namespace Rcpp;
using namespace RcppParallel;

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


// Worker for calculating AUC
struct AUCWorker : public Worker {
  const RMatrix<double> matrix;
  const RVector<int> group;
  RMatrix<double> aucResults;
  const int numGroups;

  AUCWorker(const NumericMatrix matrix, const IntegerVector group, NumericMatrix aucResults, int numGroups)
    : matrix(matrix), group(group), aucResults(aucResults), numGroups(numGroups) {}

  void operator()(std::size_t begin, std::size_t end) {
    size_t nSamples = matrix.nrow();
    std::vector<std::pair<double, int>> temp(nSamples);
    std::vector<size_t> rank_index(nSamples);
    std::vector<double> ranks(nSamples);

    for (std::size_t col = begin; col < end; ++col) {
      for (std::size_t i = 0; i < nSamples; ++i) {
        temp[i].first = matrix(i, col);
        temp[i].second = group[i];
      }
      std::iota(rank_index.begin(), rank_index.end(), 0);
      std::sort(rank_index.begin(), rank_index.end(),
                [&](size_t a,size_t b){ return temp[a].first < temp[b].first; });
      
      for (size_t r = 0; r < nSamples; ++r)
        ranks[rank_index[r]] = r + 1; 

      // Calculate AUC for each group
      for (int g = 0; g < numGroups; ++g) {
        double pos_sum = 0.0;
        int pos_count = 0;
        int neg_count = 0;

        for (size_t i = 0; i < nSamples; ++i) {
          if (temp[i].second == g) {
            pos_sum += ranks[i];
            pos_count++;
          } else {
            neg_count++;
          }
        }

        double auc = 0.0;
        if (pos_count > 0 && neg_count > 0) {
          auc = (pos_sum - (pos_count * (pos_count + 1)) / 2.0) /
                (pos_count * (double)neg_count);
        } 
        aucResults(col, g) = auc;
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix calculateAUCParallel(const NumericMatrix matrix, const IntegerVector group, int n_threads = 0) {
  int minLabel = Rcpp::min(group);
  IntegerVector adjustedGroup = group - minLabel;
  int numGroups = Rcpp::max(adjustedGroup) + 1; // Assuming adjustedGroup labels are 0-indexed
  NumericMatrix aucResults(matrix.ncol(), numGroups);
  if (n_threads <= 0)
    n_threads = std::max(1u, std::thread::hardware_concurrency()-4);

  AUCWorker worker(matrix, adjustedGroup, aucResults, numGroups);
  parallelFor(0, matrix.ncol(), worker, std::max(1, (int)(matrix.ncol() / n_threads + 1)));
  return aucResults;
}

// [[Rcpp::export]]
DataFrame calculatePairedMarkersAUC(NumericMatrix matrix,
                                    IntegerVector group,
                                    double threshold,
                                    int batch_size = 1000,
                                    int n_threads = 0) {
  int n = matrix.ncol();
  int total_pairs = (n * (n - 1)) / 2;
  int num_chunks = ceil(total_pairs / (double)batch_size);

  std::vector<int> final_col1, final_col2, final_groups;
  std::vector<double> final_aucs;

  std::vector<std::pair<int,int>> all_pairs;
  all_pairs.reserve(total_pairs);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      all_pairs.emplace_back(i, j);
    }
  }

  for (int chunk = 0; chunk < num_chunks; ++chunk) {
    int start_pair = chunk * batch_size;
    int end_pair   = std::min(total_pairs, start_pair + batch_size);
    int current_batch_size = end_pair - start_pair;

    IntegerMatrix combinations(current_batch_size, 2);
    for (int k = 0; k < current_batch_size; ++k) {
      combinations(k, 0) = all_pairs[start_pair + k].first;
      combinations(k, 1) = all_pairs[start_pair + k].second;
    }

    // Parallel calculation
    NumericMatrix ratio_matrix(matrix.nrow(), current_batch_size);
    RatioWorker ratio_worker(matrix, combinations, ratio_matrix);
    ratio_worker(0, current_batch_size);

    NumericMatrix auc_results = calculateAUCParallel(ratio_matrix, group, n_threads);

    for (int i = 0; i < current_batch_size; ++i) {
      int global_index = start_pair + i; 
      for (int g = 0; g < auc_results.ncol(); ++g) {
        double auc = auc_results(i, g);
        if (auc > threshold) {
          final_col1.push_back(all_pairs[global_index].first + 1);
          final_col2.push_back(all_pairs[global_index].second + 1);
          final_groups.push_back(g + 1);
          final_aucs.push_back(auc);
        }
      }
    }
  }

  return DataFrame::create(
    Named("col1") = wrap(final_col1),
    Named("col2") = wrap(final_col2),
    Named("group") = wrap(final_groups),
    Named("AUC") = wrap(final_aucs)
  );
}

