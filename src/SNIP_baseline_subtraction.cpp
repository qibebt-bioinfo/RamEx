#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif

#include <Rcpp.h>

#include <algorithm>
#include <thread>
#include <vector>

namespace {

inline R_xlen_t clampIterations(int iterations, R_xlen_t n) {
  if (n <= 1 || iterations <= 0) {
    return 0;
  }
  R_xlen_t maxHalfWindow = (n - 1) / 2;
  R_xlen_t iter = static_cast<R_xlen_t>(iterations);
  return std::min(iter, maxHalfWindow);
}

void snipInplace(double* data, double* buffer, R_xlen_t n,
                 int iterations, bool decreasing) {
  R_xlen_t k = clampIterations(iterations, n);
  if (n == 0 || k == 0) {
    return;
  }

  if (decreasing) {
    for (R_xlen_t i = k; i > 0; --i) {
      for (R_xlen_t j = i; j < n - i; ++j) {
        double a = data[j];
        double b = (data[j - i] + data[j + i]) / 2.0;
        buffer[j] = (b < a) ? b : a;
      }
      for (R_xlen_t j = i; j < n - i; ++j) {
        data[j] = buffer[j];
      }
    }
  } else {
    for (R_xlen_t i = 1; i <= k; ++i) {
      for (R_xlen_t j = i; j < n - i; ++j) {
        double a = data[j];
        double b = (data[j - i] + data[j + i]) / 2.0;
        buffer[j] = (b < a) ? b : a;
      }
      for (R_xlen_t j = i; j < n - i; ++j) {
        data[j] = buffer[j];
      }
    }
  }
}

void processSpectraRange(const Rcpp::NumericMatrix& spectra,
                         Rcpp::NumericMatrix& corrected,
                         int startRow,
                         int endRow,
                         int iterations,
                         bool decreasing) {
  const int n = spectra.ncol();
  if (n == 0) {
    return;
  }

  std::vector<double> baseline(n);
  std::vector<double> buffer(n);

  for (int row = startRow; row < endRow; ++row) {
    for (int col = 0; col < n; ++col) {
      baseline[col] = spectra(row, col);
    }
    snipInplace(baseline.data(), buffer.data(),
                static_cast<R_xlen_t>(n), iterations, decreasing);
    for (int col = 0; col < n; ++col) {
      corrected(row, col) = spectra(row, col) - baseline[col];
    }
  }
}

}  // namespace

extern "C" SEXP C_snip(SEXP ySEXP, SEXP iterationsSEXP, SEXP decreasingSEXP) {
  Rcpp::NumericVector y(ySEXP);
  Rcpp::NumericVector baseline = Rcpp::clone(y);

  const R_xlen_t n = baseline.size();
  const int k = Rcpp::as<int>(iterationsSEXP);
  const bool decreasing = Rcpp::as<bool>(decreasingSEXP);

  std::vector<double> buffer(static_cast<std::size_t>(n));
  snipInplace(baseline.begin(), buffer.data(), n, k, decreasing);

  return baseline;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix SNIPBaselineCpp(const Rcpp::NumericMatrix spectra,
                                    int iterations = 100,
                                    bool decreasing = false,
                                    int n_threads = 1) {
  const int numSpectra = spectra.nrow();
  const int n = spectra.ncol();

  Rcpp::NumericMatrix corrected(numSpectra, n);
  if (numSpectra == 0 || n == 0) {
    return corrected;
  }

  if (n_threads <= 0) {
    n_threads = 1;
  }
  n_threads = std::min(n_threads, numSpectra);

  auto worker = [&](int startRow, int endRow) {
    processSpectraRange(spectra, corrected, startRow, endRow,
                        iterations, decreasing);
  };

  if (n_threads == 1 || numSpectra == 1) {
    worker(0, numSpectra);
    return corrected;
  }

  std::vector<std::thread> threads;
  threads.reserve(n_threads);
  const int block = (numSpectra + n_threads - 1) / n_threads;
  for (int t = 0; t < n_threads; ++t) {
    const int startRow = t * block;
    if (startRow >= numSpectra) {
      break;
    }
    const int endRow = std::min(numSpectra, startRow + block);
    threads.emplace_back(worker, startRow, endRow);
  }

  for (auto& th : threads) {
    if (th.joinable()) {
      th.join();
    }
  }

  return corrected;
}
