// [[Rcpp::depends(RcppEigen)]]
#include <vector>
#include <stdlib.h>
#include <RcppEigen.h>

#include <cmath>
#include <memory>
#include <algorithm>
#include <thread>

using namespace std;

using namespace Eigen;

static VectorXd solveBandedSystem(const VectorXd& diagBase,
								  const VectorXd& lower1Base,
								  const VectorXd& lower2Base,
								  const VectorXd& weights,
								  const VectorXd& rhs)
{
	const int N = diagBase.size();
	VectorXd diag = diagBase + weights;

	VectorXd Ldiag(N);
	VectorXd Llower1 = VectorXd::Zero(std::max(0, N - 1));
	VectorXd Llower2 = VectorXd::Zero(std::max(0, N - 2));

	for (int i = 0; i < N; ++i)
	{
		double diagVal = diag[i];
		if (i - 1 >= 0)
		{
			diagVal -= Llower1[i - 1] * Llower1[i - 1];
		}
		if (i - 2 >= 0 && i - 2 < Llower2.size())
		{
			diagVal -= Llower2[i - 2] * Llower2[i - 2];
		}
		diagVal = std::max(diagVal, 1e-12);
		Ldiag[i] = std::sqrt(diagVal);

		if (i + 1 < N)
		{
			double val = lower1Base[i];
			if (i - 1 >= 0 && i - 1 < Llower2.size())
			{
				val -= Llower2[i - 1] * Llower1[i - 1];
			}
			Llower1[i] = val / Ldiag[i];
		}
		if (i + 2 < N)
		{
			double val = lower2Base[i];
			Llower2[i] = val / Ldiag[i];
		}
	}

	VectorXd y(N);
	for (int i = 0; i < N; ++i)
	{
		double sum = rhs[i];
		if (i - 1 >= 0)
		{
			sum -= Llower1[i - 1] * y[i - 1];
		}
		if (i - 2 >= 0 && i - 2 < Llower2.size())
		{
			sum -= Llower2[i - 2] * y[i - 2];
		}
		y[i] = sum / Ldiag[i];
	}

	VectorXd x(N);
	for (int i = N - 1; i >= 0; --i)
	{
		double sum = y[i];
		if (i + 1 < N)
		{
			sum -= Llower1[i] * x[i + 1];
		}
		if (i + 2 < N && i < Llower2.size())
		{
			sum -= Llower2[i] * x[i + 2];
		}
		x[i] = sum / Ldiag[i];
	}

	return x;
}

// Solve ALS for a single spectrum and return the estimated baseline
static VectorXd computeBaselineForSpectrum(const VectorXd& sig,
										   const VectorXd& diagBase,
										   const VectorXd& lower1Base,
										   const VectorXd& lower2Base,
										   double p,
										   int maxIter)
{
	const int N = sig.size();
	VectorXd vweights = VectorXd::Ones(N);
	VectorXd baseline = VectorXd::Zero(N);
	VectorXd rhs(N);

	for (int iter = 0; iter < maxIter; ++iter)
	{
		rhs = sig.cwiseProduct(vweights);
		baseline = solveBandedSystem(diagBase, lower1Base, lower2Base, vweights, rhs);

		bool changed = false;
		for (int jj = 0; jj < N; ++jj)
		{
			double newWeight = sig[jj] > baseline[jj] ? p : 1.0 - p;
			if (std::abs(newWeight - vweights[jj]) > 1e-6)
			{
				changed = true;
			}
			vweights[jj] = newWeight;
		}
		if (!changed)
		{
			break;
		}
	}

	return baseline;
}

static void buildBandedComponents(
    int N, double smooth,
    VectorXd& diagBase,
    VectorXd& lower1Base,
    VectorXd& lower2Base)
{
    diagBase   = VectorXd::Zero(N);
    lower1Base = VectorXd::Zero(std::max(0, N - 1));
    lower2Base = VectorXd::Zero(std::max(0, N - 2));

    if (N < 3) {
        diagBase.array() *= smooth;
        lower1Base.array() *= smooth;
        lower2Base.array() *= smooth;
        return;
    }

    for (int k = 0; k <= N - 3; ++k) {
        diagBase[k]     += 1.0;
        diagBase[k + 1] += 4.0;
        diagBase[k + 2] += 1.0;

        lower1Base[k]     -= 2.0;
        lower1Base[k + 1] -= 2.0;

        lower2Base[k] += 1.0;
    }

    diagBase   *= smooth;
    lower1Base *= smooth;
    lower2Base *= smooth;
}


static void processSpectraRange(const Rcpp::NumericMatrix& spectra,
								Rcpp::NumericMatrix& corrected,
								int startIdx,
								int endIdx,
								int N,
								const VectorXd& diagBase,
								const VectorXd& lower1Base,
								const VectorXd& lower2Base,
								double p,
								int maxIter)
{
	VectorXd sig(N);
	for (int spectrumIdx = startIdx; spectrumIdx < endIdx; ++spectrumIdx)
	{
		const int offset = spectrumIdx * N;
		for (int ii = 0; ii < N; ++ii)
		{
			sig[ii] = spectra(spectrumIdx, ii);
		}

		VectorXd baseline = computeBaselineForSpectrum(sig, diagBase, lower1Base, lower2Base, p, maxIter);
		for (int ii = 0; ii < N; ++ii)
		{
			corrected(spectrumIdx, ii) = sig[ii] - baseline[ii];
		}
	}
}

// [[Rcpp::export]]
Rcpp::NumericMatrix ALSBaselineCpp(const Rcpp::NumericMatrix spectra,
					double lambda = 6,
                    double p = 0.05,
                    int maxIter = 20,
                    int n_threads = 1)
{
	const int numSpectra = spectra.nrow();
	const int N = spectra.ncol();
	const double smooth = std::pow(10, lambda);

	Rcpp::NumericMatrix corrected(numSpectra, N);

	// Initialize banded D^T D (5 x N representation)
	VectorXd diagBase;
	VectorXd lower1Base;
	VectorXd lower2Base;
	buildBandedComponents(N, smooth, diagBase, lower1Base, lower2Base);

	if (n_threads <= 0)
	{
		n_threads = 1;
	}
	n_threads = std::min(n_threads, numSpectra);

	auto worker = [&](int startSpec, int endSpec)
	{
		processSpectraRange(spectra, corrected, startSpec, endSpec, N, diagBase, lower1Base, lower2Base, p, maxIter);
	};

	if (n_threads == 1 || numSpectra == 1)
	{
		worker(0, numSpectra);
		return corrected;
	}

	std::vector<std::thread> threads;
	threads.reserve(n_threads);
	const int block = (numSpectra + n_threads - 1) / n_threads;
	for (int t = 0; t < n_threads; ++t)
	{
		const int startSpec = t * block;
		if (startSpec >= numSpectra)
		{
			break;
		}
		const int endSpec = std::min(numSpectra, startSpec + block);
		threads.emplace_back(worker, startSpec, endSpec);
	}

	for (auto& th : threads)
	{
		if (th.joinable())
		{
			th.join();
		}
	}

	return corrected;
}
