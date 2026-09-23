#pragma once
#include <vector>
#include <complex>

using Cplx = std::complex<double>;

// FFT primitives
void fft_inplace(std::vector<Cplx>& a, bool inverse);
std::vector<Cplx> rfft(const std::vector<double>& x);
std::vector<double> irfft(const std::vector<Cplx>& spec, int n);

// ── new: precomputed fit kernel ───────────────────────────────────────────────
// const_fft_term[j][k]  (n_blocks-1 rows, n_oob_freq cols)
// populated by build_const_fft_term(); caller owns the storage
// struct FitKernel {
//     int n_params;          // n_blocks - 1
//     int n_freq;            // number of in-band frequencies
//     // flat row-major: index [j * n_freq + k]
//     std::vector<Cplx> term;
//     std::vector<Cplx> spectrum_oob;   // target spectrum (in-band bins only)
// };

// FitKernel build_fit_kernel(
//     const std::vector<Cplx>& spectrum,   // full one-sided spectrum (size n/2+1)
//     int    n,                            // trace length
//     int    block_size,
//     double sampling_rate,
//     double max_frequency);

// // chi2 and gradient w.r.t. a (length n_params)
// double chi2_and_grad(const FitKernel& K,
//                      const std::vector<double>& a,
//                      std::vector<double>& grad);

// // L-BFGS-B minimizer (pure C++, no external deps)
// // returns optimised parameters
// std::vector<double> lbfgs(const FitKernel& K,
//                            const std::vector<double>& a0,
//                            double tol      = 1e-8,
//                            int    maxiter  = 200,
//                            int    m        = 10);