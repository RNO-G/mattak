#include "mattak/fft_utils.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

void fft_inplace(std::vector<Cplx>& a, bool inverse)
{
    int n = a.size();
    for (int i = 1, j = 0; i < n; ++i) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2*M_PI/len * (inverse ? 1 : -1);
        Cplx wlen(std::cos(ang), std::sin(ang));
        for (int i = 0; i < n; i += len) {
            Cplx w(1);
            for (int j = 0; j < len/2; ++j, w *= wlen) {
                Cplx u = a[i+j], v = w * a[i+j+len/2];
                a[i+j]       = u + v;
                a[i+j+len/2] = u - v;
            }
        }
    }
    if (inverse)
        for (auto& x : a) x /= n;
}

std::vector<Cplx> rfft(const std::vector<double>& x)
{
    int n = x.size();
    std::vector<Cplx> a(n);
    for (int i = 0; i < n; ++i) a[i] = x[i];
    fft_inplace(a, false);
    
    return std::vector<Cplx>(a.begin(), a.begin() + n/2 + 1);
}

std::vector<double> irfft(const std::vector<Cplx>& spec, int n)
{
    std::vector<Cplx> a(n);
    a[0]   = spec[0];
    a[n/2] = spec[n/2];
    for (int i = 1; i < n/2; ++i) {
        a[i]   = spec[i];
        a[n-i] = std::conj(spec[i]);
    }
    fft_inplace(a, true);
    std::vector<double> out(n);
    for (int i = 0; i < n; ++i) out[i] = a[i].real();
    return out;
}

// new 

// // ─── fit kernel ──────────────────────────────────────────────────────────────
// #include <cassert>
// #include <cstring>
// #include <functional>

// FitKernel build_fit_kernel(
//     const std::vector<Cplx>& spectrum,
//     int    n,
//     int    block_size,
//     double sampling_rate,
//     double max_frequency)
// {
//     double dt     = 1.0 / sampling_rate;
//     int    n_freq = (int)spectrum.size();   // n/2 + 1
//     int    n_blocks = n / block_size;

//     // collect in-band frequency indices (f > 0 && f < max_frequency)
//     std::vector<int> oob_idx;
//     for (int k = 0; k < n_freq; ++k) {
//         double f = (double)k / (n * dt);
//         // Only include if frequency is in band AND spectrum value is non-zero
//         if (f > 0.0 && f < max_frequency && std::abs(spectrum[k]) > 1e-30)
//             oob_idx.push_back(k);
//     }


//     int n_oob = (int)oob_idx.size();
//     int n_params = n_blocks - 1;   // gauge: last block absorbed into global shift

//     FitKernel K;
//     K.n_params = n_params;
//     K.n_freq   = n_oob;
//     K.term.resize((size_t)n_params * n_oob);
//     K.spectrum_oob.resize(n_oob);

//     // NuRadio FFT normalisation factor
//     double norm = 1.0 / sampling_rate * std::sqrt(2.0);
    
//     norm /= n;  // ← Normalize spectrum by n here

//     for (int j = 0; j < n_params; ++j) {
//         for (int ki = 0; ki < n_oob; ++ki) {
//             int    k = oob_idx[ki];
//             double f = (double)k / (n * dt);

//             // exponent: -2i * pi * f * dt * ((j+0.5)*block_size - 0.5)
//             double phase = -2.0 * M_PI * f * dt
//                            * ((j + 0.5) * block_size - 0.5);
//             Cplx expo(std::cos(phase), std::sin(phase));

//             // sinc-like envelope: sin(pi*f*block_size*dt) / sin(pi*f*dt)
//             double sin_num = std::sin(M_PI * f * block_size * dt);
//             double sin_den = std::sin(M_PI * f * dt);
//             double envelope = (std::abs(sin_den) > 1e-30)
//                               ? sin_num / sin_den : (double)block_size;

//             K.term[(size_t)j * n_oob + ki] = norm * expo * envelope;
//         }
//     }

//     for (int ki = 0; ki < n_oob; ++ki)
//         // K.spectrum_oob[ki] = spectrum[oob_idx[ki]];
//         K.spectrum_oob[ki] = spectrum[oob_idx[ki]] / (double)n;  

//     return K;
// }

// // ─── chi2 + gradient ─────────────────────────────────────────────────────────
// double chi2_and_grad(const FitKernel& K,
//                      const std::vector<double>& a,
//                      std::vector<double>& grad)
// {
//     int np = K.n_params, nf = K.n_freq;
//     assert((int)a.size() == np);
//     grad.assign(np, 0.0);

//     // residual r[k] = sum_j a[j] * term[j,k]  -  spectrum_oob[k]
//     std::vector<Cplx> r(nf, {0.0, 0.0});
//     for (int j = 0; j < np; ++j)
//         for (int ki = 0; ki < nf; ++ki)
//             r[ki] += a[j] * K.term[(size_t)j * nf + ki];
//     for (int ki = 0; ki < nf; ++ki)
//         r[ki] -= K.spectrum_oob[ki];

//     // chi2 = sum |r|^2
//     double chi2 = 0.0;
//     for (int ki = 0; ki < nf; ++ki)
//         chi2 += std::norm(r[ki]);

//     // grad[j] = d(chi2)/d(a[j]) = 2 * Re( sum_k conj(term[j,k]) * r[k] )
//     for (int j = 0; j < np; ++j) {
//         Cplx g{0.0, 0.0};
//         for (int ki = 0; ki < nf; ++ki)
//             g += std::conj(K.term[(size_t)j * nf + ki]) * r[ki];
//         grad[j] = 2.0 * g.real();
//     }

//     return chi2;
// }

// // ─── minimal L-BFGS ──────────────────────────────────────────────────────────
// // Two-loop recursion, Nocedal & Wright §7.4
// std::vector<double> lbfgs(const FitKernel& K,
//                            const std::vector<double>& a0,
//                            double tol, int maxiter, int m)
// {
//     int n = (int)a0.size();
//     std::vector<double> x = a0, g(n), q(n), r(n);
//     double f = chi2_and_grad(K, x, g);

//     // circular buffers
//     std::vector<std::vector<double>> s_buf(m, std::vector<double>(n,0));
//     std::vector<std::vector<double>> y_buf(m, std::vector<double>(n,0));
//     std::vector<double> rho_buf(m, 0.0);
//     int head = 0, filled = 0;

//     auto dot = [&](const std::vector<double>& u, const std::vector<double>& v){
//         double s = 0; for (int i=0;i<n;i++) s+=u[i]*v[i]; return s;
//     };

//     for (int iter = 0; iter < maxiter; ++iter) {
//         // check convergence on gradient norm
//         double gnorm = std::sqrt(dot(g, g));

//         // new, debug
//         if (iter < 5 || iter % 10 == 0) {
//             std::cout << "Iter " << iter << ": f=" << std::setprecision(12) << f 
//                     << ", gnorm=" << gnorm << "\n";
//         }
        
//         if (gnorm < tol) {
//             std::cout << "Converged at iter " << iter << " (gnorm=" << gnorm << ")\n";
//             break;
//         }
//         //if (gnorm < tol) break;

//         // two-loop recursion to get search direction d = -H*g
//         q = g;
//         std::vector<double> alpha_buf(m, 0.0);
//         for (int i = 0; i < filled; ++i) {
//             int idx = (head - 1 - i + m) % m;
//             double a = rho_buf[idx] * dot(s_buf[idx], q);
//             alpha_buf[i] = a;
//             for (int k = 0; k < n; k++) q[k] -= a * y_buf[idx][k];
//         }
//         // initial Hessian: H0 = (s^T y)/(y^T y) * I
//         double scale = 1.0;
//         if (filled > 0) {
//             int idx = (head - 1 + m) % m;
//             scale = dot(s_buf[idx], y_buf[idx]) /
//                     (dot(y_buf[idx], y_buf[idx]) + 1e-30);
//         }
//         r = q;
//         for (int k = 0; k < n; k++) r[k] *= scale;
//         for (int i = filled-1; i >= 0; --i) {
//             int idx = (head - 1 - i + m) % m;
//             double beta = rho_buf[idx] * dot(y_buf[idx], r);
//             for (int k = 0; k < n; k++) r[k] += s_buf[idx][k] * (alpha_buf[i] - beta);
//         }
//         // d = -r
//         std::vector<double> d(n);
//         for (int k = 0; k < n; k++) d[k] = -r[k];

//         // Wolfe line search
//         double step = 1.0;
//         double f0   = f, dg0 = dot(d, g);
//         std::vector<double> x_new(n), g_new(n);
//         double f_new = f;
//         const double c1 = 1e-4, c2 = 0.9;
//         for (int ls = 0; ls < 30; ++ls) {
//             for (int k = 0; k < n; k++) x_new[k] = x[k] + step * d[k];
//             f_new = chi2_and_grad(K, x_new, g_new);
//             double dg_new = dot(d, g_new);
//             if (f_new > f0 + c1 * step * dg0) { step *= 0.5; continue; }
//             if (dg_new < c2 * dg0)             { step *= 2.1; continue; }
//             break;
//         }

//         // update L-BFGS buffers
//         int idx = head % m;
//         for (int k = 0; k < n; k++) s_buf[idx][k] = x_new[k] - x[k];
//         for (int k = 0; k < n; k++) y_buf[idx][k] = g_new[k]  - g[k];
//         double sy = dot(s_buf[idx], y_buf[idx]);
//         rho_buf[idx] = (std::abs(sy) > 1e-30) ? 1.0/sy : 0.0;
//         head = (head + 1) % m;
//         filled = std::min(filled + 1, m);

//         x = x_new; g = g_new; f = f_new;
//     }
//     return x;
// }