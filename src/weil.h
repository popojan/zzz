// weil.h -- Weil explicit-formula window fit (zzz --weil).
//
// Refines a window of consecutive zeros around a centre ordinate by
// least-squares fitting their positions to exact Riemann-Weil functionals
// of the primes p^m <= X. Effective above the saturation threshold
// gap * log X > pi (see doc/notes/band-saturation.md); below it the fit
// provably reproduces the seeds (method B is information-optimal there).
// Ported from the validated prototype doc/ghy/e5-above-threshold.wls.
//
// Functional family, centred at t_base, everything window-relative
// (x = gamma - t_base), W(x) = sinc(bw x)^(2 qw), band delta = 2 qw bw,
// gHat = cardinal B-spline of order 2 qw with support [-delta, delta]:
//
//   zero side                          prime side (computable, exact)
//   2 sum_gamma W(x) cos(w x)  =  2 log(t_base/2pi) gHat(w)
//                                 - 2 sum_{p^m<=X} A cos(t_base u) [gHat(u-w)+gHat(u+w)]
//   2 sum_gamma W(x) sin(w x)  =  - 2 sum_{p^m<=X} A sin(t_base u) [gHat(u+w)-gHat(u-w)]
//
// with u = m log p, A = (log p) p^(-m/2). The archimedean term uses the
// psi asymptotics (corrections O(1/t_base^2); fine for t_base >> 10^4);
// pole terms vanish at height. Only the phases t_base*u mod 2pi need
// high precision; data, model, Jacobian run in doubles (target ~1e-7).

#ifndef WEIL_H
#define WEIL_H

#include <flint/arb.h>
#include <flint/ulong_extras.h>

typedef struct {
    double delta;      // band of the window kernel W            [3.0]
    double sig_d;      // data floor weight                      [1e-7]
    double sig_p;      // seed prior width in mean gaps          [0.10]
    slong lm_iters;    // Levenberg-Marquardt rounds             [24]
    int verbose;       // iteration trace to stderr
} weil_opts;

void weil_opts_default(weil_opts *o);

// Refine xs[0..M) (seed offsets from t_base) in place against the Weil
// data of primes p^m <= X; far[0..nfar) are fixed background offsets.
// Caller must have run n_compute_primes() up to X. PREC governs the
// phase folding only. Returns 0 on success.
int weil_fit(double *xs, slong M, const double *far, slong nfar,
             const arb_t t_base, ulong X, const weil_opts *opt, slong PREC);

#endif
