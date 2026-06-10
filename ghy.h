// ghy.h -- Gonek-Hughes-Young hybrid product primitives.
//
// Reference: Gonek, Hughes, Young 2007,
//   "A hybrid Euler-Hadamard product for the Riemann zeta function".
//
// Shared across zghy/zhad/zhybrid. All functions zeta-evaluation-free:
// they work from primes and a list of known zero ordinates only.

#ifndef GHY_H
#define GHY_H

#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/ulong_extras.h>

// log P_X(s) = sum_{p^m <= X} 1/(m p^{ms})    (GHY eq. 6)
// Caller must have called n_compute_primes() for primes up to X.
void ghy_log_px(acb_ptr out, const acb_t s, ulong X, slong PREC);

// log Z_X(s) = -sum_j U((s - rho_j) log X)    (GHY eq. 7)
//   rho_j = 1/2 + i gammas[j]
//   U(z) = E_1(z) in the X -> infinity limit (GHY eq. 4 with u -> delta_e).
// n_gammas is the length of the gammas array. Zeros whose |theta|=|t-gamma_j|
// exceed window_units / log X are skipped (window_units=20 gives |U| < 1e-3).
// Pass window_units = 0 to sum over all provided zeros.
void ghy_log_zx(acb_ptr out,
                const acb_t s,
                const double *gammas,
                slong n_gammas,
                ulong X,
                double window_units,
                slong PREC);

// log Z_X, window-relative variant for extreme heights: zero j sits at
// rho_j = 1/2 + i (t_base + dgammas[j]). Offsets fit comfortably in doubles
// even when the absolute ordinate (e.g. ~10^36) does not. skip_j >= 0
// excludes that index from the sum (leave-one-out, used by the bootstrap);
// pass skip_j = -1 to sum over all. window_units as in ghy_log_zx.
void ghy_log_zx_rel(acb_ptr out,
                    const acb_t s,
                    const arb_t t_base,
                    const double *dgammas,
                    slong n_gammas,
                    slong skip_j,
                    ulong X,
                    double window_units,
                    slong PREC);

// N_0(T) = (T/2pi) log(T/2pi e) + 7/8   (smooth Riemann-von Mangoldt part)
void ghy_n0_smooth(arb_ptr out, arb_srcptr t, slong PREC);

#endif
