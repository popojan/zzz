// ghy.c -- implementations for ghy.h

#include "ghy.h"
#include <flint/acb_hypgeom.h>
#include <flint/arb_hypgeom.h>
#include <math.h>

void ghy_log_px(acb_ptr out, const acb_t s, ulong X, slong PREC) {
    acb_t total, p_acb, p_neg_s, psm_neg, term, neg_s;

    acb_init(total);
    acb_init(p_acb);
    acb_init(p_neg_s);
    acb_init(psm_neg);
    acb_init(term);
    acb_init(neg_s);

    acb_zero(total);
    acb_neg(neg_s, s);

    for (slong i = 1; ; ++i) {
        ulong p = n_nth_prime(i);
        if (p > X) break;

        acb_set_ui(p_acb, p);
        acb_pow(p_neg_s, p_acb, neg_s, PREC);
        acb_set(psm_neg, p_neg_s);

        ulong pm = p;
        for (ulong m = 1; ; ++m) {
            acb_div_ui(term, psm_neg, m, PREC);
            acb_add(total, total, term, PREC);

            if (pm > X / p) break;
            pm *= p;
            acb_mul(psm_neg, psm_neg, p_neg_s, PREC);
        }
    }

    acb_set(out, total);

    acb_clear(total);
    acb_clear(p_acb);
    acb_clear(p_neg_s);
    acb_clear(psm_neg);
    acb_clear(term);
    acb_clear(neg_s);
}

// E_1(z) via acb_hypgeom_expint with nu = 1.
//
// Purely imaginary z = iy goes through real Si/Ci instead:
//   E_1(iy) = -Ci(|y|) + i (Si(|y|) - pi/2),  conjugated for y < 0.
// The complex E_1 power series cancels ~|y|/log 2 bits in the mid-range
// before the asymptotic regime takes over, so the generic path silently
// returns wide balls there (and on the critical line every Z_X kernel
// argument is purely imaginary, with |y| up to window-span * log X).
static void e1(acb_ptr out, const acb_t z, slong PREC) {
    if (arb_is_zero(acb_realref(z))) {
        arb_t y, si, ci, half_pi;
        int neg = arb_is_negative(acb_imagref(z));

        arb_init(y);
        arb_init(si);
        arb_init(ci);
        arb_init(half_pi);

        arb_abs(y, acb_imagref(z));
        arb_hypgeom_si(si, y, PREC);
        arb_hypgeom_ci(ci, y, PREC);
        arb_const_pi(half_pi, PREC);
        arb_mul_2exp_si(half_pi, half_pi, -1);
        arb_sub(si, si, half_pi, PREC);
        if (neg) arb_neg(si, si);

        arb_neg(acb_realref(out), ci);
        arb_set(acb_imagref(out), si);

        arb_clear(y);
        arb_clear(si);
        arb_clear(ci);
        arb_clear(half_pi);
        return;
    }
    acb_t nu;
    acb_init(nu);
    acb_one(nu);
    acb_hypgeom_expint(out, nu, z, PREC);
    acb_clear(nu);
}

void ghy_log_zx(acb_ptr out,
                const acb_t s,
                const double *gammas,
                slong n_gammas,
                ulong X,
                double window_units,
                slong PREC)
{
    acb_t total, rho, diff, arg_z, kernel, log_X_c;
    arb_t half, log_X, t_arb, gamma_arb;
    double log_X_d = log((double)X);
    double t_d = arf_get_d(arb_midref(acb_imagref(s)), ARF_RND_NEAR);

    acb_init(total);
    acb_init(rho);
    acb_init(diff);
    acb_init(arg_z);
    acb_init(kernel);
    acb_init(log_X_c);
    arb_init(half);
    arb_init(log_X);
    arb_init(t_arb);
    arb_init(gamma_arb);

    acb_zero(total);

    arb_set_d(half, 0.5);
    arb_set_ui(log_X, X);
    arb_log(log_X, log_X, PREC);
    acb_set_arb(log_X_c, log_X);

    for (slong j = 0; j < n_gammas; ++j) {
        double theta_d = t_d - gammas[j];
        if (window_units > 0.0 &&
            fabs(theta_d) * log_X_d > window_units) {
            continue;
        }

        // rho = 1/2 + i gamma_j
        arb_set_d(gamma_arb, gammas[j]);
        acb_set_arb_arb(rho, half, gamma_arb);

        // diff = s - rho
        acb_sub(diff, s, rho, PREC);

        // arg_z = (s - rho) * log X
        acb_mul(arg_z, diff, log_X_c, PREC);

        // kernel = E_1(arg_z)
        e1(kernel, arg_z, PREC);

        // total += kernel  (so that log Z_X = -total later)
        acb_add(total, total, kernel, PREC);
    }

    // log Z_X = -sum
    acb_neg(out, total);

    acb_clear(total);
    acb_clear(rho);
    acb_clear(diff);
    acb_clear(arg_z);
    acb_clear(kernel);
    acb_clear(log_X_c);
    arb_clear(half);
    arb_clear(log_X);
    arb_clear(t_arb);
    arb_clear(gamma_arb);
}

void ghy_log_zx_rel(acb_ptr out,
                    const acb_t s,
                    const arb_t t_base,
                    const double *dgammas,
                    slong n_gammas,
                    slong skip_j,
                    ulong X,
                    double window_units,
                    slong PREC)
{
    acb_t total, diff, arg_z, kernel, log_X_c;
    arb_t log_X, t_rel, half, dg;
    double log_X_d = log((double)X);
    double t_rel_d;

    acb_init(total);
    acb_init(diff);
    acb_init(arg_z);
    acb_init(kernel);
    acb_init(log_X_c);
    arb_init(log_X);
    arb_init(t_rel);
    arb_init(half);
    arb_init(dg);

    acb_zero(total);

    arb_set_ui(log_X, X);
    arb_log(log_X, log_X, PREC);
    acb_set_arb(log_X_c, log_X);

    // t_rel = Im(s) - t_base: small in the bootstrap window, exact in arb
    arb_sub(t_rel, acb_imagref(s), t_base, PREC);
    t_rel_d = arf_get_d(arb_midref(t_rel), ARF_RND_NEAR);

    // Re(s - rho) = Re(s) - 1/2, shared by all terms
    arb_set_d(half, 0.5);
    arb_sub(acb_realref(diff), acb_realref(s), half, PREC);

    for (slong j = 0; j < n_gammas; ++j) {
        if (j == skip_j) continue;
        if (window_units > 0.0 &&
            fabs(t_rel_d - dgammas[j]) * log_X_d > window_units) {
            continue;
        }

        // Im(s - rho_j) = t_rel - dgamma_j
        arb_set_d(dg, dgammas[j]);
        arb_sub(acb_imagref(diff), t_rel, dg, PREC);

        acb_mul(arg_z, diff, log_X_c, PREC);
        e1(kernel, arg_z, PREC);
        acb_add(total, total, kernel, PREC);
    }

    acb_neg(out, total);

    acb_clear(total);
    acb_clear(diff);
    acb_clear(arg_z);
    acb_clear(kernel);
    acb_clear(log_X_c);
    arb_clear(log_X);
    arb_clear(t_rel);
    arb_clear(half);
    arb_clear(dg);
}

void ghy_n0_smooth(arb_ptr out, arb_srcptr t, slong PREC) {
    arb_t x, a, b;
    arb_init(x); arb_init(a); arb_init(b);

    arb_set(x, t);
    arb_const_pi(b, PREC);
    arb_div_ui(x, x, 2, PREC);
    arb_div(x, x, b, PREC);
    arb_set(a, x);
    arb_const_e(b, PREC);
    arb_div(a, a, b, PREC);
    arb_log(a, a, PREC);
    arb_mul(x, a, x, PREC);
    arb_set_d(a, 0.875);
    arb_add(x, x, a, PREC);
    arb_set(out, x);

    arb_clear(x); arb_clear(a); arb_clear(b);
}
