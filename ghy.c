// ghy.c -- implementations for ghy.h

#include "ghy.h"
#include <flint/acb_hypgeom.h>
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

// E_1(z) via acb_hypgeom_expint with nu = 1
static void e1(acb_ptr out, const acb_t z, slong PREC) {
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
