#include <flint/arb.h>
#include <flint/acb.h>
#include "flint/ulong_extras.h"
#include <stdlib.h>
#include <math.h>
#include "argp.h"
#include "ghy.h"

#define CUBIC 1

void coord(arb_ptr out, arb_srcptr q, arb_srcptr t, slong PREC) {
    arb_t pi;
    arb_t pi_div_log_q;
    arb_t x;
    arb_t a;
    arb_t b;

    arb_init(pi);
    arb_init(pi_div_log_q);
    arb_init(x);
    arb_init(a);
    arb_init(b);

    // PI / Log[q]
    arb_const_pi(pi, PREC);
    arb_log(pi_div_log_q, q, PREC);
    arb_inv(pi_div_log_q, pi_div_log_q, PREC);
    arb_mul(pi_div_log_q, pi_div_log_q, pi, PREC);

    // -t + PI / Log[q]
    arb_sub(x, pi_div_log_q, t, PREC);

    //Mod[-t + PI / Log[q], 2 PI / Log[q]]
    arb_add(a, pi_div_log_q, pi_div_log_q, PREC);
    arb_div(b, x, a, PREC);
    arb_floor(b, b, PREC);
    arb_mul(b, b, a, PREC);
    arb_sub(x, x, b, PREC);

    //% - PI/Log[q]
    arb_sub(x, x, pi_div_log_q, PREC);
    arb_set(out, x);

    arb_clear(pi);
    arb_clear(pi_div_log_q);
    arb_clear(x);
    arb_clear(a);
    arb_clear(b);
}

void wave2(arb_ptr out, arb_srcptr q, arb_srcptr t, slong PREC) {
    arb_t pi;
    arb_t sqrt_q;
    arb_t log_q;
    arb_t x;
    arb_t a;
    arb_t b;

    arb_init(pi);
    arb_init(x);
    arb_init(a);
    arb_init(b);
    arb_init(sqrt_q);
    arb_init(log_q);

    arb_const_pi(pi, PREC);
    arb_log(log_q, q, PREC);
    arb_sqrt(sqrt_q, q, PREC);


    //4 sqrt[q] log[q] t
    arb_mul_ui(a, t, 4,PREC);
    arb_mul(b, log_q, sqrt_q, PREC);
    arb_mul(a, a, b, PREC);

    //PI (sqrt[q]-1)^2
    arb_one(b);
    arb_sub(b, sqrt_q, b, PREC);
    arb_mul(b, b, b, PREC);
    arb_mul(b, b, pi, PREC);

    //(%1 + %2)
    arb_add(a, a, b, PREC);

    //q Log[q]^2
    arb_mul(b, q, log_q, PREC);
    arb_mul(b, b, log_q, PREC);

    //sqrt(%/%)
    arb_div(a, a, b, PREC);
    arb_sqrt(a, a, PREC);

    // (sqrt[PI] sqrt[q] sqrt(%/%) - 2 t) Log[q]
    arb_sqrt(b, pi, PREC);
    arb_mul(b, b, sqrt_q, PREC);
    arb_mul(b, b, a, PREC);
    arb_sub(b, b, t, PREC);
    arb_sub(b, b, t, PREC);
    arb_mul(a, b, log_q, PREC);

    arb_add(b, pi, pi, PREC);
    arb_div(a, a, b, PREC);

    arb_div_ui(b, sqrt_q, 2, PREC);
    arb_sub(a, a, b, PREC);
    arb_set_d(b, 0.5);
    arb_add(a, a, b, PREC);

    arb_set(out, a);

    arb_clear(pi);
    arb_clear(x);
    arb_clear(a);
    arb_clear(b);
    arb_clear(sqrt_q);
    arb_clear(log_q);
}


void wave3(arb_ptr out, arb_srcptr q, arb_srcptr t, slong PREC) {
    arb_t pi;
    arb_t sqrt_q;
    arb_t sqrt_q3;
    arb_t log_q;
    arb_t x;
    arb_t y;
    arb_t a;
    arb_t b;
    arb_t c;
    arb_t d;

    arb_init(pi);
    arb_init(x);
    arb_init(y);
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    arb_init(sqrt_q);
    arb_init(sqrt_q3);
    arb_init(log_q);

    arb_sqrt(sqrt_q, q, PREC);

    arb_mul(sqrt_q3, q, q, PREC);
    arb_mul(sqrt_q3, sqrt_q3, q, PREC);
    arb_sqrt(sqrt_q3, sqrt_q3, PREC);

    arb_log(log_q, q, PREC);
    arb_const_pi(pi, PREC);

    arb_mul(a, t, log_q, PREC);
    arb_neg(a, a);
    arb_add(a, a, pi, PREC);

    arb_add(b, sqrt_q3, sqrt_q3, PREC);
    arb_mul(a, a, b, PREC);

    arb_add(b, q, q, PREC);
    arb_add(b, b, q, PREC);
    arb_mul(b, b, pi, PREC);
    arb_neg(b, b);
    arb_add(b,b, pi, PREC);
    arb_add(a, a, b, PREC);

    arb_mul(b, t, log_q, PREC);
    arb_mul_ui(c, sqrt_q, 8, PREC);
    arb_mul(b, b, c, PREC);
    arb_mul(a, a, b, PREC);
    //a == -8 sqrt[q]...

    arb_mul(b, pi, pi, PREC);
    arb_one(c);
    arb_sub(c, sqrt_q, c, PREC);
    arb_mul(d, c, c, PREC);
    arb_mul(d, d, c, PREC);
    arb_mul(b, b, d, PREC);
    //b == pi^2 ( sqrt[q]-1)^3

    arb_mul_ui(c, sqrt_q, 5,PREC);
    arb_add_ui(c, c, 3, PREC);
    arb_mul(b, b, c, PREC);
    arb_sub(a, b, a, PREC);
    //a == pi^2..... (...Log[q]), numerator

    arb_mul(b, q, q, PREC);
    arb_div(a, a, b, PREC);
    arb_sqrt(a, a, PREC);
    //a == sqrt(numerator / q^2)

    arb_mul(c, t, log_q, PREC);
    arb_mul_ui(b, c, 4,PREC);
    arb_add(a, a, b, PREC);

    arb_mul_ui(b, q, 3, PREC);
    arb_one(c);
    arb_sub(b, b, c, PREC);
    arb_div(b, b, sqrt_q3, PREC);
    arb_sub_ui(b, b, 2, PREC);
    arb_mul(b, b, pi, PREC);
    arb_add(a, a, b, PREC);
    //a==a

    arb_root_ui(pi, pi, 3, PREC);
    arb_root_ui(b, a, 3, PREC);

    arb_one(d);
    arb_sub(d, sqrt_q, d, PREC);
    arb_mul(d, d, sqrt_q, PREC);
    arb_mul(d, d, b, PREC);
    arb_mul(d, d, pi, PREC);

    arb_one(x);
    arb_sub(x, q, x, PREC);
    arb_mul(x, x, pi, PREC);
    arb_mul(x, x, pi, PREC);
    arb_sub(x, x, d, PREC);
    arb_mul(d, q, b, PREC);
    arb_mul(d, d, b, PREC);
    arb_sub(x, x, d, PREC);

    arb_one(d);
    arb_add(d, sqrt_q, d, PREC);
    arb_mul(d, d, sqrt_q, PREC);
    arb_mul(d, d, b, PREC);
    arb_mul(d, d, pi, PREC);

    arb_one(y);
    arb_sub(y, q, y, PREC);
    arb_mul(pi, pi, pi, PREC);
    arb_mul(y, y, pi, PREC);
    arb_add(y, y, d, PREC);
    arb_mul(d, q, b, PREC);
    arb_mul(d, d, b, PREC);
    arb_sub(y, y, d, PREC);

    arb_mul(x, x, y, PREC);
    arb_set_d(y, -0.375);
    arb_mul(x, x, y, PREC);
    arb_mul(y, pi, b, PREC);
    arb_mul(y, y, b, PREC);
    //arb_set_d(b, 2.5);
    //arb_pow(b, q, b, PREC);
    arb_pow_ui(b, q, 5, PREC);
    arb_sqrt(b, b, PREC);
    arb_mul(y, y, b, PREC);
    arb_div(x, x, y, PREC);
    arb_set(out, x);

    arb_clear(pi);
    arb_clear(x);
    arb_clear(y);
    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(sqrt_q);
    arb_clear(sqrt_q3);
    arb_clear(log_q);
}

void wave_complex(arb_ptr out, arb_srcptr q, arb_srcptr t, slong PREC) {
     //-((I Log[1 - p^(-(1/2) + I x)])/\[Pi])
    acb_t x;
    acb_t y;
    acb_t z;
    arb_t a;
    arb_init(a);

    acb_init(x);
    acb_init(y);
    acb_init(z);

    arb_const_pi(a, PREC);
    arb_inv(a, a, PREC);
    acb_onei(x);
    acb_neg(x, x);
    acb_mul_arb(x, x, a, PREC); // -I/Pi

    acb_one(y);
    acb_div_ui(y, y, 2, PREC);
    acb_neg(y, y);

    acb_onei(z);
    acb_mul_arb(z, z, t, PREC);
    acb_add(z, z, y, PREC);

    acb_set_arb(y, q);
    acb_pow(y, y, z, PREC);
    acb_neg(y, y);

    acb_add_ui(y, y, 1, PREC);
    acb_log(y, y, PREC);
    acb_mul(y, y, x, PREC);

    arb_set(out, acb_real_ptr(y));

    arb_clear(a);
    acb_clear(z);
    acb_clear(y);
    acb_clear(x);
}

void wave_complex_opt(acb_ptr out, arb_srcptr q, arb_srcptr t, slong PREC) {
     // Log[1 - p^(-(1/2) + I x)]
    acb_t y;
    acb_t z;

    acb_init(y);
    acb_init(z);

    acb_one(y);
    acb_div_ui(y, y, 2, PREC);
    acb_neg(y, y);

    acb_onei(z);
    acb_mul_arb(z, z, t, PREC);
    acb_add(z, z, y, PREC);

    acb_set_arb(y, q);
    acb_pow(y, y, z, PREC);
    acb_neg(y, y);

    acb_add_ui(y, y, 1, PREC);
    acb_log(out, y, PREC);

    acb_clear(z);
    acb_clear(y);
}

void zeta(acb_ptr out, arb_srcptr t, slong ZETA_PREC) {
    acb_t s;
    acb_t z;
    acb_init(s);
    acb_init(z);
    arb_t a;
    arb_init(a);
    arb_set_d(a, 0.5);
    acb_set_arb_arb(s, a, t);
    acb_zeta(z, s, ZETA_PREC);
    acb_set(out, z);

    acb_clear(s);
    acb_clear(z);
}

void nt(arb_ptr out, arb_srcptr t, slong PREC) {
    arb_t x;
    arb_t a;
    arb_t b;
    arb_init(x);
    arb_init(a);
    arb_init(b);

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

    arb_clear(x);
    arb_clear(a);
    arb_clear(b);
}

void sawp(arb_ptr out, arb_srcptr q, arb_srcptr x, slong PREC) {
    arb_t a;
    arb_t b;
    arb_t c;
    arb_t d;
    arb_t xx;
    arb_t pi;
    arb_t sqrt_q;
    arb_t log_q;

    arb_init(pi);
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    arb_init(xx);
    arb_init(sqrt_q);
    arb_init(log_q);

    arb_sqrt(sqrt_q, q, PREC);
    arb_log(log_q, q, PREC);
    arb_const_pi(pi, PREC);

    //a == pi p - pi sqrt@p
    arb_mul(a, pi, q, PREC);
    arb_mul(b, pi, sqrt_q, PREC);
    arb_sub(a, a, b, PREC);

    //Mod[x, 2 PI / Log[q], -a/(2 p Log[p])]
    arb_mul_ui(b, pi, 2, PREC);
    arb_div(b, b, log_q, PREC);

    arb_div(c, a, log_q, PREC);
    arb_div(c, c, q, PREC);
    arb_div_ui(d, c, 2, PREC);
    arb_add(xx, x, d, PREC);

    arb_div(c, xx, b, PREC);
    arb_floor(c, c, PREC);
    arb_mul(c, c, b, PREC);
    arb_sub(c, xx, c, PREC);
    arb_sub(c, c, d, PREC);

    arb_mul(b, c, log_q, PREC);
    arb_mul(b, b, sqrt_q, PREC);
    arb_div(b, b, a, PREC);

    arb_div(a, pi, log_q, PREC);
    arb_sub(a, a, c, PREC);
    arb_mul(a, a, log_q, PREC);
    arb_div(a, a, pi, PREC);
    arb_sub_ui(c, q, 1, PREC);
    arb_div(a, a, c, PREC);
    arb_sub_ui(c, sqrt_q, 1, PREC);
    arb_mul(a, a, c, PREC);
    arb_min(a, a, b, PREC);
    arb_neg(a, a);
    arb_set(out, a);

    arb_clear(a);
    arb_clear(b);
    arb_clear(pi);
    arb_clear(c);
    arb_clear(d);
    arb_clear(xx);
    arb_clear(sqrt_q);
    arb_clear(log_q);

}

void nt_inv(arb_ptr out, arb_srcptr m, slong PREC) {
    arb_t x;
    arb_t a;
    arb_t b;
    arb_t nom;
    arb_t den;
    arb_init(x);
    arb_init(a);
    arb_init(b);
    arb_init(den);
    arb_init(nom);

    arb_mul_ui(x, m, 8, PREC);
    arb_sub_ui(x, x, 11, PREC);

    arb_set(nom, x);
    arb_const_pi(b, PREC);
    arb_mul(nom, nom, b, PREC);

    arb_set(den, x);
    arb_const_e(b, PREC);
    arb_div(den, den, b, PREC);
    arb_div_ui(den, den, 8, PREC);
    arb_lambertw( den, den, 0, PREC);
    arb_mul_ui(den, den, 4, PREC);

    arb_set(x, nom);
    arb_div(x, x, den, PREC);

    arb_set(out, x);

    arb_clear(x);
    arb_clear(a);
    arb_clear(b);
    arb_clear(den);
    arb_clear(nom);
}

void zero_count_exact(arb_ptr out,  arb_srcptr t, slong k, slong PREC) {
    arb_t a;
    acb_t x;
    acb_t u;
    arb_t q;
    acb_t z;
    arb_t att;

    arb_init(att);
    arb_init(a);
    acb_init(x);
    acb_init(u);
    arb_init(q);
    acb_init(z);

    acb_zero(u);

    for(int i = 1; i <= k; ++i) {
        arb_set_ui(q, n_nth_prime(i));
        wave_complex_opt(z, q, t, PREC);

        arb_div(att, t, q, PREC);
        arb_sqrt(att,att,PREC);
        arb_neg(att, att);
        arb_exp(att, att, PREC);
        arb_sub_ui(att, att, 1, PREC);
        arb_neg(att, att);

        acb_mul_arb(z, z, att, PREC);
        acb_add(u, u, z, PREC);
    }

    arb_const_pi(a, PREC);
    arb_inv(a, a, PREC);
    acb_onei(x);
    acb_neg(x, x);
    acb_mul_arb(x, x, a, PREC); // -I/Pi

    acb_mul(u, u, x, PREC);

    nt(a, t, PREC);
    acb_add_arb(u, u, a, PREC);

    arb_set(out, acb_real_ptr(u));

    arb_clear(a);
    acb_clear(x);
    acb_clear(u);
    arb_clear(q);
    acb_clear(z);
    arb_clear(att);
}

// GHY partial-Euler-only counting function: F_B(T) = N_0(T) + arg P_X(1/2+iT) / pi
// where P_X is the Gonek-Hughes-Young partial Euler factor. Zeta-evaluation-free
// and primes-only, like zero_count_exact, but with the heuristic damping replaced
// by the rigorous GHY smooth cutoff (Thm 1 of Gonek-Hughes-Young 2007).
// Parameter X is the prime-power cutoff; typically pass X = p_k to match the
// "k primes" specification of zero_count_exact.
void zero_count_ghy(arb_ptr out, arb_srcptr t, ulong X, slong PREC) {
    acb_t s;
    acb_t logp;
    arb_t pi;
    arb_t n0;

    acb_init(s);
    acb_init(logp);
    arb_init(pi);
    arb_init(n0);

    arb_const_pi(pi, PREC);

    // s = 1/2 + i t
    arb_one(acb_realref(s));
    arb_div_ui(acb_realref(s), acb_realref(s), 2, PREC);
    arb_set(acb_imagref(s), t);

    ghy_log_px(logp, s, X, PREC);

    ghy_n0_smooth(n0, t, PREC);
    arb_div(out, acb_imagref(logp), pi, PREC);
    arb_add(out, out, n0, PREC);

    acb_clear(s);
    acb_clear(logp);
    arb_clear(pi);
    arb_clear(n0);
}

void zero_count_approx(arb_ptr out, arb_srcptr t, slong k, slong PREC) {
    arb_t u;
    arb_t w;
    arb_t q;
    arb_t x;
    arb_t y;
    arb_t z;
    arb_t v;
    arb_t s;
    arb_t att;
    int neg;

    arb_init(u);
    arb_init(w);
    arb_init(q);
    arb_init(y);
    arb_init(z);
    arb_init(v);
    arb_init(s);
    arb_init(att);
    arb_init(x);

    nt(u, t, PREC);

#if CUBIC == 1
    //v = 1/sqrt(2)
    arb_set_d(v, 2.0);
    arb_sqrt(v, v, PREC);
    arb_inv(v, v, PREC);
#endif

    for(int i = 1; i <= k; ++i) {
        arb_set_ui(q, n_nth_prime(i));

        arb_div(att, t, q, PREC);
        arb_sqrt(att,att,PREC);
        arb_neg(att, att);
        arb_exp(att, att, PREC);
        arb_one(w);
        arb_sub(att, w, att, PREC);

        coord(x, q, t, PREC);
        arb_abs(w, x);

        wave2(z, q, w, PREC);

        neg = arb_is_negative(x);
        if(neg) {
            arb_neg(z, z);
        }
#if CUBIC == 1
        sawp(x, q, x, PREC);
        arb_abs(x, x);
        arb_sqrt(y, q, PREC);
        arb_mul(x, x, y, PREC);
        arb_set_ui(y, i);
        arb_inv(y, y, PREC);
        arb_mul_ui(y, y, 2, PREC);
        arb_div_ui(y, y, 3, PREC);
        arb_neg(y, y);
        arb_add_ui(y, y, 1, PREC);
        arb_mul(x, x, y, PREC);
        arb_one(y);
        arb_div_ui(y, y, 2, PREC);
        arb_sub(x, y, x, PREC);

        arb_set_d(y, 0.5);
        arb_sub(y, y, x, PREC);
        arb_mul(z, z, y, PREC);
        arb_mul(z, z, att, PREC);
        arb_add(u, u, z, PREC);

        wave3(z, q, w, PREC);
        if(neg) {
            arb_neg(z, z);
        }
        arb_set_d(y, 0.5);
        arb_add(y, y, x, PREC);
        arb_mul(z, z, y, PREC);
#endif
        arb_mul(z, z, att, PREC);
        arb_add(u, u, z, PREC);
    }
    arb_set(out, u);

    arb_clear(u);
    arb_clear(w);
    arb_clear(q);
    arb_clear(y);
    arb_clear(z);
    arb_clear(v);
    arb_clear(s);
    arb_clear(att);
    arb_clear(x);
}

// --- GHY hybrid bootstrap (method C with self-computed seed zeros) ----------
//
// Pass 1 locates 2W+1 consecutive zeros around the target with the primes-only
// counter F_B = N_0 + arg P_X / pi. These approximate ordinates then seed the
// GHY local Hadamard factor Z_X, and every zero is re-located with the
// leave-one-out hybrid count
//   F_C^(j)(T) = N_0(T) + (arg P_X + arg Z_X^{excl j})(1/2 + iT) / pi
// Gauss-Seidel style for a configurable number of rounds. Z_X subtracts the
// neighbours' truncation oscillation ~cos((T-gamma_j) log X)/((T-gamma_j) log X)
// -- the dominant F_B location error -- while the excluded target's own smeared
// step stays centred at its true ordinate. Zeta-evaluation-free throughout;
// the seeds come from the method itself, not from precomputed zero tables.

typedef struct {
    int method;               // 1 = B (P_X only), 2 = C leave-one-out hybrid
    slong k;                  // prime count (method A)
    ulong X;                  // prime-power cutoff for P_X / Z_X
    const arb_struct *t_base; // C: window base ordinate
    const double *dgam;       // C: zero offsets from t_base
    slong n_gam;
    slong skip;               // C: index left out of Z_X, -1 for none
} zc_cfg;

static void zero_count_eval(arb_ptr out, arb_srcptr t, const zc_cfg *cfg, slong PREC) {
    if (cfg->method == 0) {
        zero_count_exact(out, t, cfg->k, PREC);
        return;
    }
    zero_count_ghy(out, t, cfg->X, PREC);
    if (cfg->method == 2) {
        acb_t s, logz;
        arb_t pi;

        acb_init(s);
        acb_init(logz);
        arb_init(pi);

        arb_one(acb_realref(s));
        arb_div_ui(acb_realref(s), acb_realref(s), 2, PREC);
        arb_set(acb_imagref(s), t);

        ghy_log_zx_rel(logz, s, cfg->t_base, cfg->dgam, cfg->n_gam,
                       cfg->skip, cfg->X, 0.0, PREC);

        arb_const_pi(pi, PREC);
        arb_div(acb_imagref(logz), acb_imagref(logz), pi, PREC);
        arb_add(out, out, acb_imagref(logz), PREC);

        acb_clear(s);
        acb_clear(logz);
        arb_clear(pi);
    }
}

// Bisect F(T) = m_half starting from guess t0 with half-width w0, doubling the
// bracket up to widen_max times if the initial window misses. Returns 0 on
// success. Pass widen_max = 1 when a far-away crossing would be worse than no
// answer (bootstrap relocations must not hop into a neighbour's basin).
static int bisect_zero(arb_t out_t, const arb_t m_half, const arb_t t0,
                       double w0, double tol, int widen_max,
                       const zc_cfg *cfg, slong PREC)
{
    arb_t lo_t, hi_t, mid_t, lo, hi, mid, step, width, tolb;
    int bracketed = 0;

    arb_init(lo_t);
    arb_init(hi_t);
    arb_init(mid_t);
    arb_init(lo);
    arb_init(hi);
    arb_init(mid);
    arb_init(step);
    arb_init(width);
    arb_init(tolb);

    arb_set_d(step, w0);
    arb_set_d(tolb, tol);

    for (int widen = 0; widen < widen_max && !bracketed; ++widen) {
        arb_sub(lo_t, t0, step, PREC);
        arb_add(hi_t, t0, step, PREC);
        zero_count_eval(lo, lo_t, cfg, PREC);
        zero_count_eval(hi, hi_t, cfg, PREC);
        if (!(arb_gt(lo, m_half) || arb_lt(hi, m_half) || arb_lt(hi, lo)))
            bracketed = 1;
        else
            arb_mul_ui(step, step, 2, PREC);
    }

    if (bracketed) {
        while (1) {
            arb_add(mid_t, lo_t, hi_t, PREC);
            arb_mul_2exp_si(mid_t, mid_t, -1);
            zero_count_eval(mid, mid_t, cfg, PREC);
            if (arb_gt(mid, m_half))
                arb_set(hi_t, mid_t);
            else
                arb_set(lo_t, mid_t);
            arb_sub(width, hi_t, lo_t, PREC);
            if (arb_lt(width, tolb)) break;
        }
        arb_set(out_t, mid_t);
    }

    arb_clear(lo_t);
    arb_clear(hi_t);
    arb_clear(mid_t);
    arb_clear(lo);
    arb_clear(hi);
    arb_clear(mid);
    arb_clear(step);
    arb_clear(width);
    arb_clear(tolb);

    return bracketed ? 0 : 1;
}

// Locate zero #m0 via the self-consistent leave-one-out hybrid. W neighbours
// each side, `rounds` relocation rounds. Returns 0 on success.
static int boot_locate(arb_t out_t, const arb_t m0, slong W, slong rounds,
                       ulong X, double w0, double tol, slong PREC, int verbose)
{
    slong lo_off = -W;
    {
        arb_t lim;
        arb_init(lim);
        for (; lo_off < 0; ++lo_off) {
            arb_set_si(lim, 1 - lo_off);
            if (!arb_lt(m0, lim)) break;    // m0 >= 1 - lo_off: window fits
        }
        arb_clear(lim);
    }
    slong n = W - lo_off + 1;
    slong ci = -lo_off;

    double *dgam = malloc(n * sizeof(double));
    arb_struct *tj = malloc(n * sizeof(arb_struct));
    arb_t t_base, m_half, mj, guess, tmp;
    int fail = 0;

    for (slong i = 0; i < n; ++i) arb_init(tj + i);
    arb_init(t_base);
    arb_init(m_half);
    arb_init(mj);
    arb_init(guess);
    arb_init(tmp);

    zc_cfg cfg = { 0 };
    cfg.method = 1;
    cfg.X = X;
    cfg.skip = -1;

    // pass 1: primes-only F_B for every window ordinal
    for (slong i = 0; i < n && !fail; ++i) {
        arb_add_si(mj, m0, lo_off + i, PREC);
        arb_set_d(tmp, 0.5);
        arb_sub(m_half, mj, tmp, PREC);
        nt_inv(guess, mj, PREC);
        fail = bisect_zero(tj + i, m_half, guess, w0, tol, 6, &cfg, PREC);
    }

    if (!fail) {
        arb_set(t_base, tj + ci);
        for (slong i = 0; i < n; ++i) {
            arb_sub(tmp, tj + i, t_base, PREC);
            dgam[i] = arf_get_d(arb_midref(tmp), ARF_RND_NEAR);
        }
        if (verbose) {
            flint_fprintf(stderr, "boot pass 1: %wd zeros seeded around target\n", n);
            if (verbose > 1) {
                flint_fprintf(stderr, "seed offsets:");
                for (slong i = 0; i < n; ++i)
                    flint_fprintf(stderr, " %.4f", dgam[i]);
                flint_fprintf(stderr, "\n");
            }
        }

        // bracket the relocations by the local mean gap 2 pi / log(t / 2 pi)
        double t_d = arf_get_d(arb_midref(t_base), ARF_RND_NEAR);
        double wloc = 0.75 * 2.0 * M_PI / log(t_d / (2.0 * M_PI));
        if (wloc < 16.0 * tol) wloc = 16.0 * tol;

        cfg.method = 2;
        cfg.t_base = t_base;
        cfg.dgam = dgam;
        cfg.n_gam = n;

        // Only the inner core gets relocated; the outer guard ring keeps its
        // P_X seeds. Relocating window-edge zeros is biased (one-sided
        // neighbour coverage) and the bias would propagate inward with the
        // rounds. Guards still serve as Z_X subtraction terms, where their
        // seed error is second order. Jacobi (snapshot) updates keep the
        // rounds free of sweep-direction artifacts.
        slong Wc = (W + 1) / 2;
        double *dnew = malloc(n * sizeof(double));

        for (slong r = 0; r < rounds && !fail; ++r) {
            double max_shift = 0.0;
            slong rejected = 0;
            for (slong i = 0; i < n; ++i) {
                dnew[i] = dgam[i];
                if (i - ci > Wc || ci - i > Wc) continue;
                cfg.skip = i;
                arb_add_si(mj, m0, lo_off + i, PREC);
                arb_set_d(tmp, 0.5);
                arb_sub(m_half, mj, tmp, PREC);
                arb_set_d(tmp, dgam[i]);
                arb_add(guess, t_base, tmp, PREC);
                // no widening, and reject basin hops: a relocation that does
                // not converge near the seed keeps the seed (B quality)
                if (bisect_zero(tj + i, m_half, guess, wloc, tol, 1, &cfg, PREC)) {
                    ++rejected;
                    continue;
                }
                arb_sub(tmp, tj + i, t_base, PREC);
                double nd = arf_get_d(arb_midref(tmp), ARF_RND_NEAR);
                if (fabs(nd - dgam[i]) > 0.45 * wloc) {
                    ++rejected;
                    continue;
                }
                dnew[i] = nd;
                if (fabs(nd - dgam[i]) > max_shift)
                    max_shift = fabs(nd - dgam[i]);
            }
            for (slong i = 0; i < n; ++i) dgam[i] = dnew[i];
            if (verbose) {
                flint_fprintf(stderr, "boot round %wd: max shift %.3e, %wd rejected\n",
                              r + 1, max_shift, rejected);
                if (verbose > 1) {
                    flint_fprintf(stderr, "offsets:");
                    for (slong i = 0; i < n; ++i)
                        flint_fprintf(stderr, " %.4f", dgam[i]);
                    flint_fprintf(stderr, "\n");
                }
            }
        }

        free(dnew);

        // final centre relocation with the converged neighbour set; fall back
        // to the last accepted estimate if the bracket misses
        if (!fail) {
            cfg.skip = ci;
            arb_set_d(tmp, 0.5);
            arb_sub(m_half, m0, tmp, PREC);
            arb_set_d(tmp, dgam[ci]);
            arb_add(guess, t_base, tmp, PREC);
            if (bisect_zero(out_t, m_half, guess, wloc, tol, 1, &cfg, PREC)) {
                arb_set_d(tmp, dgam[ci]);
                arb_add(out_t, t_base, tmp, PREC);
            }
        }
    }

    for (slong i = 0; i < n; ++i) arb_clear(tj + i);
    free(tj);
    free(dgam);
    arb_clear(t_base);
    arb_clear(m_half);
    arb_clear(mj);
    arb_clear(guess);
    arb_clear(tmp);

    return fail;
}

const char *argp_program_version = "zzz 0";
const char *argp_program_bug_address = "<>";
static char doc[] = "fast approximation of large Riemann zeta zeros";
static char args_doc[] = "N [offset] [count]";
static struct argp_option options[] = {
        { "k", 'k', "K", 0, "use first k primes for zero counting function approximation [default 100]"},
        { "evaluate", 'e', 0, 0, "evaluate Riemann zeta function value at the approximate zero location"},
        { "tolerance", 't', "TOL", 0, "tolerance for bisection [default 1e-6]"},
        { "window", 'w', "WIN", 0, "initial span around Lambert W asymptotic zero location +- WIN [default 1.5]"},
        { "precision", 'p', "PREC", 0, "arb precision for counting function approximation [default 256]"},
        { "zeta-prec", 'z', "ZETA_PREC", 0, "arb precision for zeta evaluation [default 64]"},
        { "digits", 'd', "DIGITS", 0, "extra digits for number formatting [default 6]"},
        { "verbose", 'v', 0, 0, "verbose progress output"},
        { "debug", 'g', 0, 0, "debug counting function from <N> to <N+offset> in <count> steps"},
        { "ghy",   'G', 0, 0, "use GHY partial Euler P_X (X = p_k) instead of heuristic damping"},
        { "boot",  'B', "W", 0, "self-consistent hybrid bootstrap: seed 2W+1 zeros with P_X, then iterate leave-one-out P_X*Z_X relocation (implies --ghy)"},
        { "rounds",'R', "R", 0, "bootstrap relocation rounds [default 2]"},
        { 0 }
};

struct arguments {
    slong k;
    slong eval;
    double step0;
    double w0;
    slong PREC;
    slong ZETA_PREC;
    slong DIGITS;
    slong verbose;
    slong debug;
    slong ghy;
    slong boot;
    slong rounds;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = state->input;

    switch (key) {
        case 'k': arguments->k = atol(arg); break;
        case 'e': arguments->eval = 1; break;
        case 't': arguments->step0 = atof(arg); break;
        case 'w': arguments->w0 = atof(arg); break;
        case 'p': arguments->PREC = atol(arg); break;
        case 'z': arguments->ZETA_PREC = atol(arg); break;
        case 'd': arguments->DIGITS = atol(arg); break;
        case 'v': arguments->verbose += 1; break;
        case 'g': arguments->debug = 1; break;
        case 'G': arguments->ghy = 1; break;
        case 'B': arguments->boot = atol(arg); break;
        case 'R': arguments->rounds = atol(arg); break;
        case ARGP_KEY_ARG: return 0;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };

void iter_print(arb_srcptr lo_t, arb_srcptr lo, arb_srcptr hi_t, arb_srcptr hi, slong DIGITS, slong verbose) {

    if(verbose == 0) return;

    flint_fprintf(stderr, "lower x N(x)\t");
    arb_fprintd(stderr, lo, DIGITS);
    flint_fprintf(stderr, "\t");
    arb_fprintd(stderr, lo_t, DIGITS);
    flint_fprintf(stderr, "\n");

    flint_fprintf(stderr, "upper x N(x)\t");
    arb_fprintd(stderr, hi, DIGITS);
    flint_fprintf(stderr, "\t");
    arb_fprintd(stderr, hi_t, DIGITS);
    flint_fprintf(stderr, "\n");
}

int main(int argc, char *argv[])
{

    if(argc == 1) {
        flint_fprintf(stderr, "zzz --usage for help\n");
        return 1;
    }

    struct arguments arguments;

    arguments.k = 100;
    arguments.DIGITS = 6;
    arguments.PREC = 256;
    arguments.ZETA_PREC = 64;
    arguments.step0 = 0.000001;
    arguments.w0 = 1.5;
    arguments.eval = 0;
    arguments.verbose = 0;
    arguments.debug = 0;
    arguments.ghy = 0;
    arguments.boot = 0;
    arguments.rounds = 2;

    int arg_index = 1;
    argp_parse(&argp, argc, argv, ARGP_NO_ARGS, &arg_index, &arguments);

    // GHY mode: resolve X = p_k so the cost matches the heuristic's "k primes"
    ulong ghy_X = 0;
    if ((arguments.ghy || arguments.boot > 0) && arguments.k > 0) {
        ghy_X = n_nth_prime(arguments.k);
    }

    arb_t m;
    arb_t u;
    arb_t m0;
    arb_t mm;
    arb_t m_lo;
    arb_t m_hi;
    arb_t tt;
    arb_t step;

    arb_init(m);
    arb_init(u);
    arb_init(m0);
    arb_init(mm);
    arb_init(m_lo);
    arb_init(m_hi);
    arb_init(tt);
    arb_init(step);

    arb_one(m0);
    arb_zero(m);
    if(argc > arg_index) {
        arb_set_str(m0, argv[arg_index], arguments.PREC);
    }
    if(argc > arg_index+1) {
        arb_set_str(m, argv[arg_index+1], arguments.PREC);
    }

    slong count = 1;
    if(argc > arg_index+2) {
        count = atol(argv[arg_index+2]);
    }

    if(arguments.k > 0) {
        n_compute_primes(arguments.k);
    }

    if(arguments.debug) {
        arb_t tx;
        arb_t zc;
        arb_t stp;

        arb_init(tx);
        arb_init(zc);
        arb_init(stp);
        arb_zero(tx);
        arb_set(stp, m);
        arb_div_ui(stp, stp, count, arguments.PREC);
        arb_set(tx, m0);
        arb_add(m, m, m0, arguments.PREC);

        for (; arb_lt(tx, m); arb_add(tx, tx, stp, arguments.PREC)) {
            if (arguments.ghy) zero_count_ghy(zc, tx, ghy_X, arguments.PREC);
            else               zero_count_exact(zc, tx, arguments.k, arguments.PREC);
            arf_printd(&tx->mid, arguments.DIGITS);
            flint_printf("\t");
            arf_printd(&zc->mid, arguments.DIGITS);
            flint_printf("\n");
        }
    } else {
        arb_add(m0, m0, m, arguments.PREC);

        for (slong ord = 0; ord < count || count < 0; ++ord) {

            nt_inv(tt, m0, arguments.PREC);

            //calc required digits
            arb_const_log10(u, arguments.PREC);
            arb_log(m, tt, arguments.PREC);
            arb_div(u, m, u, arguments.PREC);
            arb_ceil(u, u, arguments.PREC);

            slong digits = arguments.DIGITS;
            digits += arf_get_si(&u->mid, 0);

            if (arguments.verbose) {
                flint_fprintf(stderr, "asymptotic zero location = ");
                arb_fprintd(stderr, tt, digits);
                flint_fprintf(stderr, "\n");
            }

            if (arguments.boot > 0) {
                arb_t bt;
                arb_init(bt);
                if (boot_locate(bt, m0, arguments.boot, arguments.rounds, ghy_X,
                                arguments.w0, arguments.step0, arguments.PREC,
                                arguments.verbose)) {
                    flint_fprintf(stderr, "bootstrap failed to bracket; increase the window\n");
                } else {
                    if (arguments.eval > 0) {
                        acb_t zb;
                        acb_init(zb);
                        flint_fprintf(stderr, "value    z = \t");
                        zeta(zb, bt, arguments.ZETA_PREC);
                        acb_fprintd(stderr, zb, digits);
                        flint_fprintf(stderr, "\n");
                        acb_clear(zb);
                    }
                    arf_fprintd(stdout, &bt->mid, digits);
                    flint_fprintf(stdout, "\n");
                    fflush(stdout);
                }
                arb_clear(bt);
                arb_one(m);
                arb_add(m0, m0, m, arguments.PREC);
                continue;
            }

            arb_set_d(m, -0.5);
            arb_add(m, m0, m, arguments.PREC);
            arb_set(m_lo, m);
            arb_set(m_hi, m);

            arb_t lo;
            arb_t hi;
            arb_t lo_t;
            arb_t hi_t;
            acb_t zz;

            arb_init(lo);
            arb_init(hi);
            arb_init(lo_t);
            arb_init(hi_t);
            acb_init(zz);

            arb_set_d(step, arguments.w0);

            arb_sub(lo_t, tt, step, arguments.PREC);
            arb_add(hi_t, tt, step, arguments.PREC);

            if (arguments.ghy) {
                zero_count_ghy(lo, lo_t, ghy_X, arguments.PREC);
                zero_count_ghy(hi, hi_t, ghy_X, arguments.PREC);
            } else {
                zero_count_exact(lo, lo_t, arguments.k, arguments.PREC);
                zero_count_exact(hi, hi_t, arguments.k, arguments.PREC);
            }

            if (arb_gt(lo, m) || arb_lt(hi, m) || arb_lt(hi, lo)) {
                iter_print(lo_t, lo, hi_t, hi, digits, arguments.verbose);
                flint_fprintf(stderr, "please increase the window\n");
            } else {

                arb_t mid_t;
                arb_t mid;
                arb_init(mid);
                arb_init(mid_t);

                while (1) {

                    iter_print(lo_t, lo, hi_t, hi, digits, arguments.verbose);

                    arb_add(mid_t, lo_t, hi_t, arguments.PREC);
                    arb_set_d(mm, 0.5);
                    arb_mul(mid_t, mid_t, mm, arguments.PREC);

                    if (arguments.ghy) zero_count_ghy(mid, mid_t, ghy_X, arguments.PREC);
                    else               zero_count_exact(mid, mid_t, arguments.k, arguments.PREC);

                    arb_zero(mm);
                    if (arb_gt(mid, m)) {
                        arb_set(hi_t, mid_t);
                        arb_set(hi, mid);
                    } else {
                        arb_set(lo_t, mid_t);
                        arb_set(lo, mid);
                    }

                    arb_sub(mm, hi_t, lo_t, arguments.PREC);
                    arb_set_d(u, arguments.step0);
                    if (arb_lt(mm, u)) break;

                }

                iter_print(lo_t, lo, hi_t, hi, digits, arguments.verbose);

                if (arguments.verbose) {
                    flint_fprintf(stderr, "argument s = \t");
                    arb_set_d(lo_t, 0.5);
                    acb_set_arb_arb(zz, lo_t, mid_t);
                    acb_fprintd(stderr, zz, digits);
                    flint_fprintf(stderr, "\n");
                }
                if (arguments.eval > 0) {
                    flint_fprintf(stderr, "value    z = \t");
                    zeta(zz, mid_t, arguments.ZETA_PREC);
                    acb_fprintd(stderr, zz, digits);
                    flint_fprintf(stderr, "\n");
                }
                arb_set_d(lo_t, 0.005);

                arf_fprintd(stdout, &mid_t->mid, digits);
                flint_fprintf(stdout, "\n");
                fflush(stdout);

                arb_one(m);
                arb_add(m0, m0, m, arguments.PREC);

                arb_clear(mid);
                arb_clear(mid_t);
            }

            arb_clear(lo);
            arb_clear(hi);
            arb_clear(lo_t);
            arb_clear(hi_t);
            acb_clear(zz);

        }
    }

    arb_clear(m);
    arb_clear(u);
    arb_clear(m0);
    arb_clear(mm);
    arb_clear(m_lo);
    arb_clear(m_hi);
    arb_clear(tt);
    arb_clear(step);

    n_cleanup_primes();
    flint_cleanup();

    return 0;
}
