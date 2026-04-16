// zghy.c — evaluate GHY partial Euler factor P_X(s) and the counting
// proxy F_P(T) = N_0(T) + arg(P_X(1/2+iT)) / pi.
//
// P_X(s) = exp( sum_{p^m <= X} 1/(m p^{ms}) )   -- GHY 2007 eq. (6)
// log P_X(s) = sum_{p <= X} sum_{m <= floor(log X / log p)} 1/(m p^{ms})
//
// Compared with zproxy's Z_k(s,T) (heuristic, boundary-admissible damping),
// P_X is the *provably admissible* partial Euler factor in Gonek-Hughes-Young.
// The full GHY rigorous proxy is P_X * Z_X where Z_X handles local zeros;
// this binary implements P_X alone so we can diff against the heuristic and
// measure what Z_X will need to contribute.
//
// Usage: zghy <T0> <T1> <N> <X>
//   prints N+1 rows (TSV): T  log|P_X|  arg(P_X)  F_P(T)

#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/ulong_extras.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// log P_X(s) = sum over primes p<=X of sum over m with p^m<=X of 1/(m p^{ms})
static void log_px(acb_ptr out, const acb_t s, ulong X, slong PREC) {
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
        acb_pow(p_neg_s, p_acb, neg_s, PREC);   // p^{-s}
        acb_set(psm_neg, p_neg_s);              // p^{-m s} with m=1

        ulong pm = p;
        for (ulong m = 1; ; ++m) {
            // accumulate 1/(m p^{ms}) = p^{-ms} / m
            acb_div_ui(term, psm_neg, m, PREC);
            acb_add(total, total, term, PREC);

            // advance to m+1: does p^{m+1} still fit in X?
            if (pm > X / p) break;       // pm * p would exceed X
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

// N_0(T) = (T/2pi) log(T/2pi e) + 7/8
static void n0_smooth(arb_ptr out, arb_srcptr t, slong PREC) {
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

int main(int argc, char *argv[]) {
    if (argc < 5) {
        fprintf(stderr,
            "usage: zghy <T0> <T1> <N> <X>\n"
            "  prints N+1 rows (TSV): T  log|P_X|  arg(P_X)  F_P(T)\n"
            "  P_X(s) = GHY partial Euler factor (eq. 6 of GHY 2007)\n"
            "  F_P(T) = N_0(T) + arg P_X(1/2+iT) / pi\n");
        return 1;
    }

    double t0 = atof(argv[1]);
    double t1 = atof(argv[2]);
    slong N = atol(argv[3]);
    ulong X = (ulong)atol(argv[4]);
    const slong PREC = 128;

    // Precompute enough primes. pi(X) <= 1.3 * X / log(X) for X >= 17,
    // use a safe over-estimate.
    ulong n_primes_guess;
    if (X < 20) n_primes_guess = 10;
    else {
        double dx = (double)X;
        n_primes_guess = (ulong)(1.3 * dx / (log(dx) - 1.1)) + 16;
    }
    n_compute_primes(n_primes_guess);

    arb_t t, step, n0, pi, arg_px, log_abs_px, fp;
    acb_t s, half, it, logp;

    arb_init(t); arb_init(step); arb_init(n0); arb_init(pi);
    arb_init(arg_px); arb_init(log_abs_px); arb_init(fp);
    acb_init(s); acb_init(half); acb_init(it); acb_init(logp);

    arb_const_pi(pi, PREC);
    arb_set_d(t, t0);
    arb_set_d(step, (t1 - t0) / (double)N);

    acb_one(half);
    acb_div_ui(half, half, 2, PREC);   // 1/2

    for (slong i = 0; i <= N; ++i) {
        // s = 1/2 + i t
        acb_set(s, half);
        arb_set(acb_imagref(s), t);

        log_px(logp, s, X, PREC);
        arb_set(log_abs_px, acb_realref(logp));
        arb_set(arg_px, acb_imagref(logp));

        n0_smooth(n0, t, PREC);
        arb_div(fp, arg_px, pi, PREC);
        arb_add(fp, fp, n0, PREC);

        arf_printd(arb_midref(t), 12);         printf("\t");
        arf_printd(arb_midref(log_abs_px), 10);printf("\t");
        arf_printd(arb_midref(arg_px), 10);    printf("\t");
        arf_printd(arb_midref(fp), 10);        printf("\n");

        arb_add(t, t, step, PREC);
    }

    arb_clear(t); arb_clear(step); arb_clear(n0); arb_clear(pi);
    arb_clear(arg_px); arb_clear(log_abs_px); arb_clear(fp);
    acb_clear(s); acb_clear(half); acb_clear(it); acb_clear(logp);

    n_cleanup_primes();
    flint_cleanup();
    return 0;
}
