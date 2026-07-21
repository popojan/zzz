// zghy.c — evaluate GHY partial Euler factor P_X(s) and the counting
// proxy F_P(T) = N_0(T) + arg(P_X(1/2+iT)) / pi.
//
// P_X(s) = exp( sum_{p^m <= X} 1/(m p^{ms}) )   -- GHY 2007 eq. (6)
//
// Usage: zghy <T0> <T1> <N> <X>
//   prints N+1 rows (TSV): T  log|P_X|  arg(P_X)  F_P(T)

#include "ghy.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

    // Precompute primes up to X. pi(X) <= 1.3 X / (log X - 1.1) for X >= 17.
    ulong n_primes_guess;
    if (X < 20) n_primes_guess = 10;
    else {
        double dx = (double)X;
        n_primes_guess = (ulong)(1.3 * dx / (log(dx) - 1.1)) + 16;
    }
    n_compute_primes(n_primes_guess);

    arb_t t, step, n0, pi, arg_px, log_abs_px, fp;
    acb_t s, half, logp;

    arb_init(t); arb_init(step); arb_init(n0); arb_init(pi);
    arb_init(arg_px); arb_init(log_abs_px); arb_init(fp);
    acb_init(s); acb_init(half); acb_init(logp);

    arb_const_pi(pi, PREC);
    arb_set_d(t, t0);
    arb_set_d(step, (t1 - t0) / (double)N);

    acb_one(half);
    acb_div_ui(half, half, 2, PREC);

    for (slong i = 0; i <= N; ++i) {
        acb_set(s, half);
        arb_set(acb_imagref(s), t);

        ghy_log_px(logp, s, X, PREC);
        arb_set(log_abs_px, acb_realref(logp));
        arb_set(arg_px, acb_imagref(logp));

        ghy_n0_smooth(n0, t, PREC);
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
    acb_clear(s); acb_clear(half); acb_clear(logp);

    n_cleanup_primes();
    flint_cleanup();
    return 0;
}
