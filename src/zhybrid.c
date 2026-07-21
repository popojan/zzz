// zhybrid.c -- full GHY hybrid proxy on the critical line.
//
// F_hybrid(T) = N_0(T) + arg(P_X(1/2+iT) * Z_X(1/2+iT)) / pi
//            = N_0(T) + (arg P_X + arg Z_X) / pi
//
// P_X: GHY partial Euler factor (primes & prime powers up to X)
// Z_X: GHY local Hadamard factor over supplied zeros within the window
//
// F_hybrid is a zeta-evaluation-free approximation to N(T) with explicit
// error bounds (GHY Thm 1). In the limit X fixed and zeros list complete,
// it is an *identity* on N modulo the explicit GHY error.
//
// Usage: zhybrid <zeros_file> <T0> <T1> <N> <X> [window_units]
//   Columns: T  log|P|  arg P  log|Z|  arg Z  F_hybrid

#include "ghy.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static slong load_zeros(const char *path, double **gammas) {
    FILE *fp = fopen(path, "r");
    if (!fp) { perror(path); return -1; }
    slong cap = 1024, count = 0;
    double *buf = malloc(cap * sizeof(double));
    char line[512];
    while (fgets(line, sizeof line, fp)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\0') continue;
        double g;
        if (sscanf(line, "%lf", &g) == 1) {
            if (count == cap) {
                cap *= 2;
                buf = realloc(buf, cap * sizeof(double));
            }
            buf[count++] = g;
        }
    }
    fclose(fp);
    *gammas = buf;
    return count;
}

int main(int argc, char *argv[]) {
    if (argc < 6) {
        fprintf(stderr,
            "usage: zhybrid <zeros_file> <T0> <T1> <N> <X> [window_units]\n"
            "  F_hybrid(T) = N_0(T) + (arg P_X + arg Z_X)/pi\n"
            "  Columns: T  log|P|  arg P  log|Z|  arg Z  F_hybrid\n");
        return 1;
    }

    const char *zeros_path = argv[1];
    double t0 = atof(argv[2]);
    double t1 = atof(argv[3]);
    slong N = atol(argv[4]);
    ulong X = (ulong)atol(argv[5]);
    double window_units = (argc > 6) ? atof(argv[6]) : 20.0;
    const slong PREC = 128;

    double *gammas = NULL;
    slong n_g = load_zeros(zeros_path, &gammas);
    if (n_g < 0) return 1;
    fprintf(stderr, "# loaded %ld zeros from %s\n", (long)n_g, zeros_path);

    ulong n_primes_guess;
    if (X < 20) n_primes_guess = 10;
    else {
        double dx = (double)X;
        n_primes_guess = (ulong)(1.3 * dx / (log(dx) - 1.1)) + 16;
    }
    n_compute_primes(n_primes_guess);

    arb_t t, step, n0, pi, f_hybrid, tmp;
    acb_t s, half, logp, logz;

    arb_init(t); arb_init(step); arb_init(n0); arb_init(pi);
    arb_init(f_hybrid); arb_init(tmp);
    acb_init(s); acb_init(half); acb_init(logp); acb_init(logz);

    arb_const_pi(pi, PREC);
    acb_one(half);
    acb_div_ui(half, half, 2, PREC);

    arb_set_d(t, t0);
    arb_set_d(step, N > 0 ? (t1 - t0) / (double)N : 0.0);

    for (slong i = 0; i <= N; ++i) {
        acb_set(s, half);
        arb_set(acb_imagref(s), t);

        ghy_log_px(logp, s, X, PREC);
        ghy_log_zx(logz, s, gammas, n_g, X, window_units, PREC);

        ghy_n0_smooth(n0, t, PREC);

        // f = N_0 + (arg P + arg Z) / pi
        arb_add(tmp, acb_imagref(logp), acb_imagref(logz), PREC);
        arb_div(tmp, tmp, pi, PREC);
        arb_add(f_hybrid, n0, tmp, PREC);

        arf_printd(arb_midref(t), 12);                       printf("\t");
        arf_printd(arb_midref(acb_realref(logp)), 10);       printf("\t");
        arf_printd(arb_midref(acb_imagref(logp)), 10);       printf("\t");
        arf_printd(arb_midref(acb_realref(logz)), 10);       printf("\t");
        arf_printd(arb_midref(acb_imagref(logz)), 10);       printf("\t");
        arf_printd(arb_midref(f_hybrid), 10);                printf("\n");

        arb_add(t, t, step, PREC);
    }

    free(gammas);
    arb_clear(t); arb_clear(step); arb_clear(n0); arb_clear(pi);
    arb_clear(f_hybrid); arb_clear(tmp);
    acb_clear(s); acb_clear(half); acb_clear(logp); acb_clear(logz);

    n_cleanup_primes();
    flint_cleanup();
    return 0;
}
