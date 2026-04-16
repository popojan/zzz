// zhad.c -- evaluate GHY local Hadamard factor Z_X(s) at s = 1/2 + iT on a grid.
//
// log Z_X(s) = -sum_j E_1((s - rho_j) log X)
//   (X -> infinity limit of the GHY kernel; U -> E_1)
// rho_j = 1/2 + i gammas[j] read from a zeros file (one gamma per line).
//
// Usage: zhad <zeros_file> <T0> <T1> <N> <X> [window_units]
//   Columns: T  log|Z_X|  arg(Z_X)

#include "ghy.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Load a simple one-gamma-per-line file; lines with '#' ignored.
static slong load_zeros(const char *path, double **gammas) {
    FILE *fp = fopen(path, "r");
    if (!fp) { perror(path); return -1; }

    slong cap = 1024, count = 0;
    double *buf = malloc(cap * sizeof(double));
    if (!buf) { fclose(fp); return -1; }

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
            "usage: zhad <zeros_file> <T0> <T1> <N> <X> [window_units]\n"
            "  log Z_X(1/2+iT) = -sum_j E_1((1/2+iT - rho_j) log X)\n"
            "  rho_j = 1/2 + i gamma_j  read from zeros_file (one gamma/line)\n"
            "  window_units: skip zeros with |t-gamma_j|*log X > this (default 20)\n"
            "  Columns: T  log|Z_X|  arg(Z_X)\n");
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

    arb_t t, step;
    acb_t s, half, logz;

    arb_init(t); arb_init(step);
    acb_init(s); acb_init(half); acb_init(logz);

    acb_one(half);
    acb_div_ui(half, half, 2, PREC);

    arb_set_d(t, t0);
    arb_set_d(step, N > 0 ? (t1 - t0) / (double)N : 0.0);

    for (slong i = 0; i <= N; ++i) {
        acb_set(s, half);
        arb_set(acb_imagref(s), t);

        ghy_log_zx(logz, s, gammas, n_g, X, window_units, PREC);

        arf_printd(arb_midref(t), 12);                        printf("\t");
        arf_printd(arb_midref(acb_realref(logz)), 10);        printf("\t");
        arf_printd(arb_midref(acb_imagref(logz)), 10);        printf("\n");

        arb_add(t, t, step, PREC);
    }

    free(gammas);
    arb_clear(t); arb_clear(step);
    acb_clear(s); acb_clear(half); acb_clear(logz);

    flint_cleanup();
    return 0;
}
