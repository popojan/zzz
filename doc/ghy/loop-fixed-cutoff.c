// loop-fixed-cutoff.c -- the forward step of `zzz --loop` with a FIXED prime list.
//
// Locates the zeros with ordinates in [tlo, thi] from the true primes <= X only, with the
// same F_B (box-smoothed method B) and marching bisection as src/loop.c, and writes them as
// raw doubles.  It answers "what do zeros computed from the primes <= X know about the
// integers beyond X?" -- see loop-fixed-cutoff.py and doc/notes/loop-margins.md.
//
//   gcc -O3 -march=native -ffast-math -fopenmp -o loop-fixed-cutoff doc/ghy/loop-fixed-cutoff.c -lm
//   ./loop-fixed-cutoff X tlo thi boxc out.bin        e.g.  100 120 200000 0.375 z100.bin
//
// (-ffast-math only vectorises the sine; the zeros agree with loop.c's to ~1e-12.)
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
static long NT; static double *LW, *AMP;   // term frequency m*log p, amplitude p^{-m/2}/m
static double g_boxc;
static double N0(double t) { double u = t / (2.0 * M_PI); return u * log(u / M_E) + 0.875; }
static double FB(double t) {
    double w = 0.0, u = t / (2.0 * M_PI);
    if (g_boxc > 0.0 && u > 1.5) w = g_boxc * 2.0 * M_PI / log(u);
    double acc = 0.0;
    if (w > 0.0) { for (long i = 0; i < NT; ++i) { double a = w * LW[i]; acc += AMP[i] * sin(t * LW[i]) * sin(a) / a; } }
    else         { for (long i = 0; i < NT; ++i) acc += AMP[i] * sin(t * LW[i]); }
    return N0(t) - acc / M_PI;
}
static double locate_march(long k, double prev) {
    double tgt = k - 0.5;
    double gap = 2.0 * M_PI / log((prev > 10.0 ? prev : 10.0) / (2.0 * M_PI));
    double lo = prev + 0.25 * gap, hi = prev + 1.9 * gap;
    double flo = FB(lo) - tgt, fhi = FB(hi) - tgt;
    if (flo * fhi > 0) { lo = prev + 0.04 * gap; hi = prev + 2.7 * gap; flo = FB(lo) - tgt; fhi = FB(hi) - tgt; }
    if (flo * fhi > 0) return prev + gap;
    for (int it = 0; it < 50; ++it) { double mid = 0.5 * (lo + hi), fm = FB(mid) - tgt;
        if (flo * fm <= 0) { hi = mid; fhi = fm; } else { lo = mid; flo = fm; } }
    return 0.5 * (lo + hi);
}
static double n0inv(double n) { double t = 2 * M_PI * n / log(n + 3) + 10; for (int i = 0; i < 80; ++i) { double u = t/(2*M_PI); t -= (N0(t) - n) / (log(u) / (2*M_PI)); } return t; }
int main(int argc, char **argv) {
    long X = atol(argv[1]); double tlo = atof(argv[2]), thi = atof(argv[3]); g_boxc = atof(argv[4]);
    char *c = calloc(X + 1, 1); LW = malloc((X + 64) * sizeof(double)); AMP = malloc((X + 64) * sizeof(double));
    long np = 0;
    for (long a = 2; a <= X; ++a) if (!c[a]) { ++np; for (long q = 2*a; q <= X; q += a) c[q] = 1;
        double lp = log((double)a), pm = (double)a;
        for (long m = 1; ; ++m) { LW[NT] = m * lp; AMP[NT] = pow((double)a, -0.5 * m) / m; ++NT; if (pm > (double)X / a) break; pm *= a; } }
    long kstart = (long)N0(tlo) + 1, kmax = (long)N0(thi), n = kmax - kstart + 1;
    double *z = malloc(n * sizeof(double));
    int NCH = 16 * omp_get_max_threads(); long ch = n / NCH + 1;
    #pragma omp parallel for schedule(dynamic,1)
    for (int cidx = 0; cidx < NCH; ++cidx) {
        long a = cidx * ch, b = a + ch; if (b > n) b = n; if (a >= b) continue;
        long k0 = kstart + a - 1; double t0 = n0inv(k0 - 0.5), gap = 2*M_PI/log(t0/(2*M_PI));
        double lo = t0 - 4*gap, hi = t0 + 4*gap, tgt = k0 - 0.5, flo = FB(lo) - tgt;
        for (int it = 0; it < 50; ++it) { double mid = 0.5*(lo+hi), fm = FB(mid) - tgt; if (flo * fm <= 0) hi = mid; else { lo = mid; flo = fm; } }
        double prev = 0.5 * (lo + hi);
        for (long i = a; i < b; ++i) { prev = locate_march(kstart + i, prev); z[i] = prev; }
    }
    FILE *f = fopen(argv[5], "wb"); fwrite(z, sizeof(double), n, f); fclose(f);
    fprintf(stderr, "X=%ld primes=%ld terms=%ld zeros %ld..%ld (%ld)  z[0]=%.4f z[last]=%.4f\n", X, np, NT, kstart, kmax, n, z[0], z[n-1]);
    return 0;
}
