// stage2-sim.c -- frame-by-frame replay of the live-loop shader (stage2-common.glsl +
// stage2-bufferA.glsl) on the CPU, in fp32 like the GPU, checked against a sieve.  Every
// texel of a frame is computed from the previous frame's state, as on Shadertoy.
//
//   gcc -O2 -fopenmp -o stage2-sim doc/ghy/shadertoy/stage2-sim.c -lm
//   ./stage2-sim [NLIM=8000] [old]
//
// "old" replays the tabs as published before 2026-11 (one full-NZ reading against the Li
// envelope at 0.30, reach = nz/2.7, plain method B, 40 seed zeros): 12 composites flagged up
// to 2962, the first at 417.  The current tabs start from no zeros, with 2 marked prime.  Keep this file in step with the two tabs when either changes.
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI     3.14159265358979f
#define TWO_PI 6.28318530717959f
static const float DTAU = 0.40f, UMIN = 3.0f, BOXC = 0.375f, KFIX = 5.0f;
static const int   BATCH = 12, BAND = 4;
static const float SEED[40] = {          // "old" only: the published tabs start from these
    14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073,
    48.005151,49.773832,52.970321,56.446248,59.347044,60.831779,65.112544,67.079811,
    69.546402,72.067158,75.704691,77.144840,79.337375,82.910381,84.735493,87.425275,
    88.809111,92.491899,94.651344,95.870634,98.831194,101.317851,103.725538,105.446623,
    107.168611,111.029536,111.874659,114.320221,116.226680,118.790783,121.370125,122.946829};
static int NLIM = 8000, OLD = 0;
static float *Z;          // .r channel: zeros
static char  *G;          // .g channel: prime indicator

static float pcos(float a) { float k = floorf(a*0.15915494f + 0.5f); float r = (a - k*6.2831855f) + k*1.7484555e-7f; return cosf(r); }
static float psin(float a) { return pcos(a - 1.5707963f); }
static float N0f(float t)  { float u = t/TWO_PI; return u*logf(u/2.718281828f) + 0.875f; }
static float peakf(float x, float gamma) { return (2.0f/PI)*gamma*logf(x)/x; }
static float envf_old(float x, int n) {                 // (20/4) Li(n-4) log x / x
    float z = logf((float)n - 4.0f), s = 0.5772156649015329f + logf(z), t = 1.0f;
    for (int k = 1; k < 64; k++) { t *= z/(float)k; s += t/(float)k; }
    return 5.0f*s*logf(x)/x;
}
static float smoothGamma(int k) {
    float tgt = (float)k - 0.5f, t = TWO_PI*tgt/logf(tgt + 2.0f) + 10.0f;
    for (int i = 0; i < 40; i++) { float f = N0f(t) - tgt, d = logf(t/TWO_PI)/TWO_PI; if (d < 1e-9f) d = 1e-9f; t -= f/d; if (t < 1.0f) t = 1.0f; }
    return t;
}
static int nFrontier(int X, float kmin) {
    if (X < 2) return 0;
    float T = TWO_PI*powf((float)X, TWO_PI/kmin), u = T/TWO_PI, n = floorf(u*(logf(u) - 1.0f));
    return n > 1e9f ? 1000000000 : (int)n;
}
static int isPP(int x) {
    for (int b = 2; b < 14; b++) {
        float a = floorf(powf((float)x, 1.0f/(float)b) + 0.5f);
        if (a < 2.0f) break;
        float p = 1.0f; for (int e = 0; e < b; e++) p *= a;
        if ((int)(p + 0.5f) == x) return 1;
    }
    return 0;
}
static float detFire(float x, int nz) {
    float lx = logf(x), s = 0.0f, base = 1.0f + 1.0f/(x - x*x*x), c4 = 4.0f/sqrtf(x), v[6];
    if (OLD) { for (int k = 0; k < NLIM && k < nz; k++) s += pcos(Z[k]*lx); return (base - c4*s)/envf_old(x, nz); }
    int b = 0, nb = nz/2 > 0 ? nz/2 : 1;
    for (int k = 0; k < NLIM; k++) {
        if (k >= nz) break;
        float zk = Z[k];
        s += pcos(zk*lx);
        while (b < 6 && k + 1 == nb) { v[b] = (base - c4*s)/peakf(x, zk); b++;
            nb = (b == 5) ? nz : (int)floorf((0.5f + 0.1f*(float)b)*(float)nz); if (nb < 1) nb = 1; }
    }
    for (int i = 1; i < 6; i++) { float t = v[i]; int j = i - 1; while (j >= 0 && v[j] > t) { v[j+1] = v[j]; j--; } v[j+1] = t; }
    return 0.5f*(v[2] + v[3]);
}
static float FB(float t, int Xk) {
    float acc = 0.0f, w = OLD ? 0.0f : BOXC*TWO_PI/logf(t/TWO_PI);
    for (int x = 2; x <= Xk; x++) if (G[x]) {
        float lp = logf((float)x), pm = (float)x;
        for (int m = 1; m < 20; m++) {
            float hf = powf((float)x, -0.5f*(float)m), a = (float)m*w*lp;
            acc += (hf/(float)m)*psin((float)m*t*lp)*(a > 0.0f ? sinf(a)/a : 1.0f);
            if (pm > (float)Xk/(float)x) break;
            pm *= (float)x;
        }
    }
    return N0f(t) - acc/PI;
}
static float locateSmooth(int k, int Xknown, float kmin) {
    float tgt = (float)k - 0.5f, g = smoothGamma(k);
    int Xk = (int)fminf(powf(g/TWO_PI, kmin/TWO_PI), (float)Xknown); if (Xk < 2) Xk = 2;
    float gap = TWO_PI/logf(g/TWO_PI);
    for (float w = 1.5f; w <= 9.0f; w += 1.5f) {
        float lo = g - w*gap, hi = g + w*gap, flo = FB(lo, Xk) - tgt, fhi = FB(hi, Xk) - tgt;
        if (flo*fhi <= 0.0f) {
            for (int it = 0; it < 50; it++) { float mid = 0.5f*(lo + hi), fm = FB(mid, Xk) - tgt; if (flo*fm <= 0.0f) { hi = mid; fhi = fm; } else { lo = mid; flo = fm; } }
            return 0.5f*(lo + hi);
        }
    }
    return g;
}

int main(int argc, char **argv) {
    if (argc > 1) NLIM = atoi(argv[1]);
    if (argc > 2) OLD = !strcmp(argv[2], "old");
    Z = calloc(NLIM + 16, sizeof(float)); G = calloc(NLIM + 16, 1);
    float *rlog = calloc(NLIM + 16, sizeof(float));
    int nz = 0, Xknown = 2, idle = 0; long frames = 0;
    if (OLD) { memcpy(Z, SEED, sizeof SEED); nz = 40; Xknown = 1; } else G[2] = 1;
    for (;; ++frames) {
        int reach = OLD ? (int)((float)nz/2.7f) : nz == 0 ? Xknown : (int)(Z[(nz/2 > 0 ? nz/2 : 1) - 1]/UMIN);
        int new_X = Xknown + BAND; if (new_X > reach) new_X = reach; if (new_X < Xknown) new_X = Xknown;
        int new_nz = nz;
        if (OLD ? new_X >= 8 : 1) {
            int nfr = OLD ? nFrontier(new_X, KFIX) : nFrontier(Xknown, KFIX), tgt = nz + BATCH;
            if (nfr < 1) nfr = 1;
            if (tgt > nfr) tgt = nfr;
            if (tgt > NLIM) tgt = NLIM;
            if (tgt > nz) new_nz = tgt;
        }
        float zn[16]; char gn[8];
        #pragma omp parallel for
        for (int i = nz; i < new_nz; i++) zn[i - nz] = locateSmooth(i + 1, Xknown, KFIX);
        for (int i = Xknown + 1; i <= new_X; i++) {
            float r = i >= 2 ? detFire((float)i, nz) : 0.0f;
            rlog[i] = r; gn[i - Xknown - 1] = i >= 2 && r > (OLD ? 0.30f : DTAU) && !isPP(i);
        }
        for (int i = nz; i < new_nz; i++) Z[i] = zn[i - nz];
        for (int i = Xknown + 1; i <= new_X; i++) G[i] = gn[i - Xknown - 1];
        if (new_nz == nz && new_X == Xknown) { if (++idle > 3) break; } else idle = 0;
        nz = new_nz; Xknown = new_X;
    }
    char *comp = calloc(Xknown + 2, 1);
    for (int a = 2; a <= Xknown; a++) if (!comp[a]) for (int q = 2*a; q <= Xknown; q += a) comp[q] = 1;
    int np = 0, bad = 0; float pmin = 9, cmax = -9;
    for (int a = 2; a <= Xknown; a++) {
        np += G[a];
        if (G[a] == comp[a]) { printf("  %s %d\n", G[a] ? "FALSE PRIME" : "missed prime", a); bad++; }
        if (a >= 100 && !isPP(a)) { if (!comp[a] && rlog[a] < pmin) pmin = rlog[a]; if (comp[a] && rlog[a] > cmax) cmax = rlog[a]; }
    }
    printf("%s tabs, NLIM=%d: %ld frames, %d zeros to gamma=%.1f, X=%d, %d primes flagged, %d wrong\n",
           OLD ? "old" : "current", NLIM, frames, nz, Z[nz-1], Xknown, np, bad);
    printf("readings for x >= 100: lowest prime %.3f, highest composite %.3f (threshold %.2f)\n", pmin, cmax, OLD ? 0.30f : DTAU);
    return bad != 0;
}
