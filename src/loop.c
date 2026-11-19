// loop.c -- implementation of loop.h (zzz --loop). See doc/notes/zeros-primes-bootstrap.md
// and doc/ghy/bootstrap-selfpave.wls (the validated Wolfram prototype this ports).
//
// Audited free of any zeta evaluation and any primality test in the loop body:
// the only inputs to "is x a prime power" are the psi' signal (from the loop's own
// zeros) and the loop's own discovered-prime list (perfect-power bookkeeping).

#include "loop.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define EULERG 0.57721566490153286061

// ---- seed data (finite, given; not part of the self-paving claim) ----------
// First 40 nontrivial zero ordinates gamma_n (Im rho).
static const double SEED[] = {
    14.134725141734694, 21.022039638771555, 25.010857580145686, 30.424876125859523,
    32.935061587739190, 37.586178158825680, 40.918719012147505, 43.327073280915015,
    48.005150881167160, 49.773832477672312, 52.970321477714471, 56.446247697063406,
    59.347044002602354, 60.831778524609819, 65.112544048081617, 67.079810529494184,
    69.546401711173988, 72.067157674481918, 75.704690699083943, 77.144840068874816,
    79.337375020249378, 82.910380854086042, 84.735492980517060, 87.425274613125231,
    88.809111207634465, 92.491899270558485, 94.651344040519889, 95.870634228245310,
    98.831194218193691, 101.31785100573139, 103.72553804047834, 105.44662305232611,
    107.16861118427641, 111.02953554316969, 111.87465917699264, 114.32022091545271,
    116.22668032085756, 118.79078286597622, 121.37012500242064, 122.94682929355259
};
#define NSEED ((long)(sizeof(SEED) / sizeof(SEED[0])))

// ---- detector / scale (matched-filter peak-height normalisation) -----------
// The truncated psi' has a Dirichlet-kernel peak at each prime power of height
// ~ (gamma_n/2pi) * Lambda(x)/x; with Li(n) ~ gamma_n/2pi this is the envelope
// below.  Detection = peak as a fraction tau of that expected full-prime height.
// (The textbook integral int psi' = Lambda fails here: under truncation the +1
// in psi' is not cancelled and leaves a baseline -- verified to misfire.)
static const double TAULO = 0.12;   // reach gate: "clear composite" level at the lowest truncation
static const double ENVK  = 20.0 / 4.0;
// loop_opts.contrast: fit-free normalisation.  The peak of the truncated sum at a prime is
// Landau's main term (2*gamma_n/pi) * log(x)/x, known from the loop's own top zero, so
// neither Li nor the empirical 20/4 is needed (5*Li(n-4) ~ 4*gamma_n/2pi to a few % here,
// drifting to 5:4 asymptotically).
static int g_fitfree = 0;

// ---- growable state --------------------------------------------------------
typedef struct {
    double *z;  long nz, zcap;       // zero ordinates (z[i] = gamma_{i+1})
    long   *p;  long np, pcap;       // discovered primes (ascending)
    long    Xknown;                  // largest x covered by detection
    long    iter;
    double  kmin;                    // current (annealed) kappa margin
} state;

static volatile sig_atomic_t g_stop = 0;
static void on_sigint(int sig) { (void)sig; g_stop = 1; }

static void z_push(state *st, double v) {
    if (st->nz == st->zcap) { st->zcap = st->zcap ? st->zcap * 2 : 1024;
        st->z = realloc(st->z, st->zcap * sizeof(double)); }
    st->z[st->nz++] = v;
}
static int  p_has(const state *st, long x) {     // membership (ascending list)
    long lo = 0, hi = st->np - 1;
    while (lo <= hi) { long m = (lo + hi) / 2;
        if (st->p[m] == x) return 1;
        if (st->p[m] < x) lo = m + 1; else hi = m - 1; }
    return 0;
}
static void p_push(state *st, long x) {           // append (x larger than all)
    if (st->np == st->pcap) { st->pcap = st->pcap ? st->pcap * 2 : 256;
        st->p = realloc(st->p, st->pcap * sizeof(long)); }
    st->p[st->np++] = x;
}

// ---- math primitives (zeta-free) -------------------------------------------
static double ei(double zz) {                     // exponential integral, zz>0
    double s = EULERG + log(zz), t = 1.0;
    for (int k = 1; k < 400; ++k) {
        t *= zz / k;
        double add = t / k;
        s += add;
        if (fabs(add) < 1e-15 * fabs(s) && k > zz) break;
    }
    return s;
}
static double li(double x) { return ei(log(x)); }                 // Li(x)=Ei(ln x)

static int cmp_d(const void *a, const void *b) {
    double d = *(const double *)a - *(const double *)b; return d < 0 ? -1 : d > 0 ? 1 : 0;
}

// is x a perfect power p^m (m>=2) of an ALREADY-DISCOVERED prime? (not a primality test)
static int is_known_pow(long x, const state *st) {
    for (long i = 0; i < st->np; ++i) {
        long p = st->p[i];
        if (p * p > x) break;
        long q = p * p;
        while (q < x) q *= p;
        if (q == x) return 1;
    }
    return 0;
}

// N0(t) = (t/2pi) log(t/2pi e) + 7/8   (smooth Riemann-von Mangoldt count)
static double N0(double t) {
    double u = t / (2.0 * M_PI);
    return u * log(u / M_E) + 0.875;
}
// forward box-smoothing coefficient (loop_opts.boxc): w = g_boxc * gap(t); 0 = plain B
static double g_boxc = 0.0;
// F_B(t) = N0(t) + Im log P_X / pi,  Im log P_X = -sum_{p^m<=X} (1/m) p^{-m/2} sin(m t log p)
// With g_boxc > 0 each term carries sinc(m*w*log p): that is exactly the box average
// (1/2w) int_{t-w}^{t+w} F_B dx.  The box average of the TRUE staircase still crosses
// k-1/2 at gamma_k when the window holds one zero, so the bisection target is unchanged;
// the weight only damps the noisy top of the prime band.  Gate-checked 2.5-3.5x location
// gain for kappa = gap*log X in (1.5pi, 3pi) -- the loop's operating band -- and null at
// kappa <= pi, consistent with doc/notes/band-saturation.md (in-band functional).
static double FB(double t, const state *st, long X) {
    double w = 0.0;
    if (g_boxc > 0.0) {
        double u = t / (2.0 * M_PI);
        if (u > 1.5) w = g_boxc * 2.0 * M_PI / log(u);
    }
    double acc = 0.0;
    for (long i = 0; i < st->np; ++i) {
        long p = st->p[i];
        if (p > X) break;
        double lp = log((double)p), pm = (double)p, half = 1.0;
        for (long m = 1; ; ++m) {
            half = pow((double)p, -0.5 * m);
            double term = (half / m) * sin(m * t * lp);
            if (w > 0.0) { double a = m * w * lp; term *= sin(a) / a; }
            acc += term;
            if (pm > (double)X / p) break;
            pm *= p;
        }
    }
    return N0(t) - acc / M_PI;
}

// locate zero #k from primes<=X, method-B bisection, marching bracket off prev
static double locate_march(long k, const state *st, long X, double prev) {
    double tgt = k - 0.5;
    double gap = 2.0 * M_PI / log((prev > 10.0 ? prev : 10.0) / (2.0 * M_PI));
    double lo = prev + 0.25 * gap, hi = prev + 1.9 * gap;
    double flo = FB(lo, st, X) - tgt, fhi = FB(hi, st, X) - tgt;
    if (flo * fhi > 0) { lo = prev + 0.04 * gap; hi = prev + 2.7 * gap;
        flo = FB(lo, st, X) - tgt; fhi = FB(hi, st, X) - tgt; }
    // a gap wider than 2.7 mean gaps (about 2 zeros in 10^4): keep widening rather than guess
    for (double w = 4.0; flo < 0 && fhi < 0 && w <= 16.0; w *= 2.0) {
        lo = hi; flo = fhi;
        hi = prev + w * gap; fhi = FB(hi, st, X) - tgt;
    }
    if (flo * fhi > 0) return prev + gap;                  // fallback: predicted
    for (int it = 0; it < 50; ++it) {
        double mid = 0.5 * (lo + hi), fm = FB(mid, st, X) - tgt;
        if (flo * fm <= 0) { hi = mid; fhi = fm; } else { lo = mid; flo = fm; }
    }
    return 0.5 * (lo + hi);
}

// backward: scan x from x_start upward (primes below x_start already known),
// detect prime powers from the psi' peak, append new primes. Returns X_known.
//
// A SINGLE-NZ reading of r has Dirichlet-sidelobe variance: r(x) oscillates with
// gamma_max, so near the reach edge a prime can dip into a null (a miss) and a
// prime-adjacent composite can ride a sidelobe above threshold (a false positive).
// This is NOT a fundamental wall -- the signal is stable across NZ.  So (mirroring
// src/vulkan/loopvk) evaluate r at K=6 NZ-truncations (gamma_max 0.5x..1x) and take
// the MEDIAN per candidate; the truncations are PREFIX SUMS of the same cos series,
// so one pass over the zeros yields all K (no extra cost).  The frontier advances
// STRICTLY CONTIGUOUSLY -- it stops at the first unresolved candidate and never
// leapfrogs one (a skipped prime would be a permanent miss: nothing downstream notices
// or repairs a wrong list).  Commits happen where even the LOWEST truncation resolves unit
// spacing, gamma(0.5*NZ) >~ pi*x; there the median is bimodal-clean (primes ~0.8,
// composites ~-0.05 +- 0.07).  The reach gate finds that point statistically; U_MIN states
// the requirement outright.  Measured floor (doc/notes/loop-margins.md): below
// gamma(0.5*NZ) ~ 1.8*x a composite between twin primes rides both main lobes and reads as
// a prime, at every height; the gate has not let the ratio under 2.8.  Below X_GATE the
// window is too prime-power-dense for the gate, and the floor alone bounds the reach.
static long backward(state *st, long x_start, long Xcap) {
    const long   GW = 30;            // reach-gate neighbourhood half-width
    const double FRAC_OK = 0.70;     // clear-composite fraction (lowest truncation) => inside reach
    const long   X_GATE = 512;       // statistical gate above this (below, the window is power-dense)
    const double TAU_C = 0.40, TAU_LO_M = 0.25;   // median commit / clear-composite thresholds
    const double U_MIN = 2.5;        // lag floor: commit x only while gamma(0.5*NZ) >= U_MIN*x
    enum { K = 6 };
    static const double FR[K] = { 0.50, 0.60, 0.70, 0.80, 0.90, 1.00 };
    long nz = st->nz;
    if (x_start < 2) x_start = 2;
    long s0 = x_start - GW; if (s0 < 2) s0 = 2;
    long s1 = Xcap + GW, NX = s1 - s0 + 1;
    if (NX < 1) return st->Xknown;
    long nb[K]; double LiK[K];                            // truncation sizes and peak-height scales
    for (int b = 0; b < K; ++b) { nb[b] = (b == K-1) ? nz : (long)floor(FR[b] * (double)nz);
        if (nb[b] < 6) nb[b] = (nz < 6 ? nz : 6);
        LiK[b] = g_fitfree ? (2.0 / M_PI) * st->z[nb[b] - 1] : ENVK * li((double)nb[b] - 4.0); }
    double *R = malloc((size_t)NX * K * sizeof(double));  // r at each truncation per candidate
    for (long i = 0; i < NX; ++i) {                       // prefix sum: all K truncations in one pass
        double x = (double)(s0 + i), lx = log(x), s = 0.0;
        double base = 1.0 + 1.0 / (x - x*x*x), c4 = 4.0 / sqrt(x);
        long prev = 0;
        for (int b = 0; b < K; ++b) {
            for (long k = prev; k < nb[b]; ++k) s += cos(st->z[k] * lx);
            R[i*K + b] = (base - c4 * s) / (LiK[b] * lx / x);
            prev = nb[b];
        }
    }
#define RAT(X,B) (((X) < s0 || (X) > s1) ? 0.0 : R[((X) - s0)*K + (B)])
    long Xk = x_start - 1;
    for (long x = x_start; x <= Xcap; ++x) {
        int kp = is_known_pow(x, st), resolved;
        int inreach = st->z[nb[0] - 1] >= U_MIN * (double)x;
        if (inreach && x >= X_GATE) {                    // reach frac from the lowest truncation
            long cnt = 0, tot = 0;
            for (long j = x - GW; j <= x + GW; ++j) { if (j < 2) continue; ++tot; if (RAT(j, 0) < TAULO) ++cnt; }
            inreach = tot && (double)cnt / (double)tot >= FRAC_OK;
        }
        if (!inreach) {
            resolved = 0;                                // beyond reach(0.5*NZ) -> stop
        } else {
            double v[K]; for (int b = 0; b < K; ++b) v[b] = RAT(x, b);
            qsort(v, K, sizeof(double), cmp_d);
            double m = 0.5 * (v[K/2 - 1] + v[K/2]);       // median over the K truncations
            if (m > TAU_C)         { if (!kp && !p_has(st, x)) { p_push(st, x); printf("%ld\n", x); fflush(stdout); } resolved = 1; }
            else if (m < TAU_LO_M) resolved = 1;         // clear composite
            else                   resolved = (kp || p_has(st, x));   // ambiguous -> only if known
        }
        if (resolved) Xk = x; else break;                // contiguous: stop at first unresolved
    }
#undef RAT
    free(R);
    return Xk;
}

// max zero index resolvable from primes<=X at margin kappa>=kmin
static long n_frontier(long X, double kmin, long nmax) {
    if (X < 2) return 0;
    double T = 2.0 * M_PI * pow((double)X, 2.0 * M_PI / kmin);
    double u = T / (2.0 * M_PI);
    long n = (long)floor(u * (log(u) - 1.0));
    if (n > nmax) n = nmax;
    return n;
}

// ---- checkpoint I/O --------------------------------------------------------
static int save_state(const state *st, const char *path) {
    char tmp[1024]; snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) { fprintf(stderr, "loop: cannot write %s\n", tmp); return 1; }
    fprintf(f, "# zzz --loop checkpoint\niter %ld\nxknown %ld\nkmin %.17g\n",
            st->iter, st->Xknown, st->kmin);
    fprintf(f, "nprimes %ld\n", st->np);
    for (long i = 0; i < st->np; ++i) fprintf(f, "%ld\n", st->p[i]);
    fprintf(f, "nzeros %ld\n", st->nz);
    for (long i = 0; i < st->nz; ++i) fprintf(f, "%.17g\n", st->z[i]);
    fclose(f);
    if (rename(tmp, path)) { fprintf(stderr, "loop: cannot rename to %s\n", path); return 1; }
    return 0;
}
static int load_state(state *st, const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "loop: cannot open %s for resume\n", path); return 1; }
    char line[256]; long count;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#') continue;
        if (sscanf(line, "iter %ld", &st->iter) == 1) continue;
        if (sscanf(line, "xknown %ld", &st->Xknown) == 1) continue;
        if (sscanf(line, "kmin %lf", &st->kmin) == 1) continue;
        if (sscanf(line, "nprimes %ld", &count) == 1) {
            for (long i = 0; i < count; ++i) { long p; if (fscanf(f, "%ld\n", &p) == 1) p_push(st, p); }
            continue;
        }
        if (sscanf(line, "nzeros %ld", &count) == 1) {
            for (long i = 0; i < count; ++i) { double z; if (fscanf(f, "%lf\n", &z) == 1) z_push(st, z); }
            continue;
        }
    }
    fclose(f);
    fprintf(stderr, "loop: resumed iter=%ld X_known=%ld primes=%ld zeros=%ld\n",
            st->iter, st->Xknown, st->np, st->nz);
    return 0;
}

// ---- driver ----------------------------------------------------------------
void loop_opts_default(loop_opts *o) {
    o->state_path = "zzz-loop.state";
    o->resume = 0;
    o->max_iters = 0;
    o->nmax_zeros = 2000000;
    o->kmin = 5.0;
    o->kmin_set = 0;
    o->kmin_floor = 3.5;
    o->anneal = 1;
    o->seed_n = 0;
    o->contrast = 0;
    o->batch = 2000;
    o->fresh = 0;
    o->boxc = 0.375;
    o->verbose = 0;
}

static int file_exists(const char *path) {
    FILE *f = fopen(path, "r");
    if (f) { fclose(f); return 1; }
    return 0;
}

int loop_run(const loop_opts *o) {
    state st; memset(&st, 0, sizeof st);
    g_boxc = o->boxc;
    g_fitfree = o->contrast;

    // auto-resume: bare --loop continues from a checkpoint if one is present.
    int resume = o->fresh ? 0 : (o->resume ? 1 : file_exists(o->state_path));
    if (o->resume && !file_exists(o->state_path)) {
        fprintf(stderr, "loop: --resume but no checkpoint at %s\n", o->state_path);
        return 1;
    }

    st.kmin = o->kmin;
    if (resume) {
        if (load_state(&st, o->state_path)) return 1;   // restores st.kmin too
        if (o->kmin_set) st.kmin = o->kmin;             // explicit --loop-kmin overrides
    } else {
        long ns = (o->seed_n > 0 && o->seed_n < NSEED) ? o->seed_n : NSEED;
        if (ns < 6) { fprintf(stderr, "loop: seed must be >= 6 zeros\n"); return 1; }
        for (long i = 0; i < ns; ++i) z_push(&st, SEED[i]);
        fprintf(stderr, "loop: seeded with %ld zeros%s; loop is zeta-free and primality-free\n",
                ns, o->contrast ? " (fit-free peak scale)" : "");
    }

    signal(SIGINT, on_sigint);
    printf("# zzz --loop: primes derived from zeros (one per line). trajectory on stderr.\n");
    fflush(stdout);

    long prev_X = st.Xknown, stalled = 0, best_X = st.Xknown, no_improve = 0;
    while (!g_stop && (o->max_iters == 0 || st.iter < o->max_iters)) {
        st.iter++;
        // re-detect from just below the known frontier (lower primes already found),
        // so frequent re-detection stays cheap with small forward batches; Xcap is a
        // generous margin past the frontier (the scan stops at the reach anyway, and a
        // shortfall just gets caught up in the next iteration as the frontier grows)
        long x0 = st.Xknown - 16; if (x0 < 2) x0 = 2;
        long Xcap = st.Xknown + 1024;
        long Xnew = backward(&st, x0, Xcap);               // streams new primes
        if (Xnew > st.Xknown) st.Xknown = Xnew;

        // forward only a batch of zeros, so primes stream smoothly (Ctrl+C-able)
        long nt = n_frontier(st.Xknown, st.kmin, o->nmax_zeros);
        long target = nt;
        if (o->batch > 0 && target > st.nz + o->batch) target = st.nz + o->batch;
        for (long k = st.nz + 1; k <= target && !g_stop; ++k)
            z_push(&st, locate_march(k, &st, st.Xknown, st.z[st.nz - 1]));

        if (o->verbose || st.Xknown != prev_X || st.iter % 16 == 0)  // +heartbeat
            fprintf(stderr, "loop: iter %ld  X_known=%ld  primes=%ld  zeros=%ld  kmin=%.3f\n",
                    st.iter, st.Xknown, st.np, st.nz, st.kmin);

        if (st.iter % 8 == 0) save_state(&st, o->state_path);   // throttle the O(nz) rewrite

        if (st.Xknown > best_X) best_X = st.Xknown;
        if (st.Xknown > prev_X) {                          // progress: reach grew
            prev_X = st.Xknown; stalled = 0;
        } else if (st.nz >= nt) {                          // forward exhausted, reach flat
            if (++stalled >= 2) {
                // fixed point at this kmin: anneal the margin down and keep climbing
                no_improve = (st.Xknown >= best_X) ? 0 : no_improve + 1;
                if (!o->anneal || st.kmin <= o->kmin_floor + 1e-9 || no_improve >= 2) {
                    fprintf(stderr, "loop: %s ceiling reached -- X*=%ld at kmin=%.3f "
                            "(go further with a %s forward step)\n",
                            g_boxc > 0.0 ? "box-B" : "plain-B", best_X, st.kmin,
                            g_boxc > 0.0 ? "Weil" : "box-B/Weil");
                    break;
                }
                st.kmin *= 0.9;
                if (st.kmin < o->kmin_floor) st.kmin = o->kmin_floor;
                fprintf(stderr, "loop: stalled at X=%ld; annealing kmin -> %.3f and continuing\n",
                        st.Xknown, st.kmin);
                stalled = 0; prev_X = -1;
            }
        }                                                  // else: still filling forward
        if (st.nz >= o->nmax_zeros) {
            fprintf(stderr, "loop: hit zero cap %ld; raise --loop-nmax to continue\n", o->nmax_zeros);
            break;
        }
    }

    if (g_stop) fprintf(stderr, "loop: interrupted; checkpointing to %s\n", o->state_path);
    save_state(&st, o->state_path);                    // always persist the final state
    fprintf(stderr, "loop: done. primes=%ld (largest %ld)  zeros=%ld  state=%s\n",
            st.np, st.np ? st.p[st.np - 1] : 0, st.nz, o->state_path);

    free(st.z); free(st.p);
    return 0;
}
