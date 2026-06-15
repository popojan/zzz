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
static const double TAU   = 0.30;   // strong-peak threshold (prime power)
static const double TAULO = 0.12;   // ambiguous band floor
static const double ENVK  = 20.0 / 4.0;

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
        if (st->p[m] == x) return 1; if (st->p[m] < x) lo = m + 1; else hi = m - 1; }
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
static double envelope(double x, long n) { return ENVK * li((double)n - 4.0) * log(x) / x; }

// psi'(x) ~ 1 + 1/(x-x^3) - (4/sqrt x) sum_k cos(gamma_k log x)
static double psi_prime(double x, const double *z, long n) {
    double lx = log(x), s = 0.0;
    for (long k = 0; k < n; ++k) s += cos(z[k] * lx);
    return 1.0 + 1.0 / (x - x * x * x) - (4.0 / sqrt(x)) * s;
}

// fit-free (CFAR) detector: psi'(x) rises >= C local MADs above the median of its
// neighbors. No envelope, no Li -- only a dimensionless SNR threshold. Needs a
// small prime warm-start (at tiny x almost every integer is a prime power, so the
// local median has no composite floor).
static const long  SEED_PRIMES[] = {2,3,5,7,11,13,17,19,23,29};
#define NSEEDP ((long)(sizeof(SEED_PRIMES)/sizeof(SEED_PRIMES[0])))
static const double CHI = 4.0;   // clear prime power (MADs above local median)
static const double CLO = 2.0;   // ambiguous band floor
static int cmp_d(const void *a, const void *b) {
    double d = *(const double *)a - *(const double *)b; return d < 0 ? -1 : d > 0 ? 1 : 0;
}
static double median_d(double *v, int m) {
    qsort(v, m, sizeof(double), cmp_d);
    return m == 0 ? 0.0 : (m % 2 ? v[m / 2] : 0.5 * (v[m / 2 - 1] + v[m / 2]));
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
// F_B(t) = N0(t) + Im log P_X / pi,  Im log P_X = -sum_{p^m<=X} (1/m) p^{-m/2} sin(m t log p)
static double FB(double t, const state *st, long X) {
    double acc = 0.0;
    for (long i = 0; i < st->np; ++i) {
        long p = st->p[i];
        if (p > X) break;
        double lp = log((double)p), pm = (double)p, half = 1.0;
        for (long m = 1; ; ++m) {
            half = pow((double)p, -0.5 * m);
            acc += (half / m) * sin(m * t * lp);
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
    if (flo * fhi > 0) return prev + gap;                  // fallback: predicted
    for (int it = 0; it < 50; ++it) {
        double mid = 0.5 * (lo + hi), fm = FB(mid, st, X) - tgt;
        if (flo * fm <= 0) { hi = mid; fhi = fm; } else { lo = mid; flo = fm; }
    }
    return 0.5 * (lo + hi);
}

// backward: scan x from x_start upward (primes below x_start already known),
// detect prime powers from the psi' peak, append new primes; conservative stop
// at the first ambiguous non-(known-power). Returns the reach X_known.
static long backward(state *st, long x_start, long Xcap) {
    long n = st->nz, Xk = x_start - 1;
    for (long x = x_start; x <= Xcap; ++x) {
        double r = psi_prime((double)x, st->z, n) / envelope((double)x, n);
        int kp = is_known_pow(x, st);
        if (r > TAU) {                                     // clear prime power
            if (!kp && !p_has(st, x)) { p_push(st, x); printf("%ld\n", x); fflush(stdout); }
            Xk = x;
        } else if (r > TAULO) {
            if (kp) { Xk = x; } else break;                // ambiguous non-power -> stop
        } else {
            Xk = x;                                        // clear composite
        }
    }
    return Xk;
}

// fit-free variant: local-contrast (CFAR) detection, no envelope/Li.
// psi' is evaluated lazily over the local window only, so cost tracks the scan
// length (which stops at the reach), not the full [x_start, Xcap] span.
static long backward_contrast(state *st, long x_start, long Xcap) {
    long n = st->nz, Xk = x_start - 1;
    const int W = 6;
    double buf[2 * W + 2], tmp[2 * W + 2];
    for (long x = x_start; x <= Xcap; ++x) {
        int m = 0;
        for (long j = x - W; j <= x + W; ++j)
            if (j >= 2 && j != x) buf[m++] = psi_prime((double)j, st->z, n);
        if (m < 3) { Xk = x; continue; }
        memcpy(tmp, buf, m * sizeof(double));
        double base = median_d(tmp, m);
        for (int i = 0; i < m; ++i) tmp[i] = fabs(buf[i] - base);
        double mad = median_d(tmp, m);
        double c = (psi_prime((double)x, st->z, n) - base) / (mad + 1e-9);
        int kp = is_known_pow(x, st);
        if (c > CHI) {
            if (!kp && !p_has(st, x)) { p_push(st, x); printf("%ld\n", x); fflush(stdout); }
            Xk = x;
        } else if (c > CLO) {
            if (kp || p_has(st, x)) Xk = x; else break;    // ambiguous & unknown -> stop
        } else {
            Xk = x;
        }
    }
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
    o->verbose = 0;
}

static int file_exists(const char *path) {
    FILE *f = fopen(path, "r");
    if (f) { fclose(f); return 1; }
    return 0;
}

int loop_run(const loop_opts *o) {
    state st; memset(&st, 0, sizeof st);

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
        long nmin = o->contrast ? 2 : 6;     // envelope's Li(n-4) needs n>5
        if (ns < nmin) { fprintf(stderr, "loop: seed must be >= %ld zeros\n", nmin); return 1; }
        for (long i = 0; i < ns; ++i) z_push(&st, SEED[i]);
        if (o->contrast)                     // warm-start the prime list (small x is power-dense)
            for (long i = 0; i < NSEEDP; ++i) { p_push(&st, SEED_PRIMES[i]); printf("%ld\n", SEED_PRIMES[i]); }
        fprintf(stderr, "loop: seeded with %ld zeros%s; loop is zeta-free and primality-free\n",
                ns, o->contrast ? " + 10 warm-start primes (contrast detector)" : "");
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
        long Xnew = o->contrast ? backward_contrast(&st, x0, Xcap)   // streams new primes
                                : backward(&st, x0, Xcap);
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
                    fprintf(stderr, "loop: plain-B ceiling reached -- X*=%ld at kmin=%.3f "
                            "(go further with a sharper/Weil forward step)\n", best_X, st.kmin);
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
