// loopvk.cpp -- the self-paving zeros<->primes loop, BOTH directions on the GPU.
//
//   loopvk [maxX] [kmin] [statefile] [kmin_floor]
//
// A Vulkan port of `loop.c` (see doc/notes/zeros-primes-bootstrap.md): from a
// finite seed of zero ordinates it detects primes (backward: psi.comp, a
// reduction over the zeros) and locates more zeros (forward: zeromb.comp, method
// B, a reduction over the discovered primes), alternating.  No zeta, no
// primality test in the loop body.  Pure double precision (modest heights), like
// loop.c; both O(.) sums run on the GPU.
//
// Feature parity with `zzz --loop`: auto-resume + Ctrl+C checkpoint to the SAME
// `zzz-loop.state` format (so a CPU run and a GPU run are interoperable -- resume
// either on the other), and auto-anneal of kmin on stall.  Verification is
// post-hoc and clearly separated.

#include <kompute/Kompute.hpp>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <csignal>
#include <fstream>
#include <vector>

#ifndef ZEROMB_SPV
#define ZEROMB_SPV "zeromb.comp.spv"
#endif
#ifndef PSI_SPV
#define PSI_SPV "psi.comp.spv"
#endif

// ---- seed: first 40 nontrivial zero ordinates (given data) -----------------
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
#define NSEED ((long)(sizeof(SEED)/sizeof(SEED[0])))

static const double TAU = 0.30, TAULO = 0.12, ENVK = 20.0 / 4.0;

struct PCfwd { uint32_t P, M; int32_t Wc; float rho, rho1, c0, gap; };
struct PCbwd { uint32_t NZ, X0, NX; float envc; };

static std::vector<uint32_t> load_spv(const char *path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) { fprintf(stderr, "loopvk: cannot open '%s'\n", path); exit(1); }
    size_t n = (size_t)f.tellg(); std::vector<uint32_t> v(n / 4);
    f.seekg(0); f.read((char *)v.data(), n); return v;
}

// ---- fp64 primitives (port of loop.c) --------------------------------------
static double ei(double zz) {
    double s = 0.57721566490153286061 + log(zz), t = 1.0;
    for (int k = 1; k < 400; ++k) { t *= zz / k; double a = t / k; s += a;
        if (fabs(a) < 1e-15 * fabs(s) && k > zz) break; } return s;
}
static double li(double x) { return ei(log(x)); }
static double N0(double t)  { double u = t / (2 * M_PI); return u * log(u / M_E) + 0.875; }
static double N0p(double t) { return log(t / (2 * M_PI)) / (2 * M_PI); }      // density
static double t_of_index(long n) {                                            // N0(t)=n-1/2
    double t = 2 * M_PI * n / log((double)(n > 2 ? n : 3));
    for (int i = 0; i < 80; ++i) { double f = N0(t) - (n - 0.5); t -= f / N0p(t); }
    return t;
}
static long n_frontier(long X, double kmin, long nmax) {
    if (X < 2) return 0;
    double T = 2 * M_PI * pow((double)X, 2 * M_PI / kmin);
    double u = T / (2 * M_PI); long n = (long)floor(u * (log(u) - 1.0));
    return n > nmax ? nmax : n;
}

// ---- loop state (zeros + primes), with loop.c-compatible checkpoint --------
static std::vector<long>   primes;       // ascending
static std::vector<double> z;            // the loop's own zero ordinates

static int p_has(long x) {
    long lo = 0, hi = (long)primes.size() - 1;
    while (lo <= hi) { long m = (lo + hi) / 2; if (primes[m] == x) return 1;
        if (primes[m] < x) lo = m + 1; else hi = m - 1; } return 0;
}
static int is_known_pow(long x) {
    for (long p : primes) { if (p * p > x) break; long q = p * p;
        while (q < x) q *= p; if (q == x) return 1; } return 0;
}

static volatile sig_atomic_t g_stop = 0;
static void on_sigint(int) { g_stop = 1; }
static int file_exists(const char *p) { FILE *f = fopen(p, "r"); if (f) { fclose(f); return 1; } return 0; }

// byte-identical to loop.c save_state -> CPU/GPU interop on the same file
static int save_state(const char *path, long iter, long Xknown, double kmin) {
    char tmp[1024]; snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) { fprintf(stderr, "loopvk: cannot write %s\n", tmp); return 1; }
    fprintf(f, "# zzz --loop checkpoint\niter %ld\nxknown %ld\nkmin %.17g\n", iter, Xknown, kmin);
    fprintf(f, "nprimes %zu\n", primes.size());
    for (long p : primes) fprintf(f, "%ld\n", p);
    fprintf(f, "nzeros %zu\n", z.size());
    for (double zz : z) fprintf(f, "%.17g\n", zz);
    fclose(f);
    if (rename(tmp, path)) { fprintf(stderr, "loopvk: cannot rename to %s\n", path); return 1; }
    return 0;
}
static int load_state(const char *path, long *iter, long *Xknown, double *kmin) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "loopvk: cannot open %s for resume\n", path); return 1; }
    char line[256]; long count;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#') continue;
        if (sscanf(line, "iter %ld", iter) == 1) continue;
        if (sscanf(line, "xknown %ld", Xknown) == 1) continue;
        if (sscanf(line, "kmin %lf", kmin) == 1) continue;
        if (sscanf(line, "nprimes %ld", &count) == 1) {
            for (long i = 0; i < count; ++i) { long p; if (fscanf(f, "%ld\n", &p) == 1) primes.push_back(p); }
            continue;
        }
        if (sscanf(line, "nzeros %ld", &count) == 1) {
            for (long i = 0; i < count; ++i) { double zz; if (fscanf(f, "%lf\n", &zz) == 1) z.push_back(zz); }
            continue;
        }
    }
    fclose(f);
    fprintf(stderr, "loopvk: resumed iter=%ld X_known=%ld primes=%zu zeros=%zu\n",
            *iter, *Xknown, primes.size(), z.size());
    return 0;
}

int main(int argc, char **argv) {
    long        maxX   = (argc > 1) ? strtol(argv[1], nullptr, 10) : 500;       // 0 = until ceiling
    double      cliKmin = (argc > 2) ? strtod(argv[2], nullptr) : 0.0;          // 0 = auto
    const char *spath  = (argc > 3) ? argv[3] : "zzz-loop.state";
    double      kmin_floor = (argc > 4) ? strtod(argv[4], nullptr) : 3.5;
    const long  W = 32, BATCH = 4000, NMAX = 2000000;

    auto spv_fwd = load_spv(ZEROMB_SPV);
    auto spv_bwd = load_spv(PSI_SPV);
    kp::Manager mgr;

    long iter = 0, Xknown = 1; double kmin = cliKmin > 0 ? cliKmin : 4.0;
    if (file_exists(spath)) {
        if (load_state(spath, &iter, &Xknown, &kmin)) return 1;
        if (cliKmin > 0) kmin = cliKmin;                  // explicit kmin overrides resumed
    } else {
        z.assign(SEED, SEED + NSEED);
        fprintf(stderr, "loopvk: seeded %ld zeros; zeta-free & primality-free; GPU both ways\n", NSEED);
    }
    signal(SIGINT, on_sigint);
    printf("# loopvk: primes from zeros on the GPU (one/line); state=%s\n", spath); fflush(stdout);

    // ---- backward: detect prime powers from z[] over [x0,Xcap] on the GPU ----
    auto backward = [&](long x0, long Xcap) -> long {
        if (x0 < 2) x0 = 2;
        long NX = Xcap - x0 + 1; if (NX < 1) return Xknown;
        auto tZ = mgr.tensor(std::vector<float>(z.begin(), z.end()));
        auto tR = mgr.tensor(std::vector<float>((size_t)NX, 0.0f));
        PCbwd pc{ (uint32_t)z.size(), (uint32_t)x0, (uint32_t)NX,
                  (float)(ENVK * li((double)z.size() - 4.0)) };
        auto algo = mgr.algorithm<float, PCbwd>({ tZ, tR }, spv_bwd,
                      kp::Workgroup({ (uint32_t)NX, 1, 1 }), {}, { pc });
        mgr.sequence()->record<kp::OpTensorSyncDevice>({ tZ })
            ->record<kp::OpAlgoDispatch>(algo)->record<kp::OpTensorSyncLocal>({ tR })->eval();
        std::vector<float> r = tR->vector();
        long Xk = x0 - 1;
        for (long x = x0; x <= Xcap; ++x) {
            double rr = r[(size_t)(x - x0)]; int kp = is_known_pow(x);
            if (rr > TAU) {
                if (!kp && !p_has(x)) {
                    long lo = 0, hi = (long)primes.size();
                    while (lo < hi) { long m = (lo + hi) / 2; if (primes[m] < x) lo = m + 1; else hi = m; }
                    primes.insert(primes.begin() + lo, x);
                    printf("%ld\n", x); fflush(stdout);
                }
                Xk = x;
            } else if (rr > TAULO) { if (kp) Xk = x; else break; }
            else Xk = x;
        }
        return Xk;
    };

    // ---- forward: locate zeros nz+1..target via method B (zeromb) on the GPU --
    auto forward = [&](long target) {
        std::vector<float>  ampD, omD; std::vector<double> mlpd;
        for (long p : primes) { if (p > Xknown) break; double lp = log((double)p);
            long q = p, m = 1; for (;;) { ampD.push_back((float)(1.0 / (m * pow((double)p, 0.5 * m))));
                omD.push_back((float)((double)m * lp)); mlpd.push_back((double)m * lp);
                if (q > Xknown / p) break; q *= p; ++m; } }
        uint32_t P = (uint32_t)ampD.size(); if (P == 0) return;
        auto tAmp = mgr.tensor(ampD); auto tOm = mgr.tensor(omD);
        std::vector<float> phiD(P);
        for (long base = (long)z.size() + 1; base <= target && !g_stop; ) {
            long M = (target - base + 1 < W) ? (target - base + 1) : W;
            long nc = base + M / 2; double tb = t_of_index(nc);
            double rho = N0p(tb), c0 = N0(tb) - (nc - 0.5);
            for (uint32_t k = 0; k < P; ++k) phiD[k] = (float)fmod(mlpd[k] * tb, 2.0 * M_PI);
            PCfwd pc{ P, (uint32_t)M, (int32_t)(M / 2), (float)rho,
                      (float)(1.0 / (2.0 * M_PI * tb)), (float)c0, (float)(1.0 / rho) };
            auto tPhi = mgr.tensor(phiD);
            auto tOut = mgr.tensor(std::vector<float>((size_t)M, 0.0f));
            auto algo = mgr.algorithm<float, PCfwd>({ tPhi, tAmp, tOm, tOut }, spv_fwd,
                          kp::Workgroup({ (uint32_t)M, 1, 1 }), {}, { pc });
            mgr.sequence()->record<kp::OpTensorSyncDevice>({ tPhi, tAmp, tOm })
                ->record<kp::OpAlgoDispatch>(algo)->record<kp::OpTensorSyncLocal>({ tOut })->eval();
            std::vector<float> d = tOut->vector();
            for (long j = 0; j < M; ++j) z.push_back(tb + (double)d[(size_t)j]);
            base += M;
        }
    };

    // ---- the loop (auto-anneal, like loop.c) ----
    long prev_X = Xknown, stalled = 0, best_X = Xknown, no_improve = 0;
    while (!g_stop && (maxX <= 0 || Xknown < maxX)) {
        ++iter;
        long Xnew = backward(Xknown - 16, Xknown + 1024);
        if (Xnew > Xknown) Xknown = Xnew;

        long nt = n_frontier(Xknown, kmin, NMAX);
        long target = ((long)z.size() + BATCH < nt) ? (long)z.size() + BATCH : nt;
        if (target > (long)z.size()) forward(target);

        fprintf(stderr, "loopvk: iter %ld  X_known=%ld  primes=%zu  zeros=%zu  kmin=%.3f\n",
                iter, Xknown, primes.size(), z.size(), kmin);
        if (iter % 8 == 0) save_state(spath, iter, Xknown, kmin);

        if (Xknown > best_X) best_X = Xknown;
        if (Xknown > prev_X) { prev_X = Xknown; stalled = 0; }
        else if ((long)z.size() >= nt) {
            if (++stalled >= 2) {
                no_improve = (Xknown >= best_X) ? 0 : no_improve + 1;
                if (kmin <= kmin_floor + 1e-9 || no_improve >= 2) {
                    fprintf(stderr, "loopvk: plain-B ceiling X*=%ld at kmin=%.3f\n", best_X, kmin);
                    break;
                }
                kmin *= 0.9; if (kmin < kmin_floor) kmin = kmin_floor;
                fprintf(stderr, "loopvk: stalled at X=%ld; annealing kmin -> %.3f\n", Xknown, kmin);
                stalled = 0; prev_X = -1;
            }
        }
        if ((long)z.size() >= NMAX) { fprintf(stderr, "loopvk: zero cap %ld hit\n", NMAX); break; }
    }
    if (g_stop) fprintf(stderr, "loopvk: interrupted; checkpointing to %s\n", spath);
    save_state(spath, iter, Xknown, kmin);

    // NO verification here: like loop.c, the binary stays primality-test-free
    // (the loop body's only "is x a prime power" inputs are the psi' signal and
    // is_known_pow, which is integer bookkeeping against the loop's own primes).
    // Check the checkpoint externally, post-hoc:
    //     python3 doc/ghy/check-loop.py <state>
    fprintf(stderr, "loopvk: done. primes=%zu (largest %ld)  zeros=%zu  state=%s\n",
            primes.size(), primes.empty() ? 0 : primes.back(), z.size(), spath);
    fprintf(stderr, "loopvk: verify externally -- python3 doc/ghy/check-loop.py %s\n", spath);
    return 0;
}
