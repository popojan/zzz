// zvk.cpp -- Vulkan-accelerated batch method-B zero computer (experimental).
//
//   zvk <n0> <count> [X]
//
// Computes the ordinates of zeros #n0 .. #n0+count-1 in the original-zzz spirit
// (method B / GHY P_X), split across CPU and GPU exactly where precision demands:
//
//   ARB (CPU), once per window:
//     * t_base   = Lambert-W asymptotic location of the block-centre zero,
//                  to full precision (a number of ~log10(n0) digits);
//     * fold the huge base phases  phi[k] = (m log p . t_base) mod 2pi  for every
//       prime power p^m <= X  (this is what genuinely needs ARB -- the fractional
//       part of a ~10^D-scale product);
//     * c0 = N0(t_base) - (n_centre - 1/2),  rho = N0',  rho1 = N0''   (all O(1)).
//   GPU (df32, Vulkan via kompute):
//     * one thread per zero bisects F_B(t_base+delta) in the small offset delta,
//       reusing the folded phases.  Returns delta[i]; the host adds t_base+delta
//       in ARB so the output carries t_base's leading digits + the gap-scale
//       refinement the primes determine (exactly the band-saturation structure).
//
// This is a constant-factor accelerator of method B at gap-scale accuracy; it is
// NOT a substitute for --boot/--weil (those stay CPU/ARB for sub-gap precision).
// See src/vulkan/README.md.  UNTESTED in this tree (needs Vulkan+glslang+kompute).

#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/ulong_extras.h>

#include <kompute/Kompute.hpp>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <fstream>
#include <vector>

#ifndef ZVK_SPV
#define ZVK_SPV "zeromb.comp.spv"
#endif

struct PC {                       // must match the std430 push_constant block
    uint32_t P, M;
    int32_t  Wc;
    float    rho, rho1, c0, gap;
};

static std::vector<uint32_t> load_spv(const char *path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    if (!f) { fprintf(stderr, "zvk: cannot open SPIR-V '%s'\n", path); exit(1); }
    size_t n = (size_t)f.tellg();
    std::vector<uint32_t> v(n / 4);
    f.seekg(0); f.read(reinterpret_cast<char *>(v.data()), n);
    return v;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <n0> <count> [X]\n", argv[0]);
        return 2;
    }
    const char *n0_str = argv[1];
    long  count = strtol(argv[2], nullptr, 10);
    ulong X     = (argc > 3) ? strtoul(argv[3], nullptr, 10) : 100003UL;
    if (count <= 0) { fprintf(stderr, "zvk: count must be > 0\n"); return 2; }

    // working precision: enough that the folded phase is good to ~fp32 abs.
    slong D = (slong)strlen(n0_str);
    slong PREC = 256 + (slong)(D * 3.34) + 128;

    arb_t nc, tbase, arg, w, twopi, pi, lp, phase, tmp, theta, c0a, rhoa, t_lo;
    acb_t z, lg;
    arb_init(nc); arb_init(tbase); arb_init(arg); arb_init(w);
    arb_init(twopi); arb_init(pi); arb_init(lp); arb_init(phase);
    arb_init(tmp); arb_init(theta); arb_init(c0a); arb_init(rhoa); arb_init(t_lo);
    acb_init(z); acb_init(lg);

    arb_const_pi(pi, PREC);
    arb_mul_2exp_si(twopi, pi, 1);                       // 2 pi

    // n_centre = n0 + count/2
    if (arb_set_str(nc, n0_str, PREC)) { fprintf(stderr, "zvk: bad n0\n"); return 2; }
    arb_add_ui(nc, nc, (ulong)(count / 2), PREC);

    // Lambert-W asymptotic:  t_base = 2 pi (nc - 11/8) / W( e^{-1} (nc - 11/8) )
    arb_set_ui(tmp, 11); arb_div_ui(tmp, tmp, 8, PREC);
    arb_sub(arg, nc, tmp, PREC);                         // nc - 11/8
    arb_set_si(tmp, -1); arb_exp(tmp, tmp, PREC);        // e^{-1}
    arb_mul(w, arg, tmp, PREC);                          // e^{-1}(nc-11/8)
    arb_lambertw(w, w, 0, PREC);                         // principal branch W_0
    arb_mul(tbase, twopi, arg, PREC);
    arb_div(tbase, tbase, w, PREC);                      // t_base

    // theta(t_base) = Im logGamma(1/4 + i t_base/2) - (t_base/2) log pi
    arb_set_d(tmp, 0.25);
    arb_mul_2exp_si(t_lo, tbase, -1);                    // t_base/2
    acb_set_arb_arb(z, tmp, t_lo);
    acb_lgamma(lg, z, PREC);
    arb_log(tmp, pi, PREC);
    arb_mul(tmp, tmp, t_lo, PREC);                       // (t_base/2) log pi
    arb_sub(theta, acb_imagref(lg), tmp, PREC);

    // c0 = theta/pi - nc + 3/2   (= N0(t_base) - (nc - 1/2))
    arb_div(c0a, theta, pi, PREC);
    arb_sub(c0a, c0a, nc, PREC);
    arb_set_d(tmp, 1.5); arb_add(c0a, c0a, tmp, PREC);

    // rho = (1/2pi) log(t_base/2pi)
    arb_div(tmp, tbase, twopi, PREC);
    arb_log(tmp, tmp, PREC);
    arb_div(rhoa, tmp, twopi, PREC);

    double rho_d  = arf_get_d(arb_midref(rhoa), ARF_RND_NEAR);
    double tbase_d_log = arf_get_d(arb_midref(tbase), ARF_RND_NEAR);  // only for rho1 ~0
    PC pc{};
    pc.M    = (uint32_t)count;
    pc.Wc   = (int32_t)(count / 2);
    pc.rho  = (float)rho_d;
    pc.rho1 = (float)(1.0 / (2.0 * M_PI * tbase_d_log));   // N0'' ~ tiny
    pc.c0   = (float)arf_get_d(arb_midref(c0a), ARF_RND_NEAR);
    pc.gap  = (float)(1.0 / rho_d);

    // fold base phases for every prime power p^m <= X
    std::vector<float> phiD, ampD, omD;
    n_primes_t it; n_primes_init(it);
    ulong p;
    while ((p = n_primes_next(it)) <= X) {
        ulong q = p, m = 1;
        for (;;) {
            arb_log_ui(lp, p, PREC);                     // log p (high precision)
            arb_mul_ui(phase, lp, m, PREC);
            arb_mul(phase, phase, tbase, PREC);          // m log p . t_base
            arb_div(tmp, phase, twopi, PREC);
            arb_floor(tmp, tmp, PREC);
            arb_submul(phase, tmp, twopi, PREC);         // mod 2pi  -> [0,2pi)
            phiD.push_back((float)arf_get_d(arb_midref(phase), ARF_RND_NEAR));
            ampD.push_back((float)(1.0 / (m * pow((double)p, 0.5 * m))));
            omD.push_back((float)(m * log((double)p)));
            if (q > X / p) break;                        // next power overflow-safe
            q *= p; ++m;
        }
    }
    n_primes_clear(it);
    pc.P = (uint32_t)phiD.size();
    fprintf(stderr, "zvk: t_base ~ %g, primes<=%lu -> %u prime-powers, window=%ld\n",
            arf_get_d(arb_midref(tbase), ARF_RND_NEAR), X, pc.P, count);

    // ---- GPU dispatch (kompute) -------------------------------------------
    std::vector<uint32_t> spirv = load_spv(ZVK_SPV);
    kp::Manager mgr;
    auto tPhi = mgr.tensor(phiD);
    auto tAmp = mgr.tensor(ampD);
    auto tOm  = mgr.tensor(omD);
    auto tOut = mgr.tensor(std::vector<float>((size_t)count, 0.0f));
    auto algo = mgr.algorithm<float, PC>(
        { tPhi, tAmp, tOm, tOut }, spirv,
        kp::Workgroup({ (uint32_t)count, 1, 1 }),   // one workgroup (256 threads) per zero
        std::vector<float>{}, std::vector<PC>{ pc });
    mgr.sequence()
        ->record<kp::OpTensorSyncDevice>({ tPhi, tAmp, tOm })
        ->record<kp::OpAlgoDispatch>(algo)
        ->record<kp::OpTensorSyncLocal>({ tOut })
        ->eval();
    std::vector<float> delta = tOut->vector();

    // ---- output: t_n = t_base + delta (in ARB) ----------------------------
    slong ndigits = D + 6;
    arb_t tn, dlt;
    arb_init(tn); arb_init(dlt);
    for (long i = 0; i < count; ++i) {
        arb_set_d(dlt, (double)delta[(size_t)i]);
        arb_add(tn, tbase, dlt, PREC);
        char *s = arb_get_str(tn, ndigits, ARB_STR_NO_RADIUS);
        printf("%s\n", s);
        flint_free(s);
    }
    arb_clear(tn); arb_clear(dlt);

    arb_clear(nc); arb_clear(tbase); arb_clear(arg); arb_clear(w);
    arb_clear(twopi); arb_clear(pi); arb_clear(lp); arb_clear(phase);
    arb_clear(tmp); arb_clear(theta); arb_clear(c0a); arb_clear(rhoa); arb_clear(t_lo);
    acb_clear(z); acb_clear(lg);
    return 0;
}
