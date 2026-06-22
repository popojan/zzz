// zvk.cpp -- Vulkan-accelerated batch method-B zero computer (experimental).
//
//   zvk <n0> <count> [X] [W]
//
// Computes ordinates of zeros #n0 .. #n0+count-1 in the original-zzz spirit
// (method B / GHY P_X), CPU/GPU split where precision demands it:
//
//   ARB (CPU): per sub-window of <=W zeros, the Lambert-W anchor t_base and the
//     folded base phases phi[k] = (m log p . t_base) mod 2pi for every p^m <= X
//     (the part that genuinely needs ARB).  amp=1/(m p^{m/2}) and om=m log p are
//     t_base-independent -> folded ONCE for the whole block.
//   GPU (kompute, df32): one workgroup per zero cooperatively reduces the O(P)
//     prime sum and bisects F_B(t_base+delta) in the small offset (zeromb.comp).
//
// Tiling: a single Taylor anchor (c0+rho*d+rho1*d^2/2) is valid only for small
// delta, so the block is split into sub-windows of <=W zeros, each its own
// anchor + folded phases; the edge zeros stay accurate.  Output = t_base+delta
// in ARB (leading digits from the smooth anchor, gap-scale refinement from
// primes -- the band-saturation structure).  See src/vulkan/README.md.

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
        fprintf(stderr, "usage: %s <n0> <count> [X] [W]\n", argv[0]);
        return 2;
    }
    const char *n0_str = argv[1];
    long  count = strtol(argv[2], nullptr, 10);
    ulong X     = (argc > 3) ? strtoul(argv[3], nullptr, 10) : 100003UL;
    long  W     = (argc > 4) ? strtol(argv[4], nullptr, 10) : 32;   // zeros per sub-window
    if (count <= 0 || W <= 0) { fprintf(stderr, "zvk: count, W must be > 0\n"); return 2; }

    slong D    = (slong)strlen(n0_str);
    slong PREC = 256 + (slong)(D * 3.34) + 128;

    arb_t pi, twopi, tmp, aa, lp, tbase, nc, arg, w, tlo, theta, c0a, rhoa;
    acb_t z, lg;
    arb_init(pi); arb_init(twopi); arb_init(tmp); arb_init(aa); arb_init(lp);
    arb_init(tbase); arb_init(nc); arb_init(arg); arb_init(w); arb_init(tlo);
    arb_init(theta); arb_init(c0a); arb_init(rhoa);
    acb_init(z); acb_init(lg);
    arb_const_pi(pi, PREC);
    arb_mul_2exp_si(twopi, pi, 1);

    // ---- prime powers p^m <= X : amp, om (GPU, static), mlp = m log p (ARB) ----
    std::vector<float>  ampD, omD;
    std::vector<double> mlpd;            // m log p in fp64 (fast phase fold)
    std::vector<arb_struct> mlp;         // m log p in ARB  (fold fallback, huge t_base)
    n_primes_t it; n_primes_init(it);
    ulong p;
    while ((p = n_primes_next(it)) <= X) {
        ulong q = p, m = 1;
        arb_log_ui(lp, p, PREC);                         // log p, high precision
        double lpd = log((double)p);
        for (;;) {
            arb_struct e; arb_init(&e);
            arb_mul_ui(&e, lp, m, PREC);                 // m log p  (kept in ARB)
            mlp.push_back(e);
            mlpd.push_back((double)m * lpd);
            ampD.push_back((float)(1.0 / (m * pow((double)p, 0.5 * m))));
            omD.push_back((float)((double)m * lpd));
            if (q > X / p) break;
            q *= p; ++m;
        }
    }
    n_primes_clear(it);
    uint32_t P = (uint32_t)ampD.size();
    fprintf(stderr, "zvk: primes<=%lu -> %u prime-powers; block %ld zeros, W=%ld\n",
            X, P, count, W);

    // ---- GPU setup (persistent Manager + static amp/om tensors) ----
    std::vector<uint32_t> spirv = load_spv(ZVK_SPV);
    kp::Manager mgr;
    auto tAmp = mgr.tensor(ampD);
    auto tOm  = mgr.tensor(omD);

    std::vector<float> phiD(P);
    slong ndigits = D + 6;
    arb_t tn, dlt;
    arb_init(tn); arb_init(dlt);

    // ---- sub-window loop ----
    for (long done = 0; done < count; ) {
        long M = (count - done < W) ? (count - done) : W;

        // n_center = n0 + done + M/2  (sub-window centre)
        arb_set_str(nc, n0_str, PREC);
        arb_add_ui(nc, nc, (ulong)(done + M / 2), PREC);

        // Lambert-W anchor: t_base = 2pi (nc-11/8) / W( e^{-1}(nc-11/8) )
        arb_set_ui(tmp, 11); arb_div_ui(tmp, tmp, 8, PREC);
        arb_sub(arg, nc, tmp, PREC);
        arb_set_si(tmp, -1); arb_exp(tmp, tmp, PREC);
        arb_mul(w, arg, tmp, PREC);
        arb_lambertw(w, w, 0, PREC);
        arb_mul(tbase, twopi, arg, PREC);
        arb_div(tbase, tbase, w, PREC);

        // theta(t_base), c0, rho, rho1
        arb_set_d(tmp, 0.25);
        arb_mul_2exp_si(tlo, tbase, -1);
        acb_set_arb_arb(z, tmp, tlo);
        acb_lgamma(lg, z, PREC);
        arb_log(tmp, pi, PREC); arb_mul(tmp, tmp, tlo, PREC);
        arb_sub(theta, acb_imagref(lg), tmp, PREC);
        arb_div(c0a, theta, pi, PREC);
        arb_sub(c0a, c0a, nc, PREC);
        arb_set_d(tmp, 1.5); arb_add(c0a, c0a, tmp, PREC);     // c0 = N0(t_base)-(nc-1/2)
        arb_div(tmp, tbase, twopi, PREC); arb_log(tmp, tmp, PREC);
        arb_div(rhoa, tmp, twopi, PREC);                       // rho = N0'(t_base)
        double rho_d   = arf_get_d(arb_midref(rhoa), ARF_RND_NEAR);
        double tbase_d = arf_get_d(arb_midref(tbase), ARF_RND_NEAR);

        // fold phases for this anchor.  fp64 fmod is enough while the product
        // m log p . t_base stays resolvable in a double (abs err ~ m log p .
        // t_base . 2^-52 < 1e-7  =>  t_base < ~1e7 / (m log p)); above that the
        // fractional part is lost, so fall back to the ARB fold (slow but exact).
        if (tbase_d < 1.0e7) {
            const double TWOPI_D = 2.0 * M_PI;
            for (uint32_t k = 0; k < P; ++k)
                phiD[k] = (float)fmod(mlpd[k] * tbase_d, TWOPI_D);   // m log p . t_base mod 2pi
        } else {
            for (uint32_t k = 0; k < P; ++k) {
                arb_mul(tmp, &mlp[k], tbase, PREC);
                arb_div(aa, tmp, twopi, PREC); arb_floor(aa, aa, PREC);
                arb_submul(tmp, aa, twopi, PREC);
                phiD[k] = (float)arf_get_d(arb_midref(tmp), ARF_RND_NEAR);
            }
        }

        PC pc{};
        pc.P = P; pc.M = (uint32_t)M; pc.Wc = (int32_t)(M / 2);
        pc.rho = (float)rho_d; pc.rho1 = (float)(1.0 / (2.0 * M_PI * tbase_d));
        pc.c0  = (float)arf_get_d(arb_midref(c0a), ARF_RND_NEAR);
        pc.gap = (float)(1.0 / rho_d);

        auto tPhi = mgr.tensor(phiD);
        auto tOut = mgr.tensor(std::vector<float>((size_t)M, 0.0f));
        auto algo = mgr.algorithm<float, PC>(
            { tPhi, tAmp, tOm, tOut }, spirv,        // binding order: 0=phi,1=amp,2=om,3=out
            kp::Workgroup({ (uint32_t)M, 1, 1 }),    // one workgroup (256 threads) per zero
            std::vector<float>{}, std::vector<PC>{ pc });
        mgr.sequence()
            ->record<kp::OpTensorSyncDevice>({ tPhi, tAmp, tOm })
            ->record<kp::OpAlgoDispatch>(algo)
            ->record<kp::OpTensorSyncLocal>({ tOut })
            ->eval();
        std::vector<float> delta = tOut->vector();

        for (long j = 0; j < M; ++j) {
            arb_set_d(dlt, (double)delta[(size_t)j]);
            arb_add(tn, tbase, dlt, PREC);
            char *s = arb_get_str(tn, ndigits, ARB_STR_NO_RADIUS);
            printf("%s\n", s);
            flint_free(s);
        }
        done += M;
    }

    for (uint32_t k = 0; k < P; ++k) arb_clear(&mlp[k]);
    arb_clear(tn); arb_clear(dlt);
    arb_clear(pi); arb_clear(twopi); arb_clear(tmp); arb_clear(aa); arb_clear(lp);
    arb_clear(tbase); arb_clear(nc); arb_clear(arg); arb_clear(w); arb_clear(tlo);
    arb_clear(theta); arb_clear(c0a); arb_clear(rhoa);
    acb_clear(z); acb_clear(lg);
    return 0;
}
