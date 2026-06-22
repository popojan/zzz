# `zvk` — Vulkan-accelerated batch method B (experimental)

Original-zzz-spirit method B (GHY `P_X`), split across CPU and GPU exactly where
precision demands it. Computes a contiguous block of zero ordinates in parallel.

```
zvk <n0> <count> [X]      # ordinates of zeros #n0 .. #n0+count-1, primes p <= X
```

## Why the split (and why it's honest)

The prime sum's phase is exactly linear in the ordinate, `m·log p·t = m·log p·t_base
+ m·log p·δ`. At height `t_base` the **base phase** `φ = (m·log p·t_base) mod 2π`
needs `t_base` and `log p` to ~`log₁₀(t_base)` digits (you keep the fractional part
of a `10^D`-scale product) — that is irreducibly **ARB**, but it is computed *once
per window* and shared by all `count` zeros. Everything the GPU touches is `O(1)`:
`φ∈[0,2π)`, `ω=m·log p`, `δ~gap`. So:

- **ARB (CPU), once:** `t_base` = Lambert-W location of the block centre; fold `φ`
  for every `p^m ≤ X`; `c0 = N0(t_base)−(n_c−½)`, `ρ=N0'`, `ρ'=N0''`.
- **GPU (df32), per zero:** bisect `F_B(t_base+δ)` in `δ` (one thread per zero,
  lockstep). Host adds `t_base+δ` in ARB.

The printed ordinate carries `t_base`'s **leading digits** (exact, prime-free,
from the smooth Lambert-W anchor) plus the **gap-scale** refinement the primes
determine — precisely the band-saturation / França–LeClair structure
([`../../../doc/notes/exact-primes-inexact-zeros.md`](../../doc/notes/exact-primes-inexact-zeros.md)).
Accuracy of the refinement is `κ = gap·log X`: above `κ≈π` it sharpens; with a
fixed feasible `X` at large height it is gap-scale (all method B ever gives).

This is a **constant-factor** accelerator of method B at gap-scale — not a
substitute for `--boot`/`--weil` (those stay CPU/ARB for sub-gap precision), and
not competitive with a sieve for *producing primes*. Its sweet spot is computing a
**block of zeros at one height**, where the ARB fold amortises and the GPU saturates.

## Build (Linux)

```bash
# prereqs: libvulkan-dev, glslang-tools, plus the FLINT/GMP/MPFR zzz already needs
cd build
cmake -DBUILD_VULKAN=ON ..      # FetchContent pulls kompute (needs network once)
make zvk
./src/vulkan/zvk 1000 64 5000
```

`BUILD_VULKAN` is **OFF by default**; the normal `cmake .. && make` is unaffected.

## Status / caveats (read before trusting a digit)

- **VERIFIED (in its working range).** Builds on Fedora (kompute v0.9.0,
  glslang 16) and runs on an Intel UHD 620. Output matches **serial
  `./zzz --ghy -k K`** *and* an independent mpmath root of `F_B=k-½` to **~1e-6**
  at zeros #1000 / #1032 (X = p_K), and Odlyzko to gap-scale. So `zeromb.comp` is
  faithful method B. (Note for users: `./zzz --ghy` needs `-k` to set the prime
  count; without it X=0 and it returns N₀-only.)
- **Architecture = one workgroup per zero, cooperative df32 prime-reduction
  (Axis A × B).** Each zero gets a 256-thread workgroup; the threads reduce the
  O(P) prime sum together (so a *single* high-`-k` zero is parallel, not one
  serial lane), and the verified bisection drives off the shared reduction
  result (so the answer is bit-identical to the thread-per-zero baseline).
  Confirmed: a single zero at X=10⁶ (P≈79k) now computes (≈1419.415, near true)
  where the thread-per-zero kernel fell back to `t_base`.
- **Window tiling (done):** a single Taylor anchor `c0+ρδ+½ρ'δ²` is only small-δ
  valid, so a block of `count` zeros is split into sub-windows of `W` zeros
  (default 32; `zvk <n0> <count> [X] [W]`), each with its **own** ARB anchor +
  folded phases (`amp`/`om` are anchor-independent, folded once). Verified: at
  `count=256` the edge zeros match the single-zero path to ~1e-5 (vs ~0.01 drift
  with one anchor). `amp`/`om` tensors are persistent; only `φ` is re-folded and
  re-uploaded per sub-window.
- **fp64 phase fold (done):** profiling a large block (count=2048, X=2×10⁵) showed
  the bottleneck was the *per-sub-window ARB phase-fold* (~1.16 s of CPU
  `arb_mul`/mod over 64 windows × 18k primes), **not** the GPU. So `φ` is folded
  in **fp64** (`fmod(m·log p · t_base, 2π)`) whenever `t_base ≲ 10⁷` (the product
  still resolves in a double; abs err ~ `m·log p · t_base · 2⁻⁵²`), with the ARB
  fold kept as a height-guarded fallback above that. Result: same digits
  (verified at #1000 and #10⁶ vs serial `--ghy -k`), large high-X block **1.36 s
  → 0.46 s (~3×)**.
- **Axis C — dropped by profiling.** A parallel δ-grid replacing the sequential
  bisection would shorten GPU critical path, but the GPU was never the
  bottleneck (the ARB fold was), so it isn't pursued. (Newton was also avoided
  on purpose: `F_B` is oscillatory, so its derivative misleads the iteration.)
- **Precision:** the sum uses a **df32 (two-float) accumulator** with a
  **range-reduced `sin` argument** (mod 2π — GPU `sin` degrades past ~2π). These
  are defensive; at moderate X plain fp32+Kahan already matched mpmath to 1e-6.
- **kompute API drifts between versions** — `tensor()/algorithm<…>/OpTensorSyncDevice/
  OpTensorSyncLocal` match ~v0.8–0.9; adjust if you pin a different tag.
- **Push-constant layout** `PC` must stay byte-identical to the shader block
  (7 × 4-byte scalars, std430 — no padding here).
- **`arb_lambertw`** is assumed available (Arb / FLINT≥3). If your FLINT lacks it,
  invert `θ(t)/π+1 = n_c` by a few Newton steps instead.
- **Feasibility, not precision, is the wall:** `X` is capped by what you can fold
  and upload; at huge height that's far below `√(T/2π)`, so you get gap-scale — by
  design, not by defect.
