# zzz live loop on Shadertoy (Stage 2)

**Live:** https://www.shadertoy.com/view/s3BSzD

Self-paving primes-from-zeros bootstrap running entirely on the GPU. From **no
baked-in data at all** (2 is prime by definition; the first zero comes from the
method-B sum over {2}) it runs the explicit formula as a closed loop — no sieve, no primality test, no ζ evaluation. It uses the decision
rule of the `zzz` project's `--loop` (`loop.c`) in a per-frame form; `stage2-sim.c`
replays the tabs frame by frame in fp32 and checks them against a sieve.

**Project:** https://github.com/popojan/zzz — `zzz` ("Zeta Zeros Zeal") computes
Riemann-zeta zeros from primes (and the reverse) via a ζ-evaluation-free
explicit-formula counter. The maths behind this demo lives in
`doc/notes/zeros-primes-bootstrap.md`, `doc/notes/exact-primes-inexact-zeros.md`
and `doc/notes/loop-margins.md` (why the decision rule is what it is).

## Setup (Common + two buffers + Image)
1. **Common** tab   ← `stage2-common.glsl`  (shared code, *not* a buffer; no `mainImage`)
2. **Buffer A** tab ← `stage2-bufferA.glsl` ; iChannel0 = **Buffer A** (self-feedback = state)
3. **Buffer B** tab ← `stage2-bufferB.glsl` ; iChannel0 = **Buffer A** (per-column trace cache)
4. **Image** tab    ← `stage2-image.glsl`   ; iChannel0 = **Buffer A**, iChannel1 = **Buffer B**

The single most common failure is forgetting that **Buffer A's iChannel0 must
point at Buffer A itself** — that self-feedback *is* the persistent state. Symptom
if it's missing: it runs but only ever shows "2". After wiring, rewind (⏮) to
restart.

Buffer A holds the loop state in one RGBA32F texture (zeros in `.r`, prime-indicator
in `.g`, control in texels 0–1). Buffer B caches the detector trace `r(x)` once per
screen column (row 0) so the Image pass samples it with 3 `texelFetch` instead of a
384-term `cos` sum per pixel — ~100× cheaper drawing, the right move on weak GPUs.

## What you should see
From the single prime 2 the number line grows on its own; primes light up green
as the loop discovers them, the current one flares, its value prints.

## Knobs
- Common: `NLIM` zero cap, shared by both buffers (8000 → primes to 1502, 16000 → 2716,
  60000 → 8585; raise carefully — deep zeros lengthen the detector sum and can trip the
  browser's shader time limit). `HOLD` seconds per prime. `DTAU`, `UMIN`, `BOXC` are the
  decision rule — rerun `stage2-sim.c` before changing them.
- Buffer A: `BATCH` zeros/frame, `BAND` frontier/frame (bootstrap speed).
- Buffer B has no knob: the trace sums `2.8·xHi + 128` zeros (capped at `nz`), enough to
  resolve the visible window at any depth.

## Status
Stage 1 (`primes-from-zeros.glsl`) is the baked-snapshot renderer (verified).
Stage 2 is the live loop. `stage2-sim.c` replays its Common + Buffer A logic on the
CPU (fp32, synchronous per-frame update) and compares with a sieve: exactly the
primes to 1502 at `NLIM = 8000` and to 8585 at `NLIM = 60000`, with every prime
reading ≥ 0.70 and every composite ≤ 0.21 against the 0.40 threshold. All four tabs
compile as GLSL ES 3.00 (`glslangValidator` with a Shadertoy-style wrapper); the
changed tabs have not been run on a GPU.

**Changed 2026-11 — the published copy (s3BSzD) still has the old tabs until all four
are pasted again** (Common lost `envf`, so Buffer B and Image must follow). The old
Buffer A decided each integer once, from a single full reading against the
$\mathrm{Li}$ envelope at 0.30, at a reach of `nz/2.7`; replayed, it flags 12
composites as primes up to 2962 (417, 1449, 1865, 2476, 2592, …), and its prime and
composite readings overlap (lowest prime 0.49, highest composite 0.50). Now:
- the reading is the **median over six truncations** of the zero list, each scaled
  by its own top zero (no fitted constant), threshold 0.40;
- an integer is decided only while `gamma(nz/2) >= 3·x` (the reach in the old rule
  crept toward the point where a composite between twin primes reads as a prime);
- the forward step is **box-smoothed** (`BOXC`, one extra `sin` per term), which is
  what opens the gap between the two clouds from 0.49/0.27 to 0.71/0.21.

The price is reach per zero: 1502 instead of 2962 at the default cap.

**No seed (2026-11-21).** The 40-ordinate table is gone: frame 0 marks 2 as prime with
`nz = 0`, and frame 1 locates γ₁ from {2} (14.20). Three small-`nz` fixes came with it: the
detector fills every truncation even when their sizes coincide, the reach reads
`Z[max(nz/2,1)-1]`, and each frame licenses only the zeros that the primes it uses
(`<= Xknown`) support. The forward step and the detector run in the same frame, so with
`nFrontier(new_X)` the first zeros were located before any prime was committed. That was
harmless with 40 table zeros, but from one zero it misses 7 and echoes every multiple of 7.
`stage2-sim.c` gives the same result as the old 40-zero start: exactly the primes to 1502
(`NLIM = 8000`) and to 8585 (`NLIM = 60000`), 3–5 frames later.

## Credits
Number printing: `PrintValue` by **@P_Malin**, CC0
(https://www.shadertoy.com/view/4sBSWW) — credit kept in the Image tab.
