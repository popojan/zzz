# zzz live loop on Shadertoy (Stage 2)

Self-paving primes-from-zeros bootstrap running entirely on the GPU. From a
**40-zero seed** (the only baked-in data) it runs the explicit formula as a
closed loop — no sieve, no primality test, no ζ evaluation. Transliterated from
the `zzz` project's `--loop` (`loop.c`), verified flawless to 615 primes in fp32.

**Project:** https://github.com/popojan/zzz — `zzz` ("Zeta Zeros Zeal") computes
Riemann-zeta zeros from primes (and the reverse) via a ζ-evaluation-free
explicit-formula counter. The maths behind this demo lives in
`doc/notes/zeros-primes-bootstrap.md` and `doc/notes/exact-primes-inexact-zeros.md`.

## Setup (Common + two buffers + Image)
1. **Common** tab   ← `stage2-common.glsl`  (shared code, *not* a buffer; no `mainImage`)
2. **Buffer A** tab ← `stage2-bufferA.glsl` ; iChannel0 = **Buffer A** (self-feedback = state)
3. **Buffer B** tab ← `stage2-bufferB.glsl` ; iChannel0 = **Buffer A** (per-column trace cache)
4. **Image** tab    ← `stage2-image.glsl`   ; iChannel0 = **Buffer A**, iChannel1 = **Buffer B**

The single most common failure is forgetting that **Buffer A's iChannel0 must
point at Buffer A itself** — that self-feedback *is* the persistent state. Symptom
if it's missing: it runs but only ever shows "2". After wiring, rewind (⏮) to
reseed.

Buffer A holds the loop state in one RGBA32F texture (zeros in `.r`, prime-indicator
in `.g`, control in texels 0–1). Buffer B caches the detector trace `r(x)` once per
screen column (row 0) so the Image pass samples it with 3 `texelFetch` instead of a
384-term `cos` sum per pixel — ~100× cheaper drawing, the right move on weak GPUs.

## What you should see
From 40 seed ordinates the number line grows on its own; primes light up green
as the loop discovers them, the current one flares, its value prints.

## Knobs
- Buffer A: `BATCH` zeros/frame, `BAND` frontier/frame (bootstrap speed); `NLIM`
  zero cap (8000 → primes to ~3000; raise carefully — deep zeros lengthen the `FB`
  scan and can trip the browser's shader time limit).
- Buffer B: `NVIS` zeros in the background trace — now cheap (computed once per
  column), so you can push it high for razor spikes, or low for gentler ones.
- Common: `HOLD` seconds per prime.

## Status
Stage 1 (`primes-from-zeros.glsl`) is the baked-snapshot renderer (verified).
Stage 2 is the live loop; the kernels are a faithful port of the C verified
flawless to 615 primes, but it was authored without a GPU on hand, so if you hit
a compile error or a stuck readout, the symptom localises it fast (see the
self-feedback note above) and the texel-encoding can be re-checked in C.

## Credits
Number printing: `PrintValue` by **@P_Malin**, CC0
(https://www.shadertoy.com/view/4sBSWW) — credit kept in the Image tab.
