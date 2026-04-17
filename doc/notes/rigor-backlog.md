# Roadmap: Rigorous Error Bounds for Approximate Zero Location

**Goal**: Prove that the code's approximate zero counting function is within 1/2 of the true N(T), which suffices for correct zero location via bisection.

**Conditional on**: Riemann Hypothesis throughout.

---

## What the code computes

The zero counting function approximation (see `main.c:479`, `zero_count_exact`):

```
F(T) = N_0(T) + (1/pi) Im sum_{i=1}^{k} att(p_i, T) * log(1 - p_i^{-1/2 + iT})
```

where:
- `N_0(T) = (T/2pi) log(T/2pi e) + 7/8` is the smooth Riemann-von Mangoldt counting function
- `att(p,T) = 1 - exp(-sqrt(T/p))` is the attenuation factor
- The sum runs over the first k primes

**Note on sign and damping.** The code evaluates `Re(-i/pi · sum att log(1-p^{-1/2+iT}))`, which equals `(1/pi) Im sum att log(1-p^{-1/2+iT})` (a **plus** sign on the sum, not minus). The exponent used is `-1/2+iT`, so by complex conjugation `Im log(1-p^{-1/2+iT}) = -Im log(1-p^{-1/2-iT})` and the formal sum `-sum sin(mT log p)/(m p^{m/2})` agrees with `pi·S(T) = Im log zeta(1/2+iT)`. The damping is `1 - exp(-sqrt(T/p))` as implemented in `main.c:500-506`, matching the README and the Wolfram reference `doc/heuristic/test-function-search.wls:64-65` (where `alpha = sqrt(T)` and `u = log p` so `alpha·exp(-u/2) = sqrt(T/p)`). The soft cutoff therefore sits at `p ~ T`, where `att = 1 - 1/e ~ 0.632`; for `p >> T`, `att ~ sqrt(T/p)`. (An earlier draft claimed `1 - exp(-2T/sqrt(p))` and cutoff at `p ~ 4T^2`; both were wrong.)

The true counting function is:

```
N(T) = N_0(T) + S(T) + small corrections
```

where `S(T) = (1/pi) arg zeta(1/2 + iT)`.

**Target**: Show |F(T) - N(T)| < 1/2 for T away from zero ordinates, given sufficient k.

---

## Key findings from Phase 1 exploration

See `doc/test-function-search.wls` for the computations.

1. **Boundary result**: The heuristic attenuation corresponds to Gamma(-2ir) in the inverse FT, with a pole at Im(r) = 1/2 — exactly the admissibility boundary for the Weil explicit formula.

2. **Universality vs admissibility**: The attenuation's universality in T/p (scale-free behavior) is equivalent to non-admissibility. Any admissible test function breaks this universality and degrades with T.

3. **Consequence**: We cannot rigorize the attenuation via the Weil explicit formula. Instead, we bound the truncation error directly.

---

## Revised approach: direct truncation bound

The error decomposes as:

```
|F(T) - N(T)| ≤ |N_0(T) - N_0^true(T)|     ... (A) smooth part error
              + |att correction|              ... (B) attenuation error for kept primes
              + |tail from omitted primes|    ... (C) truncation error
```

Since att ≈ 1 for all used primes (p_k << T), term (B) is negligible. Term (A) is a known quantity (the Riemann-Siegel theta approximation error). Term (C) is the main challenge.

---

## Phase 1: The smooth part error (A)

### Task 1.1: Theta function approximation

The true smooth part is `(1/pi) theta(T) + 1` where theta is the Riemann-Siegel theta function. The code uses the Stirling approximation `N_0(T) = (T/2pi) log(T/2pi e) + 7/8`.

The error is:

```
|(1/pi) theta(T) + 1 - N_0(T)| = O(1/T)
```

For T > 100 this is less than 0.01. Well-understood, explicit constants available (see DLMF 25.11).

**Deliverable**: Explicit bound, confirm < 0.01 for T > 1000.

**Status**: Routine. The Stirling expansion of theta is textbook.

---

## Phase 2: The attenuation error (B)

### Task 2.1: Bound the attenuation correction for kept primes

For p ≤ p_k with p_k << T, the actual attenuation `att(p,T) = 1 - exp(-sqrt(T/p))` satisfies:

```
|1 - att(p,T)| = exp(-sqrt(T/p)) ≤ exp(-sqrt(T/p_k))
```

For T = 10^12, k = 1000, p_k = 7919: `sqrt(T/p_k) = sqrt(10^12 / 7919) ≈ 11236`, so `|1 - att| ≤ exp(-11236) ≈ 0`.

The total attenuation error is:

```
|sum_{p ≤ p_k} (1-att) * f(p)| ≤ exp(-sqrt(T/p_k)) * sum |f(p)|
```

Not doubly-exponential (as mistakenly claimed before), but still vastly beyond any working precision: for `T = 10^12`, `p_k = 7919`, the factor is `~10^-4880`. Safe to treat as zero for any realistic k and T; specifically, the bound is tighter than machine zero whenever `sqrt(T/p_k) > ~50`, i.e. `T > 2500 * p_k`.

**Deliverable**: One-line bound. Essentially zero whenever `T >> p_k`.

**Status**: Trivial.

---

## Phase 3: The truncation error (C) — MAIN CHALLENGE

This is where the real work lives. We need to bound:

```
|S(T) - F_k(T)|
```

where F_k is the code's prime sum (with att ≈ 1).

Since S(T) cannot be expressed as a convergent Euler product sum on the critical line, we cannot simply "bound the omitted tail." Instead, we need an independent way to relate S(T) to a truncated prime sum.

### Strategy overview

Three complementary approaches, in order of increasing difficulty and tightness:

**Strategy 3A**: Use known bounds on S(T) itself (coarse but simple)
**Strategy 3B**: Use the explicit formula at σ > 1/2 and continuity to the critical line
**Strategy 3C**: Use oscillatory cancellation in the truncated sum via partial summation

---

### Strategy 3A: Bounds on S(T) directly

Under RH (Littlewood):

```
|S(T)| ≤ C log(T) / log(log(T))
```

with explicit C. For T = 10^12: |S(T)| ≤ ~8 (with C ~ 1).

Also, Selberg's mean-square result:

```
integral_0^T |S(t)|^2 dt ~ (1/2pi^2) T log(log(T))
```

So S(T) is *typically* O(sqrt(log log T)) ≈ 1.7, though it can occasionally be as large as O(log T / log log T).

**What this gives**: If |S(T)| < 1/2 (which happens for "most" T since S is typically small), then the code just needs F_k ≈ 0, which is true for any reasonable k. But S(T) CAN exceed 1/2 near clusters of zeros, and that's exactly where we need the prime sum to track it.

**Verdict**: Too coarse for a general bound, but useful as context. The prime sum must actually approximate S(T), not just be small.

**References**:
- Littlewood: "On the zeros of the Riemann zeta-function" (1924)
- Selberg: "Contributions to the theory of the Riemann zeta-function" (1946)
- Goldston: "On the function S(T)" (1984)

---

### Strategy 3B: σ-shifting — evaluate off the critical line

**Idea**: At σ = 1/2 + η (slightly off the critical line), the Euler product converges. Evaluate there, then bound the difference to σ = 1/2.

**Step 1**: For σ > 1/2, define:

```
S_σ(T) = (1/pi) Im log zeta(σ + iT)
       = -(1/pi) sum_p sum_m sin(mT log p) / (m p^{mσ})
```

This converges conditionally for σ > 1/2 under RH, absolutely for σ > 1.

**Step 2**: Bound the difference |S(T) - S_σ(T)|:

```
|S(T) - S_σ(T)| = (1/pi) |Im [log zeta(1/2+iT) - log zeta(σ+iT)]|
                ≤ (1/pi) integral_{1/2}^{σ} |zeta'/zeta(s+iT)| ds
```

Under RH, near the critical line:

```
zeta'/zeta(σ+iT) = sum_{|T-γ|<1} 1/(σ+iT-ρ) + O(log T)
```

where the sum is over nearby zeros ρ = 1/2 + iγ. The dominant contribution comes from the nearest zero. If T is distance δ from the nearest zero ordinate:

```
|zeta'/zeta(σ+iT)| ~ 1/(σ - 1/2 + δ) + O(log T)
```

Integrating from 1/2 to 1/2+η:

```
|S(T) - S_σ(T)| ~ (1/pi) log(1 + η/δ) + O(η log T)
```

For η = 1/log T and δ ~ 2pi/log T (typical zero spacing): this is O(1), which is not small enough.

**Step 3**: Truncate S_σ to k primes and bound the tail.

At σ = 1/2 + η, the tail Σ_{p>P} p^{-σ} converges better. By partial summation with PNT under RH:

```
|sum_{p>P} p^{-σ-iT}| ~ P^{1/2-η} / (T log P)  (from the Li term)
                       + T P^{-η} log(P) / η      (from the RH error term)
```

**The problem**: The RH error term gives T × P^{-η} × stuff, which grows with T. For η = O(1/log T), this requires P to grow with T, defeating the purpose.

**Verdict**: The σ-shifting approach gives a bound, but the constants are bad. The truncation at σ > 1/2 requires P to grow faster than any power of log T, making the prime count requirement much worse than empirical. The approach is sound but the bounds are too loose.

**What would be needed**: Explicit, tight versions of the ζ'/ζ bound near the critical line. The best results in this direction are by:
- Carneiro, Chandee, Milinovich (2013): "Bounding S(T) and S_1(T) on the Riemann hypothesis"
- Goldston, Gonek, Lee, "A zero density result for the Riemann zeta function" (2010s)

**Task 3B.1**: Look up the explicit constants in Carneiro-Chandee-Milinovich for the σ-shifted bound and evaluate whether they're tight enough. Their work gives:

```
|S(T)| ≤ (1/4 + o(1)) log T / log log T
```

with explicit error terms. The method might extend to give |S(T) - S_σ(T)| bounds.

---

### Strategy 3C: Direct partial summation on the critical line

**Idea**: Rather than going through log ζ, bound the truncated prime sum directly using the distribution of primes (PNT under RH) and oscillatory cancellation from sin(T log p).

**The sum to bound**:

```
Tail(P,T) = sum_{p > P} sin(T log p) / (pi sqrt(p))
```

(focusing on the m=1 terms; m ≥ 2 handled separately since they converge absolutely for m ≥ 3, and conditionally for m = 2 via PNT).

**Step 1**: Partial summation (Abel summation). Let A(x) = Σ_{p≤x} 1. Then:

```
sum_{P < p ≤ X} sin(T log p) / sqrt(p)
  = A(X) sin(T log X)/sqrt(X) - A(P) sin(T log P)/sqrt(P)
  - integral_P^X A(x) d/dx [sin(T log x)/sqrt(x)] dx
```

The derivative:

```
d/dx [sin(T log x)/sqrt(x)] = [T cos(T log x)/x - sin(T log x)/(2x)] / sqrt(x)
                              = T cos(T log x) / (x sqrt(x)) + lower order
```

Under RH: A(x) = π(x) = Li(x) + O(√x log x), where Li(x) ~ x/log x.

**Step 2**: The main term (from Li(x)):

```
integral_P^∞ Li(x) * T cos(T log x) / (x sqrt(x)) dx
≈ integral_P^∞ T cos(T log x) / (sqrt(x) log x) dx
```

Substitute u = log x:

```
= integral_{log P}^∞ T cos(Tu) * exp(u/2) / u * exp(-u) du     [wait, needs care]
```

Actually: x = e^u, dx = e^u du, √x = e^{u/2}, so:

```
= integral T cos(Tu) / (e^{u/2} u) * e^u du ... hmm let me redo
```

Let me be more careful. ∫ (T cos(T log x))/(x^{3/2} log x) dx. Let u = log x:

```
= integral T cos(Tu) exp(-3u/2) / u * exp(u) du = integral T cos(Tu) exp(-u/2) / u du
```

Integration by parts with the oscillating factor cos(Tu):

```
integral T cos(Tu) exp(-u/2)/u du = [sin(Tu) exp(-u/2)/u]
                                  + integral sin(Tu) d/du[exp(-u/2)/u] du
```

The boundary term at u = log P: `sin(T log P) exp(-log(P)/2) / log P = sin(T log P) / (√P log P)`.
The remaining integral has integrand ~ sin(Tu) exp(-u/2)/u^2, which is smaller by a factor 1/u.

**So the main contribution is**:

```
|boundary term| ~ 1/(√P log P)
```

This is the leading behavior of the truncation error from the Li part.

**Step 3**: The error from RH. The contribution from E(x) = π(x) - Li(x):

```
integral_P^∞ E(x) * T cos(T log x) / x^{3/2} dx
```

Under RH: |E(x)| ≤ C √x log x, so:

```
|error| ≤ CT integral_P^∞ √x log x * 1/x^{3/2} dx = CT integral_P^∞ log(x)/x dx
```

This DIVERGES. The RH error term is too large for a pointwise bound via this naive approach.

**Step 4**: The fix — exploit oscillation of E(x) · cos(T log x).

The key insight: E(x) is not just bounded, it OSCILLATES. Under RH, E(x) itself has a representation via zeros:

```
E(x) = π(x) - Li(x) = -(1/log x) sum_ρ Li(x^ρ) + lower order
```

The interaction of E(x)'s oscillation with cos(T log x) produces cancellation. Bounding this requires understanding the correlation between E(x) and the oscillating kernel.

**This is where the problem becomes research-level.** The required bound is:

```
|integral_P^∞ E(x) cos(T log x) / x^{3/2} dx| ≤ ???
```

Techniques available:
- **Large sieve inequalities**: Bound mean-square of sums over primes with oscillating weights
- **Montgomery-Vaughan estimates**: Explicit bounds on exponential sums over primes
- **Gallagher's lemma**: Relate pointwise bounds to mean-square bounds
- **Vaughan's identity**: Decompose the sum into bilinear forms for better cancellation

**Task 3C.1**: Formulate the truncation error as an exponential sum over primes and apply Montgomery-Vaughan type estimates. Specifically, bound:

```
|sum_{P < p ≤ X} p^{iT-1/2}| ≤ ???
```

Under RH, using the explicit formula for sums over primes.

**Task 3C.2**: Handle the m = 2 prime-power terms separately. These contribute:

```
sum_p sin(2T log p) / (2p)
```

By PNT (Mertens' theorem with oscillation), this converges conditionally. The partial sum up to P has error O(1/log P) from partial summation.

**Task 3C.3**: The m ≥ 3 terms converge absolutely:

```
sum_p sum_{m≥3} 1/(m p^{m/2}) ≤ sum_p p^{-3/2} / (1 - p^{-1/2}) < ∞
```

The tail for p > P is O(1/P). Negligible.

---

### Strategy 3D: Goldston's formula (literature shortcut)

There is a classical result connecting S(T) to truncated prime sums. The strongest version appears to be:

**Theorem (Goldston 1984, conditional on RH)**: For X ≥ 2 and T > 0:

```
S(T) = -(1/pi) sum_{n ≤ X} Lambda(n) sin(T log n) / (n^{1/2} log n) + R(T,X)
```

where `R(T,X)` satisfies explicit bounds depending on X and T.

The error R(T,X) has been studied by several authors:
- Goldston (1984): R = O(log T / log X) — coarse but general
- Fujii (1999): improved bounds in various ranges
- Carneiro-Chirre-Milinovich (2022): optimal bounds using extremal functions

**Task 3D.1**: Locate the tightest explicit version of this formula in the literature. The key question is: for the code's truncation point X = p_k, is the known bound on R tight enough to give |R| < 1/2?

**Task 3D.2**: If the best known R = O(log T / log X) with constant ~ 1, then for |R| < 1/2 we need log X > 2 log T, i.e., X > T^2. This would require k ~ T^2/log(T^2) primes — far too many.

But the O-constant matters. If the actual result is R ≤ (1/(2π) + ε) log T / log X (as suggested by some results), then X > T^{1/π} suffices, still too many.

**Task 3D.3**: Check whether Goldston's bound can be improved for our specific setting (where we also have the attenuation providing a smooth cutoff rather than a sharp truncation at X).

**Key references to obtain and study**:
- Goldston, D.A. "On the function S(T) in the theory of the Riemann zeta-function" (1984), J. Number Theory
- Fujii, A. "Explicit formulas and oscillations" (1999), in Number Theory volumes
- Carneiro, E., Chirre, A., Milinovich, M. "Bandlimited approximations and estimates for the Riemann zeta-function" (2022)
- Selberg, A. "Contributions to the theory of the Riemann zeta-function" (1946)

---

## Phase 4: The m = 1 exponential sum (computational exploration)

Before committing to a full proof, explore the truncation error numerically using wolframscript.

### Task 4.1: Compute S(T) at sample points

For moderate T (say T near the first few thousand zeros), compute S(T) exactly using `acb_dirichlet_hardy_z` or similar, and compare with the code's truncated sum.

```
error(T, k) = |S(T) - F_k(T)|
```

Map out error(T, k) as a function of both T and k.

### Task 4.2: Empirical scaling of the tail

Compute the partial sums:

```
tail(P, T) = sum_{p > P, p ≤ Q} sin(T log p) / (pi sqrt(p))
```

for various P, Q, T. Measure the actual cancellation. Compare with:
- 1/√P · 1/log P (the partial-summation prediction)
- log T / log P (the Goldston bound)

### Task 4.3: Explore the exponential sum

The key object is sum_{P < p ≤ X} p^{-1/2+iT}. Plot its magnitude as a function of P for various T. Under GRH, we expect it to be O(√X / log X) or better (square-root cancellation). Verify numerically.

---

## Phase 5: Assembly

### Task 5.1: Combine the bounds

```
|F(T) - N(T)| ≤ |smooth error|      ... Phase 1: O(1/T), negligible
              + |att error|          ... Phase 2: exp(-2T/√P), negligible
              + |m≥3 tail|           ... Task 3C.3: O(1/P), negligible
              + |m=2 tail|           ... Task 3C.2: O(1/log P), small
              + |m=1 tail|           ... THIS IS THE BOTTLENECK
```

The m=1 tail bound determines k(T). Three scenarios:

| Bound achieved | Required k(T) | Status |
|---|---|---|
| O(log T / log P) with const ~ 1 | k ~ π(T^2) ≈ ridiculous | Goldston raw |
| O(1/√P · 1/log P) from partial sum | k ~ 1/(4ε²) ≈ 2500 for ε=0.01 | If oscillatory cancellation works |
| O(√P/(T log P)) from exp sum + Li | k ~ (few hundred) | If GRH exp sum estimate applies |

The empirical scaling k ~ (log T)^2 suggests the actual error is MUCH smaller than any of these, indicating either:
- Significant additional cancellation not captured by these bounds
- The attenuation provides more help than the "att ≈ 1" analysis suggests
- There is structure in the specific T values (near zeros) that helps

### Task 5.2: Determine whether existing bounds suffice

If the literature bounds give |error| < 1/2 with k ~ (log T)^2, we're done. If not, identify the gap and determine whether it's closeable with known techniques or requires new ideas.

### Task 5.3: Write up

If the bound works: write as a self-contained theorem with explicit constants.
If not: document the gap precisely, identifying what new estimate would suffice.

---

## Phase 6: Validate numerically

### Task 6.1: Test against Odlyzko

Verify the proved bound holds (with room to spare) at:
- T ~ 10^12 (Odlyzko's zeros3 table)
- T ~ 10^20 (Odlyzko's large computation)
- Consecutive zeros near T = 10^6 (many available)

### Task 6.2: Stress test

Find the T values where the code's error is LARGEST (near clusters of zeros where S(T) is large). Verify the bound holds even there.

---

## Dependency graph

```
Phase 1 (smooth)  ─────────────────────────────────┐
Phase 2 (att)     ─────────────────────────────────┤
                                                    ├──> Phase 5 (assembly)
Phase 3 strategies:                                 │
  3A (S bounds)     ─ context ─────────────────────┤
  3B (σ-shift)      ─ if constants are tight ──────┤
  3C (partial sum)  ─ main technical work ─────────┤
  3D (literature)   ─ possible shortcut ───────────┘
                                                    │
Phase 4 (numerics)  ─ guides which strategy ────────┘
                                                    │
Phase 5 ──────────> Phase 6 (validation)
```

## Critical path

**Phase 4 (numerics) should come FIRST**. The numerical exploration will tell us:
1. How large is the actual error for various k and T?
2. How does it compare with the Goldston bound O(log T / log P)?
3. Is there a tighter empirical relationship we can try to prove?

This determines whether Strategy 3B, 3C, or 3D is the right path.

## Open question

The empirical success with k ~ (log T)^2 primes suggests an error of O(1/log T) or better. The best rigorous bound we've identified is O(log T / log P), which with P = p_k ~ k log k gives error ~ log T / log k. For this to be < 1/2 with k ~ (log T)^2, we'd need:

```
log T / log((log T)^2) ≈ log T / (2 log log T) ≈ same Littlewood bound
```

For T = 10^12: 27.6 / (2 × 3.3) ≈ 4.2. This is > 1/2, so the Goldston bound is NOT tight enough to explain the empirical performance.

**Something else must be happening.** Either:
1. The O-constant in Goldston is much smaller than 1 (need to check the paper)
2. The attenuation provides additional cancellation beyond "att ≈ 1"
3. The bisection process is more forgiving than |error| < 1/2 (the counting function is monotone, so even with error > 1/2, bisection may converge to the right zero if the error is consistent)
4. The code works "in practice" but not provably — the gap between theory and practice is real

Resolving this is the core intellectual challenge of the project.
