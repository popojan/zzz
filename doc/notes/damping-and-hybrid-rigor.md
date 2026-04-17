# Damping, Euler‑Product Projection, and a Path to Rigor via the Hybrid Product

**Status:** design note. Companion to `rigor-backlog.md`. Assumes RH throughout.

---

## 1. TL;DR

1. The active damping in `main.c:497–508` is

   ```
   att(p, T) = 1 − exp(−√(T/p))
   ```

   — **not** `1 − exp(−2T/√p)` as claimed in `rigor-backlog.md:22`. The shape is very different (see §2.3). This bug should be fixed; the Phase‑2 arithmetic that depends on it survives but the reasoning does not.

2. Using the complex conjugation identity `Im log(1 − p^{−1/2+iT}) = −Σ_m sin(mT log p)/(m p^{m/2})`, the code's contribution to `F(T)` is (formally) the Dirichlet sum for `π·S(T)`. Numerical check against `./zzz -g` confirms the sign and magnitude.

3. **Projection.** Define, for each `T` and cutoff `k`:

   ```
   Z_k(s, T)  :=  ∏_{p ≤ p_k} (1 − p^{−s})^{−att(p, T)}.
   ```

   Then `zzz`'s bisection finds `T_n` satisfying
   `N_0(T_n) + (1/π) arg Z_k(1/2 − iT_n, T_n) = n − 1/2`.
   **These are *zeros of `Z_k`*, not of `ζ`**, and `Z_k` is an explicit finite object. They approximate `ζ`‑zeros because `Z_k` tracks the partial Euler product of `ζ` truncated near `p ∼ T`.

4. The damping lives exactly on the Weil admissibility boundary (`Γ(−2ir)` pole at `Im r = 1/2`). It therefore cannot be rigorized through Weil's explicit formula. The cleanest rigorous replacement is the **Gonek–Hughes–Young (GHY) hybrid Euler–Hadamard product** — which is an *identity*, not an approximation.

5. **Zeta‑evaluation‑freeness is preserved** under GHY *if* the local Hadamard factor is supplied with already‑known zeros. Two viable modes: bootstrap from below (using previously found zeros) or seeded from a tabulated cache (Odlyzko). No `acb_zeta` call is needed at any step.

---

## 2. What the code computes (corrected)

### 2.1 Formulas

- Smooth part: `N_0(T) = (T/2π) log(T/2πe) + 7/8` (`nt()` in `main.c:345`).
- Exponent used: `−1/2 + iT`, so `wave_complex_opt` returns `log(1 − p^{−1/2+iT})`.
- Damping: `att(p, T) = 1 − exp(−√(T/p))` (`main.c:500–506`).
- Output:

  ```
  F(T) = N_0(T) + (1/π) Σ_{p ≤ p_k}  att(p, T) · Im log(1 − p^{−1/2+iT}).
  ```

### 2.2 Numerical sanity check

Running `./zzz -g 10 40 30 -k 100`:

| T     | F(T)    | expected N(T) |
|-------|---------|---------------|
| 15.33 | 0.9851  | 1             |
| 22.00 | 2.0171  | 2             |
| 26.00 | 2.9902  | 3             |
| 34.00 | 5.0243  | 5             |

Good agreement — signs and magnitudes match `N(T) = N_0(T) + S(T)` as intended.

### 2.3 Why the `rigor-backlog.md` formula is different

| Form                          | Value at `p = T`   | Decay for `p ≫ T` | FT strip |
|-------------------------------|--------------------|-------------------|----------|
| Code: `1 − exp(−√(T/p))`      | `1 − 1/e ≈ 0.632`  | `∼ √(T/p)`        | `|Im r| < 1/2` (boundary) |
| Doc: `1 − exp(−2T/√p)`        | `≈ 1` (at p≈T; stays near 1 until `p ∼ 4T²`) | ? | different |

The actual "soft cutoff" is at `p ∼ T`, not `p ∼ T²`. All of Phase 1's conclusions (pole at `Im r = 1/2`, boundary admissibility, universality in `T/p`) are properties of `1 − exp(−√(T/p))` and should remain. The Phase‑2 "doubly exponential" bound survives quantitatively but for a different reason; re‑derive with the corrected form.

---

## 3. Why the damping works (three reasons)

1. **Inert on evaluated terms.** For `p ≤ p_k ≪ T`, `exp(−√(T/p_k))` is doubly subnormal; to working precision `att ≡ 1`. The damping's effect on what you *compute* is null.

2. **Marginal tail taming.** For `p ≫ T`, `att ∼ √(T/p)`, giving tail terms `∼ √T · sin(mT log p)/p`. Summed over primes this is *conditionally* convergent under PNT/RH; without damping, `p^{−1/2}`‑weighted oscillation has no absolute‑sense limit either.

3. **Weil boundary.** The test‑function search (`test-function-search.wls:22–50`) proves the inverse FT of `1 − exp(−α e^{−u/2})` is `−2 α^{2ir} Γ(−2ir)`, whose first pole is at `Im r = 1/2` — *exactly* the admissibility threshold. The heuristic is the sharpest zero‑concentrator still sitting on the Weil boundary. Any genuinely admissible replacement (Gaussian, logistic, Vaaler) is strictly softer.

The *reason the method works empirically* is (1) plus the oscillatory cancellation inside the truncated bare sum. The damping's role is conceptual/theoretical: it makes the `k → ∞` limit well posed, at the cost of non‑admissibility.

---

## 4. The Euler‑product projection — what `zzz` actually finds

### 4.1 Finite proxy

`att(p, T) · log(1 − p^{−s}) = log((1 − p^{−s})^{att(p, T)})`, so

```
Z_k(s, T) := ∏_{p ≤ p_k} (1 − p^{−s})^{−att(p, T)}
```

is an explicit finite product of branch‑cut Euler factors with fractional, `T`‑dependent exponents. By construction,

```
F(T) = N_0(T) + (1/π) arg Z_k(1/2 − iT, T).
```

`zzz` locates the `T_n` for which this quantity passes `n − 1/2`, i.e. **the zeros of the explicit proxy `Z_k`**, not of `ζ` itself.

### 4.2 `k → ∞` limit

```
Z_∞(s, T) = ∏_p (1 − p^{−s})^{−att(p,T)}
          = ζ(s) · ∏_p (1 − p^{−s})^{exp(−√(T/p))}.
```

The correction factor's Euler exponent `exp(−√(T/p))`:
- is `≈ 0` for `p ≪ T` — those Euler factors *cancel* those of `ζ`, leaving only the `p ≲ T` tail from `ζ`.
- is `≈ 1` for `p ≫ T` — those factors contribute `(1 − p^{−s})^1`, i.e. *remove* the corresponding Euler factors of `ζ`.

Result: `Z_∞(s, T) ≈ ∏_{p ≲ T} (1 − p^{−s})^{−1}`, the partial Euler product of `ζ` softly cut at `p ∼ T`. **This is exactly the "partial Euler" half of the GHY hybrid product** — which is why rigorization naturally leads there.

---

## 5. Path to rigor: the GHY hybrid Euler–Hadamard product

Reference: Gonek, Hughes, Young, *A hybrid Euler–Hadamard product for the Riemann zeta function*, Duke Math. J. 136 (2007), 507–549.

### 5.1 The identity

For a smooth compactly supported `u(x) ≥ 0` with `∫ u = 1` supported on `[0,1]`, cutoff `X ≥ 2`, and `s = σ + it` with `t ≥ 2`:

```
log ζ(s) = log P_X(s) + log Z_X(s) + E(s, X)
```

with

```
log P_X(s) = Σ_{n ≥ 2}  Λ(n) / (n^s log n) · u(log n / log X)
log Z_X(s) = − Σ_ρ       U((s − ρ) log X)
```

and an explicit kernel

```
U(z) = ∫_0^∞ u(x) Γ(0, z·x) dx         (incomplete Γ; GHY eq. (2.4))
```

`E(s, X)` is bounded uniformly on vertical strips away from zeros/poles; with the `u` choice of GHY one obtains error `O(1 / log X)` at fixed `t`, unconditionally, and better under RH.

**Key point.** This is an *identity*, not a heuristic. Both `P_X` and `Z_X` are computable when `X` and the relevant `ρ` are fixed; `E` is a finite remainder with known bounds.

### 5.2 Why this is the natural rigor target for `zzz`

- `P_X(s)` is a prime‑power Dirichlet sum **with a smooth cutoff at `n ∼ X`** — structurally identical to what `zzz` computes, but with explicit `Λ(n)` weights (prime powers included) and an *admissible* mollifier `u`.
- `Z_X(s)` is supported only on a window `|γ − t| ≲ 1/log X`, because `|U(z)|` decays rapidly once `|z| ≳ 1`. Zeros far from height `t` contribute exponentially little.
- The `zzz` heuristic corresponds to setting `Z_X ≡ 0` and choosing a non‑admissible `u`. The identity above tells us precisely what is being neglected: the local zero contribution in a window of size `O(1/log T)` around `t`.

### 5.3 Cutoff `X` scaling

`zzz` today uses `k ∼ (log T)^2` primes, i.e. `X = p_k ∼ (log T)^2 log log T`. With that `X`:
- `log X ∼ 2 log log T`, so the GHY window is `|γ − t| ≲ 1 / (2 log log T)` — typically 0–few zeros.
- `E(s, X)` at `σ = 1/2` has a known bound `O(log T / log X) = O(log T / log log T)`. Coarse — same ballpark as Littlewood's `S(T)` bound.

To force `|E| < 1/2` may require `X ≳ T^ε`. But `P_X`'s cost is `O(π(X))`, so making `X = √T` is tolerable for moderate `T` and far cheaper than Riemann–Siegel's `O(T^{1/2})` *per zero*.

### 5.4 Local‑zero contribution `Z_X`

- For each nearby zero `γ_j` with `|γ_j − t| < c/log X`: include `−U((s − ρ_j) log X)`.
- Zeros outside that window: bound contribution by `O(exp(−c' |γ_j − t| log X))`, sum analytically.
- The **target zero itself** is one of the `ρ_j` — its contribution near its own location is a known log‑singularity `−U(0·log X) = +∞·δ`‑like. Practically this is where the bisection root lives: `F_corrected(T_n) − (n − 1/2) = 0` becomes a transcendental equation solvable by bisection exactly as now, but with a tighter, rigorous error band.

---

## 6. Zeta‑evaluation‑free?

**Short answer: yes, in all operational modes.**

| Ingredient                    | Requires `ζ(·)` eval? | Requires zero list? |
|-------------------------------|-----------------------|---------------------|
| `N_0(T)` smooth               | no                    | no                  |
| `P_X(s)` partial Euler sum    | no                    | no                  |
| `Z_X(s)` local Hadamard       | no                    | **yes (window only)** |
| `E(s, X)` error bound         | no                    | no                  |

The local Hadamard factor is the only term that consumes zero data. Three acceptable ways to supply it without invoking `ζ`:

1. **Bootstrap.** Zeros found in increasing order: the window for zero `n+1` needs only zeros at heights `≥ T_n − O(1/log X)`, which have already been found (and can themselves be re‑refined by plugging the new zero back in).
2. **Seeded from table.** For random‑access at arbitrary height, start from the nearest tabulated zero (Odlyzko, LMFDB) and bootstrap outward.
3. **Self‑consistent iteration.** Use heuristic `zzz` as initial guess for the window zeros; evaluate `P_X + Z_X`; refine. Fixed point is rigorous.

Modes 1 and 2 preserve the current design: `main` never calls `acb_zeta`; the `-e` flag remains a post‑hoc *verification* only.

---

## 7. Concrete task list

Ordered by dependency. Tasks are self‑contained; pick one at a time.

### Task R1 — Fix the documentation bug
- Correct `rigor-backlog.md:22` to `att = 1 − exp(−√(T/p))`.
- Recompute the Phase‑2 upper bound `|1 − att| = exp(−√(T/p_k))` and the numerical example.
- Update the rhetorical claim ("stays at ~1 for all `p` up to `~4T²`") — for the corrected form, the soft cutoff is at `p ∼ T`, not `p ∼ T²`.

### Task R2 — Export the finite proxy `Z_k(s, T)` as a first‑class object
**Goal:** make "what zzz actually computes" visible.
- New small binary `zproxy` (or a `--proxy` flag on `zzz`) that evaluates `Z_k(σ + iT, T)` on a grid of `σ` and `T`, using the same code path as `zero_count_exact` but exposing magnitude and argument.
- Plot `|Z_k(1/2 + iT, T)|` over `T ∈ [10, 100]` and mark `ζ`‑zeros from `ZetaZero[]` for comparison. Ship plot in `doc/ghy/`.
- Output: empirical answer to the user's question, "what zeros are these?"

### Task R3 — Read and extract explicit constants from GHY 2007
**Goal:** pin down the kernel `U` and the `E`‑bound constants.
- Target Theorem 1 and its proof (the explicit form of `E`).
- Choose a specific mollifier `u` — GHY's canonical `u(x) = (some C^∞_c bump)` or Vaaler's. Compute `U(z)` for that `u` to high precision; tabulate or write a series.
- Produce a tiny Wolfram notebook `doc/ghy/ghy-kernel.wls` that:
  - plots `U(z)` vs `|z|`,
  - verifies `U(0)`, `U(iπ)` etc. against the paper,
  - provides `fmt`ed numeric coefficients for C implementation.

### Task R4 — Implement `log P_X(s)` in C
**Goal:** drop‑in replacement for the current inner sum, with a provable mollifier.
- Arguments: `s`, `X`, choice of `u`.
- Sum: `Σ_{n ≤ X} Λ(n)/(n^s log n) · u(log n / log X)` using FLINT's factorization for prime powers (`n_factor_power` / prime‑sieve up to `X`).
- Compare to current `zero_count_exact`'s prime sum on the Odlyzko test set; measure sign of the systematic shift.

### Task R5 — Implement `log Z_X(s)` in C
**Goal:** local zero‑correction module.
- Inputs: list of `(γ_j)` in a window `|γ_j − Im(s)| ≤ C/log X`, cutoff `X`, mollifier choice.
- Compute `−Σ_j U((s − ρ_j) log X)` via the series/integral from Task R3.
- Tail bound: deterministic `O(exp(−c|γ_j − t| log X))` summation for zeros just outside the window; prove negligible.
- Unit test: input a single dummy `ρ_j = 1/2 + iT`, evaluate `Z_X(1/2 + it)` in a neighborhood of `T`, confirm it looks like the expected log singularity.

### Task R6 — Bootstrap harness
**Goal:** maintain a rolling zero cache.
- Wrap `zzz` in a loop that, for each new zero, feeds back into the `Z_X` computation for the next.
- Keep an in‑memory ring buffer of `O(log X)` recent zeros.
- Add a `--seed <file>` flag to warm‑start from a pre‑computed table (Odlyzko format).

### Task R7 — Empirical error measurement
**Goal:** quantify what GHY actually buys.
- Three variants: (A) current heuristic, (B) `P_X` alone (no `Z_X`), (C) `P_X + Z_X` bootstrap.
- Run on a shared Odlyzko reference list at `T ∈ {10^6, 10^12, 10^20}` with matched `k`/`X`.
- Plot error distribution for each variant. Expect: A ≈ C in typical cases; C beats A near zero clusters where `|S(T)|` is large.

### Task R8 — Write up the theorem
**Goal:** crystallize a publishable claim.
- Statement form: *Conditional on RH, the `zzz_rigor` output at height `T` has `|F(T) − N(T)| ≤ c·g(T,X)` for explicit `c` and `g`.*
- Appendix: software version of the error budget matching `phase4`'s empirical curves.

### Task R9 — (Stretch) Symmetric Selberg‑type fallback
**Goal:** a fully zero‑free but looser rigor path.
- Use Selberg's 1946 identity

  ```
  S(T) = −(1/π) Σ_{n ≤ X} Λ(n) sin(T log n) / (√n log n) + R(T, X)
  ```

  with explicit `R = O(log T / log X)` (Goldston 1984 refined by Carneiro–Chirre–Milinovich).
- Implement in parallel; compare bounds to GHY. Useful for cases where no local zeros are available (first ascent from `T = 0` with an empty cache).

---

## 8. Open questions (not blocking)

- What is the *tightest* admissible `u` for minimizing the GHY error at given `X`? Extremal‑function techniques (Vaaler, Carneiro–Vaaler) may give the answer.
- Can the `Z_X` kernel be expressed in closed form using `FLINT`'s `acb_hypgeom` (incomplete Γ)? If so, the implementation of Task R5 reduces to one library call per zero.
- Is there a variant that uses only `log p` (no prime powers) with a weighted sum, matching `main.c` line‑for‑line but with provable error? Selberg's original formulation is over prime powers; a prime‑only variant may exist with slightly worse constants.
- Does the observed `k ∼ (log T)^2` empirical regime correspond to a specific optimization inside GHY (e.g. balancing `P_X`‑tail vs `Z_X`‑window size)?

---

## 9. Suggested order of attack

```
R1  ──►  R2  ──►  R3  ──►  R4  ──►  R5  ──►  R6  ──►  R7  ──►  R8
                                            │
                                            └─── R9 (parallel, optional)
```

R1–R2 are one‑afternoon jobs and produce tangible artifacts immediately. R3 is the literature gate; R4–R5 are the substantial engineering. R6–R7 are wiring + measurement. R8 closes it out.

---

## 10. References

- Gonek, Hughes, Young, "A hybrid Euler–Hadamard product for the Riemann zeta function", *Duke Math. J.* **136** (2007), 507–549.
- Gonek, "Finite Euler products and the Riemann Hypothesis", *Trans. AMS* **364** (2012), 2157–2191.
- Selberg, "Contributions to the theory of the Riemann zeta‑function", *Arch. Math. Naturvid.* **48** (1946), 89–155.
- Goldston, "On the function `S(T)` in the theory of the Riemann zeta‑function", *J. Number Theory* **27** (1987), 149–177.
- Carneiro, Chandee, Milinovich, "Bounding `S(T)` and `S_1(T)` on the Riemann hypothesis", *Math. Ann.* **356** (2013), 939–968.
- Carneiro, Chirre, Milinovich, "Bandlimited approximations and estimates for the Riemann zeta‑function", *Publ. Mat.* **63** (2022).
- This repo: `main.c:479–530` (`zero_count_exact`), `doc/notes/rigor-backlog.md`, `doc/heuristic/test-function-search.wls`.
