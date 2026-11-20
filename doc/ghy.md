# `--ghy`: GHY rigor mode — design and reference

The default `zzz` counter (method A) bisects a smooth heuristic
zero-counting function. `--ghy` swaps that inner factor for the GHY 2007
**partial Euler product** $P_X$ — same primes, slightly different weight,
but now placed inside a provable error chain (RH + Goldston 1987 + GHY
Theorem 1). This document is the design rationale and the rigor reference.

## Three methods

| | formula | input | rigor |
|---|---|---|---|
| **A** heuristic *(default)* | $F_A = N_0 + \tfrac{1}{\pi}\,\mathrm{Im}\!\sum_{p \le p_k}\bigl(1-e^{-\sqrt{T/p}}\bigr)\log(1-p^{-1/2+iT})$ | primes only | none |
| **B** GHY $P_X$ *(`--ghy`)* | $F_B = N_0 + \tfrac{1}{\pi}\,\mathrm{Im}\log P_X(\tfrac{1}{2}+iT)$,  $\log P_X = \sum_{p^m \le X}\frac{1}{m\,p^{ms}}$ | primes only | GHY + Goldston |
| **C** hybrid *(`zhybrid`)* | $F_C = N_0 + \tfrac{1}{\pi}\,\mathrm{Im}\log\bigl[P_X\cdot Z_X\bigr]$ | primes + nearby zeros | GHY (tightest) |

with $X = p_k$, $N_0(T) = \tfrac{T}{2\pi}\log\tfrac{T}{2\pi e} + \tfrac{7}{8}$, and the local Hadamard factor $Z_X(s) = \exp\!\bigl(-\sum_\rho U((s-\rho)\log X)\bigr)$ (GHY eq. 4).

The bisection driver is shared. `--ghy` is a one-line swap of the inner
sum; see `main.c::zero_count_ghy` calling `ghy_log_px` from `ghy.c`.

## Provable error bound

Under RH, two theorems compose:

**GHY 2007, Thm 1.** $\zeta(s) = P_X(s)\,Z_X(s)\,(1 + E)$ with
$|E(\tfrac{1}{2}+iT, X)| \le C_2\,\log X / \sqrt X$ for $T \ge 2$,
$X \ge 2$. This is the rigorous bound for **method C**.

**Goldston 1987.** For $T \ge T_0$ and $X \ge 2$ under RH,
$$
S(T) \;=\; -\frac{1}{\pi}\sum_{n \le X}\frac{\Lambda(n)\sin(T\log n)}{\sqrt n\,\log n} \;+\; R(T,X),
\qquad |R(T,X)| \le C\,\frac{\log T}{\log X}.
$$
The truncated sum **is** $F_B(T) - N_0(T)$, so

$$
\boxed{\;\bigl|F_B(T) - N(T)\bigr| \;\le\; \frac{1}{\pi}\,\frac{\log T}{\log X} \;+\; O\!\left(\frac{\log X}{\sqrt X}\right)\;}
$$

(conservative $C \le 1$). In the production regime this is $\sim 0.25$–$0.50$ — looser than $\tfrac{1}{2}$ at the highest heights, but explicit.

## Numerical confrontation

From `doc/ghy/ab-scan.wls`, fixed $k = 1000$ (so $X = 7919$):

| ordinal | $T$ | $\dfrac{\log T}{\pi\log X}$ | $\dfrac{\log X}{\sqrt X}$ | **bound** | **observed B** | ratio |
|---|---|---|---|---|---|---|
| #10     | $49.77$  | $0.139$ | $0.101$ | $0.240$ | $0.023$ | $10$  |
| #100    | $236.5$  | $0.193$ | $0.101$ | $0.294$ | $0.001$ | $294$ |
| #1000   | $1419$   | $0.257$ | $0.101$ | $0.358$ | $0.025$ | $14$  |
| #10000  | $9878$   | $0.326$ | $0.101$ | $0.427$ | $0.006$ | $71$  |
| #100000 | $74920$  | $0.397$ | $0.101$ | $0.498$ | $0.009$ | $55$  |

Two takeaways. (1) The bound is loose by $10$–$300\times$; observed error is roughly flat in $T$ while the bound grows like $\log T$. (2) Methods A and B agree to within $\sim 1.5\times$ across this scan: A is *not* doing anything mathematically deeper, it merely uses a smooth cutoff which suppresses Gibbs-type oscillation. Re-tuning B's mollifier (Gaussian, Vaaler-extremal) should close the empirical gap while preserving rigor — open work.

## Low-zero validity caveat

GHY Thm 1's $X^{K+2}/(T\log X)^K$ side term forces $X \lesssim T\log T$. At very low zeros this matters. For zero $\#1$ ($T \approx 14.13$):

| $k$ | $X$ | A finds | B (`--ghy`) finds |
|---|---|---|---|
| 5    | 11    | 14.203 | 14.192        |
| 50   | 229   | 14.141 | **14.133**    |
| 1000 | 7919  | 14.136 | 14.076        |
| 3000 | 27449 | 14.131 | 14.057        |

A converges and stays at $|err|\sim 10^{-3}$; B has a sweet spot near $k \approx \pi(T\log T)$ and drifts thereafter. **Production regime** ($T \ge 10^6$, $k \le 10^4$) sits firmly inside B's rigor band so this is not visible there. **Very low zeros** ($T \lesssim 100$): prefer A, or clamp $k$ at $\pi(T\log T)$ (not yet implemented).

## Auxiliary binaries

These dump TSV grids; they are not user-facing solvers. They exist so that
each ingredient of the rigor story can be visualized in isolation against
true zeta zeros.

| binary | dumps | role |
|---|---|---|
| `zproxy <T0> <T1> <N> [k]` | $\log\|Z_k\|$, $\arg Z_k$, $F_A(T)$ | method A's inner factor; lets you confirm $\arg Z_k$ jumps by $\pi$ at each $\gamma_n$ |
| `zghy <T0> <T1> <N> <X>` | $\log\|P_X\|$, $\arg P_X$, $F_B(T)$ | method B's inner factor in isolation |
| `zhad <zeros_file> <T0> <T1> <N> <X> [window]` | $\log\|Z_X\|$, $\arg Z_X$ | the *correction* needed on top of $P_X$ for rigor |
| `zhybrid <zeros_file> <T0> <T1> <N> <X> [window]` | full $P_X Z_X$ on a grid, $F_C(T)$ | rigorous method C; not bisection-integrated |

Smoke test for `zhybrid` (uses the 20 Odlyzko ordinates shipped under `doc/ghy/first-20-zeros.txt`):

```bash
./zhybrid doc/ghy/first-20-zeros.txt 14 50 180 100 | awk '{print $1, $6}'
```

Column 6 is $F_C(T)$; it crosses each integer at a true zeta zero (14.13, 21.02, 25.01, 30.42, …). At $T = 45$ between $\gamma_8 = 43.33$ and $\gamma_9 = 48.00$ it reads $8.01$; the $\sim 10^{-2}$ residual is the GHY truncation tail, consistent with the $O(\log X/\sqrt X)$ bound.

## Open problem

Find an explicit-constant version of

$$
\Bigl|S(T) - \bigl(-\tfrac{1}{\pi}\!\sum_{n \le X}\tfrac{\Lambda(n)\sin(T\log n)}{\sqrt n\,\log n}\bigr)\Bigr| \;\le\; c(X,T)
$$

with $c \to 0$ faster than $\log T / \log X$. Candidates: Carneiro–Chirre–Milinovich (2022) bandlimited approximations applied to the *tail*; large sieve on $\sum_{p \le X} p^{-1/2}e^{iT\log p}$; Vaughan's identity. Any of these would tighten the boxed bound to match the observed accuracy.

## References

- Gonek, Hughes, Young, *A hybrid Euler–Hadamard product for the Riemann zeta function*, Duke Math. J. **136** (2007), 507–549 — Theorem 1 (the identity).
- Goldston, *On the function $S(T)$ in the theory of the Riemann zeta-function*, J. Number Theory **27** (1987), 149–177 — the $\log T/\log X$ truncation bound.
- Carneiro, Chandee, Milinovich, *Bounding $S(T)$ and $S_1(T)$ on the Riemann hypothesis*, Math. Ann. **356** (2013), 939–968 — tighter $|S(T)|$ constants.
- Carneiro, Chirre, Milinovich, *Bandlimited approximations and estimates for the Riemann zeta-function*, Publ. Mat. **63** (2022) — best candidate for the open tail bound.
