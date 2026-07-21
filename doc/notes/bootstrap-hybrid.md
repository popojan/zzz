# The self-consistent hybrid bootstrap (`zzz --boot`)

*2026-06-10. Companion scripts: `doc/ghy/cutoff-sweep.sh`, `doc/ghy/cutoff-osc-analysis.wls`,
`doc/ghy/boot-validate.sh`, `doc/ghy/boot-validate.wls`.*

## Idea

Method C (GHY hybrid $P_X \cdot Z_X$) is the most accurate of the three counters but
needs seed zeros — useless at heights where nobody has computed them. The bootstrap
makes method C **self-hosting**: it manufactures its own seeds with method B, then
iterates the hybrid count to a fixed point. The chain stays zeta-evaluation-free;
no zero tables are consumed.

## Why the dominant error of B is recoverable

Method B's counting function satisfies

$$F_B(T) = N(T) + \frac{1}{\pi}\sum_\rho \operatorname{Im} E_1\big(i(T-\gamma_\rho)L\big),
\qquad L = \log X,$$

and $\operatorname{Im} E_1(iyL) = \mathrm{Si}(yL) - \tfrac{\pi}{2}\,\mathrm{sign}(y)$.
Near a target zero $\gamma_n$ the slope of $F_B$ is $L/\pi$ (the target's own smeared
step), so solving $F_B = n - \tfrac12$ displaces the located zero by

$$\Delta(L) = -\frac{1}{L}\sum_{j\neq n}\Big[\mathrm{Si}\big((\gamma_n-\gamma_j)L\big)
  - \frac{\pi}{2}\,\mathrm{sign}(\gamma_n-\gamma_j)\Big].$$

This is **parameter-free** and was verified empirically by sweeping the cutoff
(`cutoff-osc-analysis.wls`): correlation 0.96–0.98 between measured and predicted error
at $n = 10^3 \ldots 10^5$ and $n = 10^{12}$; the model accounts for 54–80 % of the rms.
The displacement is exactly the quantity that the GHY local Hadamard factor

$$\log Z_X(s) = -\sum_j E_1\big((s-\rho_j)\log X\big)$$

subtracts — *if* one knows the neighbours $\gamma_j$. Approximate neighbours suffice:
a neighbour seeded with error $\varepsilon$ leaves only a second-order residual
$\sim \varepsilon \sin(\delta L)/\delta$, which is why B-quality seeds
(themselves wrong by the very oscillation being corrected) still help, and why
iterating contracts.

Two dead ends, established by the same sweep, are worth recording:

- plain averaging of $\gamma(X)$ over cutoffs buys only 2–4× (the nearest-neighbour
  oscillation period in $\log X$ exceeds any feasible cutoff range);
- free-frequency least-squares extrapolation of $\gamma(X)$ to $L \to \infty$ is
  ill-conditioned (the constant and $1/L$ basis functions are degenerate over a
  short $L$-range).

The smooth part of the missing prime tail is *not* the problem: it equals
$E_1((s-1)\log X)$, of magnitude $\sqrt X/(T\log X)$ — about $10^{-33}$ at
$T \sim 10^{36}$. The entire recoverable error lives in the local zeros.

## Algorithm (`main.c::boot_locate`)

1. **Seed**: locate zeros $n-W \ldots n+W$ with plain method B
   (bisection of $F_B$ from the asymptotic ordinate).
2. **Relocate**: for `--rounds` rounds (default 2), re-bisect every zero of the
   *inner core* ($\pm W/2$ around the target) with the leave-one-out hybrid count
   $$F_C^{(j)}(T) = N_0(T) + \frac{\arg P_X + \arg Z_X^{\,\mathrm{excl}\, j}}{\pi},$$
   Jacobi-style (all relocations of a round read the previous round's offsets).
3. **Report** the target's final relocation.

Design points, each of which was load-bearing in testing:

- **Leave-one-out is mandatory.** With the target included, $Z_X$ contributes the
  target's sharp step at its *seeded* position and $F_C$ simply reproduces the seed
  (self-fulfilling). With it excluded, the target's step enters only through
  $\arg P_X$, smeared symmetrically about the *true* ordinate — that is the signal.
- **Guard ring.** Window-edge zeros have one-sided neighbour coverage; relocating
  them is biased, and the bias propagates inward over rounds. Only the inner
  $\pm W/2$ core is relocated; the outer ring keeps its B seeds, serving as
  subtraction terms only (where seed error is second order).
- **No basin hopping.** Relocation bisections get a bracket of $0.75$ mean gaps,
  no auto-widening, and any relocation that moves a zero by more than $0.45$ of
  the bracket is rejected (the seed is kept). A single zero jumping to a
  neighbouring crossing corrupts the count for the whole window.
- **Jacobi, not Gauss–Seidel**: in-place updates propagate an edge artifact across
  the window within one round.

## Numerical pitfall: $E_1$ on the critical line

On the critical line the kernel argument $(s-\rho_j)\log X = iyL$ is *purely
imaginary*, with $|yL|$ up to (window span)$\times \log X$ — several hundred at low
heights. The complex power series for $E_1$ cancels $\approx |yL|/\log 2$ bits
before `acb_hypgeom_expint` switches to its asymptotic regime, so at 256-bit
precision the mid-range returns wide balls. Wide balls make `arb_gt`/`arb_lt`
return "indeterminate" (false), which a bisection loop silently misreads — the
observed symptom was zeros drifting *rightward* at the rejection cap, round after
round, only for $W \ge 16$, only at low heights, and curable by `-p 512`.

Fix (`ghy.c::e1`): for purely imaginary arguments use real Si/Ci,

$$E_1(iy) = -\,\mathrm{Ci}(|y|) + i\big(\mathrm{Si}(|y|) - \tfrac{\pi}{2}\big),
\quad \text{conjugated for } y<0,$$

which is stable at any window width and precision (and cheaper than complex $E_1$).

## Results

Validation against Odlyzko tables (`boot-validate.sh` / `boot-validate.tsv`:
20 consecutive zeros per height, `--boot 32 -R 2`, mean |error|):

| height | ordinals | k | gap·log X | B mean err | boot mean err | gain |
|---|---|---|---|---|---|---|
| low | 1000–1019 | 1000 | 10.4 | 0.0153 | **0.0008** | **19.3×** |
| mid | 99900–99919 | 1000 | 6.0 | 0.0179 | 0.0064 | 2.8× |
| e12 | 10¹²+30… | 10⁴ | 3.0 | 0.0089 | 0.0044 | 2.0× |
| e21 | 10²¹+30… | 10⁴ | 1.6 | 0.0242 | 0.0254 | 0.95× |
| e22 | 10²²+30… | 10⁴ | 1.4 | 0.0119 | 0.0128 | 0.93× |

W-scan at $n=1000$ (truth 1419.422481, B err +0.0254, k=1000):
W=4 → +0.0035, W=8 → +0.0032, W=16 → +0.0018, W=32 → **+0.0009**, W=64 → −0.0014.
Rounds converge there: at W=32 the R=1..4 errors are +0.0011, +0.0009, +0.0009, +0.0009.

Note what 19× means against the $c/\log X$ law: to match boot's 0.0008 at
$n\approx1000$, plain B would need $\log X \approx 170$ — no feasible prime count
reaches that. **Inside its validity domain the bootstrap breaks the $1/\log X$
wall**, not just its constant.

## The validity threshold: $X \gtrsim \sqrt{T/2\pi}$

The gain collapses to 1 almost exactly where $\mathrm{gap}\cdot\log X = \pi$,
i.e. $\log X = \tfrac12 \log(T/2\pi)$, i.e. $X = \sqrt{T/2\pi}$ — the
**Riemann–Siegel scale**. Interpretation: the kernel $\mathrm{Si}((T-\gamma)L)$
must *resolve individual neighbours* for leave-one-out subtraction to carry
information. Below the threshold the truncation deficit is a field coherent
across many consecutive zeros; B's seeds inherit that deficit *coherently*, so a
$Z_X$ built from them is wrong in exactly the way that cancels the correction
(a soft collective-translation mode). Two ensemble checks confirm this:

- e21/e22 at k=10⁴ (gap·L ≈ 1.4–1.6): gain 0.95×/0.93× over 20+10 zeros;
- the 10³⁶ regime reproduced at 10²² by lowering k to 600 (gap·L ≈ 1.0,
  `boot-regime-e22-k600.tsv`): B 0.0261, W32R1 0.0267, W32R2 0.0258,
  W32R3 0.0298, W64R2 0.0273 — no recipe beats B.

Multi-round iteration is *unstable* below the threshold (the neighbour coupling
$\sin(\mathrm{gap}\,L)/(\mathrm{gap}\,L) \to 1$): at $10^{36}$, `--boot 64 -R 3`
showed round shifts growing 0.009 → 0.018 → 0.027 with rejections, final error
+0.035 (worse than B's +0.021). Single-zero runs at $10^{36}$ that beat B
(`--boot 32 -R 2`: +0.0057) are consistent with luck; the ensemble expectation
there is no gain.

Showcase zero #$10^{36} + 42420637374017961984$
(published $\gamma$ = 81029194732694548890047854481676712.98790, Bober–Hiary), k=10⁴:

| method | located | error |
|---|---|---|
| A heuristic | …713.009431 | +0.0215 |
| B `--ghy` | …713.009348 | +0.0214 |
| `--boot 32 -R 1` | …713.009877 | +0.0220 |
| `--boot 64 -R 1` | …713.009577 | +0.0217 |
| `--boot 32 -R 2` | …712.993619 | +0.0057 (single draw, not expectation) |
| `--boot 64 -R 3` | …713.022412 | +0.0345 (unstable iteration) |

## Practical guidance

- $T \lesssim 10^9$ (k ≥ 10³): order-of-magnitude gains; `--boot 32 -R 2` is a
  good default and the rounds are convergent.
- $T \sim 10^{10}$–$10^{13}$ (k ~ 10⁴): 2–3× gains.
- $T \gg 2\pi X^2 \cdot O(100)$: stay with plain `--ghy`; the bootstrap is
  inert in expectation and multi-round can hurt. With k = 10⁴ the boundary
  sits near $T \sim 10^{13}$.

## Open ends

- Below the Riemann–Siegel threshold the deficit is a smooth field that
  single-cutoff data cannot see past (the seeds carry it coherently). The one
  untried lever: **multi-cutoff (multi-L) data at the same T** — the same zero
  field smeared at different kernel widths — used as a deconvolution problem
  with B seeds as priors. The cutoff-sweep experiment showed free-frequency
  fitting is ill-conditioned, but seeded/regularised multi-L inversion is open.
  *[Resolved, negative: pursued the same day as experiments E1–E4 (Si-kernel
  ML, multi-cutoff ML, exact Weil-identity ML with a truth control). Below
  gap·log X = π, method B saturates the information in the primes — see
  `band-saturation.md` for the synthesis, `superres-multi-cutoff.md` for the
  lab log.]*
- Cost: pass 1 dominates (every neighbour bisected from the asymptotic guess
  with the full default window). Marching seeds (neighbour $j+1$ from located
  $j$ + mean gap) would cut pass-1 cost several-fold. Not yet implemented.
- A marching variant could sweep a whole ordinal range, reusing each window's
  refined core as the next window's guard — amortised cost per zero close to
  plain B, with the full bootstrap gain inside the validity domain.
