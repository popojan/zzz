# Does the zeros' own rigidity break the forward wall? A decisive test

*2026-06-21. [`band-saturation.md`](band-saturation.md) left one crack open: below
$\kappa=\pi$, all prior negatives used band-limited **prime** data or continuous
$L^2$ estimators — both provably in the null space — so none could see whether
the **zeros' own point-process rigidity** (the unit point process; "integrality
of the zeros") selects the true configuration. "Prior's against it (five
negatives)" was a posture, not a result. This note settles it. Scripts:
`doc/ghy/zero-rigidity-resolution.wls` (resolution + synthetic control),
`doc/ghy/zero-rigidity-realdB.wls` (the real B-error field).*

## Design, engineered against a false verdict

We test the intrinsic question the five negatives structurally could not: *how
well does the zeros' process law pin one zero given the others, and does that
remove the error B actually makes?* Every choice maximizes power for the
hypothesis, so a negative is decisive:

- **exact** Odlyzko zeros (`zeros1`), unfolded to $u_n = n - N_0(\gamma_n) =
  S(\gamma_n)$ (the rigid log-correlated fluctuation; a position error
  $\delta\gamma = u\cdot\mathrm{gap}$, so $\sigma$ in $u$-units **is**
  $\sigma$/gap — height-independent);
- **true neighbours** (no self-consistency / no circularity);
- the **optimal** predictor: best linear predictor = conditional mean for the
  asymptotically-Gaussian $S$ (Selberg CLT) — so no power lost to a weak
  estimator — with two-sided **interpolation** (maximum information);
- **out-of-sample** residuals (train 60% / test 40%) — no overfit illusion;
- a **power control** in the same harness: the same denoiser is shown to remove a
  *smooth* displacement, so a null result on the real error is teeth-proven, not
  dead.

## Test 1 — rigidity resolution: a real (surprising) positive

Best two-sided interpolation of $u_n$ from $m$ neighbours each side
($N=60000$ zeros, $\gamma\le47531$, gap $0.704$; baseline $\sqrt{\mathrm{Var}\,u}=0.307$):

| $m$/side | 1 | 2 | 4 | 8 | 16 | 24 | 32 |
|---|---|---|---|---|---|---|---|
| $\sigma_{\rm cond}$ (gap) | 0.303 | 0.302 | 0.288 | 0.113 | 0.024 | 0.0153 | **0.0096** |

$\sigma_{\rm cond}$ keeps falling — from 64 true neighbours a zero is pinned to
**0.01 gap**, an order of magnitude below B's $\sim$0.13. This is genuine
**Ghosh–Peres rigidity** (the sine process is number- and centre-of-mass-rigid,
so individual points are strongly determined by the rest), here quantified. The
precision floor of the data ($\sim2\times10^{-9}$ gap) is far below 0.0096, and
out-of-sample protects against overfit, so the resolution is real.

## Test 2 — the denoiser has teeth (synthetic control)

Now the forward question: given B's estimate of **every** zero (a displaced
$u+d$), can the rigidity prior remove $d$? Optimal oracle denoiser (estimate
clean $u_n$ from the corrupt window incl. centre, weights fit on
$\text{corrupt}\!\to\!\text{clean}$, out of sample), injecting a displacement of
B's magnitude $\Delta=0.13$ gap at varying coherence length $L_c$:

| $L_c$ (gaps) | 0 (white) | 1 | 2 | 3 | 5 | 8 | 12 |
|---|---|---|---|---|---|---|---|
| corrected | 0.11 | 0.13 | 0.21 | **0.33** | 0.54 | 0.77 | 0.92 |

The denoiser corrects a **smooth/coherent** displacement well (and white noise
poorly) — because $S(\gamma_n)$ is itself a wiggly field, so a smooth $d$ is
spectrally **separable** from it. (This refutes the old intuition that "a smooth
coherent warp looks like signal.") Teeth proven.

## Test 3 — the real below-threshold B error: decisive negative

Replace the synthetic $d$ with B's **actual** error field
$d_B = \arg P_X(\gamma_n)/\pi - S(\gamma_n)$ (the prime-truncation error of $S$),
at a sub-threshold cutoff $X$ ($\kappa=\mathrm{gap}\cdot\log X$), block
$n\in[24001,36000]$, $\gamma\approx25755$, gap $0.755$ — same oracle denoiser:

| $X$ | $\pi(X)$ | $\kappa$ | $\mathrm{SD}(d_B)$/gap | corrected |
|---|---|---|---|---|
| 29 | 10 | 2.54 | 0.140 | **0.021** |
| 49 | 15 | 2.94 | 0.125 | **0.031** |
| 64 | 18 | 3.14 ($\approx\pi$) | 0.119 | 0.035 |
| 121 | 30 | 3.62 | 0.106 | 0.065 |
| 256 | 54 | 4.19 | 0.097 | 0.106 |

Below $\kappa=\pi$ the rigidity prior removes **2–3%** of the real B error —
nothing — even though Test 2 proved it removes 33–92% of a smooth displacement of
the *same magnitude*. So the synthetic positive was too easy: **$d_B$ is not a
smooth warp.** It is the high-frequency tail of $S$ (the oscillations of primes
$p>X$, frequencies $\log p > \log X$), which lands exactly in the
**locally-unpredictable** component that Test 1 leaves as irreducible residual.
Rigidity sharply determines the smooth/long-range part of the zeros — *but that
is the part B already gets right.* It is blind to precisely the part B gets wrong.

And the correction rises *through* the threshold (2–3% below, 6–10% above
$\kappa=\pi$): the rigidity denoiser **independently reproduces the
band-saturation wall** — surplus to extract only above $\kappa=\pi$, exactly
where `--boot` gains.

## Verdict

The integrality/rigidity of the zeros does **not** break the forward degeneracy
below $\sqrt T$. This is now a *result with a mechanism*, not a tally of
negatives: B's error is the high-frequency prime tail of $S$, which is orthogonal
to the (smooth, long-range) structure the zero process pins down. The crack found
a real positive — the zeros self-determine to 0.01 gap (Ghosh–Peres) — that is
**orthogonal to what's needed.** The one remaining hedge: the optimal denoiser is
*linear*, which is conditional-mean-optimal for Gaussian $S$ (Selberg) and in the
$\Delta=0.13$-gap linear-response regime; a fully nonlinear determinantal
estimator is untested, but turning 3% into anything useful would need the
process's higher correlations to carry order-1 information about the high-frequency
tail, which Gaussianity makes implausible. Sub-floor rigidity is, for practical
purposes, closed: the wall stands against primes **and** primes-plus-zero-rigidity
alike.

## Artifact index

`doc/ghy/zero-rigidity-resolution.wls` (+ `.tsv`, `.txt`) — Tests 1–2;
`doc/ghy/zero-rigidity-realdB.wls` (+ `.txt`) — Test 3. Companions:
[`band-saturation.md`](band-saturation.md) (the wall this fails to cross),
[`spectral-shadow.md`](spectral-shadow.md) (why B's error preserves the
*statistics* — same common-mode story, dual face),
[`snap-margin-detector.md`](snap-margin-detector.md) (the backward dual: integer
rigidity of the *primes*, which **is** decisive).
