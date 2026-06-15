# The spectral shadow: one invariant, the GUE texture, and the Hilbert–Pólya altitude

*Synthesis of the 2026-06-15 discussion that followed the self-paving loop
([`zeros-primes-bootstrap.md`](zeros-primes-bootstrap.md)). It records a
conceptual thread — can the alternating loop be collapsed to a single invariant,
and does the zeta-free output carry anything useful for Hilbert–Pólya? — and the
**data** that settled (and corrected) the speculation. Scripts:
`doc/ghy/loop-spectral-stats.py`. Includes a claim of mine that the data refuted;
that is left in deliberately, because the correction is the point.*

## One invariant behind the two steps

The loop alternates primes→zeros (forward, method B) and zeros→primes (backward,
$\psi'$ detector). Both are the **Riemann–Weil explicit formula** read in two
directions:
$$\sum_{\rho} h(\gamma_\rho)=\text{arch}(g)-2\sum_{p,m}\tfrac{\log p}{p^{m/2}}\,g(m\log p),\qquad h=\hat g.$$
A pair $(P,Z)$ is *consistent* iff the residual $W[P,Z;g]=0$ for all $g$. Forward
holds $P$ and solves for $Z$; backward holds $Z$ and solves for $P$ — **both
project onto the same manifold $\{W\equiv0\}$.** The loop is alternating
projection; the single invariant is that manifold's defining equation, and its
multiplicative form is the GHY hybrid $P_X\!\cdot\!Z_X\approx\zeta$ (`zhybrid`).

**A joint "Weil-residual verifier" is therefore near-vacuous on the loop's own
output:** the primes and zeros there were *derived from each other through $W$*,
so $W\approx0$ by construction, and the band-saturation errors live in $W$'s null
space. It detects only independently-injected faults, not the loop's own error.
Worth noting, not worth building as a validator.

## The Hilbert–Pólya altitude (where the loop cannot reach)

The explicit formula is **operator-agnostic**: provable from $\zeta$'s functional
equation, true whether or not a self-adjoint $\hat H$ with spectrum $\{\gamma\}$
exists. Knowing it amounts to knowing the spectrum plus its prime-dual — and a
spectrum does not determine its operator (*you can't hear the shape of a drum*).
So computing the shadow to any precision does not triangulate the caster; the
underdetermination is structural. HP is stuck not on *computing* zeros (now
doable many ways, including zeta-free) but on **positivity / self-adjointness**
(Weil positivity), which no computation touches, and on the absence of a
number-field analogue of the one solved case — Frobenius on $\ell$-adic
cohomology for curves over finite fields, the only place the caster is visible.

**Zeta-free + primality-free is an epistemic property of the *algorithm*, not a
structural handle on the *operator*.** What it genuinely changes is the *stance*:
it exhibits the primes↔zeros duality as a standalone constraint with a **solution
manifold** (the band-saturation null space), the true pair selected by prime
**integrality** + a few exact anchors. "What selects the true pair" is a shadow
of "what rigidity does the operator impose" — the right neighbourhood, but
characterising that rigidity *is* essentially constructing the selector. A lens,
not a lever.

## The data: how much spectral fingerprint survives the zeta-free output?

The answerable version of the question. Measured on the loop's own zeros
(292 376 of them, $\gamma$ up to $2\times10^5$, unfolded by $\theta(t)/\pi+1$):

**Short range — GUE survives (I had guessed it would wash out; wrong).**

| | observed | GUE | Poisson |
|---|---|---|---|
| nearest-nbr $P(s<0.5)$ | 0.163 | 0.112 | 0.393 |
| histogram $L^2$ | — | **0.56** | 3.72 |

The spacing distribution tracks the Wigner curve, $7\times$ closer to GUE than
Poisson. The $\sim5\%$-of-a-gap error is a *correlated* common mode that cancels
in spacings, so local level repulsion passes through intact. (The very smallest
$s$ is additionally floored by the marching bracket, so that end is
algorithm-limited, not physics.)

**Long range — saturated, and the loop reproduces it (corrected finding).**

$\Sigma^2(L)$ and $\Delta_3(L)$ for the loop *saturate* (~0.4, ~0.19) instead of
growing like GUE's $\tfrac1{\pi^2}\log L$. My first reading was "over-rigid
artefact." **That was wrong** — the *true* zeros saturate too. Comparing the loop
to `ZetaZero` over the **same index block** (#30000–31800, $\gamma\approx25755$;
`doc/ghy/loop-vs-true-sigma2.py`):

| $L$ | $\Sigma^2_{\rm TRUE}$ | $\Sigma^2_{\rm LOOP}$ | GUE $\tfrac1{\pi^2}\log$ |
|---|---|---|---|
| 2  | 0.370 | 0.397 | 0.416 |
| 10 | 0.335 | 0.355 | 0.579 |
| 20 | 0.383 | 0.407 | 0.650 |
| 40 | 0.405 | 0.428 | 0.720 |

The true zeros are *already saturated* at this height, far below the GUE log-law,
and the loop tracks them to ~0.03 (if anything a hair *less* rigid — the 5% jitter
adds variance, it doesn't remove it). This is **Berry's semiclassical saturation**:
at finite $T$ the number variance follows universal GUE only up to $L_{\max}(T)$,
then saturates at a value fixed by the **short** primes (the long periodic orbits).
The loop has those primes, so it reproduces the saturation. No missing mechanism,
no artefact — the same explicit formula, weighted by scale: large primes (band
edge) set the short-range correlations, small primes set the long-range saturation,
and the loop holds both within its feasible heights.

## Verdict

The zeta-free shadow carries the operator's spectral fingerprint **faithfully on
both scales** within the loop's feasible range — local GUE *and* the prime-governed
long-range (Berry-saturated) rigidity, matched to the true zeros. The
*orbit-dependent* part (the saturation value / $L_{\max}$, fixed by the actual
primes) is the genuinely caster-adjacent texture, and it survives the band-
saturation floor. What would eventually drift at asymptotically high $T$ — where
$L_{\max}$ outgrows the fixed band — is the long-range *statistical* fidelity, but
that is a harmless divergence in a collective statistic, not a per-zero breakdown
(detection is local; the per-zero error is height-independent in gap-units), and it
lies beyond any feasible climb. So for Hilbert–Pólya the quality is a clean,
circularity-free **lens** that faithfully shows the arithmetic spectral structure —
not a **lever** on positivity or construction, where the program is blocked.

*Methodological note: this section was wrong twice before it was right. "GUE will
wash out" (guess) was refuted by the spacing histogram; "long-range is an artefact"
(guess) was refuted by the same-block comparison to `ZetaZero`. Both corrections
came from looking. Look, don't guess — it's the throughline of this whole project.*
