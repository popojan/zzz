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

**Long range — over-rigid, an artefact (the genuinely new finding).**

| $L$ | $\Sigma^2_{\rm obs}$ | GUE $\tfrac{1}{\pi^2}\log$ | Poisson $L$ | $\Delta_3^{\rm obs}$ | GUE | Poisson $L/15$ |
|---|---|---|---|---|---|---|
| 5  | 0.44 | 0.51 | 5   | 0.143 | 0.156 | 0.33 |
| 20 | 0.42 | 0.65 | 20  | 0.190 | 0.297 | 1.33 |
| 50 | 0.42 | 0.74 | 50  | 0.193 | 0.389 | 3.33 |
| 200| 0.40 | 0.88 | 200 | — | — | — |

$\Sigma^2(L)$ and $\Delta_3(L)$ **saturate** (~0.4, ~0.19) instead of growing like
$\log L$ — the loop's zeros are *more rigid than the true zeros.* Mechanism: the
loop fixes $N_0(\gamma_k)+\tfrac1\pi\arg P_X = k-\tfrac12$, so unfolding by $N_0$
leaves only the **band-limited** fluctuation $\tfrac1\pi\arg P_X$, the truncated-
prime approximation to $S(t)=\tfrac1\pi\arg\zeta$. The long-range rigidity
$\Sigma^2\sim\log L$ lives in the fluctuation content $\arg P_X$ **cannot reach**;
the loop therefore pins zeros to the smooth count and manufactures an artificial
super-rigidity. Band-saturation, in spectral-statistics form: **the band carries
short-range correlations faithfully and long-range rigidity not at all.**

## Verdict

The zeta-free shadow carries the operator's *local* spectral fingerprint
faithfully and **distorts** its *long-range* one (too rigid, and misleadingly so —
one would be deceived trusting the loop's zeros for $\Sigma^2$/$\Delta_3$). So for
Hilbert–Pólya the quality buys a clean, circularity-free **lens** on the duality
and its selection-rigidity, and a measuring instrument for *which* textures the
prime-only data determine (answer: short-range yes, long-range no). It does not
buy leverage on positivity or construction, where the program is actually
blocked. Lens, not lever — now with the boundary measured rather than guessed.

*Methodological note: the GUE "wash-out" claim above was an educated guess stated
without data; the spacing histogram refuted it within the hour, and the
over-rigidity it surfaced is the real result. Look, don't guess.*
