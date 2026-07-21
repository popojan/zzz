# Exact primes from inexact zeros: where the exactness lives

*Capstone of the 2026-06-21 arc. One statement ties together the whole project:
**the exact prime sequence is recoverable from zeta zeros known only to a fixed
fraction of a gap — never exactly, never to sub-gap precision — ζ-free and
primality-free, from a finite seed.** This note argues that this is not a paradox
against "the primes are governed by the exact non-trivial zeros," but the
self-consistency of the prime↔zero duality, and that the exactness was never in
the zeros: it lives in $\mathbb{Z}$. Synthesises
[`band-saturation.md`](band-saturation.md),
[`zeros-primes-bootstrap.md`](zeros-primes-bootstrap.md),
[`snap-margin-detector.md`](snap-margin-detector.md),
[`zero-rigidity.md`](zero-rigidity.md), [`backward-reach.md`](backward-reach.md),
[`spectral-shadow.md`](spectral-shadow.md).*

## The statement, at full precision

`zzz --loop` starts from 40 seed ordinates and, alternating the explicit formula
in both directions, streams an unbounded sequence of **exact** primes (verified
flawless: 4450 primes $=\pi(42557)$ exactly, 0 false / 0 missed) while never
evaluating $\zeta$, never testing divisibility, and never holding a zero to better
than $\sim 5\%$ of its local gap (median error $0.036$, decades coarser than the
$\sqrt T$ precision wall). It is super-critical (it generates a prime *surplus*),
so the prime frontier climbs without internal limit.

The honest scope: this is **not "primes from nothing."** The loop needs the
finite seed and the gap-accurate zeros — which it manufactures from the primes it
already has. The content is a statement about **precision and direction**:

> exact primes ⟸ inexact zeros, and no exact zero ever appears.

## Why this is not a paradox — the exactness lives in $\mathbb{Z}$

"The primes are governed by the exact zeros" quietly conflates two things:

- the zeros **as a set / a distribution** — which governs the *statistics* of
  primes (the prime number theorem, the size of the fluctuations);
- the **exact individual positions** — which govern the *real-valued* function
  $\psi(x)=x-\sum_\rho x^\rho/\rho-\dots$ to high precision.

The prime *set* is integer-valued. So the map zeros $\to$ primes ends in a
**quantizer**: the backward detector reads $\psi'(x)$ and asks the discrete
question "is there von Mangoldt mass at the integer $x$?" Any quantized readout
returns the exact answer from an input good to within half a step. **Discreteness
is an error-correcting code** — the same mechanism that lets the Pell work snap a
floating-point regulator to the exact integer solution
($\mathrm{round}[\cosh R]$), here applied to $\Lambda(x)$.

Quantitatively, the precision the backward step needs is set by the requirement
that the phase $\gamma_k\log x$ be good to $O(1)$ over the scan window, i.e.
$\Delta\gamma_k \ll 1/\log x$ — **a fixed fraction of a gap, height-independent.**
Never sub-gap, never exact. The absolute precision shrinks with height (gaps
shrink), but in gap units it is a constant, coarse number. The "infinite
precision" the exact-zeros image suggests is an artifact of picturing $\psi$ as a
real number; once the output is "which $n$ are prime," the requirement collapses
to gap-scale.

## Why it is *required* by consistency — the duality cannot eat its own tail

This is the deeper reason, and it inverts the surprise. Band-saturation
([`band-saturation.md`](band-saturation.md)) proved that the primes
$p\le X=\sqrt{T/2\pi}$ **cannot determine** the sub-gap positions of the zeros up
to $T$: that information is not attenuated, it is *absent* (truth and the
displaced method-B configuration agree against every band-limited prime functional
to the $5\times10^{-6}$ data floor).

Now suppose the primes genuinely *needed* the exact zeros. Then the explicit
formula would be inconsistent: the primes would require information about the zeros
that the primes themselves provably cannot supply. The only coherent state is the
measured one —

> **"exact zeros are unnecessary for the primes" is the same fact as "the primes
> cannot pin the exact zeros."**

Forward **null-space** = backward **don't-care**. Both are one kernel of width
$\pi/\log X$ read two ways: what method B cannot resolve (forward) is exactly what
the quantizer discards (backward). So the duality is impedance-matched at
gap-scale, and it would be *paradoxical* if you **did** need exact zeros. The
zero-rigidity test ([`zero-rigidity.md`](zero-rigidity.md)) sealed the same point
from the other side: the zeros' own point-process rigidity, applied at oracle
strength, recovers only 2–3% of the real method-B error below $\kappa=\pi$ —
because that error is the high-frequency tail of $S(t)$, the part no one needs and
no one can reach.

## The three pillars, one statement

| direction | mechanism | precision verdict |
|---|---|---|
| primes → zeros (forward) | band saturation: counting through a kernel of width $\pi/\log X$ | **information-limited**: gap-scale is all primes $\le X$ can ever say; sub-gap is absent |
| zeros → primes (backward) | integer-snapping: a quantizer onto $\mathbb{Z}$ | **detector-saturated, super-critical** ($\alpha\approx1$): gap-scale zeros give exact primes, surplus to spare |
| the loop | alternating projection on the Weil manifold | rides $\kappa\approx6$, **above** the wall, climbing; dies of feasibility, never of band saturation |

The forward floor (what B cannot make finer) and the backward tolerance (what the
quantizer ignores) are the same gap-scale quantity. That coincidence is not luck —
it is the duality being self-consistent.

## Where the wall is, and why the loop never meets it

The band-saturation wall is a **curve**, $X=\sqrt{T/2\pi}$ ($\kappa=\pi$), not a
point. It is the cost floor for the *different* task of computing **one isolated
zero at height $T$ from the fewest primes** (the $10^{36}$ showcase zero needs
$X\approx10^{17}$, Riemann–Siegel-class). The self-paving loop is exempt: backward
super-criticality keeps it $10^2$–$10^6\times$ over-supplied with primes
($X_{\rm loop}\approx0.11\,T$ vs wall $0.40\sqrt T$; currently $X=42557$ at
$T\approx4\times10^5$ is $169\times$ above, and the margin widens as $0.27\sqrt T$).
So `--loop` never fails at the wall — it fails at **feasibility** ($T\sim10^8$,
$X\sim10^7$: tens of GB, $\sim10^{16}$ ops), and, if compute were free, at
**double-precision phase folding** ($T\sim10^{13}$). The wall is real; it is
simply not on the loop's road.

## What still needs the exact zeros (scope, so this is not oversold)

The *discrete* fact "$n$ is prime" does not need them. The *continuum* does:

- $\pi(x)-\operatorname{Li}(x)$ as an actual real number, and the precise
  oscillation amplitudes — these are real-valued and need the zeros to matching
  precision;
- RH itself and Weil positivity — a sign condition on a real form, untouched by
  any amount of zero *computation* ([`spectral-shadow.md`](spectral-shadow.md));
- the Hilbert–Pólya operator — a spectrum does not determine its operator;
  computing the shadow to any precision does not triangulate the caster.

The loop recovers the **arithmetic** (the set, exactly); it does not recover, and
does not need, the **analysis** (the continuum) exactly.

## Capstone

> **Coarse spectrum + integer rigidity = exact arithmetic.** The exactness never
> crosses the channel — it is manufactured at the receiver by the integers
> themselves. The zeros (as a gap-accurate set up to $\sim T$) are needed and
> used; their *exactness* is not, and — by band saturation — could not be, because
> the primes that would carry it cannot themselves resolve it.

The vivid demonstration (`--loop`) makes the point unmissable, but the point is
structural: "the exact zeros encode the primes to infinite precision" is a
category error. The primes were always the discrete shadow that a band-limited,
gap-resolved view of the spectrum suffices to cast exactly.

## Relation to the literature & status

Pressure-tested against the literature (web search, 2026-06). The honest verdict:
**the headline is largely folklore, and one piece is already published; one
empirical claim is genuinely this project's and is the load-bearing risk.**

**Folklore — and the headline is published.** The prime↔zero duality (the
explicit formula as a Fourier-type transform) and rebuilding the prime staircase
from finite cosine sums over zeros are textbook: Mazur & Stein, *Prime Numbers and
the Riemann Hypothesis* (Cambridge Univ. Press, 2016); survey/quantum-chaos
framings in [2204.00899](https://arxiv.org/abs/2204.00899),
[nlin/0212042](https://arxiv.org/abs/nlin/0212042). The specific statement that
**approximate (not exact) zeros already reconstruct $\pi(x)$** is in print: França
& LeClair derive the $n$-th zero from an *asymptotic* (Lambert-$W$) transcendental
equation and show those approximate zeros suffice to reconstruct the
prime-counting function — [1307.8395](https://arxiv.org/abs/1307.8395),
[1502.06003](https://arxiv.org/abs/1502.06003). Finite-precision computation of
$\pi(x)$ from a truncated, finite-accuracy set of zeros is standard practice
(Platt, [1203.5712](https://arxiv.org/abs/1203.5712); Büthe,
[1410.7008](https://arxiv.org/abs/1410.7008) and the analytic-$\pi(x)$ program —
the record $\pi(10^{24})$ used $\sim$36 billion zeros to $\pm2^{-102}$, "more than
suffices"). Recent work still treats the $n$-th prime via spectral/zero sums
directly ([2601.18816](https://arxiv.org/abs/2601.18816),
[2603.07641](https://arxiv.org/abs/2603.07641)). So the capstone's *headline* —
exact primes from inexact zeros — is **known, not new**; this note is an
exposition of it, made vivid by the ζ-free / primality-free self-paving loop.

**The resolution scale is the textbook one.** The truncated explicit formula has
error $\sim x\log^2 x/T$, so resolving *individual* primes near $x$ needs zeros up
to $T\sim x$, while $T\sim\sqrt x$ suffices only for smooth/averaged statements
(standard; e.g. Granville's course notes,
[Montréal](https://dms.umontreal.ca/~andrew/Courses/Chapter10.pdf)). Our
band-saturation scale $X=\sqrt{T/2\pi}$ is exactly that $\sqrt x$ / Riemann–Siegel
scale — corroborated, not invented.

**Ultra-high zeros independently corroborate the gap-scale, height-independent
law.** França & LeClair compute the $10^{100}$-th (even $10^{1000}$-th) zero to
$\sim$100 *significant figures* from $N=5\times10^6$ primes
([1601.00914](https://arxiv.org/abs/1601.00914)). They are careful about what this
means: the estimate is "correct up to the decimal point, i.e. to the number of
digits in the integer part" ([1307.8395](https://arxiv.org/abs/1307.8395),
Table I), so 100 significant figures of a $\sim10^{98}$-sized number is absolute
precision $\approx0.002$ — about $0.08$ of a local gap — with the leading digits
coming from the *prime-free* Lambert-$W$ smooth estimate and the primes refining
only the last few. Their gap-normalized error law (eq. 36),
$$\frac{t_n-t_{n;N}}{\text{gap}}\;\approx\;\frac{\cos(t_n\log p_N)}{\pi\sqrt{\log N}},$$
is **gap-scale and height-independent** (depends on $N$, not $t$) and improves only
as $1/\sqrt{\log N}$ — so any feasible $N$ is pinned at gap-scale. This is an
independent, published derivation of band saturation's "gap-scale,
height-independent in gap units," and it corroborates the floor we had flagged as
ours. What remains genuinely this project's is the narrower
*information-theoretic optimality* of that floor (no estimator on the same primes
beats it), not its existence or scaling.

**What is genuinely this project's — and the risk.** The search did *not* surface
the **strong** form of band saturation: that below $\kappa=\pi$ the sub-gap zero
positions are not merely expensive but **information-theoretically absent** from
primes $\le\sqrt{T/2\pi}$ (the null-space / "no estimator beats B" / $5\times10^{-6}$
control, [`band-saturation.md`](band-saturation.md)). The standard view is that
primes and zeros *mutually* determine each other — true in the **limit** (all
primes ↔ all zeros, exactly); our claim refines that to *finitely many feasible*
primes at *finite precision*, with a sharp cutoff at $\kappa=\pi$. That sharpening,
the integer-snapping-as-error-correction framing, and the super-critical
self-paving loop ($\alpha\approx1$) are this project's. The floor's *existence and
scaling* are now externally corroborated (França–LeClair eq. 36, above); only its
**information-theoretic optimality** — that *no* estimator on the same primes beats
it — remains the piece **not externally validated**.

**Implementation precedents.** The loop's *forward* half — zeros from primes only,
ζ-free — is implemented in França & LeClair
([1601.00914](https://arxiv.org/abs/1601.00914)) and in Balanzario & Cárdenas
Romero, who "compute the zeros … without actually using \[ζ\]"
([2312.00108](https://arxiv.org/abs/2312.00108)); recurrence-formula papers give
both directions and explicitly note the closure "primes → zeros → primes," but as
*formulas with spot-check numerics*, not a running bootstrap
([2009.02640](https://arxiv.org/abs/2009.02640),
[2012.06581](https://arxiv.org/abs/2012.06581)). A literature search did not
surface a public implementation that *runs* the closed loop as a seed-bootstrapped,
ζ-free *and* primality-free self-paver with the criticality ($\alpha\approx1$) and
integer-snapping analysis; that running artifact (`zzz --loop`) appears to be this
project's own — a demonstration-grade novelty, not a mathematical one (with the
caveats that absence of evidence is not proof, and that nobody had practical reason
to build it, since sieving is far cheaper).

**Calibrated status.**

| claim | confidence | basis |
|---|---|---|
| gap-scale zeros → exact prime *set* | **very high** | França–LeClair precedent + verified flawless vs a sieve |
| resolution scale $\sqrt{T/2\pi}$ | **high** | matches textbook $T\sim\sqrt x$ smooth / $T\sim x$ fine |
| loop super-critical, climbs above the wall | **high qualit., medium exponent** | trajectory real; $\alpha\approx1$ from brittle fits + doubles |
| band-saturation *gap-scale floor*, height-independent | **high** | independent derivation: França–LeClair eq. 36 ($\propto1/\sqrt{\log N}$) |
| band-saturation *optimality* (no estimator beats B) at $\kappa=\pi$ | **medium–high; weakest link** | tested 4 ways + rigidity, but ours, untested below the $5\times10^{-6}$ floor |
| the capstone *framing* as a new insight | **deflate** | mostly known facts; the unification is exposition |

So: not *terribly* wrong — the verified facts (coarse zeros → exact primes) are
sound and precedented. But this note is **synthesis plus one empirical claim**
(band saturation at $\kappa=\pi$), not a theorem; if anything here is overturned it
will be the information-absence claim, which the "required by consistency" argument
leans on. The near-trivial reading also deserves airtime: recovering an
integer-valued function from real data *always* needs only finite precision — the
non-trivial residue is the *quantitative* matching (gap-scale, height-independent)
and the sharp $\kappa=\pi$ cutoff, not the bare "exactness lives in $\mathbb{Z}$."

## Companions

[`band-saturation.md`](band-saturation.md) (forward, information-limited),
[`zeros-primes-bootstrap.md`](zeros-primes-bootstrap.md) (the loop, integer-snap),
[`snap-margin-detector.md`](snap-margin-detector.md) (backward tolerance),
[`backward-reach.md`](backward-reach.md) (super-critical reach, $\alpha\approx1$),
[`zero-rigidity.md`](zero-rigidity.md) (sub-gap info absent both ways),
[`spectral-shadow.md`](spectral-shadow.md) (the operator is untouched).
