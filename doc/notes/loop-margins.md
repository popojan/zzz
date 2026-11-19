# What limits `--loop`: measured margins, the one condition, and what was fixed

*Review of 2026-11-19, done from the code rather than from the earlier notes: every
backward decision of a fresh run was logged (`--loop-declog`), the loop's zeros were
compared with FLINT reference zeros, and the decision rule was replayed on real primes at
heights the loop cannot reach. Scripts: `doc/ghy/loop-margins.py`, `loop-replay.py`,
`loop-fixed-cutoff.{c,py}`, `loop-refzeros.c` + `loop-zero-accuracy.py`. This note
supersedes the "Sustainability" and "contraction map" sections of
[`zeros-primes-bootstrap.md`](zeros-primes-bootstrap.md), which describe the prototype
detector. Companions: [`box-forward-step.md`](box-forward-step.md),
[`band-saturation.md`](band-saturation.md), [`exact-primes-inexact-zeros.md`](exact-primes-inexact-zeros.md).*

## The question and the answer

*Is the loop predestined to work forever, or must it fail at some height?* Neither.
No mechanism was found that forces a failure at any height, and nothing proves there
is none. What can be measured is that the margin on each prime/composite decision
**grows** with height, fast enough that the probability of ever emitting a wrong prime
is summable and small. What ends a run in practice is cost.

## Where the loop actually runs

Write $a=\gamma_{\rm top}/(2\pi x)$ for the **lag** at which a candidate $x$ is committed
($\gamma_{\rm top}$ = highest zero in hand).

| quantity | measured | note |
|---|---|---|
| lag at commit | $a\approx0.85$–$1.2$, median $0.91$ | so $X\approx T/2\pi$, i.e. $X\propto n/\log n$, not $n^{0.5}$ |
| forward margin | $\kappa=2\pi\log X/\log(T/2\pi)\approx6.3$ | `kmin` (5.0) never binds; no anneal, no ceiling |
| zero error, box-B (default) | $0.025$ mean gap rms | vs FLINT, at $T=2.6\times10^4$ and $1.4\times10^5$ |
| zero error, plain B | $0.052\to0.040$ gap rms | index $3\times10^3\to2.6\times10^6$: slowly *improving* |
| new-prime peak, box-B / plain B | $\approx0.85$ / $\approx0.68$ of the envelope | a prime already in the forward sum reads $\approx1.0$–$1.1$ |

The forward step needs $X\gtrsim(T/2\pi)^{\kappa/2\pi}=(T/2\pi)^{0.8}$ for $\kappa\ge5$ and
is given $X\approx T/2\pi$. That slack is a growing power of $T$, so the forward step
never becomes the constraint and the loop does not approach the band-saturation curve
$X=\sqrt{T/2\pi}$.

## The margin on each decision

`loop-margins.py` on the decision log of a fresh default run (thresholds: median
$>0.40$ commits a prime, $<0.25$ a composite):

| $x$ range | lag $a$ | primes: mean ± sd (min) | composites: mean ± sd (max) | distance to 0.40 | var·$a$·$\ln x$ | $\alpha$ |
|---|---|---|---|---|---|---|
| 600–1372 | 1.20 | 0.869 ± 0.038 (0.776) | −0.047 ± 0.067 (0.136) | 6.7σ | 0.037 | 3.3 |
| 1372–3137 | 0.97 | 0.867 ± 0.040 (0.769) | −0.073 ± 0.075 (0.135) | 6.3σ | 0.042 | 2.6 |
| 3137–7173 | 0.95 | 0.856 ± 0.038 (0.742) | −0.067 ± 0.074 (0.150) | 6.3σ | 0.044 | 2.3 |
| 7173–16400 | 0.91 | 0.835 ± 0.036 (0.723) | −0.060 ± 0.067 (0.156) | 6.8σ | 0.039 | 2.5 |
| 16400–37500 | 0.91 | 0.820 ± 0.031 (0.714) | −0.054 ± 0.063 (0.156) | 7.2σ | 0.037 | 2.6 |

37 498 decisions, 3968 primes, none wrong, none missed; 158 scan stops, all "out of
reach", none ambiguous; lag at commit between 0.846 and 1.57 (median 0.91). Composite tail,
standardised per bin: $P(z>3)=1.5\times10^{-4}$ against a Gaussian $1.3\times10^{-3}$,
nothing beyond $3.5\sigma$ in 33 000.

The same 25-minute run in the other variants, last bin ($x\approx1.7$–$4\times10^4$):

| variant | reached | primes: mean ± sd (min) | composites: sd (max) | distance to 0.40 | $\alpha$ |
|---|---|---|---|---|---|
| default (box-B, Li scale) | 37 499 | 0.820 ± 0.031 (0.714) | 0.063 (0.156) | 7.2σ | 2.6 |
| `--loop-contrast` (fit-free scale) | 36 763 | 0.831 ± 0.031 (0.726) | 0.064 (0.159) | 7.1σ | 2.5 |
| `--loop-boxc 0` (plain B) | 41 617 | 0.692 ± 0.052 (0.522) | 0.071 (0.235) | 6.2σ | 1.8 |
| `--loop-seed 6` | 37 498 | (same detector as default) | | | |

All four lists are exactly the primes. Plain B has the thinner margin throughout (and two
ambiguous stops); the older 2.7M-zero plain-B checkpoint reads $0.68\pm0.047$ and $6.8\sigma$
at $x\approx2.5\times10^5$.

- **Scaling.** The composite variance falls like $1/(a\log x)$ (`var*a*lnx` is flat or
  falling), because each candidate is read with $\sim a\log x$ zeros' worth of signal.
  With a Gaussian tail the false-positive probability per candidate is $x^{-\alpha}$,
  $\alpha=z^2/(2\log x)\approx2.3$–$2.6$ (plain B: $\approx1.8$–$2.1$). Since $\alpha>1$ the sum over all $x>X$ converges
  and is tiny ($\sim X^{1-\alpha}$). The observed tail is *lighter* than Gaussian.
- **Composition of the noise** (`loop-replay.py` model, correlation $\approx0.8$ with the
  logged readings): mostly the truncation sidelobes of neighbouring primes, which the
  median over six truncations largely cancels, plus a residual of variance
  $\approx0.021/(a_b\log x)$ per truncation from the zero errors.
- **Misses** are further away still: a prime would have to read below $0.25$, about
  $15\sigma$.
- **Slow drift.** The prime peak in the default (Li) normalisation declines by a few
  percent per decade of $x$ (0.87 → 0.82 above; 0.86 → 0.83 with the fit-free scale), because $5\,\mathrm{Li}(n-4)$ creeps from $\approx4\gamma_n/2\pi$
  toward $5\gamma_n/2\pi$; the asymptote is still $\approx0.6$–$0.7$, above $0.40$, and
  crossing it would stall the loop, not corrupt it. `--loop-contrast` (fit-free scale
  $\tfrac2\pi\gamma_n$) has no such drift.

## The one condition: the lag floor

Replaying the rule on the real primes around $x_0$ (`loop-replay.py`), as a function of
the lag:

| $x_0$ | false primes at $a=0.5$ (per $10^6$) | at $a=0.55$ | at $a\ge0.6$ | gate still passes at $a=0.5$ |
|---|---|---|---|---|
| $2.5\times10^5$ | 2951 | 143 | 0 | 80% |
| $10^7$ | 2688 | 90 | 0 | 90% |
| $10^9$ | 2071 | 28 | 0 | 96% |
| $10^{12}$ | 1364 | 9 | 0 | 99% |
| $10^{15}$ | 844 | 4 | 0 | 99.7% |

Every false prime is a composite **between twin primes**: once the lowest truncation no
longer resolves unit spacing ($\gamma(0.5N_Z)\lesssim1.8\,x$) it rides both main lobes.
The floor $a^\*\approx0.55$ does not move with height; the loop commits at $a\approx0.9$,
where a twin-sandwiched composite reads *negative*.

What held the lag at 0.9 was the statistical reach gate (70% of the 61 neighbours clear
at the lowest truncation). The last column shows that this gate thins out as primes do:
at $10^{12}$ it passes 99% of windows at a lag where the rule is wrong. It still worked,
because the frontier is contiguous and is held back by the rare densest windows, which
never run out — but that is an indirect guarantee. `loop.c` now states the condition
itself, $\gamma(0.5N_Z)\ge U_{\min}\,x$ with $U_{\min}=2.5$ ($a\gtrsim0.75$), at every $x$.
In the logged runs the gate had never let the ratio under 2.8, so the floor changes no
output; it removes the dependence on prime-cluster statistics.

## What the forward step transmits

Zeros located from the primes $\le X$ are a function of those primes only, so whatever
they reveal about $x>X$ is implied by them. `loop-fixed-cutoff` ($X=100$, zeros to
$T=2\times10^5$), coherence by the margin $\kappa$ the zeros were located at
(1 = reads as a prime, 0 = as a composite):

| integers $x$ | $\kappa>5$ | $\kappa\in(4,5)$ | $\kappa\in(\pi,4)$ | $\kappa<\pi$ |
|---|---|---|---|---|
| primes in $(20X,X^2]$ | 0.5–0.8 | 0.9–1.0 | 0.8 | 0.6–0.7 |
| composites with a factor $\le X$ | $0\pm0.2$ | $0\pm0.2$ | $\approx0$ | $\approx0$ |
| $p\,q$, both $>X$ (true zeros: 0) | $\approx0$ | 0.2 | 0.55–0.7 | **0.7** |
| $p^2$, $p>X$ (true zeros: 0.5) | $\approx0$–0.4 | 0.2–0.7 | 0.7–0.8 | **0.8** |

At high margin the zeros are close enough to the true ones to inherit the true
$\Lambda$. As $\kappa\to\pi$ they transmit only *"$x$ has no prime factor $\le X$"*: the
sieve of Eratosthenes, carried out analytically. That identifies primes exactly while
$x<X^2$, and $x<X^2$ at the reach $x\approx T/2\pi$ is $X>\sqrt{T/2\pi}$ — the
band-saturation wall in sieve language. The loop works at $x<2X$, so "no known factor"
and "prime" coincide for all but a vanishing fraction ($\sim x^{-1/2}$) of its zeros.

## Errors do not spread, and are not noticed

Replacing 10% of the primes in the upper half of a checkpoint by composites and resuming
gave exactly the same new primes as the clean checkpoint (the top zero moved in its 7th
digit). One wrong prime near $X$ shifts each later zero by $O(X^{-1/2})$ of a gap. So a
mistake would neither avalanche nor heal: a long clean run is evidence about the margins
above, not a self-check.

## What was fixed (2026-11)

| problem | fix |
|---|---|
| `--loop-contrast` emitted $437=19\cdot23$ (threshold relative to composite noise: fixed false-alarm rate per candidate) | replaced by a fit-free peak scale $\tfrac2\pi\gamma_n\tfrac{\log x}{x}$ in the shared detector; no warm-start primes |
| `--loop-seed 6` emitted 18, 22, 24, … (the $x<512$ single-reading branch had no reach limit) | one rule for every $x$: truncation median behind the lag floor |
| the safety condition was only implicit | $U_{\min}=2.5$ lag floor |
| `locate_march` returned `prev + gap` when the true gap exceeded 2.7 mean gaps (≈2 zeros in $10^4$; one caught at index 600 051, off by 1.8 gaps) | the bracket widens to 4, 8, 16 gaps first |
| a checkpoint did not say which forward variant made its zeros (the 2.7M-zero one in circulation is plain B) | `boxc` line in the checkpoint; resume keeps it unless `--loop-boxc` is given |
| backward summed all ~1100 candidates up to `Xcap` over every zero, then committed a median of 56 | rows are computed a chunk ahead of the scan |
| forward: `pow`/`log` per term, and 50 bisection steps where ~35 exhaust a double | cached $\log p$, $p^{-1/2}$; stop when the bracket is two adjacent doubles |

Both speed-ups leave checkpoints byte-identical. Per iteration at $X=3.5\times10^4$
($3\times10^5$ zeros, one core): backward 7.2 s → 0.8 s, forward 22.3 s → 9.5 s. Note
the forward step is the larger part at these heights; the backward part only overtakes
it near $\log X\approx16$. `smoke.loop` in `ctest` checks that the default, plain-B,
fit-free and 6-zero-seed variants emit exactly the primes.

Not done, because they change low-order bits and would need re-validation: a
derivative-based root step in place of bisection (~5 evaluations instead of ~35), and an
incremental prefix-sum cache for the backward step (each pending candidate only needs
the newly added zeros, not all of them).

## Limits of this evidence

- The $x^{-\alpha}$ law extrapolates a tail observed to $\approx3\sigma$ out to
  $6$–$7\sigma$; the replay beyond $X\approx2.5\times10^5$ is a model.
- Scale invariance of the zero error and of the new-prime peak is measured over
  $X=300$–$10^5$ (fixed geometry) and to index $2.6\times10^6$ (plain B), and is what
  the local-statistics picture predicts; it is not proved.
- A provable loop would need rigorous error bars on the zeros (a Turing-type check
  instead of method B), a smooth kernel with bounded sidelobes, and interval
  arithmetic. `--loop` is none of these.
- Cost grows like $X^2\log^2X$: estimated at weeks of CPU for $X=10^6$ and years for
  $10^7$; doubles are adequate to $T\sim10^{13}$, far beyond that.

## Reproduce

```bash
./zzz --loop --loop-fresh --loop-declog run.dec > run.primes      # Ctrl+C after a while
python3 doc/ghy/loop-margins.py run.dec                           # the margins table
python3 doc/ghy/check-loop.py zzz-loop.state                      # primes vs a sieve
for x in 2.5e5 1e7 1e9 1e12 1e15; do python3 doc/ghy/loop-replay.py $x; done
gcc -O3 -march=native -ffast-math -fopenmp -o loop-fixed-cutoff doc/ghy/loop-fixed-cutoff.c -lm
./loop-fixed-cutoff 100 120 200000 0.375 z100.bin && python3 doc/ghy/loop-fixed-cutoff.py 100 z100.bin
gcc -O2 -o loop-refzeros doc/ghy/loop-refzeros.c -lflint -lm
./loop-refzeros 30000 2500 > ref.tsv && python3 doc/ghy/loop-zero-accuracy.py zzz-loop.state ref.tsv
```
