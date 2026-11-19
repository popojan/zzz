# The zeros ⇄ primes bootstrap: can the explicit formula pave its own way?

*Synthesis of the 2026-06-15 session. The Riemann–von Mangoldt explicit formula
couples primes and zeros in **both** directions; `zzz` normally runs it one way
(primes → zeros). This note asks whether the loop closes — first-$k$ zeros →
more primes → one more zero → … — and demonstrates that it does, for a finite
stretch, with **no zeta evaluation and no primality test**. All scripts and TSVs
are under `doc/ghy/bootstrap-*`; ground truth `doc/refs/odlyzko/` (gitignored).
Companion: [`band-saturation.md`](band-saturation.md), [`mixed-anchor-seeds.md`](mixed-anchor-seeds.md).*

## The two directions

$$\psi(x)=\sum_{p^m\le x}\log p \;=\; x-\sum_\rho\frac{x^\rho}{\rho}-\log2\pi-\tfrac12\log(1-x^{-2}).$$

- **Forward** (primes → zeros): `zzz`'s counter, $F_B(t)=N_0(t)+\tfrac1\pi\arg P_X(\tfrac12+it)$,
  $\log P_X=\sum_{p^m\le X}\tfrac1m p^{-ms}$.
- **Backward** (zeros → primes): under RH,
  $$\psi'(x)\approx 1+\frac{1}{x-x^3}-\frac{4}{\sqrt x}\sum_{k}\cos(\gamma_k\log x),$$
  which **peaks at prime powers**. The peak height obeys an empirical envelope
  (from `2023-10-27_psi_derivative_primepi.nb`)
  $$E(x,n)=\tfrac{20}{4}\,\mathrm{Li}(n-4)\,\frac{\log x}{x}\approx\frac{\gamma_n}{2\pi}\frac{\log x}{x},$$
  i.e. (highest frequency used) × (von Mangoldt weight density). Verified to ~3%
  at $n=330$ (prime peaks track $E$, prime squares peak at $\sim E/2$, composites $\approx0$).

## What was measured

| experiment | script → artifact | result |
|---|---|---|
| **backward reach** $x_\max(n,\varepsilon)$ | `bootstrap-backward.wls` → `.tsv` | $x_\max\sim n^{0.5}$; **noise-immune to $\varepsilon\approx10^{-2}$** (integer-snapping), collapses past $3\times10^{-2}$ |
| **forward accuracy** $\varepsilon(X,n)$ | `bootstrap-loop.wls` → `bootstrap-forward.tsv` | plain-B floor $\sim0.02$–$0.07$, barely shrinks with $X$ |
| **loop closure** | → `bootstrap-closure.tsv` | $\varepsilon_{\rm fwd}$ at $X=x_\max(n)$ is $1$–$3\times$ the backward tolerance → **marginal/critical** |
| **closed loop, B vs Weil-model** | `bootstrap-closeloop.wls` → `.tsv` | plain-B **stalls** $X{\approx}63$–$83$; a $10\times$-sharper (Weil-modelled) forward **tracks the ceiling** to $127$ |
| **self-paving demo** | `bootstrap-selfpave.wls` → `.tsv` | see below |

### The integer-snapping crack
Backward detection is a *thresholded* decision, so it tolerates substantial zero
error: the reach is unchanged up to $\varepsilon\approx10^{-2}$ absolute (B/floor
quality) and only collapses past $\sim3\times10^{-2}$. **Self-computed (noisy)
zeros suffice to read off primes** — this is what keeps the loop alive.

## The self-paving demonstration (pure, audited)

`bootstrap-selfpave.wls`: seed = first 80 zeros (given data). The **loop body
uses only** the $\psi'$ signal on its *own* zeros, $N_0=\theta/\pi+1$ (`LogGamma`,
no zeta), and a partial-Euler sum over its *own* discovered primes. Audited free
of `Zeta`, `ZetaZero`, `Prime`, `PrimeQ`, `PrimePowerQ`. A "new prime" is a
detected prime-power that is not a power of a known prime (integer bookkeeping,
not a primality test). Forward = plain method B; accuracy held under the detector
tolerance by extending zeros only while $\kappa\ge5$.

| iter | $X_{\rm known}$ | #primes | #zeros |
|---|---|---|---|
| 1 | 54 | 16 | 603 |
| 2 | 101 | 26 | 1585 |
| 3 | 149 | 35 | 2846 |
| 4 | **314** | **65** | 4500 (cap) |

**Verification (external, not used by the loop):** zero false positives, zero
missed primes, largest prime 313, median zero error 0.036 — every discovered
prime is correct. $X$ climbs monotonically; the plateaus coincide *exactly* with
hitting the zero-count cap (1500-cap → froze at 137; 4500-cap → froze at 314), so
they are compute-budget artifacts, not stalls. **Super-critical, demonstrated for
a finite stretch.** Plain-B + a $\kappa$ margin sufficed — the Weil fit was not
needed at this scale.

## Sustainability as the heights climb (the asymptotics)

> **Status (2026-11): this section and the next describe the prototype, not the
> current `--loop`.** They rest on the measured reach of the *single-truncation*
> detector, $X\sim n^{0.5}$. Since the truncation-median detector (2026-06-23) the
> loop commits at $X\approx T/2\pi$, i.e. $X\propto n/\log n$: the forward margin
> sits at $\kappa\approx6.3$, `kmin` never binds, and no anneal or ceiling occurs
> (fresh run to $X=3.5\times10^4$). What limits the loop, and what does not, is
> measured in [`loop-margins.md`](loop-margins.md). The analysis below stays as the
> record of why a reach of $n^{0.5}$ *would* stall.

The $\kappa$-margin that makes the demo work does **not** scale. To hold
$\kappa\ge\kappa_{\min}$ at the frontier zero (index $n$) needs
$X\gtrsim n^{\kappa_{\min}/2\pi}$, while the backward step *supplies* only
$X\sim n^{0.5}$:

| margin | $X$ needed | vs supply $n^{0.5}$ |
|---|---|---|
| $\kappa\ge5$ (B accurate) | $n^{0.80}$ | deficit — erodes |
| $\kappa=\pi$ (bare resolution) | $n^{0.50}$ | **exactly marginal** |

So the supply exponent equals the demand exponent *precisely at $\kappa=\pi$*.
Any comfortable margin erodes as the loop climbs, and $\kappa_{\rm frontier}\to\pi^+$.
There: plain-B degrades to its $\sim0.05$ floor $>$ detector tolerance → **B-only
loop stalls** (the closed-loop W-model result); the **Weil fit buys headroom**
(its measured $3$–$15\times$ over B pushes accuracy back under tolerance, extending
the reachable height); but below $\kappa=\pi$ **band-saturation** (E4c: no estimator
beats B) starves the backward step → **doomed**. Nothing crosses $\kappa=\pi$.

### The loop as a contraction map (what "self-accelerate" means)

Composing the two steps gives a map on the prime bound, with margin $\kappa\ge\mathrm{kmin}$:
$$X_{i+1}\sim\big(\underbrace{X_i^{\,2\pi/\mathrm{kmin}}}_{n=F(X_i)}\big)^{0.5}=X_i^{\,\alpha},\qquad \alpha=\frac{\pi}{\mathrm{kmin}}.$$
- $\mathrm{kmin}>\pi\Rightarrow\alpha<1$: **contraction** — $X$ converges to a finite ceiling
  $X^\*(\mathrm{kmin})$, then stalls. The growth *decelerates*; there is no acceleration.
- $\mathrm{kmin}=\pi\Rightarrow\alpha=1$: marginal (the $\sqrt T$ wall).
- $\mathrm{kmin}<\pi\Rightarrow\alpha>1$: *would* accelerate super-exponentially — but this is
  below threshold, where the forward step cannot resolve zeros (premise fails).

So **no fixed-margin run gives sustained acceleration**: every run at fixed $\mathrm{kmin}>\pi$
converges to its ceiling $X^\*(\mathrm{kmin})$, and resuming *at the same $\mathrm{kmin}$* just
restarts at $X^\*$ and stays there. But **feeding the state forward with a smaller margin keeps
climbing**: $X^\*(\mathrm{kmin})\uparrow\infty$ as $\mathrm{kmin}\to\pi^+$, so a *sequence* of runs
that anneals $\mathrm{kmin}$ downward is unbounded in principle. The CLI does this directly —
`zzz --loop --resume --loop-kmin <smaller>` continues from the saved primes/zeros toward the higher
ceiling. The catch is the asymptote: as $\mathrm{kmin}\to\pi$ the accuracy margin vanishes
(plain-B → Weil → nothing at $\pi$), so the annealing terminates at the $\sqrt T$ wall. Thus "not
doomed at finite height" = any target $X$ is reached by some $\mathrm{kmin}>\pi$ in finitely many
steps; "cannot self-accelerate forever" = the climb is gated by $\mathrm{kmin}\to\pi$, decelerating
into the wall rather than running away.

### The C demonstrator (`zzz --loop`)

`loop.c` ports the self-pave loop to self-contained **double precision** (adequate: the loop's
heights are modest and it cannot feasibly climb past the doubles ceiling $\gamma\sim10^{13}$ — the
prime count explodes first; arb is reserved for at-height runs). `./zzz --loop` streams discovered
primes, checkpoints for `--loop --resume`, and catches Ctrl+C. Verified: every streamed prime is
correct (zero false positives), confirming the contraction map empirically — $X^\*(\mathrm{kmin}{=}5)=104$,
$X^\*(\mathrm{kmin}{=}4)>503$ (96 primes). Lowering the margin raises the ceiling, exactly as predicted.

**Auto-annealing (default; dormant since the median detector, see the status note above).** Rather than make the user re-launch with a smaller margin at each
ceiling, `--loop` automatically lowers $\mathrm{kmin}$ (×0.9) whenever it stalls, down to a floor
(`--loop-kmin-floor`, default 3.5, kept above $\pi$ so plain-B stays correct) — so one invocation
climbs through successive ceilings on its own and stops at the plain-B wall with an honest message
("go further with a sharper/Weil forward step"). The conservative detector emits no false primes at
any margin, so annealing is safe. `--loop-no-anneal` keeps $\mathrm{kmin}$ fixed (stalls at one
ceiling, for the contraction-map analysis); the checkpoint carries the annealed margin.
This realises the "feed forward" loop directly: the *sequence* is unbounded in principle, gated by
the floor at $\pi$ (the $\sqrt T$ wall), not by manual orchestration.

**Out-of-the-box UX.** `./zzz --loop` **auto-resumes** if a checkpoint exists at the state path,
else starts fresh — so: run it, Ctrl+C anytime (it checkpoints), run `--loop` again to continue.
`--loop-batch N` (default 2000) computes only N zeros per re-detect, so primes **stream smoothly**
rather than in big silent bursts (the detect is frontier-incremental, so frequent re-detection is
cheap). `--loop-fresh` forces a fresh start; `--resume` is the explicit form (errors if no
checkpoint). Verified end-to-end: interrupt at X=1293, relaunch, auto-resume to X=1643 — 259 primes,
all correct.

**Two peak scales (A/B).** Same detector, same thresholds; only the normalisation differs:

| scale | fitted constant | min seed |
|---|---|---|
| envelope (default) | Li form + 20/4 (matched-filter peak scale) | 6 zeros, 0 primes |
| `--loop-contrast` | **none**: Landau's main term $\tfrac{2}{\pi}\gamma_n\,\tfrac{\log x}{x}$ from the loop's own top zero | 6 zeros, 0 primes |

Both give exactly the primes from every seed size 6…40 (`ctest`, `smoke.loop`). The fit-free scale
shows the loop needs **neither Li nor the 20/4** — $5\,\mathrm{Li}(n-4)\approx4\gamma_n/2\pi$ to a few
percent at these heights (the Li form drifts to $5:4$ asymptotically; the fit-free one does not).
*History:* until 2026-11 `--loop-contrast` was a local-MAD (CFAR) detector. Its threshold was
relative to the composite noise rather than to the prime peak, so its false-alarm rate per candidate
was fixed, and it emitted $437=19\cdot23$ from a fresh start; it was replaced, not tuned. Soundness
of the envelope rests not on its tuned constant but on two *provable* facts: composites have $\psi'\approx0$ structurally, and the truncation turns each
prime-power delta into a Dirichlet-kernel peak of height $\propto\gamma_n\log x/x$ (and
$\mathrm{Li}(n)\approx\gamma_n/2\pi$). The naive textbook detector $\int\psi'=\Lambda$ *fails* under
truncation (the uncancelled $+1$ baseline), so the matched-filter normalisation is doing real work.

**Minimum seed.** The pure "1 zero + {2,3}" is below the floor: one zero is a single cosine and
locates no prime; the explicit-formula sum needs a *handful* (~6) to resolve even x=2,3. A small
seed must also not be over-read: six zeros ($\gamma\le37.6$) resolve only $x\lesssim15$, and until
2026-11 the small-$x$ branch had no such limit (`--loop-seed 6` emitted 18, 22, 24, …). The lag
floor $\gamma(0.5\,N_Z)\ge2.5\,x$ now applies at every $x$. So the irreducible seed is "~6 zeros",
with the peak scale read off the zeros themselves.

**Verdict (prototype, 2026-06).** Informationally consistent and **asymptotically critical**: the loop
self-paves and the prime bound grows for as long as the zero supply grows, but it
converges to the band-saturation / $\sqrt T$ wall rather than outrunning it. Not
doomed at any finite height (Weil extends the affordable range); not a free
explosive ladder. Integer-snapping plus a $\kappa$ margin put plain-B on the
growing side for the demonstrated stretch.

**Verdict (current loop, 2026-11).** With the median detector the loop does not
approach that wall at all: it holds $X\approx T/2\pi$ against a requirement of
$X\gtrsim(T/2\pi)^{0.8}$, the per-decision margin grows like $\sqrt{\log x}$, and
nothing measured points to a failure height. It is not a theorem, and the cost
grows like $X^2\log^2X$. Evidence: [`loop-margins.md`](loop-margins.md).

## Artifact index

Scripts/TSVs under `doc/ghy/`: `bootstrap-backward.{wls,tsv}` (reach),
`bootstrap-loop.{wls,tsv}`/`bootstrap-{forward,closure}.tsv` (forward + closure),
`bootstrap-closeloop.{wls,tsv}` (B vs Weil-model), `bootstrap-selfpave.{wls,tsv}`
(the pure demo). Origin: `~/Documents/Wolfram/2023-10-27_psi_derivative_primepi.nb`
($\psi'$ kernel `dPsiCosLRH`, envelope). C port: `loop.c` (`zzz --loop`).
Margins of the current loop: `loop-margins.py`, `loop-replay.py`, `loop-fixed-cutoff.{c,py}`,
`loop-refzeros.c` + `loop-zero-accuracy.py` (see [`loop-margins.md`](loop-margins.md)).
