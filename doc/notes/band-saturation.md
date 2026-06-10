# Band saturation: what the primes up to X know about the zeta zeros

*Synthesis of the 2026-06-10 sessions. This is the explanatory record; the
chronological lab logs are [`bootstrap-hybrid.md`](bootstrap-hybrid.md)
(the `--boot` work) and [`superres-multi-cutoff.md`](superres-multi-cutoff.md)
(the sub-Rayleigh program E1–E4, kept as written, including claims later
retracted in place). All numbers below are measured; scripts and data live
in `doc/ghy/`, ground truth in `doc/refs/odlyzko/` (gitignored; download
commands in the script headers).*

## The result in one paragraph

For a zero of ordinate $T$, located from primes $p^m \le X$ only, define the
resolution parameter $\kappa = \mathrm{gap}\cdot\log X$, where
$\mathrm{gap} = 2\pi/\log(T/2\pi)$ is the local mean zero spacing. Then:

- **$\kappa \gtrsim \pi$** (equivalently $X \gtrsim \sqrt{T/2\pi}$, the
  Riemann–Siegel scale): the primes carry surplus information about
  individual zeros beyond what zzz's counting-function bisection (method B)
  extracts, and the self-consistent hybrid bootstrap (`zzz --boot`)
  collects much of it — measured gains 19.3× at $n\sim10^3$, 2.8× at
  $n\sim10^5$, 2.0× at $n\sim10^{12}$.
- **$\kappa < \pi$**: method B already **saturates** the information
  content of the primes about individual zeros. Four independent estimator
  families — crossing relocation, kernel maximum likelihood, multi-cutoff
  ML, and exact Weil-identity ML — all reproduce B's errors; the final
  experiment exhibits two zero configurations 0.017 apart (truth, and B's
  output) whose band-limited prime data agree to $5\times10^{-6}$.
  Below the threshold, B's error *is* the null space of the data.

The Riemann–Siegel wall $X \sim \sqrt{T/2\pi}$ is therefore
information-theoretic, not algorithmic: no estimator, weighting, damping,
prime-shifting or prime-selection scheme applied to the same primes can
cross it.

## How it was established

The day started from the question: *can convergence of the
zeta-evaluation-free counters be accelerated at extreme heights — e.g. by
shifting the truncated primes to compensate for the missing tail?* The
chain of experiments, each closing a loophole left by the previous one:

| # | experiment | script(s) | regime | result |
|---|---|---|---|---|
| 0 | cutoff sweep: B's location error vs $\log X$ | `cutoff-sweep.sh`, `cutoff-osc-analysis.wls` | $n=10^3..10^5$, $10^{12}$ | error follows the parameter-free neighbour model $\Delta(L) = -\tfrac1L\sum_j[\mathrm{Si}(\delta_j L)-\tfrac\pi2\mathrm{sgn}\,\delta_j]$, corr 0.96–0.98 |
| 1 | self-consistent hybrid bootstrap `--boot` | `main.c::boot_locate`, `boot-validate.sh/.wls` | 5 heights | gains 19.3× / 2.8× / 2.0× / 0.95× / 0.93×; switch-off at $\kappa\approx\pi$ |
| 2 | $10^{36}$ proxy regime ($10^{22}$, k=600, $\kappa\approx1.1$) | `boot-regime-e22-k600.tsv` | ensemble of 20 | no recipe beats B; multi-round iteration diverges |
| E1 | **oracle control**: true Odlyzko seeds in the leave-one-out relocation (`zzz --seeds`) | `e4`-era runs, `boot-validate` baselines | $10^{21}/10^{22}$ | gains **5.9× / 3.5× / 2.7×** where self-seeding gives nothing — the identity machinery is healthy; only self-obtained knowledge is blind |
| E2 | single-cutoff full-profile ML (Si-kernel model) | `e2-singleL-fit.wls` | proxy regime | gain 1.01× — degeneracy is in the likelihood, not the estimator |
| E3 | multi-cutoff ML | `e3-multiL-fit.wls` | proxy regime | gain 0.54× — prefix cutoffs add no information, worse geometry |
| E4a | Weil explicit formula as exact data source | `e4-weil-identity.wls` | low height | identity validated to $10^{-34}$ (sinc$^{2q}$/B-spline pairs, prime side cuts at $p^m\le X$ exactly) |
| E4b | ideal-information probe (Paley–Wiener Gram spectra) | `e4-spectrum-probe.wls` + `.log` | $10^{22}$ configs | eigen-plunge ≈ 1.7 $\log_{10}\!\lambda$/rank past Shannon; *oracle-conditional* GO (retracted in step 3) |
| E4c | **the decisive fit**: exact Weil functionals, 178 unknowns, GN to the data floor | `e4-fit.wls` + `.log` | $10^{22}$, k=10⁴, $\kappa=1.55$ | gain **0.994×**, fit ≡ B (corr 0.9993); control below |

The E4c control, the sharpest single fact of the day
($\|r\|$ = residual of the configuration against exact Weil-functional data):

| configuration | $\|r\|$ |
|---|---|
| B positions as bisected | 0.53 |
| fitted positions (errors ≈ B's, 0.0166 mean) | $3.3\times10^{-6}$ (floor) |
| **true positions** (Odlyzko) | $5.0\times10^{-6}$ (floor) |

The optimizer, starting from B, found a configuration matching the exact
band-limited data $10^5\times$ better than B's raw reading — and that
configuration has B's errors. Truth matches the data no better. The
"missing" information about individual zeros is not attenuated in the
band; relative to what primes $\le X$ provide, it is absent.

## Why, mechanistically

Method B's counting function is the zero-counting measure observed through
a kernel of width $\pi/\log X$ (frequencies of primes $p^m \le X$ span
$[0, \log X]$). Its location error is the truncation deficit — itself a
band-limited field — *aliased onto the zero lattice*. When
$\kappa = \mathrm{gap}\cdot\log X < \pi$ the lattice samples the band above
Nyquist density: every band-consistent reading of the data is a valid
configuration, B's among them, truth among them, all indistinguishable.
Above $\pi$, the kernel resolves individual neighbours; crossing-reading
then *wastes* information (one number per zero), which is the surplus the
bootstrap's leave-one-out $Z_X$ subtraction recovers.

Three plausible escapes were each closed by measurement:

1. *Better estimators on the same field* (E2): the single-cutoff likelihood
   has the same null space as the crossings.
2. *More data at other cutoffs* (E3): prefix sums are analytically
   determined by the full-cutoff field (band-limited ⇒ entire of
   exponential type); no new information exists.
3. *Exact identities + free functional design* (E4): the Weil formula
   removes all model error and gives access to every admissible
   band-limited functional; the data floor is then reached by wrong and
   right configurations alike. The E4b "GO" verdict assumed the far zone
   known — the evanescent channels' functionals carry their mass *outside*
   the fit window, so using them presupposes the very knowledge being
   sought. Self-contained, they are circular; this is the same coherence
   that defeats the bootstrap, in its final form.

E1 calibrates exactly what *external* information would buy: with true
neighbour positions supplied, leave-one-out relocation gains 2.7–5.9×
below the threshold. That is the value of an oracle, not of the primes.

## Consequences

**Feasibility, not impossibility.** The wall moves with $X$: any height is
reachable in principle. At the README's showcase zero
($\#10^{36}$, $T \approx 8.1\times10^{34}$) the threshold sits at
$X \ge \sqrt{T/2\pi} \approx 1.1\times10^{17}$, i.e. $\pi(X) \approx
3\times10^{15}$ primes per counting-function evaluation — the same effort
class as the Riemann–Siegel main sum ($\sqrt{T/2\pi} \approx 10^{17}$
terms; the prime route needs $\sim\!\ln X \approx 40\times$ fewer terms,
a constant-factor consolation). Zeta-evaluation-free methods cannot
asymptotically undercut direct evaluation: resolution arrives exactly when
the budgets meet.

**Prime selection.** Taking the smallest $k$ primes (plus their powers) is
essentially optimal:

- the explicit formula sums over *all* prime powers, so a usable test
  function must annihilate every frequency $m\log p$ *not* supplied —
  only a contiguous known band $[0, \log X]$ survives this;
- sparse large primes alone are useless: nulling their exponentially dense
  un-included neighbours forces test functions spread over $e^{+u}$-wide
  $T$-windows;
- trading mid-band primes for top-band ones gains band $\Delta a \sim a/e^a$
  per sacrificed degree of freedom — exponentially losing;
- within the band, "shifting" or re-weighting primes is a choice of
  in-band functional, a freedom E4 already exhausted. This closes the
  session's original idea definitively: missing larger primes are
  out-of-band spectrum, and no in-band manipulation synthesizes it.

**Practitioner's map** (what to run in zzz):

| regime | best method | measured accuracy |
|---|---|---|
| $\kappa \gtrsim 2\pi$ ($X \gtrsim T/2\pi$) | `--boot 32` | ~0.0008 at $n\sim10^3$, k=1000 — beyond any feasible plain-B budget |
| $\pi \lesssim \kappa \lesssim 2\pi$ | `--boot 32` | 2–3× over `--ghy` |
| $\kappa < \pi$ (e.g. $10^{36}$ at any feasible k) | plain `--ghy` | B is information-optimal; `--boot` is inert in expectation and multi-round can diverge |

## What could still change the picture (precisely bounded)

- **Sub-floor rigidity.** Today's indistinguishability is bounded at the
  $5\times10^{-6}$ data floor (far zone held at B seeds, machine
  arithmetic). The zero measure is a *unit-mass point process*, not a
  continuum; quantization/positivity rigidity could in principle break the
  degeneracy at some far smaller scale. Five independent negatives weigh
  against it; nothing today tests below the floor.
- **Above-threshold constants.** The exact-identity fit machinery has only
  been run below threshold (where it ties B by necessity). Above
  threshold it should beat the *bootstrap's* constants — it wastes no
  information on crossing-reading. Untested.
- **An oracle.** Any external source of approximate neighbour ordinates
  (e.g. a partial zero database) immediately buys the measured E1 gains.

## Artifact index

Code (branch `ghy-hybrid`): `--boot/-R/-S` in `main.c` (`boot_locate`,
`bisect_zero`, `load_seed_ordinates`), `ghy_log_zx_rel` and the
critical-line $E_1$ fix ($E_1(iy) = -\mathrm{Ci}|y| + i(\mathrm{Si}|y| -
\tfrac\pi2)$, conjugated for $y<0$ — the complex power series loses
$\sim|y|/\log 2$ bits and silently widens arb balls) in `ghy.c`.

Commits: `bab39e32` (bootstrap), `cf950404` (validation + validity law),
`7c179d1c` (`--seeds`, E1–E3), `c885786f` (E4a identity harness),
`8864e5c1` (E4b probe), `15c20b6c` (E4c fit + saturation).

Scripts and reference data, all under `doc/ghy/`: `cutoff-sweep.{sh,tsv}`,
`cutoff-osc-analysis.wls`, `boot-validate.{sh,wls,tsv}`,
`boot-regime-e22-k600.tsv`, `boot-regime-analysis.wls`,
`e2-singleL-fit.wls`, `e3-multiL-fit.wls`, `e4-weil-identity.wls`,
`e4-spectrum-probe.{wls,log}`, `e4-fit.{wls,log}`.

Reproducibility pitfalls worth knowing: absolute zero ordinates at
$10^{21+}$ do not fit doubles (resolution $2^{-52}\gamma \sim 10^{5}$) —
work window-relative everywhere; `zzz` output and Odlyzko tables mix
absolute values and per-table offsets (`zeros3/4/5` headers state the
bases); ridge-regularized normal equations with $\sigma_d/\sigma_p \sim
10^{-4}$ against a rank-deficient $J^\top J$ need ≥50-digit solves.
