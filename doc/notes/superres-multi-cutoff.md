# Sub-Rayleigh zero recovery: full-profile fits and multi-cutoff deconvolution

*Plan + formula reference, 2026-06-10. Follow-up to
[`bootstrap-hybrid.md`](bootstrap-hybrid.md), which established the validity
threshold $X \gtrsim \sqrt{T/2\pi}$ for crossing-based methods. This note
specifies the experiments that attack the regime **below** that threshold,
with every formula the code implements. Status: planned; results to be
appended.*

## 1. The problem, restated as information theory

The observable on the critical line is band-limited: primes $p^m \le X'$
contribute oscillations $e^{-imt\log p}$, i.e. frequency content in
$[0, L']$, $L' = \log X'$. Observed on a $T$-window of width $\Delta_T$, such
a signal carries about

$$\mathrm{DOF} \;\approx\; \frac{L'\,\Delta_T}{\pi}
\qquad \text{(Shannon number)}$$

well-conditioned degrees of freedom, while the window contains
$\Delta_T \cdot \frac{1}{2\pi}\log\frac{T}{2\pi}$ zeros. Unknowns exceed DOF
exactly when $L' < \tfrac12 \log(T/2\pi)$, i.e. $X' < \sqrt{T/2\pi}$ — the
empirical validity law of the bootstrap is the Shannon-number criterion.

**But** sub-Rayleigh information is attenuated, not destroyed: zero positions
enter the band analytically, and our data is *exact* — prime sums carry no
measurement noise, only the working precision (hundreds of bits) and the
explicit-formula identity error ($\sim \sqrt X/(T\log X)$, about $10^{-30}$
in our regimes). Parametric point-source recovery from exact band-limited
data (Prony 1795; matrix pencil; ESPRIT; Candès–Fernández-Granda) has **no
Rayleigh limit**; its price is conditioning $\sim \rho^{-(m-1)}$ ($\rho$ =
separation in Rayleigh cells, $m$ = sources per cell), payable in precision
bits. Crossing-based estimators (method B, `--boot`) are linear,
single-functional readers of the field and *are* Rayleigh-limited — that is
what dies at $X \approx \sqrt{T/2\pi}$, not the information.

## 2. Notation and data field

$$s = \tfrac12 + iT, \qquad X' = p_{k'}, \qquad L' = \log X',
\qquad \log P_{X'}(s) = \sum_{p^m \le X'} \frac{1}{m\,p^{ms}}.$$

The data field (what `zero_count_ghy` computes, as a function of both
arguments) is

$$\boxed{\;\mathcal F(T, L') \;=\; N_0(T) + \frac{1}{\pi}\arg P_{X'}\!\big(\tfrac12+iT\big)\;}
\qquad N_0(T) = \frac{T}{2\pi}\log\frac{T}{2\pi e} + \frac78 .$$

**Multi-cutoff data is free**: one pass over the primes yields
$\mathcal F(T, L')$ for *every* prefix cutoff simultaneously, via cumulative
sums of the per-prime terms (including $p^m \le X'$ power terms). For a
$T$-grid $\{T_i\}$ and cutoff grid $\{L'_l\}$ the cost is one
$|\{T_i\}| \times k$ phase matrix.

## 3. Forward model

From the explicit-formula decomposition (GHY) with
$\operatorname{Im} E_1(iyL') = \mathrm{Si}(yL') - \tfrac{\pi}{2}\mathrm{sign}(y)$,
each zero's sharp step $H(T-\gamma_j)$ (in $N$) and its smearing term combine
into one smooth profile:

$$H(y) - \tfrac12\mathrm{sign}(y) \equiv \tfrac12
\;\;\Longrightarrow\;\;
H(T-\gamma_j) + \frac{1}{\pi}\Big[\mathrm{Si}\big((T-\gamma_j)L'\big) - \frac{\pi}{2}\mathrm{sign}(T-\gamma_j)\Big]
= \frac12 + \frac{1}{\pi}\mathrm{Si}\big((T-\gamma_j)L'\big).$$

With $M$ consecutive window zeros $\gamma_1 < \cdots < \gamma_M$ of ordinals
$n_0, \ldots, n_0 + M - 1$ (zeros below the window contribute their full
unit steps, zeros above contribute $\approx 0$, both up to tails):

$$\boxed{\;\mathcal F(T, L') \;=\; (n_0 - 1) \;+\; \sum_{j=1}^{M}\Big[\frac12 + \frac{1}{\pi}\,\mathrm{Si}\big((T-\gamma_j)L'\big)\Big] \;+\; R(T, L')\;}$$

The residual field $R$ collects (i) far-zero tails and (ii) the
explicit-formula identity error:

$$R(T,L') \;=\; \frac{1}{\pi}\sum_{j \notin \text{window}}\Big[\mathrm{Si}\big((T-\gamma_j)L'\big) - \frac{\pi}{2}\mathrm{sign}(T-\gamma_j)\Big] + \varepsilon_{\rm id},
\qquad \Big|\,\text{tail term}\,\Big| \approx \frac{|\cos((T-\gamma_j)L')|}{\pi\,|T-\gamma_j|\,L'} .$$

Key structural fact: a far zero's tail oscillates in $T$ at frequency
$\approx L'$ — it lives at the **band edge**. Its natural nuisance basis on
the window is therefore

$$R(T,L') \;\approx\; \alpha(L')\cos(TL') + \beta(L')\sin(TL') + \text{low-order smooth in } T,$$

with $\alpha,\beta$ slowly varying. In practice we suppress $R$ by (a) a
buffer of real fitted zeros beyond the core, (b) a smooth $T$-window, and
(c) optionally the band-edge nuisance pair above.

## 4. Estimators

### 4a. Single-$L$ full-profile weighted least squares (experiment E2)

Replace "read one crossing per zero" by fitting the whole window profile.
With $\theta = (\gamma_1,\ldots,\gamma_M;\, c)$, $c$ = linear nuisance
coefficients (constant offset, optional band-edge pair / low-order poly):

$$\chi^2(\theta) = \sum_i w_i\,\big[\mathcal F(T_i, L) - \mathcal M(T_i, L; \theta)\big]^2,
\qquad w_i = \exp\!\Big(-\frac{(T_i - T_0)^2}{2\sigma_T^2}\Big),$$

$$\mathcal M(T, L; \theta) = (n_0-1) + \sum_{j=1}^{M}\Big[\frac12 + \frac{1}{\pi}\mathrm{Si}\big((T-\gamma_j)L\big)\Big] + c \cdot \mathrm{basis}(T).$$

Gauss–Newton / Levenberg–Marquardt with the analytic Jacobian
($\mathrm{sinc}(u) := \sin(u)/u$):

$$\frac{\partial \mathcal M}{\partial \gamma_j}
= -\frac{L}{\pi}\,\mathrm{sinc}\big((T-\gamma_j)L\big),
\qquad
\delta\theta = \big(J^\top W J + \lambda I\big)^{-1} J^\top W\, r .$$

Linear nuisances are eliminated exactly per iteration (variable projection):
for fixed $\gamma$, $\hat c = (B^\top W B)^{-1} B^\top W (\mathcal F - \mathcal M_\gamma)$.

Initialization: method-B seeds (validated to be within $\sim 0.3$ mean gaps —
well inside the basin). Working precision: the data is computed at
`PREC` $\ge 256$ bits; the normal equations may be solved in high precision
if conditioning demands it.

### 4b. Multi-$L$ fit (experiment E3)

Same residual extended over a cutoff grid, with an $L'$-taper $v_l$ (e.g.
Hann over $[L_{\min}, L]$) so cutoff endpoints don't dominate:

$$\chi^2(\theta) = \sum_{l}\sum_i v_l\, w_i\,\big[\mathcal F(T_i, L'_l) - \mathcal M(T_i, L'_l; \theta)\big]^2 .$$

Rationale for the second dimension even though $[0,L'] \subset [0,L]$ adds no
new frequencies: at fixed $T$ each zero contributes
$-\cos\big((T-\gamma_j)L'\big)/\big(\pi (T-\gamma_j) L'\big)$, an oscillation
in $L'$ whose **frequency is the distance** $|T - \gamma_j|$ (this is exactly
the empirically validated cutoff-sweep law, corr 0.96–0.98). The collective
translation mode that defeats the bootstrap — all seeds shifted together —
changes the multi-$L$ field visibly, so the soft direction is lifted.

### 4c. Algebraic variant (optional, if 4a/4b show promise)

The $T$-derivative field is a band-limited sum of identical kernels,

$$\partial_T \mathcal F(T, L) = \frac{L}{\pi}\sum_j \mathrm{sinc}\big((T-\gamma_j)L\big) + \partial_T R,$$

whose windowed Fourier transform gives (up to the known window convolution)
samples of the pure exponential sum
$\hat\mu(\omega) = \sum_j e^{-i\omega \gamma_j}$ on $|\omega| < L$. On a
uniform grid $\omega_m = m\,\Delta\omega$ this is a Prony system
$h_m = \sum_j a_j z_j^m$, $z_j = e^{-i\Delta\omega\,\gamma_j}$, solvable by
matrix pencil / ESPRIT with no resolution limit at exact data. Conditioning
$\sigma_{\min} \sim \rho^{\,m-1}$ with $\rho = \mathrm{gap}\cdot L / 2\pi$;
e.g. at $10^{36}$, $k=10^4$: $\rho \approx 0.145$, $m \approx 7$ zeros per
Rayleigh cell $\Rightarrow$ amplification $\sim 10^{5}$–$10^{8}$ — trivial at
256–512 bits *if* the systematic $R$-leakage is kept below the target
accuracy. The fits in 4a/4b are the pragmatic, better-conditioned first cut.

## 5. Experiments

**E1 — true-seed control (premise check).** Verify that leave-one-out
subtraction itself is healthy below the threshold and only *self*-seeding is
blind. Using **Odlyzko zeros as seeds**, locate targets via the
leave-one-out hybrid crossing

$$F_C^{(j)}(T) = N_0(T) + \frac{\arg P_X + \arg Z_X^{\mathrm{excl}\,j}}{\pi},
\qquad \log Z_X(s) = -\sum_{\rho} E_1\big((s-\rho)\log X\big),$$

bisecting $F_C^{(j)} = n - \tfrac12$. Implementation: `zzz --seeds FILE`
(new) — full-precision decimal ordinates parsed into arb (doubles cannot
hold absolute $\gamma$ at these heights: resolution
$2^{-52}\gamma \approx 3\times10^5$ at $10^{21}$, which also rules out
`zhybrid` here), combined with `-R 0` so that only the target is relocated
by the final leave-one-out bisection while the true neighbours stay put.
Regimes: $10^{21}$, $k=10^4$ (gap·$L \approx 1.6$, self-boot gain 0.95×),
$10^{22}$ at $k=10^4$, and the $10^{36}$ proxy $10^{22}$, $k=600$
(gap·$L \approx 1.1$). *Prediction*: true seeds restore a clear gain;
failure falsifies the coherent-seed explanation.

**E2 — single-$L$ full-profile fit** (§4a) at the proxy regime
$10^{22}$, $k=600$, offsets 50–69 (the `boot-regime-e22-k600.tsv` ensemble):
baselines B = 0.0261, boot = 0.0258 mean |err|. Success: mean |err|
$\le 0.017$ (≥1.5×); stretch: $\le 0.013$ (2×).

**E3 — multi-$L$ fit** (§4b), same ensemble. Cutoff grid: $k' \in$
$\{150, 200, 270, 360, 450, 600\}$ ($L' \approx 6.8$–8.4), Hann taper.

**E4 — scale up** only if E2/E3 succeed: $10^{21}/10^{22}$ at $k = 10^4$,
then the $10^{36}$ showcase zero (validation = the single published
Bober–Hiary value).

Grid defaults (E2/E3): $T$-grid $T_0 \pm 4\,\mathrm{gap}\cdot W_{\rm eff}$
… in practice $T_0 \pm 2.5$ with 150–250 points, $\sigma_T \approx 0.8$;
window zeros $M = 2\cdot 24 + 1$ (core ±12 reported, rest buffer);
arb precision 256 bits for data generation; fit prototyped in
`wolframscript` (data either Wolfram-native prime sums at 40+ digits or
`zghy` dumps).

## 6. Error budget and expectations

| term | size (proxy regime) | control |
|---|---|---|
| identity error $\varepsilon_{\rm id}$ | $\sim 10^{-30}$ | none needed |
| arithmetic | $2^{-256}$ scale | raise PREC |
| far-zero leakage through window | the real floor | buffer + taper + band-edge nuisance |
| conditioning amplification | $\rho^{-(m-1)} \sim 10^2$–$10^8$ | precision bits |
| basin of attraction (NLS) | seeds within 0.3 gap | B init, LM damping |

The honest uncertainty is whether the correlated far-field leakage can be
pushed below B's 0.026 at the proxy regime; conditioning and noise are
non-issues by construction. If E2 fails but E1 succeeds, the verdict is that
full-profile fitting still doesn't exploit the evanescent band and the
algebraic route (4c) gets its turn.

## 7. Results

**E1 — true-seed control: premise confirmed.** With Odlyzko seeds
(`zzz --seeds FILE -R 0`, window ±28 or ±32), the leave-one-out relocation
works *below* the threshold where self-seeding gives nothing
(mean |err| over 10–20 zeros; B baselines from `boot-validate.tsv` /
`boot-regime-e22-k600.tsv`):

| regime | gap·$L$ | self-boot gain | true-seed gain |
|---|---|---|---|
| $10^{21}$, k=10⁴ | 1.6 | 0.95× | **5.9×** (0.0242 → 0.0041) |
| $10^{22}$, k=10⁴ | 1.55 | 0.93× | **3.5×** (0.0119 → 0.0034) |
| $10^{22}$, k=600 | 1.1 | ~1× | **2.7×** (0.0261 → 0.0096) |

No bracket-failure fallbacks. The subtraction identity is healthy below the
Riemann–Siegel scale; only self-computed seeds are blind. The true-seed
numbers are the realistic ceiling for any seed-improvement scheme feeding
the relocation step.

**E2 — single-L full-profile fit: clean negative.** `e2-singleL-fit.wls`
($10^{22}$, k=600, offsets 50–69, $W=24$, band-edge nuisance pair, B-seed
init): mean |err| 0.0258 vs B 0.0261 — **gain 1.01×**, with per-zero errors
tightly correlated with B's. Conclusion: the degeneracy lives in the
single-cutoff likelihood itself, not in how B reads it — the far-zero
leakage field is absorbed along the collective soft mode, so replacing
crossing-reading by maximum likelihood changes nothing at one cutoff.

**E3 — multi-cutoff fit: negative, worse than B.** `e3-multiL-fit.wls`
(cutoffs $X' \in \{863, 1223, 1733, 2423, 3181, 4409\}$, per-cutoff band-edge
nuisance blocks, $W=24$): mean |err| 0.0483 vs B 0.0261 — **gain 0.54×**.
Diagnosis: (i) the low-cutoff slices ($L' \approx 6.8$, kernel $\approx 3$
gaps wide) are *more* coherent than the top cutoff and pull the fit along
the soft mode harder than the top slice constrains it; (ii) the 24 nuisance
parameters absorb precisely the cross-cutoff differences that were supposed
to identify the zeros; (iii) fundamentally, prefix sums are analytically
derivable from the full-cutoff window field (band-limited ⇒ entire of
exponential type), so multi-$L$ carries **no new information** — only a
different, and evidently worse, optimization geometry.

### Verdict and the one principled escape left

Scoreboard at the coherent-deficit proxy regime ($10^{22}$, k=600,
B = 0.0261): self-boot ~1×, single-$L$ ML 1.01×, multi-$L$ ML 0.54×,
**true-seed oracle 2.7×**. The oracle gain proves the regime is recoverable
*given* the local configuration; every practical estimator built on the
GHY kernel model fails to extract it. A likely co-culprit alongside the
soft mode is **forward-model error**: the $E_1$/Si kernel is the
$X\to\infty$ limit of GHY's smoothed kernel, and the day-one cutoff sweep
already showed the Si model captures only 54–80 % of the rms at finite $X$
— orders of magnitude above the evanescent amplitudes a sub-Rayleigh fit
must read.

The clean escape, if this is ever resumed: drop the kernel approximation
entirely and fit through the **Weil explicit formula** with compactly
supported test functions $\varphi$ (with $\hat\varphi$ inside the allowed
band):

$$\sum_\rho \varphi(\gamma_\rho)
= \frac{1}{2\pi}\int \varphi(t)\,\big(\log\tfrac{t}{2\pi} + O(t^{-2})\big)\,dt
\;-\; \frac{1}{\pi}\sum_{p,m} \frac{\log p}{p^{m/2}}\,\hat\varphi\!\Big(\frac{m\log p}{2\pi}\Big)
\;+\; (\text{pole}/\Gamma\ \text{terms}),$$

an **exact identity** — zero model error; the only floor left is far-zero
leakage through $\varphi$'s tails (Slepian-optimal concentration) versus
the same band constraint. That isolates the genuinely open question:
is the evanescent signal above the leakage floor at exact arithmetic?
Untested; everything else in this note is now measured.
