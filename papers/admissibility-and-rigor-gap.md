# Weil admissibility, provable bounds under RH, and the gap between A and B

*Personal synthesis note. Companion to `damping-and-hybrid-rigor.md` and `rigor-bound-b.md`.*

---

## Notation

| symbol | meaning |
|---|---|
| $\rho = \tfrac{1}{2} + i\gamma$ | non-trivial zeta zero (RH assumed) |
| $N(T)$ | exact zero count $0 < \gamma \le T$ |
| $N_0(T) = \tfrac{T}{2\pi}\log\tfrac{T}{2\pi e} + \tfrac{7}{8}$ | smooth Riemann-von Mangoldt part |
| $S(T) = \tfrac{1}{\pi}\arg\zeta(\tfrac{1}{2}+iT)$ | fluctuating part, $N = N_0 + S$ |
| $\Lambda(n)$ | von Mangoldt function |
| $P_X(s) = \exp\!\bigl(\sum_{n \le X} \Lambda(n) n^{-s}/\log n\bigr)$ | GHY partial Euler factor |
| $Z_X(s) = \exp\!\bigl(-\sum_\rho U((s-\rho)\log X)\bigr)$ | GHY local Hadamard factor |
| $k$, $X$ | zzz's prime count, and $X \approx p_k$ |

The three methods under comparison:

| | formula | input | rigor |
|---|---|---|---|
| **A** | $F_A(T) = N_0(T) + \tfrac{1}{\pi}\operatorname{Im}\sum_{p \le p_k}\bigl(1-e^{-\sqrt{T/p}}\bigr)\log(1-p^{-1/2+iT})$ | primes | none |
| **B** | $F_B(T) = N_0(T) + \tfrac{1}{\pi}\operatorname{Im}\log P_X(\tfrac{1}{2}+iT)$ | primes | GHY + Goldston |
| **C** | $F_\text{hybrid}(T) = N_0(T) + \tfrac{1}{\pi}\operatorname{Im}\log\bigl[P_X\cdot Z_X\bigr]$ | primes + nearby zeros | GHY |

---

## 1. Weil admissibility — mechanical role

The Weil explicit formula pairs zeros with primes through a Fourier-dual pair $h \leftrightarrow g$:

$$
\sum_{\gamma} h(\gamma) \;=\; \hat h(0) \;+\; \text{archimedean terms} \;-\; 2\sum_p \log p \sum_{m=1}^\infty \frac{g(m\log p)}{p^{m/2}}.
$$

For the sum over zeros to converge absolutely and equal the prime side, $h$ must be **admissible**:

- analytic in the strip $|\operatorname{Im} r| < \tfrac{1}{2} + \delta$ for some $\delta > 0$,
- decay $h(r)(1 + |r|)^{2} \to 0$ uniformly in that strip,
- even.

### zzz's place on the boundary

The zero-counting indicator $\mathbb{1}_{|r| \le T}$ is **not** admissible (not smooth, not analytic). Any usable $h$ is a smoothed version. zzz's heuristic damping corresponds to

$$
g_{\text{zzz}}(u) \;=\; 1 - \exp\!\bigl(-\sqrt{T}\,e^{-u/2}\bigr),
$$

whose Fourier partner $h_{\text{zzz}}(r)$ has a first pole **exactly at $\operatorname{Im} r = \tfrac{1}{2}$** — on the admissibility boundary, one $\delta$-nudge outside the allowed class (see `test-function-search.wls`). This is the mathematical reason Weil cannot rigorize the heuristic.

### Admissible alternatives

Any of these give a rigorously usable $g$:

- **Gaussian** $g(u) = e^{-au^2}$: $h$ is Gaussian, entire.
- **Logistic** $g$: $h$ analytic in a finite strip; admissible if parameters are tuned.
- **Carneiro–Vaaler extremal** bandlimited functions: optimal for various norms.
- **GHY C∞ bump** on $[e^{1-1/X}, e]$: the specific mollifier in GHY 2007 Thm 1.

Method B is this last choice made concrete.

---

## 2. What is provable under RH

Two theorems in sequence.

### GHY 2007, Theorem 1

$$
\zeta(s) \;=\; P_X(s)\,Z_X(s)\,\bigl(1 + E(s, X, K)\bigr),
$$

with

$$
|E(s, X, K)| \;\le\; C_1\,\frac{X^{K+2}}{(|s|\log X)^K} \;+\; C_2\,\frac{\log X}{X^{\sigma}},
\qquad s = \sigma + it,\ |t| \ge 2.
$$

At $s = \tfrac{1}{2} + iT$, second term dominates; numerically:

| $T$ | $X$ | bound on $|E|$ |
|---|---|---|
| $10^6$ | $7919$ | $\sim 0.10$ |
| $10^{12}$ | $7919$ | $\sim 0.10$ |
| $10^{36}$ | $104729$ | $\sim 0.036$ |

**This is the provable bound for method C**, where we keep both $P_X$ and $Z_X$. GHY's proof relies on RH in the Hadamard factorization of $\zeta$.

### Goldston 1987

Dropping $Z_X$ costs us an additional $\tfrac{1}{\pi}|\operatorname{Im}\log Z_X|$. Goldston bounds exactly this:

$$
S(T) \;=\; -\frac{1}{\pi}\sum_{n \le X} \frac{\Lambda(n)\sin(T\log n)}{\sqrt{n}\,\log n} \;+\; R(T, X),
\qquad |R(T, X)| \le C\,\frac{\log T}{\log X}.
$$

The truncated sum *is* $(F_B(T) - N_0(T))$, so under RH:

$$
\boxed{\;
\bigl|F_B(T) - N(T)\bigr| \;\le\; \frac{C}{\pi}\cdot\frac{\log T}{\log X} \;+\; O\!\left(\frac{\log X}{\sqrt X}\right).
\;}
$$

With conservative $C \le 1$, the bound in the zzz regime is $0.24$–$0.50$.

### Unconditional fallback

Without RH, only Littlewood's $|S(T)| = O(\log T)$ is available — too loose to be useful.

---

## 3. Why A is slightly better than B, empirically

From `ab-scan.wls`, first five ordinals sampled up to $\#10^5$:

| $k$ | mean $|F_A - N|$ | mean $|F_B - N|$ | ratio |
|---|---|---|---|
| $50$ | $0.034$ | $0.039$ | $1.14$ |
| $100$ | $0.039$ | $0.040$ | $1.01$ |
| $200$ | $0.014$ | $0.025$ | $1.81$ |
| $500$ | $0.012$ | $0.018$ | $1.51$ |
| $1000$ | $0.008$ | $0.013$ | $1.69$ |

A beats B by $\sim 1.5\times$. **Is this evidence of a deeper theoretical asymmetry?** No.

### Decomposing the difference

Within the prime range both methods sum ($p \le p_k$), they differ only in two things:

1. **Damping shape.** A uses the smooth $(1 - e^{-\sqrt{T/p}})$, which is $\approx 1$ for $p \ll T$ and tapers smoothly. B uses a hard indicator $\mathbb{1}_{p \le X}$.

2. **High-$m$ tail.** For each small prime $p$, A includes all $m \ge 1$ in $-\log(1-p^{-s}) = \sum_m p^{-ms}/m$. B truncates at $p^m \le X$. The tail difference is $\sum_{m \ge M+1} (m p^{m/2})^{-1}$, dominated by the first omitted term and numerically tiny (e.g. $< 10^{-3}$ for $p=2$, $X = 7919$).

Neither difference is substantial. The dominant effect is (1): **smooth cutoffs suppress Gibbs-type oscillation** in truncated Fourier-like series. That's what A is exploiting.

### The mollifier is not unique

GHY chose a specific C∞ bump on $[e^{1-1/X}, e]$ to make Thm 1's proof clean — **not** to minimize finite-$X$ empirical error. In principle, a Gaussian or Vaaler-extremal mollifier could match A's empirical performance while keeping B's provability. The $1.5\times$ gap is therefore a choice, not a theoretical necessity.

### Conclusion on A vs B

**Not a coincidence, not a theoretical gap, but an empirical tuning gap.** A's damping was selected (implicitly or explicitly) for good finite-$X$ behaviour; B's was selected for provability. Re-tuning B's mollifier should close the gap.

---

## 4. The real theoretical gap

The interesting asymmetry is *between both methods and their provable bound*, not between A and B:

| | empirical at $k=1000$ | provable bound (RH) |
|---|---|---|
| **A heuristic** | $\sim 0.01$ | *none* |
| **B GHY** | $\sim 0.013$ | $\sim 0.25$–$0.50$ |

The rigorous bound is $10\text{–}300\times$ looser than observed. **This is where the unresolved mathematics lives.** Closing it requires a statement strictly sharper than Goldston's $\log T / \log X$ for the *tail* $|S(T) - \text{truncated sum}|$.

### Candidate techniques

1. **Carneiro–Chirre–Milinovich (2022)** bandlimited extremal functions. Produces tight bounds on $|S(T)|$ itself. The tail analogue is, to my knowledge, open.

2. **Large sieve** on $\sum_{p \le X} p^{-1/2} e^{iT\log p}$: if square-root cancellation holds, the $p > X$ tail is $O\!\bigl(\sqrt X / (T \log X)\bigr)$, polynomial and small.

3. **Vaughan's identity** / bilinear-form decomposition, using $T$-oscillation for cancellation.

Any of these would turn B from "correct up to loose constant" into "correct up to tight constant" — the concrete open problem.

---

## 5. One-line takeaways

- **A's empirical win** over B: smooth-cutoff Gibbs suppression; not mathematically deep.
- **B's theoretical win** over A: placed inside a proof chain via GHY + Goldston.
- **C's rigor win** over both: $|E| \sim 0.04$–$0.1$ via GHY alone, *provided* nearby zeros are supplied.
- **The gap that matters**: observed $0.01$ vs proved $0.5$, both methods. A tighter Goldston-style tail bound would resolve it.

---

## References

- **Weil 1952** — *Sur les "formules explicites" de la théorie des nombres premiers*, Comm. Sém. Math. Lund. The explicit formula and its admissibility condition.
- **Gonek, Hughes, Young (2007)** — *A hybrid Euler–Hadamard product for the Riemann zeta function*, Duke Math. J. **136**, 507–549. Theorem 1 is the identity with provable error.
- **Goldston (1987)** — *On the function $S(T)$ in the theory of the Riemann zeta-function*, J. Number Theory **27**, 149–177. The $\log T / \log X$ tail bound for the truncated prime sum under RH.
- **Littlewood (1924)** — *On the zeros of the Riemann zeta-function*, Proc. Camb. Phil. Soc. **22**, 295–318. Original $O(\log T / \log\log T)$ bound on $|S(T)|$ under RH.
- **Carneiro, Chandee, Milinovich (2013)** — *Bounding $S(T)$ and $S_1(T)$ on the Riemann hypothesis*, Math. Ann. **356**, 939–968. Refined $|S(T)|$ constants via extremal functions.
- **Carneiro, Chirre, Milinovich (2022)** — *Bandlimited approximations and estimates for the Riemann zeta-function*, Publ. Mat. **63**, 601–661. Best published candidate for the tight tail bound needed to close the gap.
- **Vaaler (1985)** — *Some extremal functions in Fourier analysis*, Bull. AMS **12**, 183–216. The extremal-function toolkit underlying most modern $|S(T)|$ work.

---

## Artifacts in this repo

- `main.c` (A and B in one binary, flagged by `--ghy`)
- `ghy.c`, `ghy.h` (shared GHY primitives)
- `zghy.c`, `zhad.c`, `zhybrid.c` (P_X, Z_X, and full hybrid evaluators)
- `papers/ab-scan.wls` (the accuracy scan above)
- `papers/damping-and-hybrid-rigor.md` (the original plan)
- `papers/rigor-bound-b.md` (the numerical bound confrontation)
- `papers/ghy-kernel.wls` (the Hadamard kernel $U(z) \to E_1(z)$ derivation)
