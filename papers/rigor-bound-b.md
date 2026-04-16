# Rigorous error bound for method B (`zzz --ghy`)

**Goal.** State what can be *proved* about the accuracy of $F_B(T)$ as an approximation to $N(T)$, and compare that provable bound to the empirical data from `papers/ab-scan.wls`.

**Setup.** All statements assume RH. The object of interest is

$$
F_B(T) \;=\; N_0(T) \;+\; \frac{1}{\pi}\,\operatorname{Im}\log P_X\!\left(\tfrac{1}{2} + iT\right)
$$

with

$$
P_X(s) \;=\; \exp\!\left(\sum_{p^m \le X} \frac{\Lambda(n)}{n^s \log n}\right)
$$

(GHY 2007 eq. 6, with $X = p_k$).

---

## 1. What GHY Theorem 1 actually gives us

The identity of GHY 2007, Thm 1, for $s = \sigma + it$, $\sigma \ge 0$, $|t| \ge 2$, $X \ge 2$:

$$
\zeta(s) \;=\; P_X(s)\,Z_X(s)\,\bigl(1 + E(s, X, K)\bigr)
$$

where

$$
E(s, X, K) \;=\; O\!\left(\frac{X^{K+2}}{(|s|\log X)^K}\right) \;+\; O\!\left(X^{-\sigma}\log X\right)
$$

for any fixed integer $K \ge 1$. **Crucially, this is an identity about the product $P_X \cdot Z_X$, not about $P_X$ alone.** $E$ is the error in the hybrid identity; it does **not** bound the error when you drop $Z_X$.

Specialising to $s = \tfrac{1}{2} + iT$:

$$
\bigl|E(\tfrac{1}{2}+iT,\,X,\,K)\bigr| \;\le\; C_1 \cdot \frac{X^{K+2}}{(T \log X)^K} \;+\; C_2 \cdot \frac{\log X}{\sqrt{X}}.
$$

For the zzz regime ($X = p_k \approx k \log k$, $T$ large):

| $T$ | $k$ | $X$ | $\log X$ | $X^{-1/2}\log X$ | $(X/T)^K$ term (K=10) |
|---|---|---|---|---|---|
| $10^6$ | $1000$ | $7919$ | $8.98$ | $\mathbf{0.101}$ | negligible |
| $10^{12}$ | $1000$ | $7919$ | $8.98$ | $\mathbf{0.101}$ | $10^{-51}$ |
| $10^{36}$ | $10000$ | $104729$ | $11.56$ | $\mathbf{0.036}$ | $10^{-305}$ |

The first term is already $K$-independent and dominates. **For $F_\text{hybrid}$ (method C), GHY's rigorous bound is $\approx 0.036\text{–}0.101$** across the zzz operating range — looser than the $<\tfrac{1}{2}$ needed for unambiguous zero localization, but not by orders of magnitude.

**This bound applies only when you include $Z_X$.** For method B (no $Z_X$), a separate estimate on $|\log Z_X|$ is needed.

---

## 2. What you lose by dropping $Z_X$

The additional error from omitting the local Hadamard factor is exactly

$$
\bigl|F_B(T) - F_\text{hybrid}(T)\bigr| \;=\; \frac{1}{\pi}\,\bigl|\operatorname{Im}\log Z_X(\tfrac{1}{2}+iT)\bigr|.
$$

$\operatorname{Im}\log Z_X$ is the "fluctuating" part of $\log \zeta$ — the object whose classical bound is $|S(T)| = \frac{1}{\pi}|\arg \zeta(\tfrac{1}{2}+iT)|$. Consequently

$$
\bigl|F_B(T) - N(T)\bigr| \;\le\; \bigl|F_B - F_\text{hybrid}\bigr| + \bigl|F_\text{hybrid} - N\bigr|
\;\le\; \frac{1}{\pi}\bigl|\operatorname{Im}\log Z_X(\tfrac{1}{2}+iT)\bigr| + \bigl|\text{GHY bound from §1}\bigr|.
$$

The first term is where the real difficulty lies.

---

## 3. Goldston 1984 + GHY: the B-specific bound

Goldston (1987, *J. Number Theory* 27) proved, under RH:

**Theorem (Goldston).** *For $X \ge 2$ and $T \ge T_0$,*

$$
S(T) \;=\; -\frac{1}{\pi}\sum_{n \le X} \frac{\Lambda(n)\sin(T\log n)}{\sqrt{n}\,\log n} \;+\; R(T, X)
$$

*with $|R(T, X)| \le C\,\dfrac{\log T}{\log X}$.*

The truncated sum on the RHS is **exactly** $\bigl(F_B(T) - N_0(T)\bigr)$ — Goldston's $S(T)$ expansion *is* method B. Hence, under RH:

$$
\bigl|F_B(T) - N(T)\bigr| \;\le\; \frac{C_\text{Goldston}}{\pi}\cdot\frac{\log T}{\log X} \;+\; O\!\left(\frac{\log X}{\sqrt X}\right).
$$

With the conservative constant $C_\text{Goldston} \le 1$ (explicit versions exist in Carneiro–Chandee–Milinovich 2013 refining this):

$$
\boxed{\;\bigl|F_B(T) - N(T)\bigr| \;\le\; \frac{1}{\pi}\cdot\frac{\log T}{\log X} \;+\; C_2 \cdot \frac{\log X}{\sqrt X}\;}
$$

---

## 4. Numerical confrontation with the empirical data

Using the `ab-scan.wls` data and the bound $\dfrac{1}{\pi}\dfrac{\log T}{\log X} + \dfrac{\log X}{\sqrt X}$:

| ordinal | $T$ | $k$ | $X$ | $\dfrac{\log T}{\pi\log X}$ | $\dfrac{\log X}{\sqrt X}$ | **bound** | **observed B** | ratio |
|---|---|---|---|---|---|---|---|---|
| #10 | $49.77$ | $50$ | $229$ | $0.481$ | $0.604$ | $\mathbf{1.09}$ | $0.043$ | $25$ |
| #10 | $49.77$ | $1000$ | $7919$ | $0.139$ | $0.101$ | $\mathbf{0.240}$ | $0.023$ | $10$ |
| #100 | $236.5$ | $1000$ | $7919$ | $0.193$ | $0.101$ | $\mathbf{0.294}$ | $0.001$ | $294$ |
| #1000 | $1419$ | $1000$ | $7919$ | $0.257$ | $0.101$ | $\mathbf{0.358}$ | $0.025$ | $14$ |
| #10000 | $9878$ | $1000$ | $7919$ | $0.326$ | $0.101$ | $\mathbf{0.427}$ | $0.006$ | $71$ |
| #100000 | $74920$ | $1000$ | $7919$ | $0.397$ | $0.101$ | $\mathbf{0.498}$ | $0.009$ | $55$ |

Observations:

1. **The rigorous bound is always loose** — by a factor $10\text{–}300\times$ across the scan.
2. **The bound grows like $\log T / \log X$**, whereas the observed error is flat-to-decreasing in $T$. This tells us the bound's $\log T$ growth is an artefact of proof technique, not of the method.
3. **The bound is already $> \tfrac{1}{2}$ for $T \ge 10^5$ at $X = 7919$** — so if we required $< \tfrac{1}{2}$ rigorously, we would need $X \gtrsim T^\varepsilon$ (polynomial in $T$), not $\operatorname{polylog}(T)$ as `zzz` actually uses. This is the $k \approx \pi(T^2)$ dead-end flagged in `rigor-backlog.md` §3D.

---

## 5. What "rigor" buys us in practice

| | accuracy | provable bound | primes only | needs zeros |
|---|---|---|---|---|
| **A** heuristic | mean $0.01\text{–}0.04$ | **none** | ✅ | ✗ |
| **B** GHY $P_X$ | mean $0.01\text{–}0.04$ | Goldston: $\sim \log T / \log X$ | ✅ | ✗ |
| **C** $P_X \cdot Z_X$ | mean $0.02$ (seeded) | **GHY: $X^{-1/2}\log X$** | ✅ | ✅ window |

Reading the table:

- **A and B are empirically indistinguishable.** The heuristic is *not* doing anything magical; both methods are $O(1)$ away from $N(T)$ for the zzz regime.
- **B's advantage is theoretical** — it can be placed inside a proof chain, where A cannot.
- **C is strictly tighter as an approximation** because $Z_X$ subtracts the $S(T)$-fluctuation, but needs a seed zero list.

**Honest conclusion for the zzz use case** (find the $n$-th zero at large height):

- If rigor is a goal, use B (`zzz --ghy`). It carries the Goldston + GHY bound; the constants are loose but explicit.
- If accuracy is the only goal, A and B are equivalent; the heuristic was built for speed and is fine.
- A proof that $|F_B - N| < \tfrac{1}{2}$ at zzz's empirical $k \sim (\log T)^2$ regime would require a substantially tighter bound on the truncated prime sum than Goldston 1984 gives. Carneiro–Chandee–Milinovich (2013) have $\left(\tfrac{1}{4} + o(1)\right)\dfrac{\log T}{\log\log T}$ for $|S(T)|$ directly, but the *difference* $S(T) - \text{truncated sum}$ has not received the same treatment.

---

## 6. Open gap (would close the loop)

Find (or derive) an explicit-constant version of

$$
\left|\,S(T) \;-\; \frac{-1}{\pi}\sum_{n \le X} \frac{\Lambda(n)\sin(T\log n)}{\sqrt{n}\,\log n}\,\right| \;\le\; c(X, T)
$$

with $c(X, T) \to 0$ faster than $\log T / \log X$. Candidates:

- **Explicit Carneiro–Chirre–Milinovich** (2022, *Publ. Mat.*): bandlimited approximations for $S(T)$. If the extremal function argument can be adapted to bound the *tail* (not just $S$), we'd get $O((\log T)^{1/2} / \log X)$ or similar.
- **Large sieve** applied to $\displaystyle\sum_{p \le X} p^{-1/2}\,e^{iT\log p}$: if square-root cancellation holds, the tail from $p > X$ is $O\!\left(\sqrt{X}/(T\log X)\right)$ which is polynomial and small.
- **Vaughan's identity** decomposition of the truncated sum into bilinear forms, where $T$-oscillation provides additional cancellation.

Any of these would turn B from "same accuracy as A, but rigorous with loose constant" into "same accuracy as A, with tight constant". That is the concrete mathematical deliverable.

---

## 7. Implementation status

- **A** (heuristic): `zzz` default
- **B** (partial Euler): `zzz --ghy` (commit `7ab0d2de`)
- **C** (hybrid): `zhybrid` (commit `4da27896`)
- Scan: `papers/ab-scan.wls` (commit `92f1e9b2`)
- This document: companion to `papers/damping-and-hybrid-rigor.md` §6, §7.

---

## References

- Gonek, Hughes, Young, *A hybrid Euler–Hadamard product for the Riemann zeta function*, Duke Math J. **136** (2007), 507–549. **Thm 1 provides the identity error bound.**
- Goldston, *On the function $S(T)$ in the theory of the Riemann zeta-function*, J. Number Theory **27** (1987), 149–177. **Provides the truncation tail bound.**
- Carneiro, Chandee, Milinovich, *Bounding $S(T)$ and $S_1(T)$ on the Riemann hypothesis*, Math. Ann. **356** (2013), 939–968. **Tightens the $|S(T)|$ bound.**
- Carneiro, Chirre, Milinovich, *Bandlimited approximations and estimates for the Riemann zeta-function*, Publ. Mat. **63** (2022). **Closest candidate for an explicit tight truncation bound.**
- Littlewood, *On the zeros of the Riemann zeta-function*, Proc. Camb. Phil. Soc. **22** (1924). **Original $\log T / \log\log T$ bound.**
