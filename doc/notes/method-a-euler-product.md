# What method A actually computes — the effective Euler product

*Companion to `admissibility-and-rigor-gap.md`. Details what finite object the heuristic `zzz` default secretly evaluates, where its damping has bite, and why the theoretical "unlimited primes" advantage does not materialise in practice.*

---

## 1. The exact formula

Method A (the `zzz` default, `main.c::zero_count_exact`) evaluates

$$
F_A(T) \;=\; N_0(T) \;+\; \frac{1}{\pi}\operatorname{Im}\sum_{p \le p_k}\underbrace{\bigl(1 - e^{-\sqrt{T/p}}\bigr)}_{\text{att}(p,T)} \,\log\!\bigl(1 - p^{-1/2+iT}\bigr)
$$

with the damping

$$
\operatorname{att}(p, T) \;=\; 1 - e^{-\sqrt{T/p}}.
$$

Pulling the weight inside the logarithm turns the sum into the log of a finite product. Define

$$
\boxed{\;Z_k(s, T) \;:=\; \prod_{p \le p_k}\bigl(1 - p^{-s}\bigr)^{-\operatorname{att}(p, T)}.\;}
$$

Then

$$
F_A(T) \;=\; N_0(T) \;+\; \frac{1}{\pi}\arg Z_k\!\left(\tfrac{1}{2} - iT,\,T\right),
$$

and `zzz` bisects for $T$ such that $F_A(T) = n + \tfrac{1}{2}$. Geometrically: the argument of the finite proxy $Z_k$ jumps by $\pi$ near every $\zeta$-zero (approximately), and the bisection catches the midpoint of that jump.

`zproxy` in this repo evaluates $Z_k$ directly and confirms $\log|Z_k|$ dips sharply at $T = \gamma_n$ — see `doc/ghy/` / commit `109a59ee`.

---

## 2. What $Z_k$ converges to as $k \to \infty$

Taking the formal limit,

$$
Z_\infty(s, T) \;=\; \prod_{p}\bigl(1 - p^{-s}\bigr)^{-\operatorname{att}(p, T)}.
$$

Split the exponent $-\operatorname{att} = -1 + (1 - \operatorname{att}) = -1 + e^{-\sqrt{T/p}}$ and factor:

$$
\boxed{\;Z_\infty(s, T) \;=\; \zeta(s)\,\cdot\,\prod_{p}\bigl(1 - p^{-s}\bigr)^{\exp(-\sqrt{T/p})}\;.\;}
$$

The correction factor's exponent behaves cleanly:

| regime of $p$ | $e^{-\sqrt{T/p}}$ | correction Euler factor | net effect |
|---|---|---|---|
| $p \ll T$ | $\approx 0$ | $(1-p^{-s})^{0} = 1$ | leaves $\zeta$'s Euler factor alone |
| $p = T$ | $e^{-1} \approx 0.37$ | $(1-p^{-s})^{0.37}$ | partial cancellation |
| $p \gg T$ | $\approx 1$ | $(1-p^{-s})^{1}$ | cancels $\zeta$'s factor |

So to leading order

$$
Z_\infty(s, T) \;\approx\; \prod_{p \lesssim T}\bigl(1 - p^{-s}\bigr)^{-1},
$$

i.e. the **partial Euler product of $\zeta(s)$ softly truncated at $p \sim T$**. This is *exactly* the object at the core of the Gonek–Hughes–Young hybrid formula — method B's $P_X$ with the specific choice $X \sim T$.

A and B therefore have the same mathematical ceiling. The damping in A is a device to make the $k \to \infty$ limit sit at a definite location; it does not unlock new information.

---

## 3. Where the damping actually bites — numerical reality check

Damping becomes numerically distinguishable from $1$ (to double precision, say when $\operatorname{att} \le 1 - 10^{-16}$) only when

$$
\sqrt{T/p} \;\lesssim\; 37 \;\Longleftrightarrow\; p \;\gtrsim\; \frac{T}{1400}.
$$

In the `zzz` operating regime $k \sim (\log T)^2$, so $p_k \sim k \log k$:

| height $T$ | $k$ (from `(log T)^2`) | $p_k$ | threshold $T/1400$ | is any summed prime in the damping region? |
|---|---|---|---|---|
| $10^6$ | $200$ | $1223$ | $\sim 700$ | possibly |
| $10^{12}$ | $760$ | $5807$ | $\sim 7 \cdot 10^8$ | **no** |
| $10^{36}$ | $6900$ | $\sim 77000$ | $\sim 7 \cdot 10^{32}$ | **no** |

**At any height `zzz` is designed for ($T \ge 10^{12}$), every prime actually summed has $\operatorname{att}(p, T) = 1$ to working precision.** The damping is *operationally inert*.

---

## 4. Does A have an ultra-height advantage from unlimited primes?

The only way A could outperform B would be to keep summing primes *into the damping region*, where the weight tapers smoothly rather than cutting off hard. Count the cost:

- damping starts to matter at $p \sim T/1400$,
- $\pi(T/1400) \sim T / (1400 \log T)$ primes would need to be enumerated just to cross the threshold.

For $T = 10^{36}$ this is $\sim 10^{31}$ primes — physically impossible. For $T = 10^{12}$, $\sim 10^{10}$ — at $10^7$ primes/sec that's $\sim 30$ years. For $T = 10^6$, $\sim 4 \cdot 10^4$ primes — *that* we could do, but it's the only height where it's tractable, and it's also the height where zzz already works well.

**Conclusion**: A's "in principle $k \to \infty$" route is theoretically elegant but practically unreachable at the heights A is supposed to shine at. At every reachable $k$, A and B see the same primes with the same effective weight $1$.

---

## 5. Does the smooth damping kill corrections at $p \approx T$?

It *attenuates* them, does not kill. The damping value $\operatorname{att}(T, T) = 1 - 1/e \approx 0.632$: still 63 % of the full contribution. By $p = 4T$ it is $0.39$, and by $p = 100 T$ still $0.10$. The decay is polynomial ($\sim \sqrt{T/p}$ for $p \gg T$), never abrupt.

But this does not translate into practical correction because, as §3 shows, no prime summed by `zzz` lies anywhere near $p \approx T$. The damping curve is drawn on primes that never get evaluated.

---

## 6. A vs B, restated cleanly

| | kernel form | effective cutoff | $k \to \infty$ limit | reachable | bound under RH |
|---|---|---|---|---|---|
| **A** | $\operatorname{att}(p,T)\log(1-p^{-s})$, primes only | soft, at $p \sim T$ | partial Euler product of $\zeta$ up to $p \sim T$ | no, for $T \gg 1$ | none |
| **B** | $\Lambda(n)/(n^s \log n)$, prime powers $n \le X$ | hard, at $n = X$ | partial Euler product of $\zeta$ up to $n \sim X$ | $X \sim p_k$ only | Goldston + GHY |

The two limits agree when the cutoff scales match ($X \sim T$). Neither is reachable at ultra-heights; both are then pinned to the same primes $p \le p_k$ with $p_k \ll T$, and they differ only in (i) how sharply the non-summed tail is nominally truncated and (ii) whether a provable error bound is attached.

---

## 7. Practical takeaways

1. **The heuristic damping is not doing arithmetic work at ultra-heights.** It is $\equiv 1$ on every prime the code actually sees.
2. **A's theoretical edge (convergent infinite-prime extension) is unreachable.** The first prime where damping distinguishes from $1$ is at $p \sim T/1400$, which is astronomically far beyond any feasible $k$.
3. **The B-over-A gap is in provability, not accuracy.** Since the damping is inert, A and B sum essentially identical things; A's empirical $1.5\times$ edge comes from the *shape* of the nominal cutoff (smooth vs sharp), not from using different primes — see `admissibility-and-rigor-gap.md` §3.
4. **If one wanted the heuristic to genuinely incorporate primes near $T$**, the prime count would have to scale like $T / \log T$, which is exactly the Riemann–Siegel regime and destroys the polylog cost that justifies `zzz` in the first place.

---

## 8. Implementation pointers

- `main.c::zero_count_exact` — the actual summation
- `zproxy.c` — standalone binary that evaluates $Z_k(s, T)$ directly; `doc/ghy/` visualisations confirm $|Z_k|$ dips at $\zeta$-zeros
- The factorisation $Z_\infty = \zeta \cdot \prod (1 - p^{-s})^{e^{-\sqrt{T/p}}}$ is algebraic — a one-line rearrangement of the exponent, no analytic continuation involved

---

## References

- `doc/notes/admissibility-and-rigor-gap.md` — A/B/C overview and provable bounds
- `doc/notes/damping-and-hybrid-rigor.md` — original analysis of zzz's damping
- `doc/notes/rigor-bound-b.md` — numerical confrontation of Goldston bound vs observation
- Gonek, Hughes, Young, *A hybrid Euler–Hadamard product for the Riemann zeta function*, Duke Math. J. **136** (2007), 507–549.
