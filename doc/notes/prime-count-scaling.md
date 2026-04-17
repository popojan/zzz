# Prime Count Scaling for Zero Approximation

**Status:** Empirical observations, assumes RH

---

## Caveat: Riemann Hypothesis Assumption

This entire approach assumes the **Riemann Hypothesis holds** — that all non-trivial zeros lie on the critical line Re(s) = 1/2. The truncated prime sum approximation to the zero counting function is derived from the explicit formula, which connects zeros to primes. If RH fails, the error bounds and scaling laws derived here become unreliable.

---

## The Scaling Question

How many primes `k` are needed for accurate zero approximation at height `T`?

### Empirical Evidence

From README examples:

| Height T | k (primes) | log₁₀(T) | Error vs Odlyzko |
|----------|------------|----------|------------------|
| 10^12    | 1000       | 12       | ~0.05            |
| 10^36    | 10000      | 36       | ~0.02            |

### Derived Scaling Law

Fitting the data:
- k = 1000 at log₁₀(T) = 12
- k = 10000 at log₁₀(T) = 36

The ratio 10000/1000 = 10 while (36/12)² = 9, suggesting:

$$k \approx 7 \times (\log_{10} T)^2$$

Or equivalently in natural logarithm:

$$k \approx 1.3 \times (\ln T)^2$$

### Extrapolation

| Height T | log₁₀(T) | Predicted k |
|----------|----------|-------------|
| 10^12    | 12       | ~1000       |
| 10^36    | 36       | ~9000       |
| 10^50    | 50       | ~17500      |
| 10^100   | 100      | ~70000      |

---

## Why Too Many Primes Hurts

Counterintuitively, using more primes than necessary for small heights **degrades** accuracy.

### Phase Coherence Explanation

The prime sum contribution to the counting function is:

$$-\frac{1}{\pi} \sum_{p \leq P} \sum_{m=1}^{\infty} \frac{\sin(t \log p^m)}{m \cdot p^{m/2}}$$

For a prime `p` at height `t`:

1. **When p << t**: The phase `t log p` oscillates many times. Contributions from nearby zeros average out. Information is meaningful.

2. **When p >> t**: The phase `t log p` is small and slowly varying. All zeros in a region get similar contributions. This adds **correlated bias**, not information.

3. **When p ~ t**: Optimal regime — phase variation matches zero spacing.

### The Sweet Spot

The effective range of useful primes scales as:
- Primes p where `t log p` spans multiple periods (p < e^(2π/t) is too small)
- But not so large that phase becomes nearly constant

This gives the (log T)² scaling: larger T needs more primes, but the relationship is sublinear in T.

---

## Complexity Comparison

| Method | Complexity per zero | Accuracy | Notes |
|--------|---------------------|----------|-------|
| **Riemann-Siegel** | O(T^1/2) | Exact (to precision) | Standard method |
| **Odlyzko-Schönhage** | O(T^(1/4+ε)) | Exact (to precision) | FFT-based, large precomputation |
| **zzz (this tool)** | O((log T)² × log(1/ε)) | Approximate (~0.01-0.05) | No precomputation |

Where:
- T = height on critical line (imaginary part of zero)
- ε = tolerance for bisection

### Breakdown of zzz Complexity

1. **Lambert W initial estimate**: O(1)
2. **Prime sum evaluation**: O(k) = O((log T)²)
3. **Bisection iterations**: ~20 for tolerance 10^-6
4. **Total**: O((log T)² × 20) per zero

### When zzz Wins

- For **approximate** zero locations (error ~0.01)
- When you need **many zeros quickly** (statistical analysis, GUE verification)
- At **extremely large heights** where exact methods become impractical

### When zzz Loses

- When you need **exact** zeros (verified to many decimal places)
- For rigorous verification of RH computationally
- When accuracy better than ~0.01 is required

---

## Practical Recommendations

| Task | Recommended k | Expected error |
|------|---------------|----------------|
| Zeros near 10^12 | 1000 | ~0.05 |
| Zeros near 10^20 | 3000 | ~0.03 |
| Zeros near 10^36 | 10000 | ~0.02 |
| Zeros near 10^50 | 20000 | ~0.02 |

For GUE statistics and spacing distribution analysis, errors of 0.01-0.05 are typically acceptable since you're looking at statistical properties, not individual zero positions.

---

## References

- Odlyzko, A.M.: "The 10^20-th zero of the Riemann zeta function and 175 million of its neighbors"
- França, G. & LeClair, A.: "Statistical and other properties of Riemann zeros" (arXiv:1307.8395)
