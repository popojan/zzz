# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`zzz` ("Zeta Zeros Zeal") computes high-ordinal Riemann zeta zeros on the critical line by **bisecting a zeta-evaluation-free counting function** of the form `F(T) = N_0(T) + arg(...)/π`. No `acb_zeta` is called in the inner loop; the work is a sum over primes (and optionally over nearby known zeros).

## Build

```bash
cd build && cmake .. && make            # dynamic link against system libs
cmake -DSTATIC_LINKING=ON .. && make    # link against deps/{gmp,flint,mpfr,arb}
```

Dependencies: FLINT (provides `arb`/`acb`), GMP, MPFR. Static-link mode expects sources unpacked under `deps/` (gmp-6.2.1, flint-2.9.0, mpfr-4.2.0, arb).

There is no test suite. Empirical verification lives in Wolfram scripts under `doc/ghy/` (e.g. `wolframscript -file doc/ghy/ab-scan.wls` for the A-vs-B accuracy scan).

## Binaries and what each does

All built from `CMakeLists.txt`. `libghy.a` is a static library shared by everything except `zproxy`.

| Binary | Source | Purpose |
|---|---|---|
| `zzz` | `main.c` | User-facing CLI. Bisects for the n-th zero. Default = method A (heuristic). `--ghy` switches the inner counter to method B (rigorous, primes-only). |
| `zhybrid` | `zhybrid.c` | Method C dumper: full hybrid `P_X·Z_X` evaluator on a T-grid. Requires a seed-zeros file. **Not integrated into the `zzz` CLI.** |
| `zghy` | `zghy.c` | Diagnostic: dumps `P_X(1/2+iT)` and `F_P(T)` on a grid. |
| `zhad` | `zhad.c` | Diagnostic: dumps the local Hadamard factor `Z_X(1/2+iT)` from a zeros file. |
| `zproxy` | `zproxy.c` | Diagnostic: dumps method A's finite damped Euler proxy `Z_k(s,T)`. |

`pellreg`, `pellreg2`, `psi.c` in the tree are unrelated to the zeta-zeros work and are not in `CMakeLists.txt`.

## The three counting methods

All three compute `F(T) ≈ N(T)` without evaluating ζ. They differ only in the inner factor whose argument is taken; the bisection driver is the same.

| Method | Formula | Inputs | Rigor | Implementation |
|---|---|---|---|---|
| **A** heuristic | `F_A = N_0 + (1/π) Im Σ_{p≤p_k} (1 - exp(-√(T/p))) · log(1 - p^{-1/2+iT})` | primes only | none (Weil-non-admissible damping) | `main.c::zero_count_exact` |
| **B** GHY `P_X` | `F_B = N_0 + arg P_X / π`, `log P_X = Σ_{p^m ≤ X} 1/(m p^{ms})` | primes only | GHY 2007 Thm 1 + Goldston 1987 | `main.c::zero_count_ghy` → `ghy_log_px` |
| **C** hybrid | `F_C = N_0 + (arg P_X + arg Z_X)/π`, `log Z_X = -Σ_ρ E_1((s-ρ)log X)` | primes + nearby zeros | GHY 2007 Thm 1 (tightest) | `zhybrid.c` |

Method A and B sum the same primes and produce different but close numbers (~1.5× error gap; A's damping is operationally inert in feasible regimes). Method C is the rigorously bounded one but needs seed zeros.

`ghy.h` documents the GHY primitives (`ghy_log_px`, `ghy_log_zx`, `ghy_n0_smooth`) and is the contract between `main.c` and the helper binaries.

## Code structure

- `main.c` (961 lines) — argp-based CLI, bisection driver, both A and B counters live here. Search for `zero_count_exact` (A), `zero_count_ghy` (B), and the option table near the `argp_option options[]` declaration to add CLI flags.
- `ghy.c` / `ghy.h` — shared primitives: `P_X`, `Z_X`, smooth `N_0`. Anything that adds a new counting variant should expose its kernel here so that all binaries can share it.
- `z*.c` (`zhybrid`, `zghy`, `zhad`, `zproxy`) — small standalone grid-dumpers. Each has a self-contained `main()` and prints TSV; they exist for visualization and benchmarking, not as the user-facing path.

When adding a new method, follow the pattern of method B: implement the kernel in `ghy.{c,h}`, write a thin `zero_count_<name>` wrapper in `main.c`, and add an argp flag — this preserves the original `zzz` CLI rather than spawning yet another binary.

## Documentation layout

- `doc/notes/` — synthesis markdowns (KaTeX): `admissibility-and-rigor-gap.md`, `rigor-bound-b.md`, `method-a-euler-product.md`, `damping-spike-shape.md`, etc. Read these before changing the damping or making rigor claims.
- `doc/ghy/` — Wolfram scripts and reference TSVs (committed): kernel reference, A-vs-B scans, FWHM verification.
- `doc/heuristic/` — earlier method-A development notes.
- `doc/refs/`, `doc/pell/` — gitignored, external copyright; do not commit.

## Conventions

- Precision flows through `slong PREC` (default 256 bits in `zzz`); never hardcode a bit count inside `ghy.c`.
- All complex arithmetic uses `acb_t`; reals use `arb_t`. Init/clear pairs are mandatory — there is no automatic cleanup.
- Counting functions return `arb_t` (real); only the inner factor is `acb_t`.
- `n_compute_primes(k)` must be called before any prime sum; method B's `X = n_nth_prime(k)` ties the cutoff to the same `k` the user passes.
