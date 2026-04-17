# `doc/` — index

Documentation, reference material, scripts and plots live under six subdirectories.

## `notes/` — synthesis markdowns (KaTeX)

The analytical writeups produced during the current rigor-via-GHY work. Read these first if you want the story.

| file | one-line summary |
|---|---|
| `admissibility-and-rigor-gap.md` | Weil admissibility, what is provable under RH, A vs B vs C in plain language |
| `damping-and-hybrid-rigor.md` | original plan document — damping analysis, Euler-product projection, path to rigor |
| `rigor-backlog.md` | phase-by-phase roadmap for bounding \|F − N\| < ½ |
| `rigor-bound-b.md` | method-B error bound: GHY Thm 1 + Goldston 1987, empirical confrontation |
| `prime-count-scaling.md` | heuristic $k \sim (\log T)^2$ complexity note |

## `refs/` — external papers (PDFs by other authors) — **gitignored**

Local-only mirror of the papers cited in `notes/`. Not committed (copyright of third parties). Populate by hand from the citations in the references section of each note.

| local filename | citation |
|---|---|
| `hybridformula.pdf` | Gonek, Hughes, Young, *A hybrid Euler–Hadamard product for the Riemann zeta function*, Duke Math. J. **136** (2007) |
| `09e_guinand_explicit_fml.pdf` | Guinand's explicit-formula derivation (Weil admissibility context) |
| `rnoti-p36.pdf` | AMS Notices short article |
| `pdfs.jsp.pdf` | miscellaneous reference |

## `heuristic/` — original zzz-A development (2023 and earlier)

Scripts, plots and PDFs that underpin the *current* shipped `zzz` default mode (method A, heuristic damping $1 - e^{-\sqrt{T/p}}$).

- `2023-03-20.wls`, `2023-03-27_convergence.wls` — original Wolfram notebooks
- `test-function-search.wls` — admissibility / attenuation FT analysis
- `attenuation-{comparison,universality,transition,zoom}.pdf` — damping-shape plots
- `truncation-error-vs-k.pdf`, `error-vs-ratio.pdf`, `fitted.{pdf,png}` — convergence empirics
- `phase4{,b,c}-numerics/scaling/density.wls` — Phase-4 error-scan scripts
- `waves.png`, `convergence.png`, `counting.png`, `errors.png` — README figures
- `psi-*.png` — Chebyshev ψ reconstruction plots

## `ghy/` — current GHY rigor work (this branch)

| file | role |
|---|---|
| `ghy-kernel.wls` | derive & plot $U(z) \to E_1(z)$ (GHY mollifier limit) |
| `ghy-kernel-ref.tsv` | reference values for C unit tests of `acb_hypgeom_expint` |
| `ghy-hadamard-{window,jump}.pdf` | kernel visualisation near a zero |
| `ab-scan.wls` | A-vs-B accuracy scan on first five ordinal decades |
| `ab-scan-{A,B,ratio}.pdf` | results of that scan |
| `abc-compare.wls` | A/B/C boundary-behaviour benchmark (needs zero list) |
| `first-20-zeros.txt` | Odlyzko ordinates for seeding `zhybrid` |

All `.wls` scripts assume you run `wolframscript -file doc/ghy/<script>.wls` from the project root.

## `pell/` — unrelated Pell-branch material — **gitignored**

Documents and reference PDFs from the `quadregulator` branch work (compact representations of Pell solutions). Not committed on this branch.

- `pell-compact-reconstruction.md`
- `pell_{convergence,doubling,nearsquare}.pdf`
- `R-Williams.pdf`, `Solving the Pell Equation.pdf`, `real_naf_final.pdf`
- `S0025-5718-03-01518-7.pdf` (Van der Poorten NUCOMP-with-Distance)

None of this is used by `zzz`.

## `tbd/` — historical kitchen sink (untouched)

Older ψ-reconstruction and zero-counting-function plots kept as-is for reference.

---

## Build artifacts

Plots that ship *into* this directory are produced by running the `.wls` script next to them. TSV reference files likewise. Everything else is either input (external PDFs) or hand-written markdown.
