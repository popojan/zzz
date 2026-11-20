// loop.h -- self-paving zeros<->primes bootstrap (zzz --loop).
//
// Demonstrates the explicit formula running in a closed loop with NO zeta
// evaluation and NO primality test: from a finite seed of zero ordinates it
// detects primes (Chebyshev psi' signal) and locates more zeros (method-B
// partial-Euler sum over the discovered primes), alternating until it is
// interrupted.  It commits x once the zeros reach ~2*pi*x, so the forward
// margin stays near kappa = 6.3, well above kmin.
// Pure double precision (heights are modest).  See doc/notes/loop-margins.md
// (measured margins) and doc/notes/zeros-primes-bootstrap.md (mechanism).
//
// Streams newly discovered primes to stdout, a trajectory to stderr, and
// checkpoints full state to a file so --resume can continue (Ctrl+C is caught
// and checkpoints before exit).

#ifndef LOOP_H
#define LOOP_H

typedef struct {
    const char *state_path;  // checkpoint file [zzz-loop.state]
    int resume;              // resume from state_path instead of the seed
    long max_iters;          // stop after this many more iterations (0 = until stalled or interrupted)
    long nmax_zeros;         // safety cap on the zero list [2000000]
    double kmin;             // kappa margin the forward step must keep [5.0]
    int kmin_set;            // user passed --loop-kmin (overrides resumed value)
    long seed_n;             // use only the first seed_n seed zeros (0 = all)
    long batch;              // zeros per re-detect (smooth streaming) [2000; 0=all]
    int fresh;               // force a fresh start, ignoring any checkpoint
    double boxc;             // forward box-smoothing window w = boxc*gap(t) [0.375; 0 = plain B]
    int boxc_set;            // user passed --loop-boxc (overrides the value a checkpoint records)
    const char *declog;      // log every backward decision to this TSV (see doc/ghy/loop-margins.py)
    int verbose;
} loop_opts;

void loop_opts_default(loop_opts *o);
int  loop_run(const loop_opts *o);

#endif
