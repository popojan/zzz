// loop.h -- self-paving zeros<->primes bootstrap (zzz --loop).
//
// Demonstrates the explicit formula running in a closed loop with NO zeta
// evaluation and NO primality test: from a finite seed of zero ordinates it
// detects primes (Chebyshev psi' signal) and locates more zeros (method-B
// partial-Euler sum over the discovered primes), alternating until the kappa
// margin erodes (it stalls) or it is interrupted.  Pure double precision
// (heights are modest).  See doc/notes/zeros-primes-bootstrap.md.
//
// Streams newly discovered primes to stdout, a trajectory to stderr, and
// checkpoints full state to a file so --resume can continue (Ctrl+C is caught
// and checkpoints before exit).

#ifndef LOOP_H
#define LOOP_H

typedef struct {
    const char *state_path;  // checkpoint file [zzz-loop.state]
    int resume;              // resume from state_path instead of the seed
    long max_iters;          // 0 = until stalled or interrupted
    long nmax_zeros;         // safety cap on the zero list [200000]
    double kmin;             // (starting) kappa margin for the forward step [5.0]
    int kmin_set;            // user passed --loop-kmin (overrides resumed value)
    double kmin_floor;       // anneal stops here (>~pi keeps plain-B correct) [3.5]
    int anneal;              // on stall, lower kmin and keep climbing [1]
    long seed_n;             // use only the first seed_n seed zeros (0 = all)
    int contrast;            // use the fit-free local-contrast detector (no Li)
    int verbose;
} loop_opts;

void loop_opts_default(loop_opts *o);
int  loop_run(const loop_opts *o);

#endif
