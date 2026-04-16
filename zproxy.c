// zproxy.c — dump the finite damped Euler proxy Z_k(s, T) on the critical line.
//
// Z_k(s, T) := prod_{p<=p_k} (1 - p^{-s})^{-att(p, T)}
// att(p, T) =  1 - exp(-sqrt(T/p))   (matches main.c::zero_count_exact)
//
// At s = 1/2 + iT, F(T) = N_0(T) + arg Z_k(1/2+iT, T) / pi is exactly the
// counting function main.c bisects. Use this binary to visualize the
// magnitude/argument of Z_k and to compare its behavior against true zeta
// zeros (zeros of ζ should correspond to unit jumps of F).
//
// Usage: zproxy <T0> <T1> <N> [k]
//   prints N+1 rows: T  log|Z_k|  arg(Z_k)  F(T)

#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/ulong_extras.h>
#include <stdio.h>
#include <stdlib.h>

static void compute_att(arb_ptr att, arb_srcptr t, arb_srcptr q, slong PREC) {
    // att = 1 - exp(-sqrt(t/q))
    arb_div(att, t, q, PREC);
    arb_sqrt(att, att, PREC);
    arb_neg(att, att);
    arb_exp(att, att, PREC);
    arb_sub_ui(att, att, 1, PREC);
    arb_neg(att, att);
}

// N_0(T) = (T/2pi) log(T/2pi e) + 7/8
static void n0_smooth(arb_ptr out, arb_srcptr t, slong PREC) {
    arb_t x, a, b;
    arb_init(x); arb_init(a); arb_init(b);

    arb_set(x, t);
    arb_const_pi(b, PREC);
    arb_div_ui(x, x, 2, PREC);
    arb_div(x, x, b, PREC);       // x = t / (2 pi)
    arb_set(a, x);
    arb_const_e(b, PREC);
    arb_div(a, a, b, PREC);       // a = t / (2 pi e)
    arb_log(a, a, PREC);
    arb_mul(x, a, x, PREC);       // x = (t/2pi) log(t/2pi e)
    arb_set_d(a, 0.875);
    arb_add(x, x, a, PREC);       // + 7/8
    arb_set(out, x);

    arb_clear(x); arb_clear(a); arb_clear(b);
}

// logz = log Z_k(1/2 + iT, T) = -sum_{p<=p_k} att(p,T) * log(1 - p^{-1/2 - iT})
static void log_zk(acb_ptr logz, arb_srcptr t, slong k, slong PREC) {
    acb_t term, one_minus_ps, ps, minus_s;
    arb_t half_neg, neg_t, q, att;

    acb_init(term);
    acb_init(one_minus_ps);
    acb_init(ps);
    acb_init(minus_s);
    arb_init(half_neg);
    arb_init(neg_t);
    arb_init(q);
    arb_init(att);

    // minus_s = -1/2 - iT (the exponent for p^{-s} when s = 1/2 + iT)
    arb_set_d(half_neg, -0.5);
    arb_neg(neg_t, t);
    acb_set_arb_arb(minus_s, half_neg, neg_t);

    acb_zero(logz);

    for (slong i = 1; i <= k; ++i) {
        arb_set_ui(q, n_nth_prime(i));
        compute_att(att, t, q, PREC);

        // p^{-s}
        acb_set_arb(ps, q);
        acb_pow(ps, ps, minus_s, PREC);

        // 1 - p^{-s}
        acb_neg(one_minus_ps, ps);
        acb_add_ui(one_minus_ps, one_minus_ps, 1, PREC);

        // att * log(1 - p^{-s})
        acb_log(term, one_minus_ps, PREC);
        acb_mul_arb(term, term, att, PREC);

        // log Z_k -= att * log(1 - p^{-s})
        acb_sub(logz, logz, term, PREC);
    }

    acb_clear(term);
    acb_clear(one_minus_ps);
    acb_clear(ps);
    acb_clear(minus_s);
    arb_clear(half_neg);
    arb_clear(neg_t);
    arb_clear(q);
    arb_clear(att);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr,
            "usage: zproxy <T0> <T1> <N> [k]\n"
            "  prints N+1 rows (TSV): T  log|Z_k|  arg(Z_k)  F(T)\n"
            "  Z_k(s,T) = prod_{p<=p_k}(1-p^{-s})^{-att(p,T)} at s = 1/2 + iT\n"
            "  F(T) = N_0(T) + arg Z_k(1/2+iT, T) / pi\n"
            "  (F matches zzz's counting-function approximation)\n");
        return 1;
    }

    double t0 = atof(argv[1]);
    double t1 = atof(argv[2]);
    slong N = atol(argv[3]);
    slong k = (argc > 4) ? atol(argv[4]) : 100;
    const slong PREC = 128;

    n_compute_primes(k);

    arb_t t, step, n0, pi, arg_zk, log_abs, f;
    acb_t logz;

    arb_init(t); arb_init(step); arb_init(n0); arb_init(pi);
    arb_init(arg_zk); arb_init(log_abs); arb_init(f);
    acb_init(logz);

    arb_const_pi(pi, PREC);
    arb_set_d(t, t0);
    arb_set_d(step, (t1 - t0) / (double)N);

    for (slong i = 0; i <= N; ++i) {
        log_zk(logz, t, k, PREC);

        arb_set(log_abs, acb_realref(logz));
        arb_set(arg_zk, acb_imagref(logz));

        n0_smooth(n0, t, PREC);
        arb_div(f, arg_zk, pi, PREC);
        arb_add(f, f, n0, PREC);

        arf_printd(arb_midref(t), 12);       printf("\t");
        arf_printd(arb_midref(log_abs), 10); printf("\t");
        arf_printd(arb_midref(arg_zk), 10);  printf("\t");
        arf_printd(arb_midref(f), 10);       printf("\n");

        arb_add(t, t, step, PREC);
    }

    arb_clear(t); arb_clear(step); arb_clear(n0); arb_clear(pi);
    arb_clear(arg_zk); arb_clear(log_abs); arb_clear(f);
    acb_clear(logz);

    n_cleanup_primes();
    flint_cleanup();
    return 0;
}
