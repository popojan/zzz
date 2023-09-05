#include "arb/acb.h"
#include <argp.h>

#define CUBIC 1

void coord(arb_ptr out, arb_srcptr q, arb_srcptr t, slong PREC) {
    arb_t pi;
    arb_t pi_div_log_q;
    arb_t x;
    arb_t a;
    arb_t b;

    arb_init(pi);
    arb_init(pi_div_log_q);
    arb_init(x);
    arb_init(a);
    arb_init(b);

    // PI / Log[q]
    arb_const_pi(pi, PREC);
    arb_log(pi_div_log_q, q, PREC);
    arb_inv(pi_div_log_q, pi_div_log_q, PREC);
    arb_mul(pi_div_log_q, pi_div_log_q, pi, PREC);

    // -t + PI / Log[q]
    arb_sub(x, pi_div_log_q, t, PREC);

    //Mod[-t + PI / Log[q], 2 PI / Log[q]]
    arb_add(a, pi_div_log_q, pi_div_log_q, PREC);
    arb_div(b, x, a, PREC);
    arb_floor(b, b, PREC);
    arb_mul(b, b, a, PREC);
    arb_sub(x, x, b, PREC);

    //% - PI/Log[q]
    arb_sub(x, x, pi_div_log_q, PREC);
    arb_set(out, x);

    arb_clear(pi);
    arb_clear(pi_div_log_q);
    arb_clear(x);
    arb_clear(a);
    arb_clear(b);
}

void wave2(arb_ptr out, arb_srcptr q, arb_srcptr t, slong PREC) {
    arb_t pi;
    arb_t sqrt_q;
    arb_t log_q;
    arb_t x;
    arb_t a;
    arb_t b;

    arb_init(pi);
    arb_init(x);
    arb_init(a);
    arb_init(b);
    arb_init(sqrt_q);
    arb_init(log_q);

    arb_const_pi(pi, PREC);
    arb_log(log_q, q, PREC);
    arb_sqrt(sqrt_q, q, PREC);


    //4 sqrt[q] log[q] t
    arb_mul_ui(a, t, 4,PREC);
    arb_mul(b, log_q, sqrt_q, PREC);
    arb_mul(a, a, b, PREC);

    //PI (sqrt[q]-1)^2
    arb_one(b);
    arb_sub(b, sqrt_q, b, PREC);
    arb_mul(b, b, b, PREC);
    arb_mul(b, b, pi, PREC);

    //(%1 + %2)
    arb_add(a, a, b, PREC);

    //q Log[q]^2
    arb_mul(b, q, log_q, PREC);
    arb_mul(b, b, log_q, PREC);

    //sqrt(%/%)
    arb_div(a, a, b, PREC);
    arb_sqrt(a, a, PREC);

    // (sqrt[PI] sqrt[q] sqrt(%/%) - 2 t) Log[q]
    arb_sqrt(b, pi, PREC);
    arb_mul(b, b, sqrt_q, PREC);
    arb_mul(b, b, a, PREC);
    arb_sub(b, b, t, PREC);
    arb_sub(b, b, t, PREC);
    arb_mul(a, b, log_q, PREC);

    arb_add(b, pi, pi, PREC);
    arb_div(a, a, b, PREC);

    arb_div_ui(b, sqrt_q, 2, PREC);
    arb_sub(a, a, b, PREC);
    arb_set_d(b, 0.5);
    arb_add(a, a, b, PREC);

    arb_set(out, a);

    arb_clear(pi);
    arb_clear(x);
    arb_clear(a);
    arb_clear(b);
    arb_clear(sqrt_q);
    arb_clear(log_q);
}


void wave3(arb_ptr out, arb_srcptr q, arb_srcptr t, slong PREC) {
    arb_t pi;
    arb_t sqrt_q;
    arb_t sqrt_q3;
    arb_t log_q;
    arb_t x;
    arb_t y;
    arb_t a;
    arb_t b;
    arb_t c;
    arb_t d;

    arb_init(pi);
    arb_init(x);
    arb_init(y);
    arb_init(a);
    arb_init(b);
    arb_init(c);
    arb_init(d);
    arb_init(sqrt_q);
    arb_init(sqrt_q3);
    arb_init(log_q);

    arb_sqrt(sqrt_q, q, PREC);

    arb_mul(sqrt_q3, q, q, PREC);
    arb_mul(sqrt_q3, sqrt_q3, q, PREC);
    arb_sqrt(sqrt_q3, sqrt_q3, PREC);

    arb_log(log_q, q, PREC);
    arb_const_pi(pi, PREC);

    arb_mul(a, t, log_q, PREC);
    arb_neg(a, a);
    arb_add(a, a, pi, PREC);

    arb_add(b, sqrt_q3, sqrt_q3, PREC);
    arb_mul(a, a, b, PREC);

    arb_add(b, q, q, PREC);
    arb_add(b, b, q, PREC);
    arb_mul(b, b, pi, PREC);
    arb_neg(b, b);
    arb_add(b,b, pi, PREC);
    arb_add(a, a, b, PREC);

    arb_mul(b, t, log_q, PREC);
    arb_mul_ui(c, sqrt_q, 8, PREC);
    arb_mul(b, b, c, PREC);
    arb_mul(a, a, b, PREC);
    //a == -8 sqrt[q]...

    arb_mul(b, pi, pi, PREC);
    arb_one(c);
    arb_sub(c, sqrt_q, c, PREC);
    arb_mul(d, c, c, PREC);
    arb_mul(d, d, c, PREC);
    arb_mul(b, b, d, PREC);
    //b == pi^2 ( sqrt[q]-1)^3

    arb_mul_ui(c, sqrt_q, 5,PREC);
    arb_add_ui(c, c, 3, PREC);
    arb_mul(b, b, c, PREC);
    arb_sub(a, b, a, PREC);
    //a == pi^2..... (...Log[q]), numerator

    arb_mul(b, q, q, PREC);
    arb_div(a, a, b, PREC);
    arb_sqrt(a, a, PREC);
    //a == sqrt(numerator / q^2)

    arb_mul(c, t, log_q, PREC);
    arb_mul_ui(b, c, 4,PREC);
    arb_add(a, a, b, PREC);

    arb_mul_ui(b, q, 3, PREC);
    arb_one(c);
    arb_sub(b, b, c, PREC);
    arb_div(b, b, sqrt_q3, PREC);
    arb_sub_ui(b, b, 2, PREC);
    arb_mul(b, b, pi, PREC);
    arb_add(a, a, b, PREC);
    //a==a

    arb_root_ui(pi, pi, 3, PREC);
    arb_root_ui(b, a, 3, PREC);

    arb_one(d);
    arb_sub(d, sqrt_q, d, PREC);
    arb_mul(d, d, sqrt_q, PREC);
    arb_mul(d, d, b, PREC);
    arb_mul(d, d, pi, PREC);

    arb_one(x);
    arb_sub(x, q, x, PREC);
    arb_mul(x, x, pi, PREC);
    arb_mul(x, x, pi, PREC);
    arb_sub(x, x, d, PREC);
    arb_mul(d, q, b, PREC);
    arb_mul(d, d, b, PREC);
    arb_sub(x, x, d, PREC);

    arb_one(d);
    arb_add(d, sqrt_q, d, PREC);
    arb_mul(d, d, sqrt_q, PREC);
    arb_mul(d, d, b, PREC);
    arb_mul(d, d, pi, PREC);

    arb_one(y);
    arb_sub(y, q, y, PREC);
    arb_mul(pi, pi, pi, PREC);
    arb_mul(y, y, pi, PREC);
    arb_add(y, y, d, PREC);
    arb_mul(d, q, b, PREC);
    arb_mul(d, d, b, PREC);
    arb_sub(y, y, d, PREC);

    arb_mul(x, x, y, PREC);
    arb_set_d(y, -0.375);
    arb_mul(x, x, y, PREC);
    arb_mul(y, pi, b, PREC);
    arb_mul(y, y, b, PREC);
    //arb_set_d(b, 2.5);
    //arb_pow(b, q, b, PREC);
    arb_pow_ui(b, q, 5, PREC);
    arb_sqrt(b, b, PREC);
    arb_mul(y, y, b, PREC);
    arb_div(x, x, y, PREC);
    arb_set(out, x);

    arb_clear(pi);
    arb_clear(x);
    arb_clear(y);
    arb_clear(a);
    arb_clear(b);
    arb_clear(c);
    arb_clear(d);
    arb_clear(sqrt_q);
    arb_clear(sqrt_q3);
    arb_clear(log_q);
}

void zeta(acb_ptr out, arb_srcptr t, slong ZETA_PREC) {
    acb_t s;
    acb_t z;
    acb_init(s);
    acb_init(z);
    arb_t a;
    arb_init(a);
    arb_set_d(a, 0.5);
    acb_set_arb_arb(s, a, t);
    acb_zeta(z, s, ZETA_PREC);
    acb_set(out, z);

    acb_clear(s);
    acb_clear(z);
}

void nt(arb_ptr out, arb_srcptr t, slong PREC) {
    arb_t x;
    arb_t a;
    arb_t b;
    arb_init(x);
    arb_init(a);
    arb_init(b);

    arb_set(x, t);
    arb_const_pi(b, PREC);
    arb_div_ui(x, x, 2, PREC);
    arb_div(x, x, b, PREC);
    arb_set(a, x);
    arb_const_e(b, PREC);
    arb_div(a, a, b, PREC);
    arb_log(a, a, PREC);
    arb_mul(x, a, x, PREC);
    arb_set_d(a, 0.875);
    arb_add(x, x, a, PREC);
    arb_set(out, x);

    arb_clear(x);
    arb_clear(a);
    arb_clear(b);
}

void nt_inv(arb_ptr out, arb_srcptr m, slong PREC) {
    arb_t x;
    arb_t a;
    arb_t b;
    arb_t nom;
    arb_t den;
    arb_init(x);
    arb_init(a);
    arb_init(b);
    arb_init(den);
    arb_init(nom);

    arb_mul_ui(x, m, 8, PREC);
    arb_sub_ui(x, x, 11, PREC);

    arb_set(nom, x);
    arb_const_pi(b, PREC);
    arb_mul(nom, nom, b, PREC);

    arb_set(den, x);
    arb_const_e(b, PREC);
    arb_div(den, den, b, PREC);
    arb_div_ui(den, den, 8, PREC);
    arb_lambertw( den, den, 0, PREC);
    arb_mul_ui(den, den, 4, PREC);

    arb_set(x, nom);
    arb_div(x, x, den, PREC);

    arb_set(out, x);

    arb_clear(x);
    arb_clear(a);
    arb_clear(b);
    arb_clear(den);
    arb_clear(nom);
}

void zero_count_approx(arb_ptr out, arb_srcptr t, slong k, slong PREC) {
    arb_t u;
    arb_t w;
    arb_t q;
    arb_t x;
    arb_t y;
    arb_t z;
    arb_t v;
    arb_t s;
    arb_t att;
    int neg;

    arb_init(u);
    arb_init(w);
    arb_init(q);
    arb_init(y);
    arb_init(z);
    arb_init(v);
    arb_init(s);
    arb_init(att);
    arb_init(x);

    nt(u, t, PREC);

#if CUBIC == 1
    //v = 1/sqrt(2)
    arb_set_d(v, 2.0);
    arb_sqrt(v, v, PREC);
    arb_inv(v, v, PREC);
#endif

    for(int i = 1; i <= k; ++i) {
        arb_set_d(q, n_nth_prime(i));

        arb_div(att, t, q, PREC);
        arb_sqrt(att,att,PREC);
        arb_neg(att, att);
        arb_exp(att, att, PREC);;
        arb_one(w);
        arb_sub(att, w, att, PREC);

        coord(x, q, t, PREC);
        arb_abs(w, x);

        wave2(z, q, w, PREC);

        neg = arb_is_negative(x);
        if(neg) {
            arb_neg(z, z);
        }
#if CUBIC == 1
        arb_set_d(y, i);
        arb_pow(y, y, v, PREC);
        arb_mul_ui(y, y, 3, PREC);
        arb_inv(x, y, PREC);

        arb_set_d(y, 0.5);
        arb_sub(y, y, x, PREC);
        arb_mul(z, z, y, PREC);
        arb_mul(z, z, att, PREC);
        arb_add(u, u, z, PREC);

        wave3(z, q, w, PREC);
        if(neg) {
            arb_neg(z, z);
        }
        arb_set_d(y, 0.5);
        arb_add(y, y, x, PREC);
        arb_mul(z, z, y, PREC);
#endif
        arb_mul(z, z, att, PREC);
        arb_add(u, u, z, PREC);
    }
    arb_set(out, u);

    arb_clear(u);
    arb_clear(w);
    arb_clear(q);
    arb_clear(y);
    arb_clear(z);
    arb_clear(v);
    arb_clear(s);
    arb_clear(att);
    arb_clear(x);
}

const char *argp_program_version = "zzz 0";
const char *argp_program_bug_address = "<>";
static char doc[] = "fast approximation of large Riemann zeta zeros";
static char args_doc[] = "N [offset] [count]";
static struct argp_option options[] = {
        { "k", 'k', "K", 0, "use first k primes for zero counting function approximation [default 100]"},
        { "evaluate", 'e', 0, 0, "evaluate Riemann zeta function value at the approximate zero location"},
        { "tolerance", 't', "TOL", 0, "tolerance for bisection [default 1e-6]"},
        { "window", 'w', "WIN", 0, "initial span around Lambert W asymptotic zero location +- WIN [default 1.5]"},
        { "precision", 'p', "PREC", 0, "arb precision for counting function approximation [default 256]"},
        { "zeta-prec", 'z', "ZETA_PREC", 0, "arb precision for zeta evaluation [default 64]"},
        { "digits", 'd', "DIGITS", 0, "extra digits for number formatting [default 6]"},
        { "verbose", 'v', 0, 0, "verbose progress output"},
        { "debug", 'g', 0, 0, "debug counting function from <N> to <N+offset> in <count> steps"},
        { 0 }
};

struct arguments {
    slong k;
    slong eval;
    double step0;
    double w0;
    slong PREC;
    slong ZETA_PREC;
    slong DIGITS;
    slong verbose;
    slong debug;
};

static error_t parse_opt(int key, char *arg, struct argp_state *state) {
    struct arguments *arguments = state->input;

    switch (key) {
        case 'k': arguments->k = atol(arg); break;
        case 'e': arguments->eval = 1; break;
        case 't': arguments->step0 = atof(arg); break;
        case 'w': arguments->w0 = atof(arg); break;
        case 'p': arguments->PREC = atol(arg); break;
        case 'z': arguments->ZETA_PREC = atol(arg); break;
        case 'd': arguments->DIGITS = atol(arg); break;
        case 'v': arguments->verbose = 1; break;
        case 'g': arguments->debug = 1; break;
        case ARGP_KEY_ARG: return 0;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };

void iter_print(arb_srcptr lo_t, arb_srcptr lo, arb_srcptr hi_t, arb_srcptr hi, slong DIGITS, slong verbose) {

    if(verbose == 0) return;

    flint_fprintf(stderr, "lower x N(x)\t");
    arb_fprintd(stderr, lo, DIGITS);
    flint_fprintf(stderr, "\t");
    arb_fprintd(stderr, lo_t, DIGITS);
    flint_fprintf(stderr, "\n");

    flint_fprintf(stderr, "upper x N(x)\t");
    arb_fprintd(stderr, hi, DIGITS);
    flint_fprintf(stderr, "\t");
    arb_fprintd(stderr, hi_t, DIGITS);
    flint_fprintf(stderr, "\n");
}

int main(int argc, char *argv[])
{

    if(argc == 1) {
        flint_fprintf(stderr, "zzz --usage for help\n");
        return 1;
    }

    struct arguments arguments;

    arguments.k = 100;
    arguments.DIGITS = 6;
    arguments.PREC = 256;
    arguments.ZETA_PREC = 64;
    arguments.step0 = 0.000001;
    arguments.w0 = 1.5;
    arguments.eval = 0;
    arguments.verbose = 0;
    arguments.debug = 0;

    int arg_index = 1;
    argp_parse(&argp, argc, argv, ARGP_NO_ARGS, &arg_index, &arguments);

    arb_t m;
    arb_t u;
    arb_t m0;
    arb_t mm;
    arb_t m_lo;
    arb_t m_hi;
    arb_t tt;
    arb_t step;

    arb_init(m);
    arb_init(u);
    arb_init(m0);
    arb_init(mm);
    arb_init(m_lo);
    arb_init(m_hi);
    arb_init(tt);
    arb_init(step);

    arb_one(m0);
    arb_zero(m);
    if(argc > arg_index) {
        arb_set_str(m0, argv[arg_index], arguments.PREC);
    }
    if(argc > arg_index+1) {
        arb_set_str(m, argv[arg_index+1], arguments.PREC);
    }

    slong count = 1;
    if(argc > arg_index+2) {
        count = atol(argv[arg_index+2]);
    }

    if(arguments.k > 0) {
        n_compute_primes(arguments.k);
    }

    if(arguments.debug) {
        arb_t tx;
        arb_t zc;
        arb_t stp;

        arb_init(tx);
        arb_init(zc);
        arb_init(stp);
        arb_zero(tx);
        arb_set(stp, m);
        arb_div_ui(stp, stp, count, arguments.PREC);
        arb_set(tx, m0);
        arb_add(m, m, m0, arguments.PREC);

        for (; arb_lt(tx, m); arb_add(tx, tx, stp, arguments.PREC)) {
            zero_count_approx(zc, tx, arguments.k, arguments.PREC);
            arf_printd(&tx->mid, arguments.DIGITS);
            flint_printf("\t");
            arf_printd(&zc->mid, arguments.DIGITS);
            flint_printf("\n");
        }
    } else {
        arb_add(m0, m0, m, arguments.PREC);

        for (slong ord = 0; ord < count || count < 0; ++ord) {

            nt_inv(tt, m0, arguments.PREC);

            //calc required digits
            arb_const_log10(u, arguments.PREC);
            arb_log(m, tt, arguments.PREC);
            arb_div(u, m, u, arguments.PREC);
            arb_ceil(u, u, arguments.PREC);

            slong digits = arguments.DIGITS;
            digits += arf_get_si(&u->mid, 0);

            if (arguments.verbose) {
                flint_fprintf(stderr, "asymptotic zero location = ");
                arb_fprintd(stderr, tt, digits);
                flint_fprintf(stderr, "\n");
            }

            arb_set_d(m, -0.5);
            arb_add(m, m0, m, arguments.PREC);
            arb_set(m_lo, m);
            arb_set(m_hi, m);

            arb_t lo;
            arb_t hi;
            arb_t lo_t;
            arb_t hi_t;

            arb_init(lo);
            arb_init(hi);
            arb_init(lo_t);
            arb_init(hi_t);

            // tt / 2 PI
            arb_const_pi(lo_t, arguments.PREC);
            arb_div(hi, tt, lo_t, arguments.PREC);
            arb_set_d(lo_t, 0.5);
            arb_mul(hi, hi, lo_t, arguments.PREC);

            acb_t zz;
            acb_init(zz);

            arb_pos_inf(lo_t);
            arb_neg_inf(hi_t);

            arb_pos_inf(lo);
            arb_neg_inf(hi);

            arb_set_d(step, arguments.w0);

            arb_sub(lo_t, tt, step, arguments.PREC);
            arb_add(hi_t, tt, step, arguments.PREC);

            zero_count_approx(lo, lo_t, arguments.k, arguments.PREC);
            zero_count_approx(hi, hi_t, arguments.k, arguments.PREC);

            if (arb_gt(lo, m) || arb_lt(hi, m) || arb_lt(hi, lo)) {
                iter_print(lo_t, lo, hi_t, hi, digits, arguments.verbose);
                flint_fprintf(stderr, "please increase the window\n");
            } else {

                arb_t mid_t;
                arb_t mid;
                arb_init(mid);
                arb_init(mid_t);

                while (1) {

                    iter_print(lo_t, lo, hi_t, hi, digits, arguments.verbose);

                    arb_add(mid_t, lo_t, hi_t, arguments.PREC);
                    arb_set_d(mm, 0.5);
                    arb_mul(mid_t, mid_t, mm, arguments.PREC);

                    zero_count_approx(mid, mid_t, arguments.k, arguments.PREC);

                    arb_zero(mm);
                    if (arb_gt(mid, m)) {
                        arb_set(hi_t, mid_t);
                        arb_set(hi, mid);
                    } else {
                        arb_set(lo_t, mid_t);
                        arb_set(lo, mid);
                    }

                    arb_sub(mm, hi_t, lo_t, arguments.PREC);
                    arb_set_d(u, arguments.step0);
                    if (arb_lt(mm, u)) break;

                }

                iter_print(lo_t, lo, hi_t, hi, digits, arguments.verbose);

                if (arguments.verbose) {
                    flint_fprintf(stderr, "argument s = \t");
                    arb_set_d(lo_t, 0.5);
                    acb_set_arb_arb(zz, lo_t, mid_t);
                    acb_fprintd(stderr, zz, digits);
                    flint_fprintf(stderr, "\n");
                }
                if (arguments.eval > 0) {
                    flint_fprintf(stderr, "value    z = \t");
                    zeta(zz, mid_t, arguments.ZETA_PREC);
                    acb_fprintd(stderr, zz, digits);
                    flint_fprintf(stderr, "\n");
                }
                arb_set_d(lo_t, 0.005);

                arf_fprintd(stdout, &mid_t->mid, digits);
                flint_fprintf(stdout, "\n");
                fflush(stdout);

                arb_one(m);
                arb_add(m0, m0, m, arguments.PREC);

                arb_clear(mid);
                arb_clear(mid_t);
            }

            arb_clear(lo);
            arb_clear(hi);
            arb_clear(lo_t);
            arb_clear(hi_t);
            acb_clear(zz);

        }
    }

    arb_clear(m);
    arb_clear(u);
    arb_clear(m0);
    arb_clear(mm);
    arb_clear(m_lo);
    arb_clear(m_hi);
    arb_clear(tt);
    arb_clear(step);

    n_cleanup_primes();
    flint_cleanup();

    return 0;
}
