#include "acb.h"
#include <argp.h>

#define CUBIC 1

void coord(arb_ptr out, arb_t q, arb_t t, slong PREC) {
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
    arb_set(x, pi_div_log_q);
    arb_sub(x, x, t, PREC);

    //Mod[-t + PI / Log[q], 2 PI / Log[q]]
    arb_set_d(b, 2.0);
    arb_set(a, pi_div_log_q);
    arb_mul(a, a, b, PREC);
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

void wave2(arb_ptr out, arb_t q, arb_t t, slong PREC) {
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
    arb_set_d(a, 4.0);
    arb_mul(a, a, t, PREC);
    arb_mul(a, a, log_q, PREC);
    arb_mul(a, a, sqrt_q, PREC);

    //PI (sqrt[q]-1)^2
    arb_set_d(b, -1.0);
    arb_add(b, b, sqrt_q, PREC);
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

    arb_set_d(b, 0.5);
    arb_mul(b, b, sqrt_q, PREC);
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


void wave3(arb_ptr out, arb_t q, arb_t t, slong PREC) {
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
    arb_set_d(c, 8.0);
    arb_mul(c, c, sqrt_q, PREC);
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

    arb_set_d(c, 5.0);
    arb_mul(c, c, sqrt_q, PREC);
    arb_set_d(d, 3.0);
    arb_add(c, c, d, PREC);
    arb_mul(b, b, c, PREC);
    arb_sub(a, b, a, PREC);
    //a == pi^2..... (...Log[q]), numerator

    arb_mul(b, q, q, PREC);
    arb_div(a, a, b, PREC);
    arb_sqrt(a, a, PREC);
    //a == sqrt(numerator / q^2)

    arb_set_d(b, 4.0);
    arb_mul(c, t, log_q, PREC);
    arb_mul(b, b, c, PREC);
    arb_add(a, a, b, PREC);

    arb_set_d(b, 3.0);
    arb_mul(b, b, q, PREC);
    arb_one(c);
    arb_sub(b, b, c, PREC);
    arb_div(b, b, sqrt_q3, PREC);
    arb_set_d(c, 2.0);
    arb_sub(b, b, c, PREC);
    arb_mul(b, b, pi, PREC);
    arb_add(a, a, b, PREC);
    //a==a

    arb_set_d(c, 1.0);
    arb_set_d(d, 3.0);
    arb_div(c, c, d, PREC);

    arb_pow(pi, pi, c, PREC);
    arb_pow(b, a, c, PREC);

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
    arb_mul(y, y, pi, PREC);
    arb_mul(y, y, pi, PREC);
    arb_add(y, y, d, PREC);
    arb_mul(d, q, b, PREC);
    arb_mul(d, d, b, PREC);
    arb_sub(y, y, d, PREC);

    arb_mul(x, x, y, PREC);
    arb_set_d(y, -3.0);
    arb_mul(x, x, y, PREC);
    arb_set_d(y, 8.0);
    arb_mul(y, y, pi, PREC);
    arb_mul(y, y, pi, PREC);
    arb_mul(y, y, b, PREC);
    arb_mul(y, y, b, PREC);
    arb_set_d(b, 2.5);
    arb_pow(b, q, b, PREC);
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

void zeta(acb_ptr out, arb_t t, slong ZETA_PREC) {
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

void nt(arb_ptr out, arb_t t, slong PREC) {
    arb_t x;
    arb_t a;
    arb_t b;
    arb_init(x);
    arb_init(a);
    arb_init(b);

    arb_set(x, t);
    arb_set_d(a, 2.0);
    arb_const_pi(b, PREC);
    arb_div(x, x, a, PREC);
    arb_div(x, x, b, PREC);
    arb_set(a, x);
    arb_const_e(b, PREC);
    arb_div(a, a, b, PREC);
    arb_log(a, a, PREC);
    arb_mul(x, a, x, PREC);
    arb_set_d(a, 7.0);
    arb_set_d(b, 8.0);
    arb_div(a, a, b, PREC);
    arb_add(x, x, a, PREC);
    arb_set(out, x);

    arb_clear(x);
    arb_clear(a);
    arb_clear(b);
}

void nt_inv(arb_ptr out, arb_t m, slong PREC) {
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

    arb_set_d(x, 8.0);
    arb_mul(x, x, m, PREC);
    arb_set_d(a, 11.0);
    arb_sub(x, x, a, PREC);

    arb_set(nom, x);
    arb_const_pi(b, PREC);
    arb_mul(nom, nom, b, PREC);

    arb_set(den, x);
    arb_const_e(b, PREC);
    arb_div(den, den, b, PREC);
    arb_set_d(b, 8.0);
    arb_div(den, den, b, PREC);
    arb_lambertw( den, den, 0, PREC);
    arb_set_d(b, 4.0);
    arb_mul(den, den, b, PREC);

    arb_set(x, nom);
    arb_div(x, x, den, PREC);

    arb_set(out, x);

    arb_clear(x);
    arb_clear(a);
    arb_clear(b);
    arb_clear(den);
    arb_clear(nom);
}

void zero_count_approx(arb_ptr out, arb_t t, slong k, slong PREC) {
    arb_t u;
    arb_t w;
    arb_t q;
    arb_t x;
    arb_t y;
    arb_t z;
    arb_t v;
    arb_t s;

    arb_init(u);
    arb_init(w);
    arb_init(q);
    arb_init(y);
    arb_init(z);
    arb_init(v);
    arb_init(s);
    arb_init(x);

    nt(u, t, PREC);

    for(int i = 1; i <= k; ++i) {
        arb_set_d(q, n_nth_prime(i));
        coord(x, q, t, PREC);
        arb_abs(w, x);

        wave2(z, q, w, PREC);
        if(arb_is_negative(x)) {
            arb_neg(z, z);
        }
#if CUBIC == 1
        arb_set_d(v, 2.0);
        arb_sqrt(v, v, PREC);
        arb_inv(v, v, PREC);
        arb_set_d(y, i);
        arb_pow(y, y, v, PREC);
        arb_set_d(v, 3.0);
        arb_mul(y, y, v, PREC);
        arb_inv(v, y, PREC);

        arb_set_d(y, 0.5);
        arb_sub(y, y, v, PREC);
        arb_mul(z, z, y, PREC);
        arb_add(u, u, z, PREC);

        wave3(z, q, w, PREC);
        if(arb_is_negative(x)) {
            arb_neg(z, z);
        }
        arb_set_d(y, 0.5);
        arb_add(y, y, v, PREC);
        arb_mul(z, z, y, PREC);
#endif
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
    arb_clear(x);
}

const char *argp_program_version = "zzz 0";
const char *argp_program_bug_address = "<>";
static char doc[] = "fast approximation of large Riemann zeta zeros";
static char args_doc[] = "N [offset]";
static struct argp_option options[] = {
        { "k", 'k', "K", 0, "use first k primes for zero counting function approximation [default 100]"},
        { "evaluate", 'e', 0, 0, "evaluate Riemann zeta function value at the approximate zero location"},
        { "tolerance", 't', "TOL", 0, "tolerance for bisection [default 0.5]"},
        { "window", 'w', "WIN", 0, "initial span around Lambert W asymptotic zero location +- WIN [default 1.5]"},
        { "precision", 'p', "PREC", 0, "arb precision for counting function approximation [default 256]"},
        { "zeta-prec", 'z', "ZETA_PREC", 0, "arb precision for zeta evaluation [default 64]"},
        { "digits", 'd', "DIGITS", 0, "extra digits for number formatting [default 4]"},
        { "verbose", 'v', 0, 0, "verbose progress output"},
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
        case ARGP_KEY_ARG: return 0;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc, 0, 0, 0 };

void iter_print(arb_ptr lo_t, arb_ptr lo, arb_ptr hi_t, arb_ptr hi, slong DIGITS, slong verbose) {

    if(verbose == 0) return;

    flint_printf("lower x N(x)\t");
    arb_printd(lo, DIGITS);
    flint_printf("\t");
    arb_printd(lo_t, DIGITS);
    flint_printf("\n");

    flint_printf("upper x N(x)\t");
    arb_printd(hi, DIGITS);
    flint_printf("\t");
    arb_printd(hi_t, DIGITS);
    flint_printf("\n");
}
int main(int argc, char *argv[])
{

    if(argc == 1) {
        flint_printf("zzz --usage for help\n");
        return 1;
    }

    struct arguments arguments;

    arguments.k = 100;
    arguments.DIGITS = 4;
    arguments.PREC = 256;
    arguments.ZETA_PREC = 64;
    arguments.step0 = 0.01;
    arguments.w0 = 1.5;
    arguments.eval = 0;
    arguments.verbose = 0;

    int arg_index = 1;
    argp_parse(&argp, argc, argv, ARGP_NO_ARGS, &arg_index, &arguments);

    arb_t m;
    arb_init(m);

    if(arguments.k > 0)
        n_compute_primes(arguments.k);

    arb_t u;
    arb_t m0;
    arb_t m_lo;
    arb_t m_hi;
    arb_t tt;
    arb_t step;

    arb_init(u);
    arb_init(m0);
    arb_init(m_lo);
    arb_init(m_hi);
    arb_init(tt);
    arb_init(step);

    arb_set_str(m0, argv[arg_index], arguments.PREC);
    if(argc > arg_index+1) {
        arb_set_str(m, argv[arg_index+1], arguments.PREC);
        arb_add(m0, m0, m, arguments.PREC);
    } else {
        arb_zero(m);
    }

    arb_const_log10(u, arguments.PREC);
    arb_inv(u, u, arguments.PREC);
    arb_log(m, m0, arguments.PREC);
    arb_mul(u, m, u, arguments.PREC);
    arb_floor(u, u, arguments.PREC);
    arguments.DIGITS += arf_get_si(&u->mid, 0);

    arb_set_d(m, -0.5);
    nt_inv(tt, m0, arguments.PREC);
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

    if( arb_gt(lo, m) || arb_lt(hi, m) || arb_lt(hi, lo)) {

        iter_print(lo_t, lo, hi_t, hi, arguments.DIGITS, arguments.verbose);

        flint_printf("please increase the window\n");
        return 0;
    }

    arb_t mid_t;
    arb_t mid;
    arb_init(mid);
    arb_init(mid_t);

    while(1) {

        iter_print(lo_t, lo, hi_t, hi, arguments.DIGITS, arguments.verbose);

        arb_add(mid_t, lo_t, hi_t, arguments.PREC);
        arb_set_d(m0, 0.5);
        arb_mul(mid_t, mid_t, m0, arguments.PREC);

        zero_count_approx(mid, mid_t, arguments.k, arguments.PREC);

        arb_sub(m0, hi, lo, arguments.PREC);
        arb_set_d(u, arguments.step0);
        if(arb_lt(m0, u)) break;

        arb_zero(m0);
        if(arb_gt(mid, m)) {
            arb_set(hi_t, mid_t);
            arb_set(hi, mid);
        }
        else {
            arb_set(lo_t, mid_t);
            arb_set(lo, mid);
        }
    }

    iter_print(lo_t, lo, hi_t, hi, arguments.DIGITS, arguments.verbose);

    flint_printf("argument s = \t");
    arb_set_d(lo_t, 0.5);
    acb_set_arb_arb(zz, lo_t, mid_t);
    acb_printd(zz, arguments.DIGITS);
    flint_printf("\n");

    if(arguments.eval > 0) {
        flint_printf("value    z = \t");
        zeta(zz, mid_t, arguments.ZETA_PREC);
        acb_printd(zz, arguments.DIGITS);
        flint_printf("\n");
    }
    arb_set_d(lo_t, 0.005);

    arf_printd(&mid_t->mid, arguments.DIGITS);
    flint_printf("\n");

    arb_clear(mid);
    arb_clear(mid_t);
    arb_clear(m);
    arb_clear(u);
    arb_clear(m0);
    arb_clear(m_lo);
    arb_clear(m_hi);
    arb_clear(tt);
    arb_clear(step);
    arb_clear(lo);
    arb_clear(hi);
    arb_clear(lo_t);
    arb_clear(hi_t);
    acb_clear(zz);
    return 0;
}