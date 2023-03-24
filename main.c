#include "acb.h"

slong PREC = 256;
slong ZETA_PREC = 64;
slong DIGITS = 4;

#define CUBIC 1

void coord(arb_ptr out, arb_t q, arb_t t) {
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

void wave2(arb_ptr out, arb_t q, arb_t t) {
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


void wave3(arb_ptr out, arb_t q, arb_t t) {
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

void zeta(acb_ptr out, arb_t t) {
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

void nt(arb_ptr out, arb_t t) {
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

void nt_inv(arb_ptr out, arb_t m) {
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

void zero_count_approx(arb_ptr out, arb_t t, slong k) {
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

    nt(u, t);

    for(int i = 1; i <= k; ++i) {
        arb_set_d(q, n_nth_prime(i));
        coord(x, q, t);
        arb_abs(w, x);

        wave2(z, q, w);
        if(arb_is_negative(x)) {
            arb_neg(z, z);
        }
#if CUBIC == 1
        arb_one(y);
        arb_set_d(v, 3.0);
        arb_div(y, y, v, PREC);
        arb_pow(y, q, y, PREC);
        arb_set_d(v, 0.4);
        arb_div(v, v, y, PREC);
        arb_set_d(y, 0.5);
        arb_sub(y, y, v, PREC);
        arb_mul(z, z, y, PREC);
        arb_add(u, u, z, PREC);

        wave3(z, q, w);
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
int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        flint_printf("zzz n [off [eval [step [k [PREC [ZETA_PREC [DIGITS]]]]]]]\n");
        flint_printf("    n         ... (n+off)-th zeta zero on critical line calculation [obligatory]\n");
        flint_printf("    off       ... (n+off)-th zeta zero on critical line calculation [default 0]\n");
        flint_printf("    eval      ... evaluate Riemann zeta function at the approximate zero location [default 1]\n");
        flint_printf("    step      ... heuristic width (0.0, 0.25> [default 0.1]\n");
        flint_printf("    k         ... use first k primes for zero counting function approximation [default 2]\n");
        flint_printf("    PREC      ... arb precision for counting function approximation [default 256]\n");
        flint_printf("    ZETA_PREC ... arb precision for zeta evaluation [default 64]\n");
        flint_printf("    DIGITS    ... extra decimal digits for number formatting [default 4]\n");
        return 1;
    }

    slong k = 2;
    double step0 = 0.1;
    slong ai = 2;
    slong eval = 1;
    arb_t m;
    arb_init(m);

    if (argc > ai) arb_set_str(m, argv[ai++], PREC);
    if (argc > ai) eval = atol(argv[ai++]);
    if (argc > ai) step0 = atof(argv[ai++]);
    if (argc > ai) k = atol(argv[ai++]);
    if (argc > ai) PREC = atol(argv[ai++]);
    if (argc > ai) ZETA_PREC = atol(argv[ai++]);
    if (argc > ai) DIGITS = atol(argv[ai++]);

    n_compute_primes(k);

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

    arb_set_str(m0, argv[1], PREC);
    arb_add(m0, m0, m, PREC);

    arb_const_log10(u, PREC);
    arb_inv(u, u, PREC);
    arb_log(m, m0, PREC);
    arb_mul(u, m, u, PREC);
    arb_floor(u, u, PREC);
    DIGITS += arf_get_si(&u->mid, 0);

    arb_set_d(m, -0.5);
    nt_inv(tt, m0);
    arb_add(m, m0, m, PREC);
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

    arb_set_d(lo, step0);
    arb_set_d(hi, 0.5);
    arb_mul(lo, lo, hi, PREC);
    arb_sub(tt, tt, lo, PREC);

    arb_set_d(step, 1.0);
    arb_set_d(hi, step0);
    arb_mul(step, step, hi, PREC);

    arb_pos_inf(lo);
    arb_neg_inf(hi);

    acb_t zz;
    acb_init(zz);

    for(int i = 0; i <= 1; ++i) {
        zero_count_approx(u, tt, k);
        if(arb_lt(u, lo)) {
            arb_set(lo, u);
            arb_set(lo_t, tt);
        }
        if(arb_gt(u, hi)) {
            arb_set(hi, u);
            arb_set(hi_t, tt);
        }
        if(i<1) arb_add(tt, tt, step, PREC);
    }
    flint_printf("zero counting function lower bound lo_t = ");
    arb_printd(lo, DIGITS);
    flint_printf("\n");
    flint_printf("zero counting function upper bound hi_t = ");
    arb_printd(hi, DIGITS);
    flint_printf("\n\n");

    flint_printf("zeta zero imaginary part lower approximation lo = ");
    arb_printd(lo_t, DIGITS);
    flint_printf("\n");
    flint_printf("zeta zero imaginary part upper approximation hi = ");
    arb_printd(hi_t, DIGITS);
    flint_printf("\n\n");

    arb_t aa;
    arb_t bb;
    arb_init(aa);
    arb_init(bb);

    arb_set(aa, lo_t);
    arb_set(bb, hi_t);

    arb_sub(tt, hi, lo, PREC);
    arb_inv(tt, tt, PREC);
    arb_sub(m, m, lo, PREC);
    arb_mul(tt, m, tt, PREC);
    arb_sub(m, bb, aa, PREC);
    arb_mul(tt, tt, m, PREC);
    arb_add(tt, tt, aa, PREC);


    flint_printf("Riemann zeta argument s = ");
    arb_set_d(m, 0.5);
    acb_set_arb_arb(zz, m, tt);
    acb_printd(zz, DIGITS);

    if(eval > 0) {
        flint_printf("\nRiemann zeta value    z = ");
        zeta(zz, tt);
        acb_printd(zz, DIGITS);
    }
    flint_printf("\n\nTo be refined.\n\n");
    arf_printd(&tt->mid, DIGITS);
    flint_printf("\n");

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
    arb_clear(aa);
    arb_clear(bb);
    acb_clear(zz);

    return 0;
}