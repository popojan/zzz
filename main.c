#include "acb.h"

slong PREC = 256;
slong ZETA_PREC = 64;
slong DIGITS = 4;

void zeta(acb_ptr out, arb_t t) {
    acb_t s;
    acb_t z;
    acb_init(s);
    acb_init(z);
    arb_t a;
    arb_init(a);
    arb_set_str(a, "0.5", PREC);
    acb_set_arb_arb(s, a, t);
    acb_zeta(z, s, ZETA_PREC);
    acb_set(out, z);
}

void nt(arb_ptr out, arb_t t) {
    arb_t x;
    arb_t a;
    arb_t b;
    arb_init(x);
    arb_init(a);
    arb_init(b);

    arb_set(x, t);
    arb_set_str(a, "2.0", PREC);
    arb_const_pi(b, PREC);
    arb_div(x, x, a, PREC);
    arb_div(x, x, b, PREC);
    arb_set(a, x);
    arb_const_e(b, PREC);
    arb_div(a, a, b, PREC);
    arb_log(a, a, PREC);
    arb_mul(x, a, x, PREC);
    arb_set_str(a, "7.0", PREC);
    arb_set_str(b, "8.0", PREC);
    arb_div(a, a, b, PREC);
    arb_add(x, x, a, PREC);
    arb_set(out, x);

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

    arb_set_str(x, "8.0", PREC);
    arb_mul(x, x, m, PREC);
    arb_set_str(a, "11.0", PREC);
    arb_sub(x, x, a, PREC);

    arb_set(nom, x);
    arb_const_pi(b, PREC);
    arb_mul(nom, nom, b, PREC);

    arb_set(den, x);
    arb_const_e(b, PREC);
    arb_div(den, den, b, PREC);
    arb_set_str(b, "8.0", PREC);
    arb_div(den, den, b, PREC);
    arb_lambertw( den, den, 0, PREC);
    arb_set_str(b, "4.0", PREC);
    arb_mul(den, den, b, PREC);

    arb_set(x, nom);
    arb_div(x, x, den, PREC);

    arb_set(out, x);
}
void pw(arb_ptr out, arb_t q, arb_t t) {

    // PI / Log[q]
    arb_t pi;
    arb_t pi_div_log_q;
    arb_t x;
    arb_t a;
    arb_t b;
    arb_t sqrt_q;
    arb_t xx;

    arb_init(xx);
    arb_init(pi);
    arb_init(pi_div_log_q);
    arb_init(x);
    arb_init(sqrt_q);
    arb_init(a);
    arb_init(b);

    arb_const_pi(pi, PREC);
    arb_log(pi_div_log_q, q, PREC);
    arb_inv(pi_div_log_q, pi_div_log_q, PREC);
    arb_mul(pi_div_log_q, pi_div_log_q, pi, PREC);

    // -t + PI / Log[q]
    arb_set(x, pi_div_log_q);
    arb_sub(x, x, t, PREC);

    //Mod[-t + PI / Log[q], 2 PI / Log[q]]
    arb_set_str(b, "2.0", PREC);
    arb_set(a, pi_div_log_q);
    arb_mul(a, a, b, PREC);
    arb_div(b, x, a, PREC);
    arb_floor(b, b, PREC);
    arb_mul(b, b, a, PREC);
    arb_sub(x, x, b, PREC);
    arb_set(xx, x);

    //Abs[% - PI/Log[q]]
    arb_sub(x, x, pi_div_log_q, PREC);
    arb_abs(x, x);

    // 4/PI Sqrt(q) Log(q) Abs[...]
    arb_set(sqrt_q, q);
    arb_sqrt(sqrt_q, sqrt_q, PREC);
    arb_set(b, q);
    arb_log(b, b, PREC);
    arb_mul(b, b, sqrt_q, PREC);
    arb_div(b, b, pi, PREC);
    arb_set_str(a, "4.0", PREC);
    arb_mul(b, b, a, PREC);
    arb_mul(x, x, b, PREC);

    // Sqrt((Sqrt(q)-1)^2 + 4/PI) ...
    arb_set(b, sqrt_q);
    arb_set_str(a, "-1.0", PREC);
    arb_add(b, b, a, PREC);
    arb_mul(b, b, b, PREC);
    arb_add(x, x, b, PREC);
    arb_sqrt(x, x, PREC);

    //(sqrt(q)-1-#)((sqrt(q)+1-#))&@ ...
    arb_set_str(a, "-1.0", PREC);
    arb_set_str(b, "1.0", PREC);
    arb_add(a, a, sqrt_q, PREC);
    arb_add(b, b, sqrt_q, PREC);
    arb_sub(a, a, x, PREC);
    arb_sub(b, b, x, PREC);
    arb_mul(x, a, b, PREC);

    arb_set_d(a, 1.0);
    if(arb_lt(xx, pi_div_log_q)) {
        arb_neg(a, a);
    }
    arb_set_d(b, 4.0);
    arb_mul(b, b, sqrt_q, PREC);
    arb_div(a, a, b, PREC);
    arb_mul(x, x, a, PREC);
    arb_neg(x, x);

    arb_set(out, x);
}

void zest(arb_ptr out, arb_t t, slong k) {
    arb_t u;
    arb_t w;
    arb_t q;

    arb_init(u);
    arb_init(w);
    arb_init(q);

    arb_one(u);
    arb_sub(u, u, u, PREC);
    for(int i = 1; i <= k; ++i) {
        arb_set_d(q, n_nth_prime(i));
        pw(w, q, t);
        arb_add(u, u, w, PREC);
    }

    nt(w, t);
    arb_add(u, u, w, PREC);
    arb_set(out, u);
}

int main(int argc, char *argv[])
{

    if (argc < 3)
    {
        flint_printf("zzz n off [[[[[[eval] step] k] PREC] ZETA_PREC] DIGITS]\n");
        flint_printf("    n         ... (n+off)-th zeta zero on critical line calculation [obligatory]\n");
        flint_printf("    off       ... (n+off)-th zeta zero on critical line calculation [obligatory]\n");
        flint_printf("    eval      ... evaluate Riemann zeta function value at the approximate zero location [default 1]\n");
        flint_printf("    step      ... heuristic width (0.0, 0.25> [default 0.1]\n");
        flint_printf("    PREC      ... arb precision for counting functiom approximation [default 256]\n");
        flint_printf("    ZETA_PREC ... arb precision for zeta evaluation [default 64]\n");
        flint_printf("    DIGITS    ... extra decimal digits for number formating [default 4]\n");
        return 1;
    }

    slong k = 2;
    double step0 = 0.1;
    slong ai = 3;
    slong eval = 1;
    if (argc > ai) eval = atol(argv[ai++]);
    if (argc > ai) step0 = atof(argv[ai++]);
    if (argc > ai) k = atol(argv[ai++]);
    if (argc > ai) PREC = atol(argv[ai++]);
    if (argc > ai) ZETA_PREC = atol(argv[ai++]);
    if (argc > ai) DIGITS = atol(argv[ai++]);

    n_compute_primes(k);

    arb_t u;
    arb_t m;
    arb_t m0;
    arb_t m_lo;
    arb_t m_hi;
    arb_t tt;
    arb_t step;

    arb_init(u);
    arb_init(m);
    arb_init(m0);
    arb_init(m_lo);
    arb_init(m_hi);
    arb_init(tt);
    arb_init(step);

    arb_set_str(m0, argv[1], PREC);
    arb_set_str(m, argv[2], PREC);
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
    arb_neg(lo, lo);
    arb_set_d(hi, 0.5);
    arb_mul(lo, lo, hi, PREC);
    arb_sub(tt, tt, lo, PREC);

    arb_pos_inf(lo);
    arb_neg_inf(hi);
    arb_set_d(step, step0);

    acb_t zz;
    acb_init(zz);

    for(int i = 0; i <= 1; ++i) {
        zest(u, tt, k);
        if(arb_lt(u, lo)) {
            arb_set(lo, u);
            arb_set(lo_t, tt);
        }
        if(arb_gt(u, hi)) {
            arb_set(hi, u);
            arb_set(hi_t, tt);
        }
        if(i==0) arb_add(tt, tt, step, PREC);
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
    return 0;
}