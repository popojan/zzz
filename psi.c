#include "acb.h"

arb_t pi;
arb_t pi2;
arb_t pi2inv;
arb_t div5_pi4;
arb_t div4_pipi;
arb_t mx;
arb_t my;

void fsin(arb_ptr out, arb_srcptr x, slong PREC) {
    //mx= x % 2pi
    arb_mul(mx, x, pi2inv, PREC);
    arb_floor(mx, mx, PREC);
    arb_mul(mx, mx, pi2, PREC);
    arb_sub(mx, x, mx, PREC);

    arb_sub(mx, pi, mx, PREC);
    arb_abs(my, mx);
    arb_sub(my, pi, my, PREC);
    arb_mul(mx, mx, my, PREC);
    arb_mul(out, mx, div4_pipi, PREC);

}

void psi_base(arb_ptr out, arb_srcptr x, slong PREC) {
    // x - 1/2 Log[(x^2-1)/x^2]-Log[2 PI]
    arb_t ret;
    arb_t w;
    arb_init(ret);
    arb_init(w);

    arb_one(ret);
    arb_neg(ret, ret);
    arb_mul(w, x, x, PREC);
    arb_add(ret, ret, w, PREC);
    arb_div(ret, ret, w, PREC);
    arb_log(ret, ret, PREC);
    arb_set_d(w, 0.5);
    arb_mul(ret, ret, w, PREC);
    arb_neg(ret, ret);
    arb_add(ret, ret, x, PREC);

    arb_log(w, pi2, PREC);
    arb_sub(out, ret, w, PREC);

    arb_clear(ret);
    arb_clear(w);
}
void psi_coord(arb_ptr out, arb_srcptr x, arb_srcptr t, slong PREC) {
    // t Log[x]
    arb_t ret;
    arb_init(ret);
    arb_log(ret, x, PREC);
    arb_mul(out, ret, t, PREC);
    arb_clear(ret);
}

void psi_amp(arb_ptr out, arb_srcptr x, arb_srcptr t, slong PREC) {
    // - 8 Sqrt[x] / (4 t^2 + 1)
    arb_t ret;
    arb_t w;

    arb_init(ret);
    arb_init(w);

    arb_mul(ret, t, t, PREC);
    arb_set_d(w, 4.0);
    arb_mul(ret, ret, w, PREC);
    arb_one(w);
    arb_add(ret, ret, w, PREC);
    arb_sqrt(w, x, PREC);
    arb_div(w, w, ret, PREC);
    arb_set_d(ret, -8.0);
    arb_mul(out, w, ret, PREC);

    arb_clear(ret);
    arb_clear(w);
}

void psi_wave(arb_ptr out, arb_srcptr x, arb_srcptr t, slong PREC) {
    //1/2 Cos [x] + t Sin[x]
    arb_t ret;
    arb_t w;

    arb_init(ret);
    arb_init(w);

    arb_set_d(w, 0.5);
    arb_mul(w, pi, w, PREC);
    arb_add(w, x, w, PREC);
    fsin(ret, w, PREC);
    arb_set_d(w, 0.5);
    arb_mul(ret, ret, w, PREC);

    fsin(w, x, PREC);
    arb_mul(w, w, t, PREC);
    arb_add(out, ret, w, PREC);

    arb_clear(ret);
    arb_clear(w);
}
void psi(arb_ptr out, arb_srcptr x, arb_srcptr zeros, slong k, slong PREC) {
    arb_t ret;
    arb_t w;
    arb_t u;
    arb_init(ret);
    arb_init(w);
    arb_init(u);
    psi_base(ret, x, PREC);
    for(slong i = 0; i < k; ++i) {
        psi_coord(u, x, zeros + i, PREC);
        psi_wave(u, u, zeros + i, PREC);
        psi_amp(w, x, zeros + i, PREC);
        arb_mul(w, w, u, PREC);
        arb_add(ret, ret, w, PREC);
    }
    arb_set(out, ret);
    arb_clear(ret);
    arb_clear(w);
    arb_clear(u);
}

int main(int argc, char *argv[])
{
    const slong PREC = 256;

    //init fast inexact sin
    arb_init(pi);
    arb_init(pi2);
    arb_init(pi2inv);
    arb_init(div5_pi4);
    arb_init(div4_pipi);
    arb_init(mx);
    arb_init(my);

    arb_const_pi(pi, PREC);
    arb_set_d(pi2, 4.0);
    arb_div(div4_pipi, pi2, pi, PREC);
    arb_div(div4_pipi, div4_pipi, pi, PREC);
    arb_set_d(div5_pi4, 5.0);
    arb_div(div5_pi4, div5_pi4, pi2, PREC);
    arb_div(div5_pi4, div5_pi4, pi, PREC);
    arb_add(pi2, pi, pi, PREC);
    arb_inv(pi2inv, pi2, PREC);

    const slong ZERO_COUNT = 2;
    arb_ptr zeros = _arb_vec_init(ZERO_COUNT);
    arb_set_d(zeros + 0, 14.11); //default approx
    arb_set_d(zeros + 1, 20.99); //default approx

    arb_t mm;
    arb_init(mm);

    for(double x = 2.0; x < 10.0; x += 0.05) {
        arb_set_d(mm, x);
        psi(my, mm, zeros, ZERO_COUNT, PREC);

        arb_printn(mm, 8, ARB_STR_NO_RADIUS);
        flint_printf("\t");
        arb_printn(my, 8, ARB_STR_NO_RADIUS);
        flint_printf("\n");
    }

    arb_clear(pi);
    arb_clear(pi2);
    arb_clear(pi2inv);
    arb_clear(div5_pi4);
    arb_clear(div4_pipi);
    arb_clear(mx);
    arb_clear(my);
    arb_clear(mm);

    _arb_vec_clear(zeros, ZERO_COUNT);
    flint_cleanup();

    return 0;
}
