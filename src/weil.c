// weil.c -- implementation of weil.h; see also doc/ghy/e5-above-threshold.wls
// (the Wolfram prototype validated against Odlyzko zeros3) and
// doc/notes/band-saturation.md for the theory and the validity regime.

#include "weil.h"
#include <flint/arb_mat.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void weil_opts_default(weil_opts *o) {
    o->delta = 3.0;
    o->sig_d = 1e-7;
    // tight prior (in mean gaps): also the only thing regularizing the
    // window-edge zeros, whose data support W(x) is ~0
    o->sig_p = 0.03;
    o->lm_iters = 48;
    o->verbose = 0;
}

// cardinal B-spline N_16 (window kernel order 2 qw = 16, fixed), support
// |y| < 8; ~5 digits cancel in the alternating sum, leaving ~1e-11
static double bspline16(double y) {
    static const double binom[17] = {1, 16, 120, 560, 1820, 4368, 8008,
        11440, 12870, 11440, 8008, 4368, 1820, 560, 120, 16, 1};
    static const double fact15 = 1307674368000.0;
    if (y < 0) y = -y;
    if (y >= 8.0) return 0.0;
    double s = 0.0, sign = 1.0;
    for (int j = 0; j <= 16; ++j, sign = -sign) {
        double t = y + 8.0 - j;
        if (t > 0) {
            double t3 = t * t * t, t5 = t3 * t * t;
            s += sign * binom[j] * t5 * t5 * t5;
        }
    }
    return s / fact15;
}

static double ww16(double bw, double x) {
    if (x == 0) return 1.0;
    double s = sin(bw * x) / (bw * x);
    double s2 = s * s, s4 = s2 * s2, s8 = s4 * s4;
    return s8 * s8;
}

static double dww16(double bw, double x) {
    if (x == 0) return 0.0;
    double bx = bw * x;
    double s = sin(bx) / bx;
    double ds = bw * cos(bx) / bx - bw * sin(bx) / (bx * bx);
    double s2 = s * s, s4 = s2 * s2, s8 = s4 * s4;
    return 16.0 * s8 * s4 * s2 * s * ds;
}

static double phase_mod_2pi(const arb_t t, const arb_t u, slong PREC) {
    arb_t ph, twopi, q;
    double r;
    arb_init(ph);
    arb_init(twopi);
    arb_init(q);
    arb_mul(ph, t, u, PREC);
    arb_const_pi(twopi, PREC);
    arb_mul_2exp_si(twopi, twopi, 1);
    arb_div(q, ph, twopi, PREC);
    arb_floor(q, q, PREC);
    arb_submul(ph, q, twopi, PREC);
    r = arf_get_d(arb_midref(ph), ARF_RND_NEAR);
    arb_clear(ph);
    arb_clear(twopi);
    arb_clear(q);
    return r;
}

// prime side + archimedean term on the frequency grid
static void weil_data(double *dcos, double *dsin, const double *ws, slong nw,
                      const arb_t t_base, ulong X, const weil_opts *opt,
                      slong PREC)
{
    double bw = opt->delta / 16.0;
    double inv2bw = 1.0 / (2.0 * bw);
    double dw = ws[1] - ws[0];
    double tc_d = arf_get_d(arb_midref(t_base), ARF_RND_NEAR);
    arb_t lp, u;
    slong i, j;

    for (i = 0; i < nw; ++i) {
        dcos[i] = 2.0 * log(tc_d / (2.0 * M_PI)) * bspline16(ws[i] * inv2bw) * inv2bw;
        dsin[i] = 0.0;
    }

    arb_init(lp);
    arb_init(u);
    for (i = 1; ; ++i) {
        ulong p = n_nth_prime(i);
        if (p > X) break;
        arb_set_ui(lp, p);
        arb_log(lp, lp, PREC);
        double lp_d = log((double) p);
        ulong pm = p;
        for (ulong m = 1; ; ++m) {
            double u_d = m * lp_d;
            double amp = lp_d * pow((double) p, -0.5 * (double) m);
            arb_mul_ui(u, lp, m, PREC);
            double ph = phase_mod_2pi(t_base, u, PREC);
            double c = cos(ph), s = sin(ph);
            slong j0 = (slong) ceil((u_d - opt->delta - ws[0]) / dw);
            slong j1 = (slong) floor((u_d + opt->delta - ws[0]) / dw);
            if (j0 < 0) j0 = 0;
            if (j1 >= nw) j1 = nw - 1;
            for (j = j0; j <= j1; ++j) {
                double ghm = bspline16((u_d - ws[j]) * inv2bw) * inv2bw;
                double ghp = bspline16((u_d + ws[j]) * inv2bw) * inv2bw;
                dcos[j] -= 2.0 * amp * c * (ghm + ghp);
                dsin[j] -= 2.0 * amp * s * (ghp - ghm);
            }
            if (pm > X / p) break;
            pm *= p;
        }
    }
    arb_clear(lp);
    arb_clear(u);
}

// model: out[2i] = cos row, out[2i+1] = sin row
static void weil_model(double *out, const double *xs, slong M,
                       const double *far, slong nfar,
                       const double *ws, slong nw, double bw)
{
    slong i, j;
    memset(out, 0, sizeof(double) * 2 * nw);
    for (j = 0; j < M + nfar; ++j) {
        double x = (j < M) ? xs[j] : far[j - M];
        double W = ww16(bw, x);
        if (fabs(W) < 1e-14) continue;
        for (i = 0; i < nw; ++i) {
            out[2 * i]     += 2.0 * W * cos(ws[i] * x);
            out[2 * i + 1] += 2.0 * W * sin(ws[i] * x);
        }
    }
}

static double norm2(const double *v, slong n) {
    double s = 0;
    for (slong i = 0; i < n; ++i) s += v[i] * v[i];
    return sqrt(s);
}

// solve (M x M) A dx = b at 192 bits (the ridge against rank-deficient
// J^T J exceeds double precision conditioning)
static int solve_hp(double *dx, const double *A, const double *b, slong M) {
    arb_mat_t Am, bm, xm;
    int ok;
    slong i, j;
    arb_mat_init(Am, M, M);
    arb_mat_init(bm, M, 1);
    arb_mat_init(xm, M, 1);
    for (i = 0; i < M; ++i) {
        for (j = 0; j < M; ++j)
            arb_set_d(arb_mat_entry(Am, i, j), A[i * M + j]);
        arb_set_d(arb_mat_entry(bm, i, 0), b[i]);
    }
    ok = arb_mat_solve(xm, Am, bm, 192);
    if (ok)
        for (i = 0; i < M; ++i)
            dx[i] = arf_get_d(arb_midref(arb_mat_entry(xm, i, 0)), ARF_RND_NEAR);
    arb_mat_clear(Am);
    arb_mat_clear(bm);
    arb_mat_clear(xm);
    return ok ? 0 : 1;
}

int weil_fit(double *xs, slong M, const double *far, slong nfar,
             const arb_t t_base, ulong X, const weil_opts *opt, slong PREC)
{
    double bw = opt->delta / 16.0;
    double a = log((double) X);
    double wmax = a - opt->delta - 0.05;
    double tc_d = arf_get_d(arb_midref(t_base), ARF_RND_NEAR);
    double gap = 2.0 * M_PI / log(tc_d / (2.0 * M_PI));
    double step_max = 0.3 * gap;
    double sigp_abs = opt->sig_p * gap;
    double ridge = (opt->sig_d / sigp_abs) * (opt->sig_d / sigp_abs);
    slong i, j, l, it;
    int fail = 0;

    double xmax = 0;
    for (j = 0; j < M; ++j) if (fabs(xs[j]) > xmax) xmax = fabs(xs[j]);
    slong nw = (slong) ceil(wmax * xmax / M_PI) + 8;
    if (nw < 40) nw = 40;
    if (nw > 160) nw = 160;
    slong nd = 2 * nw;

    double *ws = malloc(nw * sizeof(double));
    for (i = 0; i < nw; ++i) ws[i] = wmax * i / (nw - 1);

    double *data = malloc(nd * sizeof(double));
    double *modl = malloc(nd * sizeof(double));
    double *r    = malloc(nd * sizeof(double));
    double *J    = malloc(nd * M * sizeof(double));
    double *A    = malloc(M * M * sizeof(double));
    double *rhs  = malloc(M * sizeof(double));
    double *dx   = malloc(M * sizeof(double));
    double *x0   = malloc(M * sizeof(double));
    double *xtry = malloc(M * sizeof(double));
    memcpy(x0, xs, M * sizeof(double));

    weil_data(data, data + nw, ws, nw, t_base, X, opt, PREC);
    // weil_data fills two contiguous halves; interleave into (cos,sin) rows
    {
        double *tmp = malloc(nd * sizeof(double));
        for (i = 0; i < nw; ++i) {
            tmp[2 * i] = data[i];
            tmp[2 * i + 1] = data[nw + i];
        }
        memcpy(data, tmp, nd * sizeof(double));
        free(tmp);
    }

    double lam = 1e-3;
    for (it = 0; it < opt->lm_iters && !fail; ++it) {
        weil_model(modl, xs, M, far, nfar, ws, nw, bw);
        for (i = 0; i < nd; ++i) r[i] = data[i] - modl[i];
        double rn = norm2(r, nd);

        for (j = 0; j < M; ++j) {
            double W = ww16(bw, xs[j]), dW = dww16(bw, xs[j]);
            for (i = 0; i < nw; ++i) {
                double cw = cos(ws[i] * xs[j]), sw = sin(ws[i] * xs[j]);
                J[(2 * i) * M + j]     = 2.0 * (dW * cw - W * ws[i] * sw);
                J[(2 * i + 1) * M + j] = 2.0 * (dW * sw + W * ws[i] * cw);
            }
        }

        // normal equations (data part); prior added below per try
        for (j = 0; j < M; ++j) {
            for (l = j; l < M; ++l) {
                double s = 0;
                for (i = 0; i < nd; ++i) s += J[i * M + j] * J[i * M + l];
                A[j * M + l] = A[l * M + j] = s;
            }
            double s = 0;
            for (i = 0; i < nd; ++i) s += J[i * M + j] * r[i];
            rhs[j] = s - ridge * (xs[j] - x0[j]);
        }

        int accepted = 0;
        for (int tries = 0; tries < 8 && !accepted; ++tries) {
            double *At = malloc(M * M * sizeof(double));
            memcpy(At, A, M * M * sizeof(double));
            for (j = 0; j < M; ++j)
                At[j * M + j] += ridge + lam * (A[j * M + j] + 1e-12);
            if (solve_hp(dx, At, rhs, M)) { free(At); fail = 1; break; }
            free(At);
            // per-component clip: window-edge zeros (W ~ 0) produce runaway
            // transient steps from cross-talk noise; a global rescale would
            // freeze the well-supported zeros instead
            for (j = 0; j < M; ++j)
                if (fabs(dx[j]) > step_max)
                    dx[j] = dx[j] > 0 ? step_max : -step_max;
            for (j = 0; j < M; ++j) xtry[j] = xs[j] + dx[j];
            weil_model(modl, xtry, M, far, nfar, ws, nw, bw);
            for (i = 0; i < nd; ++i) r[i] = data[i] - modl[i];
            if (norm2(r, nd) < rn) {
                memcpy(xs, xtry, M * sizeof(double));
                lam = lam / 3 > 1e-7 ? lam / 3 : 1e-7;
                accepted = 1;
            } else {
                lam *= 10;
            }
        }
        if (opt->verbose) {
            double mx = 0;
            for (j = 0; j < M; ++j) if (fabs(dx[j]) > mx) mx = fabs(dx[j]);
            flint_fprintf(stderr, "weil LM %wd: |r| = %.3e, max|dx| = %.3e, lam = %.1e%s\n",
                          it + 1, norm2(r, nd), mx, lam, accepted ? "" : " (stalled)");
        }
        if (!accepted) break;
        {
            double mx = 0;
            for (j = 0; j < M; ++j) if (fabs(dx[j]) > mx) mx = fabs(dx[j]);
            if (mx < 1e-9) break;
        }
    }

    free(ws); free(data); free(modl); free(r); free(J);
    free(A); free(rhs); free(dx); free(x0); free(xtry);
    return fail;
}
