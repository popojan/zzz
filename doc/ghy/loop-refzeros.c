// loop-refzeros.c -- reference zeta zeros n0 .. n0+len-1 from FLINT, as "index<TAB>gamma".
// Ground truth for loop-zero-accuracy.py only; the loop never evaluates zeta.
//
//   gcc -O2 -o loop-refzeros doc/ghy/loop-refzeros.c -lflint -lm
//   ./loop-refzeros 30000 2500 > ref_30000.tsv          (~20 ms per zero near index 10^6)
#include <flint/arb.h>
#include <flint/acb_dirichlet.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: loop-refzeros n0 len\n"); return 1; }
    slong n0 = atol(argv[1]), len = atol(argv[2]);
    arb_ptr res = _arb_vec_init(len);
    fmpz_t n; fmpz_init_set_si(n, n0);
    acb_dirichlet_hardy_z_zeros(res, n, len, 80);
    for (slong i = 0; i < len; ++i)
        printf("%ld\t%.15f\n", n0 + i, arf_get_d(arb_midref(res + i), ARF_RND_NEAR));
    _arb_vec_clear(res, len); fmpz_clear(n);
    return 0;
}
