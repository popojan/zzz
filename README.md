# Zeta Zeros Zeal
fast approximation of large Riemann zeta zeros on the critical line

```text
Usage: zzz [OPTION...] N [offset] [count]
fast approximation of large Riemann zeta zeros

  -d, --digits=DIGITS        extra digits for number formatting [default 6]
  -e, --evaluate             evaluate Riemann zeta function value at the
                             approximate zero location
  -g, --debug                debug counting function from <N> to <N+offset> in
                             <count> steps
  -k, --k=K                  use first k primes for zero counting function
                             approximation [default 100]
  -p, --precision=PREC       arb precision for counting function approximation
                             [default 256]
  -t, --tolerance=TOL        tolerance for bisection [default 1e-6]
  -v, --verbose              verbose progress output
  -w, --window=WIN           initial span around Lambert W asymptotic zero
                             location +- WIN [default 1.5]
  -z, --zeta-prec=ZETA_PREC  arb precision for zeta evaluation [default 64]
  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version
```

## Zero counting function approximation

Note: obsolete inner sum approximation, not used any more.

Combines quadratic and cubic spline with correct frequency and tangents to match the amplitude.

![waves](doc/waves.png)

## Towards convergence

![waves](doc/convergence.png)

## Error distribution

In comparison with k=-∞ (basic Lambert W approximation).

![errors](doc/errors.png)

# Approximate n-th zero locations

## Zero # 10^12 + 1

see [~odlyzko/zeta_tables/zeros3](https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros3)

```bash
$ time ./zzz -ve -k 1000 1e12 +1
```

```
argument s =    (0.500000000000000000 + 267653395648.844684j)  +/-  (0, 1.05e-65j)
value    z =    (0.355290959100415380 + 0.132397324302526229j)  +/-  (3.70e-20, 2.34e-20j)
267653395648.844684

real    0m0.470s
user    0m0.456s
sys     0m0.007s
```

## Zero # 10^36 + 42420637374017961984

```bash
  $ time ./zzz -k 10000 1e36 42420637374017961984
```

```text
81029194732694548890047854481676713.009431

real    0m2.381s
user    0m2.366s
sys     0m0.005s
```

```
81029194732694548890047854481676712.94002   prev approximate     #10^36+42420637374017961983
81029194732694548890047854481676712.98790          published     #10^36+42420637374017961984
81029194732694548890047854481676713.00943        approximate     #10^36+42420637374017961984
81029194732694548890047854481676713.08748   next approximate     #10^36+42420637374017961985
```

## Chebyshev Psi Exact Formula

Using zeros approximated by `zzz -k 1000`.

range 0 to 20  (50 zeros)          |    range 541 to 661 (1,000 zeros)     | range 7920-8020 (10,000 zeros)  
:---------------------------------:|:-------------------------------------:|:--------------------------------------:
![](doc/psi-50-zeros-k1000-p1.png) | ![](doc/psi-10k-zeros-k1000-p100.png) | ![](doc/psi-10k-zeros-k1000-p1000.png)


# Literature

* Bernhard Riemann: *On the Number of Prime Numbers less than a Given Quantity*.
  * https://www.claymath.org/sites/default/files/ezeta.pdf
* Jonathan W. Bober, Ghaith A. Hiary: *New computations of the Riemann zeta function on the critical line*
  * https://arxiv.org/abs/1607.00709
* M. V. Berry, J. P. Keating: *The Riemann Zeros and Eigenvalue Asymptotics*
  * https://empslocal.ex.ac.uk/people/staff/mrwatkin/zeta/berry-keating1.pdf
* Guilherme França, André LeClair: *Statistical and other properties of Riemann zeros based on an explicit equation for the n-th zero on the critical line*
  * https://arxiv.org/abs/1307.8395
