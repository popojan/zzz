# Zeta Zeros Zeal
fast approximate location of large Riemann zeta zeros on the critical line

```text
zzz n [off [eval [step [k [PREC [ZETA_PREC [DIGITS]]]]]]]
    n         ... (n+off)-th zeta zero on critical line calculation [obligatory]
    off       ... (n+off)-th zeta zero on critical line calculation [default 0]
    eval      ... evaluate Riemann zeta function value at the approximate zero location [default 1]
    step      ... heuristic width (0.0, 0.25> [default 0.1]
    PREC      ... arb precision for counting function approximation [default 256]
    ZETA_PREC ... arb precision for zeta evaluation [default 64]
    DIGITS    ... extra decimal digits for number formatting [default 4]
```


## Quadratic approximation

matching periods amplitudes and tangents 

![waves](doc/waves.png)

## Approximate zero counting function

using first `k` primes

![counting](doc/counting.png)

## Error distribution

in comparison with k=0 (basic ProductLog approximation)

![errors](doc/errors.png)

# Approximate n-th zero locations

## Zero # 10^22 + 1

see [~odlyzko/zeta_tables/zeros5](https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros5)

```bash
$ time ./zzz 1e22 +1 0 0.1 10000
```

```text
zero counting function lower bound lo_t = 9999999999999999999999.9469 +/- 3.3660e-53
zero counting function upper bound hi_t = 10000000000000000000000.892 +/- 3.3660e-53

zeta zero imaginary part lower approximation lo = 1370919909931995308226.6531 +/- 5.9472e-56
zeta zero imaginary part upper approximation hi = 1370919909931995308226.7531 +/- 5.9472e-56

Riemann zeta argument s = (0.50000000000000000000000000 + 1370919909931995308226.7117j)  +/-  (0, 7.88e-54j)

To be refined.

1370919909931995308226.7117

real	0m0.083s
user	0m0.076s
sys	0m0.006s

```

## Zero # 10^100

```bash
$ time ./zzz 1e100 +0 0 0.1 10000 1024 64 | tail -n1
```

```text
280690383842894069903195445838256400084548030162846045192360059224930922349073043060335653109252473.2488

real	0m0.255s
user	0m0.247s
sys	0m0.008s
```

## Zeros # 10^1000000 +1 +2

```text
$ ./zzz 1e1000000 +1 0 0.001 2 8000000 64 24 | tail -c35
[...]4655.00980610558528085362029025640
$ time ./zzz 1e1000000 +2 0 0.001 2 8000000 64 24 | tail -c35
[...]4655.00980883435558616339614196111

real	0m36.620s
user	0m36.109s
sys	0m0.375s
```

# Literature

* Bernhard Riemann: *On the Number of Prime Numbers less than a Given Quantity*.
  * https://www.claymath.org/sites/default/files/ezeta.pdf
* Jonathan W. Bober, Ghaith A. Hiary: *New computations of the Riemann zeta function on the critical line*
  * https://arxiv.org/abs/1607.00709
* M. V. Berry, J. P. Keating: *The Riemann Zeros and Eigenvalue Asymptotics*
  * https://empslocal.ex.ac.uk/people/staff/mrwatkin/zeta/berry-keating1.pdf
