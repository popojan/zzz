# Zeta Zeros Zeal
fast approximate location of large zeta zeros on the critical line

## Quadratic approximation

matching periods amplitudes and tangents 

![waves](doc/waves.png)

## Approximate zero counting function

using first `k` primes

![counting](doc/counting.png)

## Error distribution

in comparison with k=0 (basic ProductLog approximation)

![errors](doc/errors.png)

```bash
$ time ./zzz 10000000000000000000001 0 0.1 10000
```

# Approximate n-th zero locations

## Zero # 10^22 + 1

see [~odlyzko/zeta_tables/zeros5](https://www-users.cse.umn.edu/~odlyzko/zeta_tables/zeros5)

```text
zero counting function lower bound lo_t = 10000000000000000000001.492 +/- 3.3661e-53
zero counting function upper bound hi_t = 10000000000000000000002.237 +/- 3.3675e-53

zeta zero imaginary part lower approximation lo = 1370919909931995308226.8202 +/- 5.9440e-56
zeta zero imaginary part upper approximation hi = 1370919909931995308226.9202 +/- 5.9440e-56

Riemann zeta argument s = (0.50000000000000000000000000 + 1370919909931995308226.6871j)  +/-  (0, 1.68e-53j)

To be refined.

1370919909931995308226.6871

real	0m0.127s
user	0m0.123s
sys	0m0.003s
```

## Zero # 10^100

```$ time ./zzz 1e100 0 0.1 10000 1024 64 | tail -n1```

```text
280690383842894069903195445838256400084548030162846045192360059224930922349073043060335653109252473.2500

real	0m0.344s
user	0m0.329s
sys	0m0.014s
```