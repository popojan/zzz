# Dependencies

* gmp https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz
* mpfr https://www.mpfr.org/mpfr-current/mpfr-4.2.0.tar.gz
* flint http://www.flintlib.org/flint-2.9.0.tar.gz
* arb https://github.com/fredrik-johansson/arb.git

```bash
cd arb
./configure --disable-shared \
  --with-mpfr=../mpfr-4.2.0/ \
  --with-gmp=../gmp-6.2.1/ --with-flint=../flint-2.9.0/
```
