# Dependencies

## Distro packages (recommended, flint >= 3)

Flint 3 absorbed arb, so the default (dynamic) CMake build needs only:

```bash
sudo apt install cmake gcc libgmp-dev libmpfr-dev libflint-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Verified on Ubuntu 24.04 (flint 3.0.1), including under WSL2.

## Vendored static build (legacy, flint 2.9 + external arb)

For `-DSTATIC_LINKING=ON`, unpack/build the following here:

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
