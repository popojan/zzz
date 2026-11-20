# Dependencies

## Distro packages (recommended, flint >= 3)

Flint 3 absorbed arb, so the default (dynamic) CMake build needs only:

```bash
sudo apt install cmake gcc libgmp-dev libmpfr-dev libflint-dev     # Debian/Ubuntu
sudo dnf install cmake gcc gmp-devel mpfr-devel flint-devel        # Fedora
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build            # ~6 s smoke tests (tests/smoke.sh)
```

Verified on Ubuntu 24.04 (flint 3.0.1, including under WSL2) and Fedora 44
(flint 3.4.0, gcc 16).

## Vendored static build (`-DSTATIC_LINKING=ON`)

Three libraries, no separate arb. Install each into one prefix; CMake takes
headers from `<prefix>/include` and the `.a` archives from `<prefix>/lib`.
The default prefix is `deps/local` (override with `-DZZZ_DEPS_PREFIX=...`);
everything under `deps/` except this file is gitignored.

* gmp   https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz
* mpfr  https://www.mpfr.org/mpfr-4.2.0/mpfr-4.2.0.tar.gz
* flint https://github.com/flintlib/flint/releases/download/v3.1.3/flint-3.1.3.tar.gz

```bash
cd deps && P=$PWD/local
(cd gmp-6.2.1   && ./configure --prefix=$P --disable-shared && make -j && make install)
(cd mpfr-4.2.0  && ./configure --prefix=$P --disable-shared --with-gmp=$P && make -j && make install)
(cd flint-3.1.3 && ./configure --prefix=$P --disable-shared --with-gmp=$P --with-mpfr=$P && make -j && make install)
cd .. && cmake -S . -B build -DSTATIC_LINKING=ON && cmake --build build -j
```
