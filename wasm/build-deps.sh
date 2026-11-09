#!/usr/bin/env bash
# Cross-compile GMP, MPFR and FLINT to WebAssembly with Emscripten.
#
# The native build links the system FLINT/GMP; those are x86 shared objects
# and cannot be linked into a WASM module, so the WASM build carries its own
# static copies, built here into wasm/deps/prefix. Run once (it takes a few
# minutes); `make -f Makefile.wasm` then picks them up.
#
#   wasm/build-deps.sh            # download + build everything
#   JOBS=8 wasm/build-deps.sh     # limit parallelism
#
# Requires emcc/emconfigure/emmake on PATH (source emsdk_env.sh), curl, tar,
# m4 (GMP's configure) and the usual autotools runtime (sh, make).
set -euo pipefail

GMP_VERSION=${GMP_VERSION:-6.3.0}
MPFR_VERSION=${MPFR_VERSION:-4.2.2}
FLINT_VERSION=${FLINT_VERSION:-3.6.0}

HERE=$(cd "$(dirname "$0")" && pwd)
DEPS=$HERE/deps
SRC=$DEPS/src
PREFIX=$DEPS/prefix
JOBS=${JOBS:-$(nproc 2>/dev/null || echo 4)}

# -O3 for speed; no -msimd128 / -pthread so the libraries run on every
# WASM engine without SIMD or SharedArrayBuffer.
export CFLAGS="-O3"
export CXXFLAGS="-O3"

for tool in emcc emconfigure emmake curl tar m4; do
    command -v "$tool" >/dev/null || { echo "build-deps: '$tool' not found on PATH" >&2; exit 1; }
done

mkdir -p "$SRC" "$PREFIX"

fetch() { # url file
    if [ ! -f "$SRC/$2" ]; then
        curl -fL --retry 3 -o "$SRC/$2.part" "$1"
        mv "$SRC/$2.part" "$SRC/$2"
    fi
}

# --- GMP -------------------------------------------------------------------
# --host=none selects the generic C limb code (no assembly); on wasm32 a limb
# is 32 bits. The build-time helpers (gen-fac, gen-bases, ...) must run on the
# build machine, so CC_FOR_BUILD points at the host compiler.
if [ ! -f "$PREFIX/lib/libgmp.a" ]; then
    fetch "https://gmplib.org/download/gmp/gmp-$GMP_VERSION.tar.xz" "gmp-$GMP_VERSION.tar.xz"
    rm -rf "$SRC/gmp-$GMP_VERSION"
    tar -C "$SRC" -xf "$SRC/gmp-$GMP_VERSION.tar.xz"
    (
        cd "$SRC/gmp-$GMP_VERSION"
        emconfigure ./configure --prefix="$PREFIX" --host=none --build="$(./config.guess)" \
            --disable-assembly --disable-shared --enable-static --disable-cxx \
            CC_FOR_BUILD="${CC_FOR_BUILD:-cc}"
        emmake make -j"$JOBS"
        emmake make install
    )
fi

# --- MPFR (FLINT >= 3 requires it) ------------------------------------------
if [ ! -f "$PREFIX/lib/libmpfr.a" ]; then
    fetch "https://www.mpfr.org/mpfr-$MPFR_VERSION/mpfr-$MPFR_VERSION.tar.xz" "mpfr-$MPFR_VERSION.tar.xz"
    rm -rf "$SRC/mpfr-$MPFR_VERSION"
    tar -C "$SRC" -xf "$SRC/mpfr-$MPFR_VERSION.tar.xz"
    (
        cd "$SRC/mpfr-$MPFR_VERSION"
        emconfigure ./configure --prefix="$PREFIX" --host=wasm32-unknown-emscripten \
            --disable-shared --enable-static --disable-thread-safe \
            --with-gmp="$PREFIX"
        emmake make -j"$JOBS"
        emmake make install
    )
fi

# --- FLINT -------------------------------------------------------------------
if [ ! -f "$PREFIX/lib/libflint.a" ]; then
    fetch "https://github.com/flintlib/flint/releases/download/v$FLINT_VERSION/flint-$FLINT_VERSION.tar.gz" "flint-$FLINT_VERSION.tar.gz"
    rm -rf "$SRC/flint-$FLINT_VERSION"
    tar -C "$SRC" -xf "$SRC/flint-$FLINT_VERSION.tar.gz"
    (
        cd "$SRC/flint-$FLINT_VERSION"
        emconfigure ./configure --prefix="$PREFIX" --host=wasm32-unknown-emscripten \
            --disable-shared --enable-static --disable-assembly --disable-pthread \
            --without-blas --without-ntl \
            --with-gmp="$PREFIX" --with-mpfr="$PREFIX"
        emmake make -j"$JOBS" library
        emmake make install
    )
fi

echo "build-deps: done, static libraries in $PREFIX/lib"
