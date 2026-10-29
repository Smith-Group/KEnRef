#!/usr/bin/env bash
# Install the Eigen release $KENREF_EIGEN_VERSION from its release archive, verified against
# $KENREF_EIGEN_SHA256, into ~/.cache/kenref/eigen-<version> (or the prefix given as $1). Used by the
# release/1.x CI, whose CMake has no fetch fallback.
#
# The prefix and its content are what cmake/Dependencies.cmake stages on release/2.x and master with
# KENREF_EIGEN_SHARED_INSTALL=ON, because every branch's CI restores the SAME cache entry. That includes
# share/eigen3/kenref-source.txt: keep its line identical to the one Dependencies.cmake writes, or the
# CMake side refuses to reuse a prefix this script made.
set -euo pipefail
ver="${KENREF_EIGEN_VERSION:?}"
sha="${KENREF_EIGEN_SHA256:?}"
prefix="${1:-$HOME/.cache/kenref/eigen-$ver}"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT

# The second URL is an unmodified copy of the first, attached to a KEnRef release, for when gitlab.com is
# down or overloaded (on 2026-10-26 it failed ~30 CI jobs at this download).
urls=(
    "https://gitlab.com/libeigen/eigen/-/archive/$ver/eigen-$ver.tar.gz"
    "https://github.com/Smith-Group/KEnRef/releases/download/deps%2Feigen-$ver/eigen-$ver.tar.gz"
)
ok=
for round in 1 2 3; do
    for url in "${urls[@]}"; do
        echo "downloading $url"
        curl -fsSL --max-time 300 -o "$work/eigen.tar.gz" "$url" || continue
        got="$(sha256sum "$work/eigen.tar.gz" | cut -d' ' -f1)"
        if [ "$got" = "$sha" ]; then ok=1; break 2; fi
        echo "::warning::$url served a file with SHA-256 $got, not $sha"
    done
    [ "$round" = 3 ] || sleep $((15 * round))
done
[ -n "$ok" ] || { echo "::error::no URL gave eigen-$ver.tar.gz with SHA-256 $sha"; exit 1; }

mkdir "$work/src"
tar xzf "$work/eigen.tar.gz" -C "$work/src" --strip-components=1
# Header-only install: Eigen's optional reference BLAS/LAPACK take minutes to compile, KEnRef uses Eigen
# header-only, and distro packages do not ship them either.
cmake -S "$work/src" -B "$work/build" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$prefix" \
    -DEIGEN_BUILD_BLAS=OFF -DEIGEN_BUILD_LAPACK=OFF -DEIGEN_BUILD_TESTING=OFF \
    -DBUILD_TESTING=OFF -DEIGEN_BUILD_DOC=OFF -DEIGEN_BUILD_DEMOS=OFF \
    -DCMAKE_EXPORT_PACKAGE_REGISTRY=OFF
cmake --install "$work/build"
printf 'Eigen %s, archive SHA-256 %s\n' "$ver" "$sha" > "$prefix/share/eigen3/kenref-source.txt"
