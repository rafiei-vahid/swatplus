#!/usr/bin/env bash
# Build one engine variant from the wt_race worktree.
#   usage: build_variant.sh <name> "<extra ifx flags>"
# Produces build_<name>/swatplus. Never touches production lib/.
set -eo pipefail
NAME="$1"; EXTRA="${2:-}"
ROOT="/data/SWATGenXApp/codes/swatplus_perf/wt_race"
BUILD="$ROOT/build_$NAME"

export SETVARS_COMPLETED=1
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1 || true
export PKG_CONFIG_PATH=/data/SWATGenXApp/codes/swatplus_perf/deps/netcdf-ifx/lib/pkgconfig:${PKG_CONFIG_PATH:-}

rm -rf "$BUILD"; mkdir -p "$BUILD"; cd "$BUILD"
cmake -G "Unix Makefiles" \
  -DCMAKE_Fortran_COMPILER=/opt/intel/oneapi/compiler/2026.0/bin/ifx \
  -DSWATPLUS_OPENMP=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_FLAGS_RELEASE="-O3 $EXTRA" \
  "$ROOT" > cmake_cfg.log 2>&1
grep -E "Found netcdf|FFLAGS|RELEASE" cmake_cfg.log | head -5
make -j"$(nproc)" > build.log 2>&1
EXE=$(find . -maxdepth 1 -type f -executable -name 'swatplus*' | head -1)
echo "BUILT: $BUILD/$(basename "$EXE")"
"$BUILD/$(basename "$EXE")" --version 2>/dev/null | head -2 || true
