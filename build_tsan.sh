#!/usr/bin/env bash
set -eo pipefail
ROOT=/data/SWATGenXApp/codes/swatplus_perf/wt_race
BUILD=$ROOT/build_tsan
export SETVARS_COMPLETED=1
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1 || true
export PKG_CONFIG_PATH=/data/SWATGenXApp/codes/swatplus_perf/deps/netcdf-ifx/lib/pkgconfig:${PKG_CONFIG_PATH:-}
rm -rf "$BUILD"; mkdir -p "$BUILD"; cd "$BUILD"
cmake -G "Unix Makefiles" \
  -DCMAKE_Fortran_COMPILER=/opt/intel/oneapi/compiler/2026.0/bin/ifx \
  -DSWATPLUS_OPENMP=ON -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_FLAGS_RELEASE="-O1 -g -fsanitize=thread -fpe3" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=thread" \
  "$ROOT" > cmake_cfg.log 2>&1
make -j"$(nproc)" > build.log 2>&1 || { tail -20 build.log; exit 1; }
find . -maxdepth 1 -type f -executable -name 'swatplus*' | head -1
