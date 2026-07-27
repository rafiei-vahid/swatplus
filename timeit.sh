#!/usr/bin/env bash
ENG="$1"; THREADS="$2"; TAG="$3"
SRC=/data/SWATGenXApp/codes/swatplus_perf/wt_race/scratch/model_02270000
D=/data/SWATGenXApp/codes/swatplus_perf/wt_race/scratch/timing_$TAG
rm -rf "$D"; cp -r "$SRC" "$D"; cd "$D" || exit 1
export LD_LIBRARY_PATH=/data/SWATGenXApp/codes/lib/netcdf-ifx:/opt/intel/oneapi/compiler/2026.0/lib
export OMP_NUM_THREADS=$THREADS SWATPLUS_ROUTING_SERIAL=1
S=$(python3 -c 'import time;print(time.time())')
"$ENG" > engine.log 2>&1
E=$(python3 -c 'import time;print(time.time())')
python3 -c "print('%-26s threads=%-2s %7.1f s' % ('$TAG','$THREADS',$E-$S))"
