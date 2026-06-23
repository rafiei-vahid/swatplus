# SWAT+ — SWATGenX accelerated fork

> **This is a divergent fork of [`swat-model/swatplus`](https://github.com/swat-model/swatplus),
> maintained for the SWATGenX platform.** It is an **independent engine line**: it adds
> shared-memory (OpenMP) parallelism, NetCDF output, and an in-stream + land-phase PFAS
> fate-and-transport module on top of the upstream science. It **fetches** upstream changes
> selectively but does **not** track toward upstream contribution (see *Relationship to upstream*).
> Because of the reentrancy refactor (see below) the engine **requires the Intel `ifx` compiler**;
> `gfortran` is not supported for production-scale models.

The **Soil and Water Assessment Tool Plus** [SWAT+](https://swatplus.gitbook.io/docs) is an open-source
model jointly developed by the USDA Agricultural Research Service ([USDA-ARS](http://ars.usda.gov)) and
Texas A&M AgriLife Research. SWAT+ simulates the quantity and quality of surface and ground water at
small-watershed to river-basin scale. This fork preserves that science unchanged and accelerates the
runtime.

## What this fork adds over upstream

1. **Shared-memory (OpenMP) parallelism** — a routing-aware wavefront over the daily object DAG with
   per-thread "current-object" state, giving multi-core speedup without changing results. Build option
   `SWATPLUS_OPENMP=ON`.
2. **NetCDF output backend** — per-stream NetCDF-4 output (`*_day.nc`, …) when `cdfout = y` in
   `print.prt`. Build option `SWATPLUS_NETCDF=ON`.
3. **PFAS fate-and-transport** — land-phase (`pfas_lch`, `pfas_sed`) and in-stream (`pfas_cha`)
   modules (PFOS-like Freundlich sorption, point-source injection, daily channel concentrations),
   activated only when a model carries PFAS inputs (`pfas.dat`).
4. **Reentrancy refactor** — engine-wide removal of implicit-`SAVE` locals (the Fortran
   "`var = 0` initializer ⇒ static" hazard) and per-thread scratch, so every routine is thread-safe.
   This is the change that makes the fork `ifx`-only.

## Building this fork (Intel `ifx` + NetCDF + OpenMP)

Requirements: Intel oneAPI (`ifx`), CMake, an `ifx`-built NetCDF-Fortran (`libnetcdff`), `git`.

```bash
# 1. one-time: build NetCDF-Fortran against ifx (see swatplus_perf/scripts/build_netcdf_ifx.sh)
# 2. source the Intel runtime and point pkg-config at the ifx NetCDF
source /opt/intel/oneapi/setvars.sh
export PKG_CONFIG_PATH=/path/to/netcdf-ifx/lib/pkgconfig:$PKG_CONFIG_PATH

# 3. configure + build (the two options below are what make this the production engine)
cmake -S . -B build/ifx -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_Fortran_FLAGS="-O3" \
      -DSWATPLUS_NETCDF=ON -DSWATPLUS_OPENMP=ON
cmake --build build/ifx -j"$(nproc)"
```

Runtime: `OMP_NUM_THREADS=<n>` sets the core count; `SWATPLUS_ROUTING_SERIAL=1` selects the
byte-identical (serial-routing) mode.

## Reproducibility / correctness standard

The acceleration preserves the science to a documented standard (methodology paper: *"Accelerating a
Regional Hyper-Resolution SWAT+ Model Without Changing the Science"*):

- **Byte-parity** — bit-for-bit identical output is the goal, and **thread-count invariance**
  (output independent of `OMP_NUM_THREADS`) doubles as an automatable data-race detector. `N=1`
  reproduces the original serial order exactly.
- **Documented model-equivalence** — where a parallel reduction reorders a summation, aggregates may
  differ at the last ULP (≤ 1e-3 absolute, ~1e-7 relative); the core science (flow, sediment, water
  balance, PFAS) stays byte-identical.
- **Standing gate** — `swatplus_perf/scripts/byteid_rogue_pfas.sh` runs the full SW+GW PFAS Rogue model
  at `N=1` vs `N=4` and asserts this standard before any engine is promoted to production.

## Relationship to upstream

This fork **diverges permanently** from `swat-model/swatplus`. The reentrancy refactor and `ifx`
dependency make the changes a poor fit for general SWAT+ users, so they are not mainlined. The fork
**pulls** future upstream engine changes and accepts/rejects them per-change; it does not push back.

## Directory structure & upstream docs

The CMake layout, scenario tests, and coding conventions follow upstream:

- [Configuring, Building, Installing SWAT+ using cmake](doc/Building.md)
- [Scenario Testing](doc/Testing.md) · [Tagging and Versioning](doc/Tagging.md)
- [Developing in Visual Studio](doc/VS-Win.md) · [VS Code Codespaces](doc/VSCode_Codespace.md)
- [SWAT+ Source Documentation](https://swat-model.github.io/swatplus) ·
  [SWAT+ I/O Documentation](https://swatplus.gitbook.io/docs) · [SWAT at TAMU](https://swat.tamu.edu)
