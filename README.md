# SWAT+ — SWATGenX engine (parallel · NetCDF · PFAS · MODFLOW 6 coupling)

> An extended build of [SWAT+](https://github.com/swat-model/swatplus) that **preserves the upstream
> science unchanged** and adds, on top of it, four production-grade capabilities developed for the
> [SWATGenX](https://swatgenx.com) platform: shared-memory (OpenMP) parallelism, a NetCDF output
> backend, land-phase and in-stream **PFAS fate-and-transport**, and a daily two-way **SWAT+ ↔
> MODFLOW 6** surface-water/groundwater coupling. Each addition is inert unless explicitly enabled, so
> a stock SWAT+ model runs here unchanged. Results track upstream except where two
> order-dependence defects were fixed (`gra` in `ch_watqual4`, `enratio` in `varinit`), which
> shift in-stream CBOD and chlorophyll-a; both are reported upstream as swat-model/swatplus#242. The reentrancy work that enables safe
> parallelism makes the engine **Intel `ifx`-only** for production-scale models (`gfortran` is fine for
> small serial runs).

The **Soil and Water Assessment Tool Plus** — [SWAT+](https://swatplus.gitbook.io/docs) — is an
open-source watershed model jointly developed by the USDA Agricultural Research Service
([USDA-ARS](http://ars.usda.gov)) and Texas A&M AgriLife Research. It simulates the quantity and
quality of surface and ground water from small-watershed to river-basin scale. This engine builds
directly on SWAT+: the hydrology, water-quality, and management science are upstream's; the work here
is an **engineering and process extension** that lets that science run at high resolution, write
analysis-ready output, and carry contaminant transport across the surface-water/groundwater interface.

---

## Has anyone parallelized SWAT+ — for both HRUs and streams? Yes: this engine.

This fork is the **SWATGenX parallel SWAT+ engine**. It parallelizes **both** phases of the daily
simulation loop, which is the part usually considered impossible to split:

- **HRU land phase** — OpenMP across HRUs with per-thread object state. Results are
  **byte-identical** to serial SWAT+ at any thread count.
- **Channel / stream routing** — an opt-in **wavefront over the full daily object dependency graph**
  (channels, reservoirs, aquifers, recall points), so independent reaches route concurrently while
  upstream→downstream order is preserved.

Measured on a dedicated AWS c8a (32 physical cores) node with a production river-basin model
(one-simulated-year benchmark, daily channel output; numbers from the manuscript below):

| Metric | Value |
|---|---|
| Peak thread speedup, full routing wavefront | **5.33× at 24 threads** (still climbing) |
| Byte-identical mode (HRU-parallel, routing serial) | 2.81× at 16 threads |
| Single-thread serial engineering gain vs stock SWAT+ | 1.67–2.47× (machine-dependent) |
| **End-to-end vs stock serial SWAT+** | **7.14× at 24 threads** (≈27 s per simulated year) |

**Honest caveats.** Production-scale builds require **Intel `ifx`** (the reentrancy refactor is
`ifx`-only at scale). The parallel routing wavefront reorders in-stream summations, so **N>1 routing
carries a documented floating-point round-off** (~1e-7 relative on channel aggregates) and is
opt-in; **HRU-parallel mode is byte-identical** and, together with fully-serial, is the production
default on [SWATGenX](https://swatgenx.com).

Benchmarks, figures, and the full write-up:
**[swatgenx.com/swat-plus-parallel-engine](https://swatgenx.com/swat-plus-parallel-engine)**.
The acceleration methodology and validation standard are described in a manuscript submitted to
*Geoscientific Model Development* (2026, under review).

---

## What this engine adds over upstream SWAT+

1. **Shared-memory (OpenMP) parallelism.** A wavefront over the daily object dependency graph, with
   per-thread "current-object" state, gives multi-core speedup on a single model. **Both parallel modes are
   byte-identical** to a serial run of the same binary -- the HRU land phase alone, and the full
   wavefront that parallelizes channel routing as well (see *Running*). This required an engine-wide reentrancy
   refactor (below). Build option `SWATPLUS_OPENMP=ON`; threads via `OMP_NUM_THREADS`.

2. **NetCDF output backend.** Per-stream NetCDF-4 output (`*_day.nc`, …) when `cdfout = y` in
   `print.prt`, in place of the wide fixed-width text files — far smaller and directly readable by
   xarray/Python for large model archives. Build option `SWATPLUS_NETCDF=ON`.

3. **PFAS fate-and-transport.** Land-phase (`pfas_lch`, `pfas_sed`) and in-stream (`pfas_cha`) modules:
   three-phase soil partitioning (aqueous, Freundlich solid-phase, and air–water-interface sorption),
   point-source injection, sediment-bound transport, and daily channel concentrations. Active only when
   a model supplies PFAS inputs (`pfas.dat` / `pfas_calib.dat`); otherwise compiled-in but dormant.

4. **SWAT+ ↔ MODFLOW 6 two-way coupling.** A daily exchange through the MODFLOW 6 library
   (`mf6_coupler.f90`): SWAT+ passes recharge **down** to a MODFLOW 6 groundwater flow-and-transport
   model and receives groundwater discharge — and groundwater-borne **PFAS** — back **up** into the
   channel network, closing a continuous surface-water + groundwater contaminant mass balance. Enabled
   per-model via the coupling control file (`mf6.con`); a recharge-multiplier knob is exposed for
   calibration. With no coupling file present the engine runs as standard SWAT+.

5. **Reentrancy refactor (enabling change).** Engine-wide removal of implicit-`SAVE` locals (the
   Fortran "`var = 0` initializer ⇒ static storage" hazard) plus per-thread scratch, making every
   routine thread-safe. This is what allows (1) and what makes the production engine `ifx`-only.

---

## Building (Intel `ifx` + NetCDF + OpenMP)

Requirements: Intel oneAPI (`ifx`), CMake, an `ifx`-built NetCDF-Fortran (`libnetcdff`), `git`.

```bash
# 1. one-time: build NetCDF-Fortran against ifx (see swatplus_perf/scripts/build_netcdf_ifx.sh)
# 2. source the Intel runtime and point pkg-config at the ifx NetCDF
source /opt/intel/oneapi/setvars.sh
export PKG_CONFIG_PATH=/path/to/netcdf-ifx/lib/pkgconfig:$PKG_CONFIG_PATH

# 3. configure + build (these two options define the production engine)
cmake -S . -B build/ifx -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_Fortran_FLAGS="-O3" \
      -DSWATPLUS_NETCDF=ON -DSWATPLUS_OPENMP=ON
cmake --build build/ifx -j"$(nproc)"
```

A stock SWAT+ build (`gfortran`, no options) still compiles and runs small serial models; the
`ifx`/OpenMP/NetCDF build is what SWATGenX runs in production.

## Running

Run the engine from inside a SWAT+ model directory. Two environment variables select the parallelism
mode; **the default is fully serial and byte-identical**:

| Mode | Environment | Result |
|---|---|---|
| Serial (default) | `OMP_NUM_THREADS=1` | byte-identical reference results |
| HRU-parallel | `OMP_NUM_THREADS=N` + `SWATPLUS_ROUTING_SERIAL=1` | **byte-identical**, multi-core speedup |
| Fully parallel | `OMP_NUM_THREADS=N` + `SWATPLUS_ROUTING_SERIAL=0` | **byte-identical**, fastest |

- `OMP_NUM_THREADS` sets the **HRU land-phase** thread count; the HRU phase is byte-identical in parallel.
- `SWATPLUS_ROUTING_SERIAL=1` keeps channel routing in serial command order; `=0` enables the parallel
  routing wavefront. Both are byte-identical to serial; the flag is a speed and cross-check choice only.
- **PFAS** activates when the model directory carries the PFAS inputs (`pfas.dat` / `pfas_calib.dat`); otherwise dormant.
- **MODFLOW 6 coupling** activates when `mf6.con` is present (and the MODFLOW 6 shared library is on the
  library path); otherwise the run is plain SWAT+.

In the [SWATGenX](https://swatgenx.com) deployment the engine is wrapped in a launcher that exposes the
same modes as a convenience CLI (serial by default):

```bash
swatplus                                              # serial, byte-identical (default)
swatplus -n 4 -hru-parallel on                        # HRU parallel, routing serial (byte-identical)
swatplus -n 8 -hru-parallel on -routing-parallel on   # both parallel (byte-identical, fastest)
```

## Correctness & reproducibility standard

The acceleration and coupling preserve the science to a documented, automatable standard:

- **Byte-parity & thread-count invariance.** Bit-for-bit identical output is the target, and output
  independence from `OMP_NUM_THREADS` doubles as a data-race detector; `N=1` reproduces the original
  serial order exactly.
- **Bitwise equivalence, not tolerance.** Every variable of every output file matches the serial run
  exactly, in both parallel modes -- verified across a five-basin fleet, 1.63e10 compared values,
  none differing. There is no ULP allowance to negotiate.
- **Standing gate.** `swatplus_perf/scripts/byteid_rogue_pfas.sh` runs the full coupled SW+GW PFAS
  Rogue River model at `N=1` vs `N=4` and asserts this standard before any engine is promoted to
  production.

## Relationship to upstream SWAT+

This engine is a **respectful extension of, not a replacement for, SWAT+**. It is built on
`swat-model/swatplus`, keeps the upstream science intact, and tracks upstream so that its scientific
updates can be incorporated. Several of the underlying improvements are general-purpose, and we have
contributed fixes back to the SWAT+ project. The larger capabilities here — the reentrancy/OpenMP
refactor, the NetCDF backend, and the PFAS and MODFLOW 6 modules — are maintained as a research engine
line that advances SWAT+ toward high-resolution and coupled contaminant-transport applications, and we
welcome collaboration with the SWAT+ developer community on bringing these advances to the wider model.

## Branches — where each piece lives

**`main` is the consolidated engine line and already contains everything described above** (OpenMP
HRU land phase + wavefront routing, NetCDF backend, PFAS fate-and-transport, MODFLOW 6 coupling).
The other branches are the development history and upstream-contribution lines:

| Branch | What it holds |
|---|---|
| `main` | **Consolidated production engine** — parallel + NetCDF + PFAS + MF6 coupling |
| `prod/engine-consolidated-20260623` | the consolidation merge point of the parallel / PFAS / NetCDF lines |
| `exp/openmp-hru-20260616` | OpenMP development line — reentrancy refactor + DAG wavefront as it evolved |
| `fix/n1-byte-identity-20260618` | N=1 byte-identity restoration (`ch_watqual4` / `et_pot` threadprivate scratch) |
| `feat/pfas-surface-water` | PFAS land-phase + in-stream development line |
| `feature/netcdf-cdfout` | NetCDF output backend development line |
| `perf/hru-read-name-index` | O(1) HRU name-index startup fix |
| `upstream-pr/*` | changes prepared as contributions to upstream `swat-model/swatplus` |

## Citing

If you use this engine, please cite SWAT+ (USDA-ARS / Texas A&M) together with the SWATGenX
publications describing the acceleration, the PFAS fate-and-transport implementation, and the SWAT+ ↔
MODFLOW 6 coupling — the acceleration manuscript is under review at *Geoscientific Model
Development* (2026). See [swatgenx.com](https://swatgenx.com) for the current reference list.

## Directory structure & upstream docs

The CMake layout, scenario tests, and coding conventions follow upstream:

- [Configuring, Building, Installing SWAT+ using cmake](doc/Building.md)
- [Scenario Testing](doc/Testing.md) · [Tagging and Versioning](doc/Tagging.md)
- [Developing in Visual Studio](doc/VS-Win.md) · [VS Code Codespaces](doc/VSCode_Codespace.md)
- [SWAT+ Source Documentation](https://swat-model.github.io/swatplus) ·
  [SWAT+ I/O Documentation](https://swatplus.gitbook.io/docs) · [SWAT at TAMU](https://swat.tamu.edu)
