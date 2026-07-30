# SWAT+ — SWATGenX engine (parallel · NetCDF · PFAS · MODFLOW 6 coupling)

> An extended build of [SWAT+](https://github.com/swat-model/swatplus) that **preserves the upstream
> science unchanged** and adds, on top of it, four production-grade capabilities developed for the
> [SWATGenX](https://swatgenx.com) platform: shared-memory (OpenMP) parallelism, a NetCDF output
> backend, land-phase and in-stream **PFAS fate-and-transport**, and a daily two-way **SWAT+ ↔
> MODFLOW 6** surface-water/groundwater coupling. Each addition is inert unless explicitly enabled, so
> a stock SWAT+ model runs here unchanged. Results track upstream except where two
> order-dependence defects were fixed -- `gra` (`ch_watqual4`) and `enratio` (`varinit`), which shift
> in-stream CBOD and chlorophyll-a and are reported upstream as swat-model/swatplus#242, and
> `pl_nut_demand` (a stale `ipl` index), which shifts six phosphorus variables. Production-scale runs use **Intel `ifx`**;
> `gfortran` builds the engine but is not the production toolchain.

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
  **byte-identical** to a serial run of the same binary, at any thread count.
- **Channel / stream routing** — an opt-in **wavefront over the full daily object dependency graph**
  (channels, reservoirs, aquifers, recall points), so independent reaches route concurrently while
  upstream→downstream order is preserved.

Measured on a dedicated AWS c8a (32 physical cores) node with a production river-basin model
(one-simulated-year benchmark, daily channel output; source data in
`publication/engine-acceleration/repro/results/crosshw_c8a.csv`):

| Metric | Value |
|---|---|
| Peak thread speedup, full routing wavefront | **5.33× at 24 threads** (still climbing) |
| HRU-parallel mode, routing serial | 2.81× at 16 threads |
| Single-thread serial engineering gain vs stock SWAT+ | 1.67–2.47× (machine-dependent) |
| **End-to-end vs stock serial SWAT+** | **7.14× at 24 threads** (≈27 s per simulated year) |

**Honest caveats.** Production builds use **Intel `ifx`**; `gfortran` compiles the engine but is
not what SWATGenX runs at scale. **Both parallel modes are byte-identical** to serial -- the earlier
statement that the routing wavefront carried an unavoidable round-off was wrong and is retracted; the
disagreement was three order-dependence defects, since fixed. Carbon (`cswat = 2`) is **not yet
supported in parallel**: upstream's carbon routines are not reentrant, so enable carbon only at
`OMP_NUM_THREADS=1`. Results differ from upstream SWAT+ by design -- see *Relationship to upstream*.

Benchmarks, figures, and the full write-up:
**[swatgenx.com/swat-plus-parallel-engine](https://swatgenx.com/swat-plus-parallel-engine)**.
The acceleration methodology and the bitwise validation standard are written up at that page,
with a walkthrough of the routing wavefront and the measured scaling on YouTube:
**[youtu.be/ltFvXS6ISGY](https://youtu.be/ltFvXS6ISGY)**.

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

`gfortran` builds the full engine (verified on Linux and Apple Silicon); the `ifx`/OpenMP/NetCDF
build is what SWATGenX runs in production. Note upstream's CMakeLists passes `-finit-local-zero` on
the `gfortran` path, which zeroes uninitialized locals -- useful for stability, but it also masks
read-before-write defects that the `ifx` build will still hit.

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
- **What that fleet does and does not cover.** Those basins were chosen for reservoir and wetland
  density, to stress the defect class known at the time. They are all SWATGenX-generated with carbon
  disabled, so they do not exercise carbon (`cswat = 2`), which is not reentrant and is unsupported
  in parallel. A fourth order-dependence defect (`pl_nut_demand`, stale `ipl`) was later found by an
  independent build on a second toolchain, not by this gate -- the gate is necessary, not sufficient.
- **Standing gate.** `swatplus_perf/scripts/byteid_rogue_pfas.sh` runs the full coupled SW+GW PFAS
  Rogue River model at `N=1` vs `N=4` and asserts this standard before any engine is promoted to
  production.

## Relationship to upstream SWAT+

This engine is a **respectful extension of, not a replacement for, SWAT+**. It is built on
`swat-model/swatplus`, keeps the upstream science intact, and tracks upstream so that its scientific
updates can be incorporated -- `main` currently carries all upstream commits.

**It is not bug-for-bug identical to upstream, by design.** Requiring a parallel run to match a
serial run bit for bit forces order-dependence defects to be fixed, and fixing them changes results:
`gra` (`ch_watqual4`) and `enratio` (`varinit`) shift in-stream chlorophyll-a and CBOD, and
`pl_nut_demand` shifts six phosphorus variables. Those two goals -- bitwise parallel/serial agreement
and bit-for-bit agreement with the sequential original -- are mutually exclusive. Calibrated
parameters, phosphorus especially, will not transfer unchanged. Several of the underlying improvements are general-purpose and have been
reported upstream (swat-model/swatplus#242 documents three order-dependence defects found by
requiring bitwise agreement); none have been merged upstream yet. The larger capabilities here — the reentrancy/OpenMP
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

If you use this engine, please cite SWAT+ (USDA-ARS / Texas A&M) alongside the SWATGenX write-ups
covering the acceleration, the PFAS fate-and-transport implementation, and the SWAT+ ↔ MODFLOW 6
coupling. The method and benchmarks are documented at
**[swatgenx.com/swat-plus-parallel-engine](https://swatgenx.com/swat-plus-parallel-engine)** and
explained in [this video](https://youtu.be/ltFvXS6ISGY). See
[swatgenx.com](https://swatgenx.com) for the current reference list.

## Directory structure & upstream docs

The CMake layout, scenario tests, and coding conventions follow upstream:

- [Configuring, Building, Installing SWAT+ using cmake](doc/Building.md)
- [Scenario Testing](doc/Testing.md) · [Tagging and Versioning](doc/Tagging.md)
- [Developing in Visual Studio](doc/VS-Win.md) · [VS Code Codespaces](doc/VSCode_Codespace.md)
- [SWAT+ Source Documentation](https://swat-model.github.io/swatplus) ·
  [SWAT+ I/O Documentation](https://swatplus.gitbook.io/docs) · [SWAT at TAMU](https://swat.tamu.edu)
