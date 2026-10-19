# SWAT+ fork: shared-memory parallel engine, NetCDF output and MODFLOW 6 coupling

This repository is a fork of [swat-model/swatplus](https://github.com/swat-model/swatplus), the Soil and Water
Assessment Tool Plus (SWAT+). It keeps upstream's licence (GNU LGPL 2.1, see `LICENSE`). The additions and fixes are
listed below, each with the commit or branch that carries it.

`main` holds this fork's production engine release, v2.2.0 (tag `ship/v1.0-51.5`, merged into `main` as `5debb34` on
2026-10-16). It is based on upstream SWAT+ 62.0.0 plus two commits (merge base `cb442f7`). On 2026-10-19 `main` was
181 commits ahead of upstream `main` and 86 behind.

## What the fork adds

| Addition | Source | Introduced in (all on `main`) | How it is switched on |
|---|---|---|---|
| Shared-memory parallel engine. HRUs, channels and the other spatial objects are scheduled as one dependency graph and run level by level with OpenMP. | `src/command.f90`, `src/command_wave_module.f90` | `82b7558` (HRU wave index, 2026-06-16), `25ef8cf` (wavefront over all objects, 2026-06-17), `857b9c1` (serial-routing mode, 2026-06-23) | build with `-DSWATPLUS_OPENMP=ON`; run with `OMP_NUM_THREADS` and `SWATPLUS_ROUTING_SERIAL` |
| Reentrancy changes the parallel engine needs: implicit-`SAVE` locals removed across the engine, per-thread scratch | files across `src/` | `cc7cccd` (2026-07-27), `8a3bf09` (2026-08-03) and the other `fix(omp)` commits | always on |
| NetCDF-4 output | `src/netcdf_output_module.f90`, `src/netcdf_writer_module.f90` | `82d0d86` and `2948849` (2026-05-28); later fixes include `9482c07` and `2f881ab` | build with `-DSWATPLUS_NETCDF=ON` (the default); run with `cdfout = y` in `print.prt` |
| Channel print filter: `channel_sd` output for a listed subset of channels | `src/output_filter_module.f90` | `56e1f46` (2026-05-28) | a `print_filter.prt` file in the model folder |
| Daily two-way SWAT+ ↔ MODFLOW 6 coupling: recharge passed down, groundwater discharge passed back to the channels | `src/mf6_coupler.f90` | `18a7722` (2026-06-25), `6a4b7f1` (recharge multiplier) | an `mf6.con` file in the model folder starts MODFLOW 6; `mf6_recharge.map` and `mf6_baseflow.map` switch on the two directions (below) |
| PFAS fate and transport: land phase, channels, reservoirs, and transport through the MODFLOW 6 coupling | `src/pfas_*.f90`, `src/cha_pfas_output.f90` | `e421811`, `373a688`, `cb6e353` (2026-06-21), `adae327` (2026-06-25) | a `pfas.dat` file in the model folder; into MODFLOW 6 with `pfas_leach.map` (below) |
| Fixes to upstream's water-temperature module | `src/ch_temp.f90` | `237348b` (the reported channel temperature carries the energy balance), `60bf61c`, `ae82b51` (channel outflow carries the computed temperature downstream; tributary mixing is flow-weighted over all inflows), `407730d` | always on |

Not yet on `main`:

- Later MODFLOW 6 coupler fixes are on `engine/q1137-sfr-drive`, which contains `engine/q1144-lockstep` and
  `engine/q1040c-coupler-sync`, and on `ship/q1130-q1139-combined` (tag `ship-v1.0-51.2-q1130-q1139`).
- Four salt-module fixes are prepared for upstream on `pr/salt-start`, `pr/salt-hru-gwup`, `pr/salt-ion-units` and
  `pr/salt-recall-read` (upstream pull requests #313, #317, #318 and #319).
- `main` does not read salt recall loads: as on upstream, `recall_read_salt` is never called. `pr/salt-recall-read`
  (#319) adds the call. Two fixes to that routine, `8f3d382` and `87ea2fc`, are already on `main` and take effect once
  it is called.

## Build

Requirements:

- CMake 3.22 or later.
- A Fortran compiler. Intel `ifx` is the reference toolchain; `gfortran` also builds the engine.
- For the parallel build, OpenMP.
- For NetCDF output, NetCDF-Fortran built with the same compiler and visible to `pkg-config` as `netcdf-fortran`.

```bash
# parallel build (ifx, Release, NetCDF, OpenMP)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=ifx \
      -DSWATPLUS_NETCDF=ON -DSWATPLUS_OPENMP=ON
cmake --build build --parallel

# serial build (SWATPLUS_OPENMP defaults to OFF)
cmake -S . -B build-serial -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=ifx
cmake --build build-serial --parallel
```

Add `-DSWATPLUS_NETCDF=OFF` to build without NetCDF. The executable is named
`swatplus-<git describe>-<ifx|gnu|ifo>-<platform>-Rel`, for example `swatplus-ship-v1.0-51.5-…-ifx-lin_x86_64-Rel`.

Notes:

- The timings below used `ifx` 2026.0.0 in the CMake Release configuration (`-O2`, `-fpe0`, no interprocedural
  optimisation), with OpenMP through `-fiopenmp`.
- On a single core, use the serial build. The OpenMP build run at one thread took 1.18 to 1.33 times the serial
  build's wall time across seven watersheds.
- With `gfortran`, a large model stops at run time with an end-of-record error while writing `area_calc.out`. The fix
  (#301) is merged into upstream's `dev` branch; neither upstream `main` nor this fork's `main` carries it yet.
- Testing is as upstream: see `doc/Testing.md`.

## Run

Run the executable from the model's `TxtInOut` folder. The thread count and one environment variable choose the mode:

| Mode | Environment | Agreement with a one-thread run of the same executable |
|---|---|---|
| One thread | `OMP_NUM_THREADS=1` | the reference; at one thread the engine takes the original serial code path |
| Default | `OMP_NUM_THREADS=N` (N ≥ 2) | bit for bit in every NetCDF output file the release gate read (five fixtures); two outputs the gate does not test, the pesticide channel outputs and the PFAS balance diagnostic, carry known data races |
| Serial routing | `OMP_NUM_THREADS=N` and `SWATPLUS_ROUTING_SERIAL=1` | bit for bit in every output file read in each test so far; the HRU phase runs in parallel and channel routing in command order |

- If `OMP_NUM_THREADS` is not set, the OpenMP runtime chooses the thread count.
- `SWATPLUS_ROUTING_SERIAL=1` selects serial routing; any other value, or none, leaves the default mode.
- The all-object wavefront is not used when water allocation (`wallo`/`mallo`) or `gwflow` is active
  (`src/command.f90:116-117`). With water allocation, the HRU phase is not parallel either, so the run is serial
  at any thread count; with `gwflow`, the HRU phase runs in parallel and routing is serial
  (`src/command.f90:194-195`).
- `SWATPLUS_NC_DEFLATE=<0-9>` sets the NetCDF deflate level.
- `SWATPLUS_SCHED_AUDIT=1` (or `y`) writes a census of scheduling conflicts to `conflict_audit.out` and stops before
  the first simulated day; it is a diagnostic, not a run mode.
- The carbon model (`cswat = 2`) has not been tested in parallel; run it at one thread.

### MODFLOW 6 coupling (`mf6.con` and the map files)

MODFLOW 6 starts beside SWAT+ when `mf6.con` is present in `TxtInOut`. One value per line, `!` starts a comment:

1. MODFLOW 6 workspace: the folder holding `mfsim.nam`, relative to `TxtInOut` (for example `./mf6`)
2. flow step in days (default 1)
3. transport step in days (default 30)
4. recharge multiplier (optional, default 1.0)
5. path to `libmf6.so` (optional; set it, because the built-in default is an absolute path on this fork's own servers)

If `libmf6.so` cannot be loaded or MODFLOW 6 does not initialise, the run prints one `MF6 COUPLER:` line and continues
without coupling. Check that standard output carries `MF6 COUPLER: active`.

`mf6.con` alone passes nothing between the models. Each exchange needs its own map file in `TxtInOut`; without it,
that exchange is off and the run goes on:

| File | Exchange | Without it |
|---|---|---|
| `mf6_recharge.map` | SWAT+ soil percolation to MODFLOW 6 recharge, area-weighted from HRUs to cells | MODFLOW 6 uses its own recharge |
| `mf6_baseflow.map` | groundwater–stream exchange of MODFLOW 6's SFR reaches to SWAT+ channels | nothing is passed back to the channels, PFAS included |
| `pfas_leach.map` | PFAS leached from HRUs to MODFLOW 6 transport source cells | no PFAS is passed from SWAT+ to the groundwater model |

Each file is a header line of counts followed by one link per line; the readers are `mf6_recharge_setup`,
`mf6_baseflow_setup` and `mf6_pfas_src_setup` in `src/mf6_coupler.f90`. Standard output names each map that was read
(`recharge map:`, `baseflow map:`, `pfas leach map:`).

`libmf6.so` is loaded at run time with `dlopen` and `RTLD_LOCAL`, not linked at build time, so that its bundled Fortran
runtime does not replace the engine's own.

## Correctness standard

A parallel run must reproduce a one-thread run of the same executable bit for bit, in every output file it writes
(NetCDF compared by value, text files without their date lines).

Where the current engine stands:

- **Serial-routing mode** met the standard in every file each test read: seven watersheds at four threads over one
  year, and two models at two threads.
- **Default mode** met it in every NetCDF output file the release gate read. Two outputs are known exceptions, the
  pesticide channel outputs and the PFAS balance diagnostic, where two data races remain. They are named exceptions,
  not tolerances.
- **Release gate.** Each parallel run is compared with the same binary's one-thread run in every NetCDF output file. For
  this release, all 27 parallel runs in each mode, at two, four and eight threads on five gate fixtures, matched with
  no differing cell (three thread counts times nine trials: three on each of two canary fixtures, one on each of the
  other three). The gate does not test the two exceptions above.

Working to this standard found four data races and four order dependences (results that depended on the order in
which spatial objects are visited). Two of the order dependences are defects of released sequential SWAT+; their
fixes are merged into upstream's `dev` branch (#294, #295).

Results from this fork also differ from upstream where its corrections change serial results. Two corrections in
the current release change simulated flow: `fd4f15a` (a zero-height row in `weir.res` no longer sets every wetland's
principal volume to zero) and `8ed18f9` (restores the use limit on the plant decision-table action, so a crop is not
replanted after harvest).

## Performance

Measured in October 2026 on tag `ship/v1.0-51.2` (`ae82b51`): watershed 03100101 (57,998 HRUs, 8,181 channels), one
simulated year, daily `channel_sd` output for all channels, threads pinned to cores. Each speedup is the instance's
own one-thread to N-thread ratio in the default mode. Two paired tests on one server put the one-to-eight-thread
speedup of the release on `main` within 0.7 % of that revision's. These timings are not in the archive cited under
Citing, which holds an earlier engine (`ship/v1.0-49`).

| Instance | Cores | Peak speedup |
|---|---|---|
| AWS c8a.8xlarge | 32, no SMT | 5.53× at 24 threads |
| AWS c8i.8xlarge | 16, with SMT | 3.99× at 16 threads (lower beyond) |

- Serial-routing mode was 56 % (c8a) and 61 % (c8i) slower than the default mode at eight threads.
- Over multi-year runs the speedup at 24 threads on c8a approaches 7.0 to 7.7× on three watersheds; on each of them a
  one-year run reports 15 to 22 % less.
- The serial build of the same revision ran in 129.8 s on c8a and 127.4 s on c8i, against 197.3 s and 276.7 s for
  R0 (1.52× and 2.17×). R0 is the first build of this fork, which already had the NetCDF backend and the print
  filter; it is not unmodified upstream SWAT+.

## Branches and tags

Read with `git ls-remote` on 2026-10-19: 118 branches and 43 tags.

| Branches | Count | Contents |
|---|---|---|
| `main` | 1 | the current engine release |
| `pr/*` | 52 | one topic each, prepared as an upstream pull request; 44 are a single commit on an upstream `main` commit |
| `fix/*` | 26 | fixes |
| `release/*` | 11 | the arms of the current release; `release/r1-all` is the merged release, and `release/arm-g3-no2out` and `release/arm-g5ctl-usle` are test arms that were never released |
| `engine/*` | 8 | engine work, including the MODFLOW 6 coupler fixes not yet on `main` |
| `ship/*` | 4 | deployed or combined engine builds |
| `upstream-pr/*` | 3 | earlier upstream pull-request branches |
| `feat/*` | 3 | features |
| `claude/*` | 2 | upstream contribution candidates |
| `test/*`, `sync/*`, `research/*`, `perf/*` | 1 each | a test fixture (never merged), an upstream sync, research, performance work |
| `Development`, `jga_develop`, `gfortranchanges041924`, `CMakePreset_json` | 4 | from 2024, before this fork's work |

Tags: `ship/v1.0-*` mark deployed engine builds (`ship/v1.0-51.5` is engine v2.2.0, `3a07e3b`; `ship/v1.0-51.3` is
v2.1.0, `fde3d8e`); `candidate/*` mark
release candidates; the 14 `61.0*` tags are upstream's.

## Relationship to upstream

Upstream changes reach `main` by merge; the latest is `bcb882c` (2026-07-30). A later sync,
`sync/upstream-2026-09-24`, is not yet merged into `main`. Fixes go back to upstream as small pull requests, one topic
each.

Merged upstream (10). Upstream merges pull requests into its `dev` branch; two of the ten are also on its `main`.

| # | Title | On upstream |
|---|---|---|
| [#220](https://github.com/swat-model/swatplus/pull/220) | perf(varinit): reset only the active HRU's hhsedy row in the daily init | `dev`, `main` |
| [#236](https://github.com/swat-model/swatplus/pull/236) | Fix infinite loop in mgt_operatn when a schedule's operations all share one date | `dev`, `main` |
| [#294](https://github.com/swat-model/swatplus/pull/294) | varinit: reset the hru_module enratio instead of a shadowing local | `dev` |
| [#295](https://github.com/swat-model/swatplus/pull/295) | ch_watqual4: reset the algal growth rate gra before the growth-option select | `dev` |
| [#297](https://github.com/swat-model/swatplus/pull/297) | res_sediment: release sediment at or below nsed | `dev` |
| [#301](https://github.com/swat-model/swatplus/pull/301) | proc_bsn: open area_calc.out without a record length | `dev` |
| [#302](https://github.com/swat-model/swatplus/pull/302) | ch_temp: do not read hin_d at index 0 for the first channel inflow | `dev` |
| [#326](https://github.com/swat-model/swatplus/pull/326) | Free the tillage work arrays also when the tillage depth is 0 | `dev` |
| [#336](https://github.com/swat-model/swatplus/pull/336) | Report an irrigation operation that is missing from irr.ops | `dev` |
| [#349](https://github.com/swat-model/swatplus/pull/349) | Report a fertilizer that is missing from fertilizer.frt | `dev` |

On 2026-10-19, 42 further pull requests from this fork were open (#307 to #357), and 7 had been closed without
merging, among them the NetCDF backend (#213) and the print filter (#214).

## Citing

- **SWAT+:** cite upstream as the `CITATION.cff` on its `main` asks: SWAT+ 62.0.0, the release this fork is based on
  (Zenodo, https://doi.org/10.5281/zenodo.22919722).
- **Parallel engine:** Rafiei, V. (2026). *Reproducibility archive: Bitwise-Exact Parallelisation of Semi-Distributed
  Watershed Models Reveals Order Dependence in Sequential Simulation* (2.0.0). Zenodo.
  https://doi.org/10.5281/zenodo.21524100. This version archives the engine at tag `ship/v1.0-49`.
- **MODFLOW 6 coupling:** Rafiei, V. (2026). *Coupled SWAT+ ↔ MODFLOW 6 PFAS fate & transport in surface water and
  groundwater — reproducibility archive (Rogue River, Michigan)*. Zenodo. https://doi.org/10.5281/zenodo.20838388 (all
  versions, as the archive asks).

## Licence

GNU Lesser General Public License, version 2.1, as upstream (`LICENSE`).

Maintained by SWATGenX.
