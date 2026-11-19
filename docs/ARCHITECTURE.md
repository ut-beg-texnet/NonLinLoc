# NLLoc source architecture

This document describes how the NLLoc location code is organized after
`NLLocLib.c` was split into smaller files. It is meant as an orientation map
for anyone reading the code, writing documentation, or building wrappers
around it — not as a full API reference.

## Background

Historically almost all of the NLLoc location logic lived in a single
~15,500-line file, `NLLocLib.c`. That file has been split, along its existing
sections, into nine files. The split is **pure code movement**: functions were
relocated verbatim (verified byte-for-byte), all declarations were left in
`NLLocLib.h`, and no calling code changed. The set of functions is the same as
before the split, and the output of the Alaska, global and Parkfield runs is
byte-identical (regression tests in `tests/`). In other words, this is a
reorganization, not a rewrite — results are unchanged.

`NLLocLib.c` is now ~2,700 lines and holds the core driver plus the helpers
that genuinely tie the pipeline together.

## Where the location library sits in the NLLoc toolchain

NLLoc is a suite of command-line tools. A typical run is a short pipeline:

```
Vel2Grid    velocity model     ->  3D velocity grid
Grid2Time   velocity grid      ->  per-station travel-time (and angle) grids
NLLoc       picks + TT grids   ->  located hypocentres
```

(`Time2EQ` generates synthetic arrival times for testing; `Loc2ssst` computes
source-specific station-term (SSST) corrections; `NLDiffLoc` does
double-difference relocation.)

The modules in this document make up the **location library** — the code behind
the `NLLoc` step. It is compiled once (`NLLOC_LIB_OBJS`) and linked into both
the `NLLoc` tool and the `NLDiffLoc` double-difference tool, which therefore
share the same observation-reading, search, quality and output code.

**Entry points.** `main()` lives in `NLLoc_main.c`, a thin wrapper that calls
`NLLoc()` in `NLLoc1.c`. `NLLoc()` reads the control file and, for each location
grid, calls `Locate()` in `NLLocLib.c`, which drives everything else. So the
top-to-bottom path is:

```
main (NLLoc_main.c) -> NLLoc (NLLoc1.c) -> Locate (NLLocLib.c) -> the stage modules
```

## The location library modules

All of the files below are compiled into the `NLLOC_LIB_OBJS` object library
(`src/CMakeLists.txt`) and linked into `NLLoc`, `NLDiffLoc`, `NLLoc_func_test`
and `mag_func_test`.

| File | ~lines | Responsibility | Representative functions |
|------|-------:|----------------|--------------------------|
| `NLLocLib.c` | 2,700 | Core driver and shared helpers: orchestrates a location, weighting, travel-time lookup, date/time handling, summary I/O | `Locate`, `SaveLocation`, `ConstWeightMatrix`, `CalcCenteredTimesObs/Pred`, `getTravelTimes`, `applyCrustElevCorrection`, `isAboveTopo`, `StdDateTime`, `OpenSummaryFiles`/`CloseSummaryFiles`, `setStationDistributionWeights` |
| `NLLocObs.c` | 5,500 | Reading phase observations and preparing arrivals | `GetObservations`, `GetNextObs` (the multi-format phase reader), `checkObs`, `InitializeArrivalFields`, `EvaluateArrivalAlias`, `ApplyTimeDelays`, `CalcSimpleElevCorr` |
| `NLLocInput.c` | 1,900 | Parsing the NLLoc control file | `ReadNLLoc_Input`, `is_nll_control_json`, the `GetNLLoc_*` / `GetLoc*` / `GetTimeDelays` / `GetTopoSurface` statement parsers |
| `NLLocSearch.c` | 820 | Grid-search and Metropolis (simulated-annealing) search | `LocGridSearch`, `LocMetropolis`, `GetNextMetropolisSample`, `MetropolisTest`, `SaveBestLocation` |
| `NLLocOctree.c` | 890 | Oct-tree search (the default search) | `InitializeOcttree`, `LocOctree`, `LocOctree_core`, `getOctTreeStationDensityWeight`, `GenEventScatterOcttree` |
| `NLLocQuality.c` | 2,000 | Solution-quality / likelihood evaluation — how good a trial hypocentre is | `CalcSolutionQuality` and the `_EDT` / `_OT_STACK` / `_ML_OT` / `_L1_NORM` / `_GAU_TEST` / `_GAU_ANALYTIC` variants, the ML/likelihood/variance origin-time helpers, `CalcMaxLikeOriginTime`, `CalcConfidenceIntrvl` |
| `NLLocDerived.c` | 560 | Quantities derived from a finished location | `CalcAzimuthGap`, `CalcArrivalDistances`, `CalcArrivalCounts`, `CalculateVpVsEstimate`, `CalculateMagnitude` (+ `Calc_ML_HuttonBoore`, `Calc_MD_FMAG`) |
| `NLLocOutput.c` | 1,050 | Writing hypocentres in the various output formats | `WriteHypo71`, `WriteHypoEll`, `WriteHypoInverseArchive`, `WriteHypoAlberto4`, `WriteHypoFmamp` |
| `NLLocStaStat.c` | 330 | Station-residual statistics hash table (for LOCDELAY / station corrections) | `InstallStaStatInTable`, `UpdateStaStat`, `WriteStaStatTable`, `FreeStaStatTable` |

(`WriteLocation`, the dispatcher that calls the `WriteHypo*` writers, lives in
`GridLib.c` alongside the other grid/geometry code.)

## How a location run flows through the modules

```
                 control file
                      |
                      v
             NLLocInput.c   (ReadNLLoc_Input, GetNLLoc_*)
                      |   populates configuration globals
                      v
   picks --->  NLLocObs.c   (GetObservations / GetNextObs)
                      |   arrival list
                      v
             NLLocLib.c  Locate()
              |- weights / travel times (ConstWeightMatrix, getTravelTimes)
              |- search --->  NLLocOctree.c  (LocOctree, the default)
              |               NLLocSearch.c  (grid search / Metropolis)
              |                     |   each trial hypocentre scored by
              |                     '->  NLLocQuality.c (CalcSolutionQuality*)
              |- derived quantities -->  NLLocDerived.c (gap, dist, Vp/Vs, mag)
              |- station residuals  -->  NLLocStaStat.c
              '- SaveLocation()     -->  NLLocOutput.c (+ WriteLocation, GridLib.c)
                      |
                      v
              hypocentre + output files
```

At a high level, for each event:

1. **Startup** — the control file is parsed once by `NLLocInput.c`
   (`ReadNLLoc_Input` dispatching to the `GetNLLoc_*` parsers), populating the
   configuration globals.
2. **Read observations** — `NLLocObs.c` (`GetObservations` → `GetNextObs`)
   reads the phase picks and builds the arrival list, applying aliases, time
   delays and exclusions.
3. **Set up the problem** — back in `NLLocLib.c`, `Locate` sets up the weight
   matrix (`ConstWeightMatrix`), centered times and travel-time access
   (`getTravelTimes`).
4. **Search** — depending on `LOCSEARCH`, `Locate` runs the oct-tree search
   (`NLLocOctree.c`, the default), grid search or Metropolis
   (`NLLocSearch.c`). Each trial hypocentre is scored by
   `NLLocQuality.c` (`CalcSolutionQuality*`).
5. **Finish the location** — derived quantities (azimuth gap, distances,
   Vp/Vs, magnitude) are computed by `NLLocDerived.c`; station residuals are
   accumulated into `NLLocStaStat.c`.
6. **Write results** — `SaveLocation` (in `NLLocLib.c`) drives output via
   `WriteLocation` (`GridLib.c`) and the format writers in `NLLocOutput.c`.

`Locate` and `SaveLocation` in `NLLocLib.c` are the orchestration points; the
rest of the modules are the stages they call into.

## Shared state and headers

- **`NLLocLib.h`** is the shared header: it holds the common `struct` and
  `typedef` definitions, the function prototypes for everything above, and the
  `extern` declarations for the file-scope globals. Because the prototypes live
  here, moving a function between translation units does not affect its callers.
- **Globals** are *defined* in `NLLocLib.c` and *declared* `extern` in
  `NLLocLib.h`, so any module that includes the header can reach them. NLLoc is
  still single-threaded / non-reentrant because of this shared global state;
  parallelism is achieved by running separate OS processes.
- **`NLLocLib.h` is not self-contained.** Sources must include a fixed set of
  headers *before* it, in this order:

  ```c
  #include "GridLib.h"
  #include "ran1/ran1.h"
  #include "velmod.h"
  #include "GridMemLib.h"
  #include "calc_crust_corr.h"
  #include "phaseloclist.h"
  #include "otime_limit.h"
  #include "NLLocLib.h"
  ```

  Every `NLLoc*.c` module begins with exactly this block (plus `json_io.h` and
  a `CUSTOM_ETH` guard). If you add a new module, copy that preamble.

## Adding a new module (the extraction pattern)

The eight modules above were all created the same way, and the same recipe
applies to any further split:

1. Move the function *definitions* to a new `.c` file; leave the *declarations*
   in `NLLocLib.h` so callers are untouched.
2. Start the new file with the include preamble shown above.
3. Add the file to `NLLOC_LIB_OBJS` in `src/CMakeLists.txt` and re-run CMake.
4. If the compiler reports an undeclared identifier, it is a global that was
   only reachable inside the monolith: add an `extern` for it in `NLLocLib.h`
   (if it is shared) or move its definition/`#define` into the new file (if
   only the moved code uses it).
5. Verify the move is byte-identical (diff the removed block against the new
   file) and run the regression suite before committing.
