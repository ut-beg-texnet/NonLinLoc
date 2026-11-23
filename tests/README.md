# NonLinLoc regression tests

## `run_regression.sh`

A dependency-light regression test that guards against unintended changes to
location results. It is the safety net for the modernization work.

### What it does

1. Builds the core tools if they are not already present in `src/bin/`.
2. Runs the core computational pipeline on the bundled `alaska` sample dataset,
   in an isolated temporary directory (with any committed `loc/` artifacts
   removed first, so only freshly produced output is checked):
   `Vel2Grid` → `Grid2Time` → `Time2EQ` → `NLLoc`.
3. Compares each freshly produced per-event hypocentre against the frozen
   reference in `nlloc_sample_test_frozen_20260625/original_output/`, within a
   tolerance (default H = 250 m, Z = 500 m).

A passing run means every located hypocentre agrees with the reference within
tolerance — i.e. behavior is preserved. The comparison is on the
maximum-likelihood (`GEOGRAPHIC`) hypocentre, which for these well-constrained
regional events reproduces across platforms to within a few tens of metres. See
[`lib_compare.sh`](lib_compare.sh) for why the estimator (ML vs. expectation) is
chosen per event class, and the global section below.

> Note: an earlier version of this test compared a single committed
> `loc/alaska.hyp` file that `NLLoc` does not actually regenerate, so it silently
> validated stale output. The test now deletes `loc/` before running and compares
> only freshly produced per-event files. (Thanks to A. Lomax for catching this.)

### Why a separate script (vs. `run_tests.bash`)

The existing `run_nlloc_sample.sh` interleaves the computation with GMT (PostScript
plots) and Java (Seismicity Viewer) visualization steps. Those require external
software and do not affect the location results. This harness runs only the
computational steps, so it can run unattended (e.g. in CI) on a machine with just a
C toolchain and CMake.

### Usage

```bash
tests/run_regression.sh
```

Use a custom tools directory:

```bash
NLL_BIN=/path/to/bin tests/run_regression.sh
```

Exit status is `0` on match and non-zero on regression (or build/run failure),
so it can be used directly as a CI gate.

## `run_regression_global.sh`

GLOBAL-mode (teleseismic) counterpart. Runs `NLLoc` on the bundled `neic` sample
(5 well-known teleseismic events: Fiji, S. East Pacific Rise, Arkansas, Morocco
2004, S. Java 2006) using the committed ak135 spherical-earth travel-time grids
in `nlloc_global_sample/taup/ak135/`, and compares the per-event output against
`tests/reference/global/`.

Those grids were generated once with the legacy `TauP_NLL` tool (TauP is Java).
They are committed so the test stays dependency-light — at test time it needs
only a C toolchain + CMake, no Java/TauP. To regenerate them, see
`nlloc_global_sample/taup/TauP_Table_NLL.sh`.

It compares the **expectation** hypocentre (the `STAT_GEOG` line) within a wider
tolerance (default H = 5 km, Z = 10 km). These teleseismic events are deep and
weakly constrained, so the maximum-likelihood point can jump between near-equal
local PDF maxima across platforms/compilers, while the expectation (the PDF's
first moment) stays stable. This is the opposite trade-off from the regional
tests above, and follows A. Lomax's guidance.

```bash
tests/run_regression_global.sh
```

## `run_regression_ssst.sh`

NLL-SSST test on a real, complex dataset: the **Parkfield 2004** example (383
events, 38 stations) from the NLL-SSST-coherence procedure. Runs
`Vel2Grid` → `Grid2Time` → `NLLoc` on the committed `nlloc_ssst_sample/` fixture
and checks every located epicentre and depth against `tests/reference/ssst/`
within a tolerance (default H = 250 m, Z = 500 m), on the ML (`GEOGRAPHIC`)
hypocentre.

NLLoc's global search is stochastic and the reference was generated on macOS, so
identical source on another platform gives numerically equivalent but not
bit-identical locations (observed max ~78 m on the ML hypocentre, vs the 500 m
grid spacing). The example's own README expects output "identical or numerically
similar". This is a platform property, not a code property — modernized and
unmodified NonLinLoc produce byte-identical output on the same machine; both
differ from the macOS reference by the same ~78 m. The tolerance sits well above
that noise and far below any real regression (which moves locations by grid
cells / km).

It takes a few minutes (383 events). The fixture is only the picks + control
files + reference (~3 MB); the full example (waveforms, QuakeML, coherence stage)
is not needed for the C-code regression.

```bash
tests/run_regression_ssst.sh                       # default H=250 m, Z=500 m
HTOL_M=150 ZTOL_M=300 tests/run_regression_ssst.sh
```

The four scripts share the location-comparison helpers in
[`lib_compare.sh`](lib_compare.sh), and each builds its parallel `make` job count
portably (`nproc` on Linux, `sysctl` on macOS). The first three run in CI via
`.github/workflows/regression.yml`, together with `run_test_qml_ellipsoid.sh`;
the full-SSST test below is opt-in (too slow for CI).

CI also builds a second time with `-funsigned-char` and runs the Alaska,
global and QML ellipsoid tests on that build. On Linux ARM (aarch64) `char` is
unsigned, so this catches code that only works when `char` is signed. An
example is storing the result of `fgetc()` in a `char` and comparing it with
`EOF`, which made `NLLoc` hang at the end of the observation file on ARM. To do
the same locally:

```bash
mkdir -p /tmp/bin-uchar
MYBIN=/tmp/bin-uchar cmake -S src -B /tmp/build-uchar -DCMAKE_C_FLAGS="-funsigned-char"
make -C /tmp/build-uchar -j4
NLL_BIN=/tmp/bin-uchar tests/run_regression.sh
```

## `run_regression_ssst_full.sh`

The **full** NLL-SSST test: the complete iterative Source-Specific Station Term
procedure, not just the initial single-pass location. It runs
`Vel2Grid` → `Grid2Time` for the initial grids, then drives A. Lomax's
`run_ssst_relocations.bash` (which calls `NLLoc` + `Loc2ssst` over several
iterations), and compares the final `loc_ssst_corr4` locations against the frozen
macOS reference in `tests/reference/ssst_full/` within a tolerance (default
H = 250 m, Z = 700 m; observed max ~141 m / ~508 m — the depth drift accumulates
over the four iterations).

Two properties make this test different:

- **Slow (~20 min)** — 383 events × ~4 iterations. It is therefore **opt-in**,
  not part of the CI gate.
- **Stochastic** — `run_ssst.bash` shuffles the station order with `sort -R`, so
  each run differs slightly. Events are matched by **origin time**, not file
  order (the parallel `NLLoc` workers append to the summary as they finish, so
  the two files list events in different orders — see `pair_by_time` in
  `lib_compare.sh`). Both facts are why the comparison is tolerance-based on the
  ML hypocentre rather than a byte diff.

The Lomax scripts are driven **unchanged** except for two environment
adaptations (not science): the parallel core count, and neutralising the Java
SeismicityViewer launch so the test runs headless.

```bash
tests/run_regression_ssst_full.sh                  # default H=250 m, Z=700 m
NLL_BIN=/path/to/bin tests/run_regression_ssst_full.sh
```

## `run_test_qml_ellipsoid.sh`

Checks the `QML_ConfidenceEllipsoid` line of the `.hyp` output against the
ellipsoid on the `STATISTICS` line of the same event. Runs the Alaska sample
(`Vel2Grid` → `Grid2Time` → `Time2EQ` → `NLLoc`) and then `LocSum`, which reads
the per-event `.hyp` files and writes them again into the summary
`loc/alaska.hyp`. For every event, in the per-event files and in the summary,
it rebuilds the major and minor axes from `majorAxisAzimuth`,
`majorAxisPlunge` and `majorAxisRotation` and compares them with the
`STATISTICS` axes (major = axis 1 × axis 2, minor = axis 1), within 1 degree,
ignoring the sign of an axis. It also checks the semi-minor and semi-major
lengths.

It was added with the fixes for two bugs it detects: a sign error in the
rotation matrix that gave a wrong `majorAxisRotation` unless the major axis was
near horizontal, and `majorAxisAzimuth` / `majorAxisPlunge` = 0 in files
written by `LocSum` (the major axis direction was not set when a `.hyp` file
was read). It takes a few seconds.

```bash
tests/run_test_qml_ellipsoid.sh
TOL_DEG=0.5 tests/run_test_qml_ellipsoid.sh
```

### Roadmap

Planned extensions:
- NLL-coherence stage (needs the Python/obspy env) — run + validate + visualise
  via the companion [NLL-SSST Studio](https://github.com/comoglu/nll-ssst-studio).
- Additional sample datasets and a larger event set.
