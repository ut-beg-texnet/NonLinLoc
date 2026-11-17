#!/usr/bin/env bash
#
# FULL NLL-SSST regression test (Parkfield 2004, iterative relocation).
#
# Unlike run_regression_ssst.sh (which does only the initial single-pass
# location), this runs the complete Source-Specific Station Term procedure:
#   Vel2Grid -> Grid2Time  (initial travel-time grids)
#   run_ssst_relocations.bash  (ITERATION_MAX iterations of NLLoc + Loc2ssst)
# and compares the final loc_ssst_corr<N+1> locations against A. Lomax's frozen
# macOS reference in tests/reference/ssst_full/, within a tolerance.
#
# It drives Lomax's own run_ssst_relocations.bash / run_ssst.bash unchanged
# except for two environment adaptations (not science): the parallel core count,
# and neutralising the Java SeismicityViewer launch so it runs headless.
#
# NOTE: this test is SLOW (~20 min: 383 events x 4 iterations) and STOCHASTIC
# (run_ssst.bash shuffles station order with `sort -R`), so it is an opt-in test
# -- not part of the fast CI gate. Comparison is on the maximum-likelihood
# hypocentre, which for these well-constrained regional events reproduces the
# reference far better than the expectation (see tests/lib_compare.sh).
#
# Usage:
#   tests/run_regression_ssst_full.sh
#   NLL_BIN=/path/to/bin  HTOL_M=250 ZTOL_M=700  tests/run_regression_ssst_full.sh
#
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BIN_DIR="${NLL_BIN:-${REPO_ROOT}/src/bin}"
REF="${REPO_ROOT}/tests/reference/ssst_full/Parkfield_2004.sum.grid0.loc.hyp"
source "${SCRIPT_DIR}/lib_compare.sh"
HTOL_M="${HTOL_M:-250}"   # horizontal tolerance, metres
ZTOL_M="${ZTOL_M:-700}"   # depth tolerance, metres (4 iterations accumulate more drift)
CTRL="Parkfield_2004_Oppenheimer_1993.in"

# Build the tools if any are missing (Loc2ssst is needed in addition to the
# regional set).
for t in Vel2Grid Grid2Time NLLoc Loc2ssst; do
    if [ ! -x "${BIN_DIR}/${t}" ]; then
        echo "[regression-ssst-full] building tools..."
        ( cd "${REPO_ROOT}/src" && rm -f CMakeCache.txt && cmake . >/dev/null && make -j"$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)" >/dev/null ) \
            || { echo "[regression-ssst-full] BUILD FAILED"; exit 2; }
        break
    fi
done
export PATH="${BIN_DIR}:${PATH}"

WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT
cp -r "${REPO_ROOT}/nlloc_ssst_sample/." "${WORK}/"
cd "${WORK}"
mkdir -p out/Oppenheimer1993/model out/Oppenheimer1993/time

# --- stage 1: initial travel-time grids ---------------------------------------
echo "[regression-ssst-full] building initial grids (Vel2Grid, Grid2Time)..."
for step in Vel2Grid Grid2Time; do
    if ! ${step} "${CTRL}" > "${WORK}/${step}.log" 2>&1; then
        echo "[regression-ssst-full] FAIL: ${step} exited non-zero:"; tail -8 "${WORK}/${step}.log"; exit 1
    fi
done

# --- stage 2: iterative SSST relocation ---------------------------------------
# Environment adaptations only: core count and the headless viewer.
NCORES="$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)"
sed -i -E "s/^NUM_CORES=[0-9]+/NUM_CORES=${NCORES}/" run_ssst_relocations.bash
sed -i -E 's/^SV_CMD=.*/SV_CMD="true"/'             run_ssst_relocations.bash

echo "[regression-ssst-full] running iterative SSST (${NCORES} cores; this takes ~20 min)..."
if ! bash run_ssst_relocations.bash > "${WORK}/ssst_run.log" 2>&1; then
    echo "[regression-ssst-full] FAIL: run_ssst_relocations.bash exited non-zero:"
    tail -15 "${WORK}/ssst_run.log"; exit 1
fi

# --- compare the final SSST locations against the reference -------------------
# Pick the highest-numbered loc_ssst_corr* directory (the final iteration).
FINAL_DIR="$(ls -d out/*_SSST/*/loc_ssst_corr* 2>/dev/null | sort -V | tail -1)"
OUT="$(ls "${FINAL_DIR}"/*/Parkfield_2004.sum.grid0.loc.hyp 2>/dev/null | head -1)"
if [ -z "${OUT}" ] || [ ! -f "${OUT}" ]; then
    echo "[regression-ssst-full] FAIL: no final loc_ssst_corr* summary produced"
    tail -15 "${WORK}/ssst_run.log"; exit 1
fi
echo "[regression-ssst-full] comparing $(basename "$(dirname "$(dirname "${OUT}")")")"

nref=$(count_events "${REF}"); nout=$(count_events "${OUT}")
if [ "${nref}" -ne "${nout}" ]; then
    echo "[regression-ssst-full] FAIL: event count ${nout} != reference ${nref}"; exit 1
fi

read -r n maxh maxz nbad < <(pair_by_time "${REF}" "${OUT}" ml | score_pairs "${HTOL_M}" "${ZTOL_M}")

echo "----------------------------------------------------------------------"
echo "[regression-ssst-full] events: ${n}   tolerance: H=${HTOL_M} m  Z=${ZTOL_M} m"
echo "  max horizontal=${maxh} m   max depth=${maxz} m   (ML hypocentre)"
if [ "${n}" -eq "${nref}" ] && [ "${nbad}" -eq 0 ]; then
    echo "[regression-ssst-full] PASS: all ${n} SSST locations within tolerance of reference"
    exit 0
else
    echo "[regression-ssst-full] FAIL: ${nbad} event(s) exceed tolerance (matched ${n}/${nref})"
    exit 1
fi
