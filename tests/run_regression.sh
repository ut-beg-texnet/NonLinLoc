#!/usr/bin/env bash
#
# Dependency-light regression test for NonLinLoc.
#
# Runs the core computational pipeline on the bundled "alaska" sample dataset
#   Vel2Grid -> Grid2Time -> Time2EQ -> NLLoc
# and compares the resulting location output against the frozen reference in
#   nlloc_sample_test_frozen_20260625/original_output/
#
# Unlike nlloc_sample/run_nlloc_sample.sh, this does NOT require GMT or Java:
# the visualization steps are intentionally omitted because they do not affect
# the location results being validated.
#
# Exit status: 0 = outputs match reference (behavior preserved); non-zero = regression.
#
# Usage:
#   tests/run_regression.sh            # uses src/bin tools, builds if missing
#   NLL_BIN=/path/to/bin tests/run_regression.sh
#
set -u

# --- locate repo root and tools ------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BIN_DIR="${NLL_BIN:-${REPO_ROOT}/src/bin}"

REQUIRED_TOOLS=(Vel2Grid Grid2Time Time2EQ NLLoc)

# Build the tools if any are missing.
missing=0
for t in "${REQUIRED_TOOLS[@]}"; do
    [ -x "${BIN_DIR}/${t}" ] || missing=1
done
if [ "${missing}" -ne 0 ]; then
    echo "[regression] tools missing in ${BIN_DIR}, building..."
    mkdir -p "${REPO_ROOT}/src/bin"
    ( cd "${REPO_ROOT}/src" && rm -f CMakeCache.txt && cmake . >/dev/null && make -j"$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)" >/dev/null ) \
        || { echo "[regression] BUILD FAILED"; exit 2; }
fi
export PATH="${BIN_DIR}:${PATH}"

source "${SCRIPT_DIR}/lib_compare.sh"
HTOL_M="${HTOL_M:-250}"   # horizontal tolerance, metres
ZTOL_M="${ZTOL_M:-500}"   # depth tolerance, metres (depth is less constrained)

# --- run the pipeline in an isolated work dir ---------------------------------
# Copy the sample, then DELETE any committed loc/ artifacts before running:
# NLLoc writes per-event files (loc/alaska.<ot>.grid0.loc.hyp), so a stale
# committed file must never be left in place to shadow freshly produced output.
WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT
cp -r "${REPO_ROOT}/nlloc_sample/." "${WORK}/"
cd "${WORK}"
rm -rf loc
mkdir -p model time obs_synth loc gmt

CTRL="run/nlloc_sample.in"
echo "[regression] running pipeline in ${WORK}"
for step in "Vel2Grid ${CTRL}" "Grid2Time ${CTRL}" "Time2EQ ${CTRL}" "NLLoc ${CTRL}"; do
    name="${step%% *}"
    if ! ${step} > "${WORK}/${name}.log" 2>&1; then
        echo "[regression] FAIL: ${name} exited non-zero. Last lines:"
        tail -10 "${WORK}/${name}.log"
        exit 1
    fi
done

# --- compare expectation hypocentres against the frozen reference -------------
# One per-event file per located event; match by NLLoc's deterministic filename.
REF_DIR="${REPO_ROOT}/nlloc_sample_test_frozen_20260625/original_output"
nref=0; missing=0
: > "${WORK}/_ref_loc"; : > "${WORK}/_out_loc"
for ref in "${REF_DIR}"/alaska.[0-9]*.grid0.loc.hyp; do
    nref=$((nref + 1))
    out="${WORK}/loc/$(basename "${ref}")"
    if [ ! -f "${out}" ]; then
        echo "[regression] FAIL: expected output $(basename "${ref}") not produced"
        missing=1
        continue
    fi
    _loc_lines "${ref}" ml >> "${WORK}/_ref_loc"
    _loc_lines "${out}" ml >> "${WORK}/_out_loc"
done
[ "${missing}" -eq 0 ] || exit 1

read -r n maxh maxz nbad < <(paste "${WORK}/_ref_loc" "${WORK}/_out_loc" | score_pairs "${HTOL_M}" "${ZTOL_M}")

echo "----------------------------------------------------------------------"
echo "[regression] events: ${n}/${nref}   tolerance: H=${HTOL_M} m  Z=${ZTOL_M} m"
echo "  max horizontal=${maxh} m   max depth=${maxz} m   (ML hypocentre)"
if [ "${n}" -eq "${nref}" ] && [ "${nbad}" -eq 0 ]; then
    echo "[regression] PASS: all ${n} locations within tolerance of frozen reference"
    exit 0
else
    echo "[regression] FAIL: ${nbad} event(s) exceed tolerance (matched ${n}/${nref})"
    exit 1
fi
