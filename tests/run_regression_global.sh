#!/usr/bin/env bash
#
# GLOBAL-mode (teleseismic) regression test for NonLinLoc.
#
# Runs NLLoc on the bundled global "neic" sample (5 well-known teleseismic
# events) using the committed ak135 spherical-earth travel-time grids, and
# compares the per-event location output against the frozen reference in
#   tests/reference/global/
#
# The ak135 grids in nlloc_global_sample/taup/ak135/ were generated once with
# the legacy TauP_NLL tool; committing them keeps this test dependency-light
# (no Java / TauP needed at test time, only a C toolchain + CMake).
#
# Exit status: 0 = outputs match reference (behavior preserved); non-zero = regression.
#
# Usage:
#   tests/run_regression_global.sh
#   NLL_BIN=/path/to/bin tests/run_regression_global.sh
#
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BIN_DIR="${NLL_BIN:-${REPO_ROOT}/src/bin}"
REF_DIR="${REPO_ROOT}/tests/reference/global"

# Build NLLoc if missing.
if [ ! -x "${BIN_DIR}/NLLoc" ]; then
    echo "[regression-global] NLLoc missing in ${BIN_DIR}, building..."
    mkdir -p "${REPO_ROOT}/src/bin"
    ( cd "${REPO_ROOT}/src" && rm -f CMakeCache.txt && cmake . >/dev/null && make -j"$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)" >/dev/null ) \
        || { echo "[regression-global] BUILD FAILED"; exit 2; }
fi
export PATH="${BIN_DIR}:${PATH}"

source "${SCRIPT_DIR}/lib_compare.sh"
# Teleseismic locations are coarse-grid and depth is weakly constrained, so the
# tolerances are wider than the regional tests. Compared on the expectation
# hypocentre, which is stable where the maximum-likelihood point is not.
HTOL_M="${HTOL_M:-5000}"    # horizontal tolerance, metres
ZTOL_M="${ZTOL_M:-10000}"   # depth tolerance, metres

# Run in an isolated work dir (drop any committed loc/ so nothing shadows output).
WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT
cp -r "${REPO_ROOT}/nlloc_global_sample/." "${WORK}/"
cd "${WORK}"
rm -rf loc
mkdir -p loc

echo "[regression-global] running NLLoc in ${WORK}"
if ! NLLoc run/neic_global.in > "${WORK}/NLLoc.log" 2>&1; then
    echo "[regression-global] FAIL: NLLoc exited non-zero. Last lines:"
    tail -10 "${WORK}/NLLoc.log"
    exit 1
fi

# Compare each frozen reference event's expectation hypocentre against output.
nref=0; missing=0
: > "${WORK}/_ref_loc"; : > "${WORK}/_out_loc"
for ref in "${REF_DIR}"/global.[0-9]*.grid0.loc.hyp; do
    nref=$((nref + 1))
    name="$(basename "${ref}")"
    out="${WORK}/loc/${name}"
    if [ ! -f "${out}" ]; then
        echo "[regression-global] FAIL: expected output ${name} not produced"
        missing=1
        continue
    fi
    _loc_lines "${ref}" expect >> "${WORK}/_ref_loc"
    _loc_lines "${out}" expect >> "${WORK}/_out_loc"
done
[ "${missing}" -eq 0 ] || exit 1

read -r n maxh maxz nbad < <(paste "${WORK}/_ref_loc" "${WORK}/_out_loc" | score_pairs "${HTOL_M}" "${ZTOL_M}")

echo "----------------------------------------------------------------------"
echo "[regression-global] events: ${n}/${nref}   tolerance: H=${HTOL_M} m  Z=${ZTOL_M} m"
echo "  max horizontal=${maxh} m   max depth=${maxz} m   (expectation hypocentre)"
# (teleseismic: expectation is more stable than the ML peak across platforms)
if [ "${n}" -eq "${nref}" ] && [ "${nbad}" -eq 0 ]; then
    echo "[regression-global] PASS: all ${n} teleseismic locations within tolerance of reference"
    exit 0
else
    echo "[regression-global] FAIL: ${nbad} event(s) exceed tolerance (matched ${n}/${nref})"
    exit 1
fi
