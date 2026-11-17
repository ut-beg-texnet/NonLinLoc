#!/usr/bin/env bash
#
# NLL-SSST regression test (Parkfield 2004 example, initial-location stage).
#
# Runs the NonLinLoc location pipeline (Vel2Grid -> Grid2Time -> NLLoc) on the
# Parkfield 2004 dataset (383 events, 38 stations) from the NLL-SSST-coherence
# example, and checks that every located epicentre and depth agrees with the
# reference in tests/reference/ssst/ within a tolerance.
#
# Why tolerance (not byte-exact): NLLoc's global search is stochastic, and the
# reference output was generated on macOS. Identical source on a different
# platform yields locations that are numerically equivalent but not bit-identical
# (observed max difference ~80 m, far below the 500 m grid spacing). The example's
# own README states output should be "identical or numerically similar". This
# behaviour is a property of the platform/algorithm, not of any one change:
# modernized and unmodified code produce byte-identical output on the same
# machine. The tolerance below sits well above platform noise and far below any
# real regression (which shifts locations by grid cells / kilometres).
#
# Exit status: 0 = within tolerance; non-zero = regression (or build/run failure).
#
# Usage:
#   tests/run_regression_ssst.sh
#   NLL_BIN=/path/to/bin TOL_M=250 tests/run_regression_ssst.sh
#
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BIN_DIR="${NLL_BIN:-${REPO_ROOT}/src/bin}"
REF="${REPO_ROOT}/tests/reference/ssst/Parkfield_2004.sum.grid0.loc.hyp"
source "${SCRIPT_DIR}/lib_compare.sh"
HTOL_M="${HTOL_M:-${TOL_M:-250}}"   # horizontal tolerance, metres
ZTOL_M="${ZTOL_M:-500}"             # depth tolerance, metres
CTRL="Parkfield_2004_Oppenheimer_1993.in"

for t in Vel2Grid Grid2Time NLLoc; do
    if [ ! -x "${BIN_DIR}/${t}" ]; then
        echo "[regression-ssst] building tools..."
        ( cd "${REPO_ROOT}/src" && rm -f CMakeCache.txt && cmake . >/dev/null && make -j"$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)" >/dev/null ) \
            || { echo "[regression-ssst] BUILD FAILED"; exit 2; }
        break
    fi
done
export PATH="${BIN_DIR}:${PATH}"

WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT
cp -r "${REPO_ROOT}/nlloc_ssst_sample/." "${WORK}/"
cd "${WORK}"
mkdir -p out/Oppenheimer1993/model out/Oppenheimer1993/time out/Oppenheimer1993/loc

echo "[regression-ssst] running pipeline (383 events; this takes a few minutes)..."
for step in Vel2Grid Grid2Time NLLoc; do
    if ! ${step} "${CTRL}" > "${WORK}/${step}.log" 2>&1; then
        echo "[regression-ssst] FAIL: ${step} exited non-zero:"; tail -8 "${WORK}/${step}.log"; exit 1
    fi
done

OUT="${WORK}/out/Oppenheimer1993/loc/Parkfield_2004.sum.grid0.loc.hyp"
[ -f "${OUT}" ] || { echo "[regression-ssst] FAIL: no summary output produced"; exit 1; }

nref=$(count_events "${REF}")
nout=$(count_events "${OUT}")
if [ "${nref}" -ne "${nout}" ]; then
    echo "[regression-ssst] FAIL: event count ${nout} != reference ${nref}"; exit 1
fi

# Match events by origin time (the sum file's event order is not guaranteed) and
# compare the ML hypocentre within tolerance. These Parkfield events are well
# constrained, so the ML peak reproduces across platforms (<80 m) far better
# than the expectation (see tests/lib_compare.sh).
read -r n maxh maxz nbad < <(pair_by_time "${REF}" "${OUT}" ml | score_pairs "${HTOL_M}" "${ZTOL_M}")

echo "----------------------------------------------------------------------"
echo "[regression-ssst] events: ${n}   tolerance: H=${HTOL_M} m  Z=${ZTOL_M} m"
echo "  max horizontal=${maxh} m   max depth=${maxz} m   (ML hypocentre)"
if [ "${n}" -eq "${nref}" ] && [ "${nbad}" -eq 0 ]; then
    echo "[regression-ssst] PASS: all ${n} Parkfield locations within tolerance of reference"
    exit 0
else
    echo "[regression-ssst] FAIL: ${nbad} event(s) exceed tolerance (matched ${n}/${nref})"
    exit 1
fi
