#!/usr/bin/env bash
#
# Round-trip test of the NLL map projections: lat,lon -> x,y -> lat,lon
#
# Compiles test_projection_roundtrip.c with the NLL source files it needs
# (those of GRID_LIB_OBJS in src/CMakeLists.txt) into a temporary directory and
# runs it. Needs only a C compiler; the NLL programs are not used.
# See the comments in test_projection_roundtrip.c and tests/README.md.
#
# Exit status: 0 = all projections agree; non-zero = mismatch or build failure.
#
# Usage:
#   tests/run_test_projection_roundtrip.sh
#   tests/run_test_projection_roundtrip.sh -v          # print every case
#   TOL_M=0.001 tests/run_test_projection_roundtrip.sh
#   CC=clang CFLAGS="-O0 -g" tests/run_test_projection_roundtrip.sh
#
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$(cd "${SCRIPT_DIR}/../src" && pwd)"
CC="${CC:-cc}"
CFLAGS="${CFLAGS:--O3 -DNDEBUG}"    # as CMake Release

WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT

echo "[projection] compiling"
# shellcheck disable=SC2086
if ! ${CC} -std=gnu11 ${CFLAGS} -w -DGMT_VER_5 -DLOC2SSST_CLUGE \
        -I"${SRC}/include" -I"${SRC}/alomax_matrix" -I"${SRC}/geometry" \
        -I"${SRC}/matrix_statistics" -I"${SRC}/ran1" -I"${SRC}/octtree" \
        -I"${SRC}/io" -I"${SRC}/io/jReadWrite" -I"${SRC}" \
        -o "${WORK}/test_projection_roundtrip" \
        "${SCRIPT_DIR}/test_projection_roundtrip.c" \
        "${SRC}/GridLib.c" "${SRC}/util.c" "${SRC}/geo.c" "${SRC}/octtree/octtree.c" \
        "${SRC}/io/json_io.c" "${SRC}/io/jReadWrite/source/jRead.c" \
        "${SRC}/io/jReadWrite/source/jWrite.c" "${SRC}/alomax_matrix/alomax_matrix.c" \
        "${SRC}/alomax_matrix/eigv.c" "${SRC}/alomax_matrix/alomax_matrix_svd.c" \
        "${SRC}/matrix_statistics/matrix_statistics.c" "${SRC}/vector/vector.c" \
        "${SRC}/ran1/ran1.c" "${SRC}/map_project.c" -lm > "${WORK}/cc.log" 2>&1; then
    echo "[projection] BUILD FAILED"
    tail -20 "${WORK}/cc.log"
    exit 2
fi

if "${WORK}/test_projection_roundtrip" "$@"; then
    echo "[projection] PASS: lat,lon -> x,y -> lat,lon agrees for all projections"
    exit 0
else
    echo "[projection] FAIL: lat,lon -> x,y -> lat,lon does not agree"
    exit 1
fi
