#!/usr/bin/env bash
#
# Check the QML_ConfidenceEllipsoid line of NLLoc .hyp output against the
# ellipsoid on the STATISTICS line of the same event.
#
# Runs the Alaska sample (Vel2Grid -> Grid2Time -> Time2EQ -> NLLoc) and then
# LocSum, which reads the per-event .hyp files and writes them again into the
# summary file loc/alaska.hyp. For every event, in the per-event files and in
# the summary file, the major and minor axes are rebuilt from the QuakeML
# angles (majorAxisAzimuth, majorAxisPlunge, majorAxisRotation; z-y'-x''
# Tait-Bryan rotation) and compared with the STATISTICS axes:
#   major axis vs. Ell axis 1 x axis 2
#   minor axis vs. Ell axis 1 (EllAz1, Dip1)
# Axes are compared without sign (an axis and its opposite are the same).
# Also checks semiMinorAxisLength = Len1 and semiMajorAxisLength = Len3.
#
# Exit status: 0 = all ellipsoids agree; non-zero = mismatch.
#
# Usage:
#   tests/run_test_qml_ellipsoid.sh
#   NLL_BIN=/path/to/bin tests/run_test_qml_ellipsoid.sh
#
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BIN_DIR="${NLL_BIN:-${REPO_ROOT}/src/bin}"
TOL_DEG="${TOL_DEG:-1.0}"   # max angle between axes, degrees

REQUIRED_TOOLS=(Vel2Grid Grid2Time Time2EQ NLLoc LocSum)

missing=0
for t in "${REQUIRED_TOOLS[@]}"; do
    [ -x "${BIN_DIR}/${t}" ] || missing=1
done
if [ "${missing}" -ne 0 ]; then
    echo "[qml-ellipsoid] tools missing in ${BIN_DIR}, building..."
    mkdir -p "${REPO_ROOT}/src/bin"
    ( cd "${REPO_ROOT}/src" && rm -f CMakeCache.txt && cmake . >/dev/null && make -j"$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)" >/dev/null ) \
        || { echo "[qml-ellipsoid] BUILD FAILED"; exit 2; }
fi
export PATH="${BIN_DIR}:${PATH}"

WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT
cp -r "${REPO_ROOT}/nlloc_sample/." "${WORK}/"
cd "${WORK}"
rm -rf loc
mkdir -p model time obs_synth loc gmt

CTRL="run/nlloc_sample.in"
echo "[qml-ellipsoid] running pipeline in ${WORK}"
for step in "Vel2Grid ${CTRL}" "Grid2Time ${CTRL}" "Time2EQ ${CTRL}" "NLLoc ${CTRL}" \
            "LocSum loc/alaska.sum.grid0.loc 1 loc/alaska loc/alaska.*.*.grid0.loc"; do
    name="${step%% *}"
    if ! ${step} > "${WORK}/${name}.log" 2>&1; then
        echo "[qml-ellipsoid] FAIL: ${name} exited non-zero. Last lines:"
        tail -10 "${WORK}/${name}.log"
        exit 1
    fi
done

# check one .hyp file (one or more events); prints one line per event and a final "SUMMARY n nbad maxang"
check_hyp() {
    awk -v tol="${TOL_DEG}" -v label="$1" '
    function rad(d) { return d * 3.14159265358979323846 / 180.0 }
    function deg(r) { return r * 180.0 / 3.14159265358979323846 }
    # angle between two axes, ignoring sign
    function axang(ax, ay, az, bx, by, bz,   c) {
        c = (ax*bx + ay*by + az*bz) / sqrt((ax*ax + ay*ay + az*az) * (bx*bx + by*by + bz*bz))
        if (c < 0) c = -c
        if (c > 1) c = 1
        return deg(atan2(sqrt(1 - c*c), c))
    }
    function val(key,   i) { for (i = 1; i < NF; i++) if ($i == key) return $(i + 1); return "" }
    /^NLLOC / { ev = $2; gsub(/"/, "", ev); sub(/.*\//, "", ev) }
    /^STATISTICS/ {
        a1 = rad(val("EllAz1")); d1 = rad(val("Dip1")); len1 = val("Len1")
        a2 = rad(val("Az2"));    d2 = rad(val("Dip2")); len3 = val("Len3")
        # axes in x east, y north, z down
        x1 = sin(a1)*cos(d1); y1 = cos(a1)*cos(d1); z1 = sin(d1)
        x2 = sin(a2)*cos(d2); y2 = cos(a2)*cos(d2); z2 = sin(d2)
        x3 = y1*z2 - z1*y2; y3 = z1*x2 - x1*z2; z3 = x1*y2 - y1*x2
        have_stat = 1
    }
    /^QML_ConfidenceEllipsoid/ {
        if (!have_stat) next
        S = rad(val("majorAxisAzimuth")); P = rad(val("majorAxisPlunge")); T = rad(val("majorAxisRotation"))
        # major axis from azimuth and plunge
        mx = sin(S)*cos(P); my = cos(S)*cos(P); mz = sin(P)
        # minor axis: (0, cos T, sin T) rotated back by plunge and azimuth
        nx = -sin(T)*sin(P)*sin(S) + cos(T)*cos(S)
        ny = -sin(T)*sin(P)*cos(S) - cos(T)*sin(S)
        nz =  sin(T)*cos(P)
        emaj = axang(mx, my, mz, x3, y3, z3)
        emin = axang(nx, ny, nz, x1, y1, z1)
        dlen = val("semiMinorAxisLength") - len1; if (dlen < 0) dlen = -dlen
        dlen3 = val("semiMajorAxisLength") - len3; if (dlen3 < 0) dlen3 = -dlen3
        bad = (emaj > tol || emin > tol || dlen > 0.001 * len1 + 1e-6 || dlen3 > 0.001 * len3 + 1e-6)
        n++; nbad += bad
        if (emaj > maxang) maxang = emaj
        if (emin > maxang) maxang = emin
        printf("  %-8s %-40s major axis %6.2f deg  minor axis %6.2f deg%s\n", label, ev, emaj, emin, bad ? "  <-- MISMATCH" : "")
        have_stat = 0
    }
    END { printf("SUMMARY %d %d %.2f\n", n, nbad, maxang) }
    ' "$2"
}

total=0; totbad=0; maxang=0
for f in loc/alaska.[0-9]*.grid0.loc.hyp; do
    check_hyp event "${f}" > "${WORK}/_check"
    grep -v '^SUMMARY' "${WORK}/_check"
    read -r _ n nbad m < <(grep '^SUMMARY' "${WORK}/_check")
    total=$((total + n)); totbad=$((totbad + nbad))
    maxang=$(awk -v a="${maxang}" -v b="${m}" 'BEGIN { print (b > a) ? b : a }')
done
nevents=${total}
check_hyp summary loc/alaska.hyp > "${WORK}/_check"
grep -v '^SUMMARY' "${WORK}/_check"
read -r _ n nbad m < <(grep '^SUMMARY' "${WORK}/_check")
total=$((total + n)); totbad=$((totbad + nbad))
maxang=$(awk -v a="${maxang}" -v b="${m}" 'BEGIN { print (b > a) ? b : a }')

echo "----------------------------------------------------------------------"
echo "[qml-ellipsoid] ellipsoids checked: ${total} (${nevents} per-event, ${n} in LocSum summary)   tolerance: ${TOL_DEG} deg"
echo "  max axis difference=${maxang} deg"
if [ "${nevents}" -gt 0 ] && [ "${n}" -eq "${nevents}" ] && [ "${totbad}" -eq 0 ]; then
    echo "[qml-ellipsoid] PASS: QML_ConfidenceEllipsoid agrees with STATISTICS ellipsoid"
    exit 0
else
    echo "[qml-ellipsoid] FAIL: ${totbad} ellipsoid(s) disagree (per-event ${nevents}, summary ${n})"
    exit 1
fi
