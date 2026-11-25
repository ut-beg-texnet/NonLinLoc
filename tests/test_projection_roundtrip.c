/*
 * Round-trip test of the NLL map projections: lat,lon -> x,y -> lat,lon
 *
 * Uses the NLL functions get_transform(), latlon2rect() and rect2latlon()
 * (GridLib.c, map_project.c). Built and run by run_test_projection_roundtrip.sh.
 *
 * Each case is one TRANS line, as in a control file, and an area around its
 * origin. The position returned by rect2latlon() must agree with the input
 * position within the tolerance, and the returned longitude must be between
 * -180 and 180 deg (except for NONE, which does not change the values).
 *
 * LIMITS:
 * - An error made the same way in both directions cancels, so this test does
 *   not show that x,y are correct, only that the two directions agree.
 * - Each case is loaded and checked on its own (with n_proj 0 to
 *   NUM_PROJ_MAX - 1 in turn), so mixing of two transforms is not checked.
 *
 * Usage: test_projection_roundtrip [-v]      (-v: print every case)
 * Environment: TOL_M, tolerance in m (default 0.01)
 *
 * Exit status: 0 = all cases pass; 1 = failure.
 */

#include "GridLib.h"

#define HALF_WIDTH_DEG 2.0      // half width of area, deg of latitude; same width in km in longitude
#define NUM_HALF 10             // number of points on each side of origin

#define AREA_REGIONAL 0     // area around origin, longitude as given and +/- 360 deg
#define AREA_GLOBE 1        // whole globe, including longitude > 180 deg
#define AREA_CARTESIAN 2    // values in km, not lat/lon

#define MAX_CASES 1000

typedef struct {
    char trans[MAXLINE];
    int area;
    double lat0, lon0;
    int npts, nbad;
    double max_err; // m; km for AREA_CARTESIAN
    double bad_lat, bad_lon, bad_lat_out, bad_lon_out;
} Case;

static Case cases[MAX_CASES];
static int num_cases = 0;

static void add_case(int area, double lat0, double lon0, const char *trans) {

    if (num_cases >= MAX_CASES) {
        fprintf(stderr, "ERROR: too many cases, increase MAX_CASES\n");
        exit(2);
    }
    Case *pc = cases + num_cases++;
    snprintf(pc->trans, sizeof (pc->trans), "%s", trans);
    pc->area = area;
    pc->lat0 = lat0;
    pc->lon0 = lon0;
}

/* check one position */

static void check_point(int n_proj, Case *pc, double lat, double lon, double tol) {

    double x = NAN, y = NAN, lat_out = NAN, lon_out = NAN, err;
    int bad = 0;

    if (latlon2rect(n_proj, lat, lon, &x, &y) != 0)
        bad = 1;
    if (rect2latlon(n_proj, x, y, &lat_out, &lon_out) != 0)
        bad = 1;

    if (pc->area == AREA_CARTESIAN) {
        err = hypot(lat_out - lat, lon_out - lon);
    } else {
        double dlon = remainder(lon_out - lon, 360.0);
        err = 1000.0 * c111 * hypot(lat_out - lat, dlon * cos(cRPD * lat));
        if (!(lon_out >= -180.0 && lon_out <= 180.0))
            bad = 1;
    }
    if (!(err <= tol)) // also true if err is NaN
        bad = 1;

    pc->npts++;
    if (bad && pc->nbad++ == 0) {
        pc->bad_lat = lat;
        pc->bad_lon = lon;
        pc->bad_lat_out = lat_out;
        pc->bad_lon_out = lon_out;
    }
    if (!(err <= pc->max_err))
        pc->max_err = err;
}

/* check all positions of a case */

static void check_case(int n_proj, Case *pc, double tol) {

    int i, j;

    if (pc->area == AREA_GLOBE) {
        for (i = -90; i <= 90; i += 10)
            for (j = -180; j <= 360; j += 10)
                check_point(n_proj, pc, i, j, tol);
    } else if (pc->area == AREA_CARTESIAN) {
        for (i = -1000; i <= 1000; i += 100)
            for (j = -1000; j <= 1000; j += 100)
                check_point(n_proj, pc, i, j, tol);
    } else {
        double coslat = cos(cRPD * pc->lat0);
        for (i = -NUM_HALF; i <= NUM_HALF; i++) {
            for (j = -NUM_HALF; j <= NUM_HALF; j++) {
                double lat = pc->lat0 + HALF_WIDTH_DEG * i / NUM_HALF;
                double lon = pc->lon0 + HALF_WIDTH_DEG * j / NUM_HALF / coslat;
                check_point(n_proj, pc, lat, lon, tol);
                // same position with longitude in other convention
                check_point(n_proj, pc, lat, lon < 0.0 ? lon + 360.0 : lon - 360.0, tol);
            }
        }
    }
}

int main(int argc, char **argv) {

    int verbose = argc > 1 && strcmp(argv[1], "-v") == 0;
    double tol = getenv("TOL_M") != NULL ? atof(getenv("TOL_M")) : 0.01;
    char trans[MAXLINE];
    int i, j, nerr = 0;

    // Alaska sample, Parkfield SSST sample, Mediterranean sample, SE Australia,
    // S of equator, across 180 deg longitude, across 0 deg longitude
    static const double orig[][2] = {
        {61.0, -150.0}, {35.96, -120.5047}, {43.75, 7.5}, {-33.0, 151.0},
        {-6.0, 106.0}, {-18.0, 179.5}, {51.5, -0.1}
    };
    static const double rot[] = {0.0, 30.0, -45.0, 135.0};
    // all ellipsoids in map_project.c
    static const char *ellipsoid[] = {
        "WGS-84", "GRS-80", "WGS-72", "Australian", "Krasovsky", "International",
        "Hayford-1909", "Clarke-1880", "Clarke-1866", "Airy", "Bessel",
        "Hayford-1830", "Sphere", "NAD-27", "NAD-83"
    };
    static const char *types[] = {
        "GLOBAL", "NONE", "SIMPLE", "SDC", "LAMBERT", "TRANS_MERC", "AZIMUTHAL_EQUIDIST"
    };
    const int num_orig = sizeof (orig) / sizeof (orig[0]);
    const int num_rot = sizeof (rot) / sizeof (rot[0]);
    const int num_ellipsoid = sizeof (ellipsoid) / sizeof (ellipsoid[0]);
    const int num_types = sizeof (types) / sizeof (types[0]);

    SetConstants();
    message_flag = 0;

    // set up cases

    add_case(AREA_GLOBE, 0.0, 0.0, "GLOBAL");
    add_case(AREA_CARTESIAN, 0.0, 0.0, "NONE");
    for (i = 0; i < num_orig; i++) {
        double lat = orig[i][0], lon = orig[i][1];
        for (j = 0; j < num_rot; j++) {
            snprintf(trans, sizeof (trans), "SIMPLE %f %f %f", lat, lon, rot[j]);
            add_case(AREA_REGIONAL, lat, lon, trans);
            snprintf(trans, sizeof (trans), "SDC %f %f %f", lat, lon, rot[j]);
            add_case(AREA_REGIONAL, lat, lon, trans);
            snprintf(trans, sizeof (trans), "LAMBERT WGS-84 %f %f %f %f %f", lat, lon, lat - 1.0, lat + 1.0, rot[j]);
            add_case(AREA_REGIONAL, lat, lon, trans);
            snprintf(trans, sizeof (trans), "TRANS_MERC WGS-84 %f %f %f", lat, lon, rot[j]);
            add_case(AREA_REGIONAL, lat, lon, trans);
            // UseFalseEasting
            snprintf(trans, sizeof (trans), "TRANS_MERC WGS-84 %f %f %f 1", lat, lon, rot[j]);
            add_case(AREA_REGIONAL, lat, lon, trans);
            // UseFalseEasting, FalseEasting
            snprintf(trans, sizeof (trans), "TRANS_MERC WGS-84 %f %f %f 1 200000", lat, lon, rot[j]);
            add_case(AREA_REGIONAL, lat, lon, trans);
            // UseFalseEasting, FalseEasting, ScaleFactor
            snprintf(trans, sizeof (trans), "TRANS_MERC WGS-84 %f %f %f 1 500000 0.9996", lat, lon, rot[j]);
            add_case(AREA_REGIONAL, lat, lon, trans);
            snprintf(trans, sizeof (trans), "AZIMUTHAL_EQUIDIST WGS-84 %f %f %f", lat, lon, rot[j]);
            add_case(AREA_REGIONAL, lat, lon, trans);
        }
    }
    for (i = 0; i < num_ellipsoid; i++) {
        snprintf(trans, sizeof (trans), "LAMBERT %s 43.75 7.5 43.1993 44.9961 30.0", ellipsoid[i]);
        add_case(AREA_REGIONAL, 43.75, 7.5, trans);
        snprintf(trans, sizeof (trans), "TRANS_MERC %s 43.75 7.5 30.0", ellipsoid[i]);
        add_case(AREA_REGIONAL, 43.75, 7.5, trans);
        snprintf(trans, sizeof (trans), "AZIMUTHAL_EQUIDIST %s 43.75 7.5 30.0", ellipsoid[i]);
        add_case(AREA_REGIONAL, 43.75, 7.5, trans);
    }
    // LAMBERT with one standard parallel
    add_case(AREA_REGIONAL, 61.0, -150.0, "LAMBERT WGS-84 61.0 -150.0 61.0 61.0 30.0");
    add_case(AREA_REGIONAL, -33.0, 151.0, "LAMBERT WGS-84 -33.0 151.0 -33.0 -33.0 30.0");

    // check cases

    for (i = 0; i < num_cases; i++) {
        Case *pc = cases + i;
        int n_proj = i % NUM_PROJ_MAX;
        if (get_transform(n_proj, pc->trans) < 0) {
            printf("  FAIL: get_transform() error: %s\n", pc->trans);
            nerr++;
            continue;
        }
        check_case(n_proj, pc, pc->area == AREA_CARTESIAN ? tol / 1000.0 : tol);
        if (pc->nbad > 0) {
            printf("  FAIL: %s: %d of %d positions\n", pc->trans, pc->nbad, pc->npts);
            printf("        first: lat %.8f lon %.8f -> lat %.8f lon %.8f\n",
                    pc->bad_lat, pc->bad_lon, pc->bad_lat_out, pc->bad_lon_out);
            nerr++;
        } else if (verbose) {
            printf("  ok    %-72s max error %.1e m\n", pc->trans,
                    pc->area == AREA_CARTESIAN ? 1000.0 * pc->max_err : pc->max_err);
        }
    }

    // summary for each type

    printf("  %-20s %6s %10s %16s\n", "type", "cases", "positions", "max error (m)");
    for (j = 0; j < num_types; j++) {
        int nc = 0, np = 0;
        double max_err = 0.0;
        size_t len = strlen(types[j]);
        for (i = 0; i < num_cases; i++) {
            if (strncmp(cases[i].trans, types[j], len) != 0
                    || (cases[i].trans[len] != ' ' && cases[i].trans[len] != '\0'))
                continue;
            nc++;
            np += cases[i].npts;
            double err = cases[i].area == AREA_CARTESIAN ? 1000.0 * cases[i].max_err : cases[i].max_err;
            if (!(err <= max_err))
                max_err = err;
        }
        printf("  %-20s %6d %10d %16.1e\n", types[j], nc, np, max_err);
    }
    printf("  cases: %d   failed: %d   tolerance: %g m\n", num_cases, nerr, tol);

    return (nerr > 0 ? 1 : 0);
}
