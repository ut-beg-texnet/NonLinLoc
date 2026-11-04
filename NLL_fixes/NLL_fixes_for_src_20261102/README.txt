NonLinLoc - "real fixes" applied to your src_20261102
=====================================================

These are your own src files (the zip you sent on 2 Nov 2026) with the fixes
from commit 1472e61 added, and nothing else. All your local changes are kept
as they were:
  - NLLocLib.c: the SP-TEST PATCH blocks and the 20261026 Vp/Vs PS-pair loop fix
  - Loc2ssst.c: the 20260924 negative SSST travel-time check
  - vector/vector.h: not touched

HOW TO USE
----------
Copy the 12 files under src/ over the same files in your src/ (same
sub-directories). Then review them in NetBeans against your version.
If you prefer, fixes_vs_your_src_20261102.patch contains the same changes as one
diff against your files (run "patch -p1 < fixes_vs_your_src_20261102.patch"
from the directory that contains src/).

The patch applies cleanly to your tree, with no conflicts, and the line
numbers below are those in your files.


WHAT I CHECKED
--------------
I built your original src and the patched src separately (Linux, gcc, cmake),
then ran both on the same data:

  - Build: 0 errors with both. The 31 "-Wunused-result" warnings in your
    version are gone; no new warnings.
  - Alaska sample (Vel2Grid, Grid2Time, Time2EQ, NLLoc), global teleseismic
    sample, and Parkfield 2004 (383 events, Vel2Grid, Grid2Time, NLLoc):
    every output file (model and time grids, .hyp, scatter, summary files) is
    byte-for-byte identical between your version and the patched version.
    The only lines that differ are the SIGNATURE run date/time.
  - Full iterative SSST on Parkfield (NLLoc + Loc2ssst, 4 iterations):
    all 383 events within tolerance of the reference for both versions, with
    identical results (max difference 55 m horizontal, 78 m depth for both;
    this run includes a random station shuffle, so it is compared within a
    tolerance rather than byte-for-byte).

As expected, the fixes change nothing when files read correctly and memory
allocation succeeds. They only change what happens on a failure (truncated or
missing file, failed malloc, failed shell command): instead of silently
continuing with garbage or crashing, the code reports an error and stops or
skips.


THE CHANGES
-----------

THE TWO MOST IMPORTANT

1. NLDiffLoc.c, ReadNLDiffLoc_Input(), line 775
   Stack buffer overflow risk: fgets is allowed to read 4*MAXLINE characters
   into `line`, which is smaller (MAXLINE_LONG).
     OLD:  while ((fgets_return = fgets(line, 4 * MAXLINE, fp_input)) != NULL
     NEW:  while ((fgets_return = fgets(line, sizeof (line), fp_input)) != NULL

2. alomax_matrix/alomax_matrix.c, two places
   Allocation failure is never detected: a pointer is tested with "< 0"
   (malloc returns NULL on failure, never a negative value).
   a) matrix_double_inverse(), line 214
     OLD:  if ((augmented_mtx = matrix_double(num_rows, 2 * num_cols)) < 0) {
     NEW:  if ((augmented_mtx = matrix_double(num_rows, 2 * num_cols)) == NULL) {
   b) matrix_double_check_diagonal_non_zero_inverse(), line 289
     OLD:  if ((mtx_checked = matrix_double(i_checked_size, i_checked_size)) < 0) {
     NEW:  if ((mtx_checked = matrix_double(i_checked_size, i_checked_size)) == NULL) {


CRASH / GARBAGE ON A BAD INPUT FILE

3. velmod.c, read_fdiff_vel(), lines 2158-2190
   - If the velocity grid file cannot be opened, the error is printed but the
     function continues and reads from a NULL file pointer (crash).
     NEW: add `return;` after the error message.
   - Header reads (fdgrid_numx, fdgrid_numz) and each row read are now
     checked; on a short read: print error, fclose, return.
   - If malloc of fdgrid_array fails, the error is printed but the function
     continues and writes into NULL. NEW: fclose and return.

4. Vel2Grid3D.c, ReadVelModel(), lines 341-485
   - Binary model: each fread of a velocity value is checked; on a short read:
     "ERROR: unexpected end of velocity model file", fclose, return -2.
   - ASCII model: each fscanf of a coordinate / velocity value is checked;
     on failure the loop stops (break) instead of continuing with an
     undefined value. This makes the existing nread short-read check work.
   - fscanf("%*f %*f %*f") skip line: warning if end of file.

5. interface2fmm.c, ReadVelModel(), lines 258-302
   Same fscanf checks as in Vel2Grid3D.c (item 4, ASCII part).

6. NLLocLib.c, GetNextObs(), lines 2768 and 2780
   fgets at end of observation file is now checked:
     OLD:  fgets(line, MAXLINE_LONG, fp_obs);
     NEW:  if (fgets(line, MAXLINE_LONG, fp_obs) == NULL)
               return (OBS_FILE_END_OF_INPUT);
   (both the plain case and the SKIP_NEXT_LINE case)

7. GridLib.c, ReadGrid3dHdr() line 2226 and OpenGrid3dFile() line 2393
   fscanf of the source line (label x y z) in time/angle grid headers is now
   checked; if it does not read 4 values:
     "WARNING: Cannot read source location from grid header file: <name>"

8. Grid2GMT.c, Scat2GMT(), lines 2549 and 2554
   fread of the scatter file header (tot_npoints) and of each point is
   checked; on a short read: error message, close both files, return -1.

9. LocSum.c, SumLocations(), line 438
   fread of scatter file header checked; on failure: error message, write
   "SCATTER Nsamples 0 / END_SCATTER", close file, skip to next event.

10. scat2latlon.c, ConvertScat(), line 253
    fread of scatter file header checked; on failure: error message, close
    file, skip to next file.

11. fmm2grid.c, ReadFmmTimes(), line 300
    fgets of the (ignored) source line checked; on failure the line is set
    empty, so the existing sscanf check then reports the error.


SHELL COMMANDS - ONLY A WARNING ADDED

12. system() return value checked, warning printed if non-zero
    (no change in behaviour otherwise):
    - Grid2GMT.c, GenGMTCommands() line 969: "chmod a+x" and running the GMT script
    - GridLib.c, ExpandWildCards_OLD() line 5094: "ls ... > list_file"
    - Loc2ssst.c, DoLoc2ssst() line 1043: GMT makecpt
    - NLDiffLoc.c, main() lines 326 and 328: the two "cp -p" of the control file
    Pattern:
      OLD:  system(sys_command);
      NEW:  if (system(sys_command) != 0)
                fprintf(stderr, "WARNING: command returned non-zero status: %s\n", sys_command);


NOTE ON TEST COVERAGE
---------------------
The test data exercise NLLoc, Vel2Grid, Grid2Time, Time2EQ and Loc2ssst.
The changes in NLDiffLoc.c, Grid2GMT.c, LocSum.c, scat2latlon.c, fmm2grid.c,
interface2fmm.c and Vel2Grid3D.c are in programs the tests do not run, and
some changed functions in the tested programs (e.g. read_fdiff_vel, the two
matrix allocation checks) are not reached by the test data. For those the
check is the build plus reading the diff. Each of those changes is one to six
lines and affects only the error path.
