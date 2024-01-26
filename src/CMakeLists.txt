# Listing the version is nice here since it sets lots of useful variables
CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

# Good tutorial:
# https://www.dealii.org/9.1.1/users/cmakelists.html

PROJECT(NonLinLoc)
set(CMAKE_BUILD_TYPE Release)
###DEBUG#set(CMAKE_BUILD_TYPE Debug)

# If you want to build the library in ${CMAKE_CURRENT_SOURCE_DIR}/lib
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# If you want to move the binaries in ${CMAKE_CURRENT_SOURCE_DIR}/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
###DEBUG#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin_debug)

# Bring the include directory with its headers into the project
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${CMAKE_CURRENT_SOURCE_DIR}/alomax_matrix
	${CMAKE_CURRENT_SOURCE_DIR}/geometry
	${CMAKE_CURRENT_SOURCE_DIR}/matrix_statistics
	${CMAKE_CURRENT_SOURCE_DIR}/ran1
	${CMAKE_CURRENT_SOURCE_DIR}/octtree
	${CMAKE_CURRENT_SOURCE_DIR}/io
	${CMAKE_CURRENT_SOURCE_DIR}/io/jReadWrite
	${CMAKE_CURRENT_SOURCE_DIR}
)

# The -fcommon places uninitialized global variables in a common block. This allows the linker to resolve all tentative definitions of the same variable in different compilation units to the same object, or to a non-tentative definition. This behavior is inconsistent with C++, and on many targets implies a speed and code size penalty on global variable references. It is mainly useful to enable legacy code to link without errors.
add_compile_options(-fcommon)

# see all warnings
add_compile_options(-Wall)

# Uncomment the following if the compiler you are using supports GNU extensions such as:
#   fmemopen, open_memstream -  open memory as stream (required to pass observations file lines as parameters
#   to a function call to NLLoc() and to run NLLoc_func_test)
#add_compile_options(-D _GNU_SOURCE)
#
# Alternatively, in the make environment, set:
# export CFLAGS='-D _GNU_SOURCE'
# before running: rm CMakeCache.txt; cmake .

# 20230318 Anthony Lomax - Uncomment the following to force reading of *DEFAULT* named travel-time grids with byte swapping.
#    Needed as a cluge when running Loc2ssst and NLLoc when doing NLL-SSST (e.g. with run_ssst_relocations.bash) where time grids with and
#    without byte swapping are used together.
#add_compile_definitions(LOC2SSST_CLUGE)

# 20220921 Sean Ho - Uncomment the following if you want Grid2GMT to generate scripts compatible with gmt 5.x or 6.x
add_compile_options(-D GMT_VER_5)
#
# Alternatively, in the make environment, set:
# export CFLAGS='-D GMT_VER_5'
# before running: rm CMakeCache.txt; cmake .

# 20211110 AJL - New fix implemented, see https://github.com/alomax/NonLinLoc/issues/7
#   This fix probably supersedes the following:
# 20210511 AJL - Possible bug fix for Loc2ddct compile errors on older(?) compilers "relocation truncated to fit: R_X86_64_PC32 against symbol"
# https://stackoverflow.com/questions/57331990/c-compiling-relocation-truncated-to-fit-r-x86-64-pc32-against-symbol
# These compile/link options not needed if following define used in Loc2ddct.c
# #define  MAX_NUM_DIFF_HYPOCENTERS 50000
#
# Adds options to the COMPILE_OPTIONS directory property. These options are used when compiling targets from the current directory and below.
#add_compile_options(-mcmodel=medium)
#
# Add options to the link step for executable, shared library or module library targets in the current directory and below that are added after this command is invoked.
#add_link_options("-Wl,-no_pie")


## Create the .o object files with add_library()
### Simplify by just creating the GRID_LIB_OBJS .o object file
add_library(GRID_LIB_OBJS OBJECT GridLib.c util.c geo.c octtree/octtree.c io/json_io.c io/jReadWrite/source/jRead.c io/jReadWrite/source/jWrite.c alomax_matrix/alomax_matrix.c alomax_matrix/eigv.c alomax_matrix/alomax_matrix_svd.c matrix_statistics/matrix_statistics.c vector/vector.c ran1/ran1.c map_project.c)

### Simplify by just creating the NLLOC_LIB_OBJS .o object file
add_library(NLLOC_LIB_OBJS OBJECT calc_crust_corr.c velmod.c NLLocLib.c GridMemLib.c phaselist.c loclist.c otime_limit.c)

#
add_library(LOC_PHS_LIST OBJECT phaselist.c loclist.c)


# Create NLLoc1.o
add_library(NLLoc1 OBJECT NLLoc1.c)

# Combine object .o files to OBJS0
set(OBJS0 NLLoc1 GRID_LIB_OBJS NLLOC_LIB_OBJS)

# --------------------------------------------------------------------------
# NLLoc
# Important: set target_link_libraries after add_executable()
# Also: add the math lib with: m
add_executable(NLLoc NLLoc_main.c)
target_link_libraries(NLLoc ${OBJS0} m)

# --------------------------------------------------------------------------
# Vel2Grid
#
add_library(velmod OBJECT velmod.c)
add_executable(Vel2Grid Vel2Grid1.c)
target_link_libraries(Vel2Grid velmod GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# Grid2Time
#
add_library(Time_3d_NLL OBJECT Time_3d_NLL.c)
target_compile_options(Time_3d_NLL PRIVATE "-DNO_IEEE_PROTOCOL")
add_executable(Grid2Time Grid2Time1.c)
target_link_libraries(Grid2Time Time_3d_NLL GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# Time2Angles
#
add_executable(Time2Angles Time2Angles1.c)
target_link_libraries(Time2Angles GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# Grid2GMT
#
add_executable(Grid2GMT Grid2GMT.c)
add_library(GridGraphLib OBJECT GridGraphLib.c)
target_link_libraries(Grid2GMT GridGraphLib GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# LocSum
#
add_executable(LocSum LocSum.c)
target_link_libraries(LocSum GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# Time2EQ
#
add_executable(Time2EQ Time2EQ1.c)
target_link_libraries(Time2EQ GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# NLLoc_func_test
#
add_executable(NLLoc_func_test NLLoc_func_test.c)
target_link_libraries(NLLoc_func_test ${OBJS0} m)

# --------------------------------------------------------------------------
# hypoe2hyp
#
add_executable(hypoe2hyp hypoe2hyp.c)
target_link_libraries(hypoe2hyp GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# fpfit2hyp
#
add_executable(fpfit2hyp fpfit2hyp.c)
target_link_libraries(fpfit2hyp GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# PhsAssoc
#
add_library(calc_crust_corr OBJECT calc_crust_corr.c)
add_executable(PhsAssoc PhsAssoc.c)
target_link_libraries(PhsAssoc calc_crust_corr GRID_LIB_OBJS m)


# --------------------------------------------------------------------------
# oct2grid
#
add_executable(oct2grid oct2grid.c)
target_link_libraries(oct2grid GRID_LIB_OBJS m)


# --------------------------------------------------------------------------
# ttime_func_test
#
add_executable(ttime_func_test ttime_func_test.c)
target_link_libraries(ttime_func_test GRID_LIB_OBJS m)


# --------------------------------------------------------------------------
# mag_func_test
#
add_executable(mag_func_test mag_func_test.c)
target_link_libraries(mag_func_test GRID_LIB_OBJS NLLOC_LIB_OBJS m)


# --------------------------------------------------------------------------
# Vel2Grid3D
#
add_executable(Vel2Grid3D Vel2Grid3D.c)
target_link_libraries(Vel2Grid3D velmod GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# interface2fmm
#
add_executable(interface2fmm interface2fmm.c)
target_link_libraries(interface2fmm GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# fmm2grid
#
add_executable(fmm2grid fmm2grid.c)
target_link_libraries(fmm2grid GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# scat2latlon
#
add_executable(scat2latlon scat2latlon.c)
target_link_libraries(scat2latlon GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# GridCascadingDecimate
#
add_executable(GridCascadingDecimate GridCascadingDecimate.c)
target_link_libraries(GridCascadingDecimate GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# sphfd_SWR_NLL
#
add_executable(sphfd_SWR_NLL sphfd_SWR_NLL.c)
target_link_libraries(sphfd_SWR_NLL GRID_LIB_OBJS m)
target_compile_options(sphfd_SWR_NLL PRIVATE "-w")

# --------------------------------------------------------------------------
# grid2scat
#
add_executable(grid2scat grid2scat.c)
target_link_libraries(grid2scat GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# Loc2ssst
#
add_executable(Loc2ssst Loc2ssst.c)
target_link_libraries(Loc2ssst GRID_LIB_OBJS LOC_PHS_LIST m)

# --------------------------------------------------------------------------
# NLDiffLoc
#
add_executable(NLDiffLoc NLDiffLoc.c)
target_link_libraries(NLDiffLoc GRID_LIB_OBJS NLLOC_LIB_OBJS m)

# --------------------------------------------------------------------------
# Loc2ddct
#
add_executable(Loc2ddct Loc2ddct.c)
target_link_libraries(Loc2ddct GRID_LIB_OBJS m)

# --------------------------------------------------------------------------
# loc_combine
#
add_executable(loc_combine loc_combine.c)
target_link_libraries(loc_combine GRID_LIB_OBJS LOC_PHS_LIST m)
