# Listing the version is nice here since it sets lots of useful variables
CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

# Good tutorial:
# https://www.dealii.org/9.1.1/users/cmakelists.html

PROJECT(NonLinLoc)
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

# If you want to build the library in ${CMAKE_CURRENT_SOURCE_DIR}/lib
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# If you want to move the binaries in ${CMAKE_CURRENT_SOURCE_DIR}/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# Bring the include directory with its headers into the project
include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${CMAKE_CURRENT_SOURCE_DIR}/alomax_matrix
	${CMAKE_CURRENT_SOURCE_DIR}/geometry
	${CMAKE_CURRENT_SOURCE_DIR}/matrix_statistics
	${CMAKE_CURRENT_SOURCE_DIR}/ran1
	${CMAKE_CURRENT_SOURCE_DIR}/octtree
	${CMAKE_CURRENT_SOURCE_DIR}
)

## Create the .o object files with add_library()
### Simplify by just creating the GRID_LIB_OBJS .o object file
add_library(GRID_LIB_OBJS OBJECT GridLib.c util.c geo.c octtree/octtree.c alomax_matrix/alomax_matrix.c alomax_matrix/eigv.c alomax_matrix/alomax_matrix_svd.c matrix_statistics/matrix_statistics.c vector/vector.c ran1/ran1.c map_project.c)

### Simplify by just creating the NLLOC_LIB_OBJS .o object file
add_library(NLLOC_LIB_OBJS OBJECT calc_crust_corr.c velmod.c NLLocLib.c GridMemLib.c phaselist.c loclist.c otime_limit.c)

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
add_library(phaselist OBJECT phaselist.c)
add_library(loclist OBJECT loclist.c)
add_executable(Loc2ssst Loc2ssst.c)
target_link_libraries(Loc2ssst GRID_LIB_OBJS phaselist loclist m)

