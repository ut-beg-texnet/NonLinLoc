cmake_minimum_required (VERSION 3.15)

project (matrix_statistics)

#set(CMAKE_INSTALL_LIBDIR "lib")
#set(CMAKE_INSTALL_INCLUDEDIR "include")

SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ../lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)

INCLUDE_DIRECTORIES(../alomax_matrix)

set (SOURCES
	${PROJECT_SOURCE_DIR}/matrix_statistics_DEV.c
	${PROJECT_SOURCE_DIR}/matrix_statistics.c
	${PROJECT_SOURCE_DIR}/test_ellipsoid.c
)

add_library(matrix_statistics_DEV SHARED ${PROJECT_SOURCE_DIR}/matrix_statistics_DEV.c)
add_library(matrix_statistics SHARED ${PROJECT_SOURCE_DIR}/matrix_statistics.c)
add_library(test_ellipsoid SHARED ${PROJECT_SOURCE_DIR}/test_ellipsoid.c)

#add_library(matrix_statistics OBJECT ${SOURCES})

#target_link_libraries(
#    test_real_symmetric_eigen
 #   ${CMAKE_CURRENT_SOURCE_DIR}/..vector/vector.o
#  ${CMAKE_CURRENT_SOURCE_DIR}../alomax_matrix/alomax_matrix.o
 #   ${CMAKE_CURRENT_SOURCE_DIR} ../alomax_matrix/alomax_matrix_svd.o
#)

#ADD_EXECUTABLE(test_real_symmetric_eigen ${SOURCES})
