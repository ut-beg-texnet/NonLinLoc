cmake_minimum_required (VERSION 3.1)

project (alomax_matrix)

#set(CMAKE_INSTALL_LIBDIR "lib")
#set(CMAKE_INSTALL_INCLUDEDIR "include")

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)

SET(SOURCES
	${PROJECT_SOURCE_DIR}/alomax_matrix.c
	${PROJECT_SOURCE_DIR}/alomax_matrix_svd.c
	${PROJECT_SOURCE_DIR}/eigv.c
	${PROJECT_SOURCE_DIR}/test_real_symmetric_eigen.c
	${PROJECT_SOURCE_DIR}/polarization.c
)

ADD_EXECUTABLE(test_real_symmetric_eigen ${SOURCES})
