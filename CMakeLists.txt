cmake_minimum_required (VERSION 3.2)

project (demo)

#option(MPI_USE "option for MPI" ON)
#if(MPI_USE)
#    add_compile_definitions(MPI_USE)
#endif()

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

aux_source_directory(src SRC_LIST)

include_directories(include)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

add_executable(main ${SRC_LIST})

set(LAPACKE_PATH /root/software/lapack-3.11.0)

target_include_directories(main PRIVATE include ${LAPACKE_PATH}/LAPACKE/include ${LAPACKE_PATH}/CBLAS/include ${MPI_INCLUDE_PATH})
target_link_directories(main PRIVATE ${LAPACKE_PATH}/build)
target_link_libraries(main OpenMP::OpenMP_CXX liblapacke.a liblapack.a libcblas.a librefblas.a gfortran ${MPI_LIBRARIES} m libscalapack-openmpi.so)
