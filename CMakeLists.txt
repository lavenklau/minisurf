cmake_minimum_required(VERSION 3.22)

# Set your project name
project(MinsurfProject LANGUAGES C CXX)

# disable all warning for clean output
add_definitions(-w)

set(USE_CONDA OFF CACHE BOOL "in conda environment?")
if(USE_CONDA)
set(PACKAGE_PREFIX $ENV{CONDA_PREFIX})
else()
message(STATUS OUTSIDE CONDA ENVIRONMENT)
endif()
set(CMAKE_PREFIX_PATH  ${PACKAGE_PREFIX})
list(APPEND CMAKE_PREFIX_PATH ${PACKAGE_PREFIX}/mkl-devel/lib/cmake/mkl)
message(STATUS prefix_path=${CMAKE_PREFIX_PATH}) 
list(APPEND CMAKE_MODULE_PATH ${PACKAGE_PREFIX}/lib/cmake/OpenVDB)

# config MKL
set(MKL_INTERFACE lp64)
option(ENABLE_BLAS95      OFF ) #Enables BLAS Fortran95 API
option(ENABLE_LAPACK95    OFF ) #Enables LAPACK Fortran95 API
option(ENABLE_BLACS       OFF ) #Enables cluster BLAS library
option(ENABLE_CDFT        OFF ) #Enables cluster DFT library
option(ENABLE_CPARDISO    OFF ) #Enables cluster PARDISO functionality
option(ENABLE_SCALAPACK   OFF ) #Enables cluster LAPACK library
option(ENABLE_OMP_OFFLOAD OFF ) #Enables OpenMP Offload functionality
find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
message(STATUS "found openmp")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
message(STATUS "openmp option c: ${OpenMP_C_FLAGS} | cpp : ${OpenMP_CXX_FLAGS}")
endif()
find_package(MKL CONFIG REQUIRED)

#find_package (SuiteSparse REQUIRED)


# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


set(SOURCE_DIR  ${CMAKE_SOURCE_DIR})

# Include directories
include_directories(${SOURCE_DIR})
include_directories(${SOURCE_DIR}/geometry)

# Add source files
file(GLOB_RECURSE  SOURCES
    "${SOURCE_DIR}/cgal/*.cpp"
    "${SOURCE_DIR}/geometry/*.cpp"
    "${SOURCE_DIR}/figure/*.cpp"
    "${SOURCE_DIR}/flow/*.cpp"
    "${SOURCE_DIR}/matlab/*.cpp"
    "${SOURCE_DIR}/material/*.cpp"
    "${SOURCE_DIR}/mesher/*.cpp"
    "${SOURCE_DIR}/vdb/*.cpp"
    "${SOURCE_DIR}/triangle++/*.cpp"
)
file(GLOB SOURCES1 "${SOURCE_DIR}/*.cpp")
list(APPEND SOURCES ${SOURCES1})

message(STATUS "src = ${SOURCES}")

find_package(Eigen3 3.4 REQUIRED)

# fetch package
include(FetchContent)

FetchContent_Declare(
  libigl
  SOURCE_DIR "${SOURCE_DIR}/extern/libigl-2.4.0"
)

FetchContent_Declare(
  OpenMesh
  SOURCE_DIR "${SOURCE_DIR}/extern/openmesh9"
)

FetchContent_Declare(
  Spectra
  SOURCE_DIR "${SOURCE_DIR}/extern/spectra-1.0.1"
)

FetchContent_MakeAvailable(libigl)
FetchContent_MakeAvailable(OpenMesh)
FetchContent_MakeAvailable(Spectra)

message(STATUS "spectra header = ${Spectra_SOURCE_DIR}/include")

add_definitions(-D_CRT_SECURE_NO_WARNINGS)

find_package(fmt REQUIRED)
find_package(CGAL REQUIRED)
find_package(CLI11 CONFIG REQUIRED)
find_package(TBB COMPONENTS tbbmalloc tbbmalloc_proxy tbb_preview)
find_package(OpenMesh REQUIRED)
find_package(OpenVDB  REQUIRED)


# Add executable
add_executable(minsurf ${SOURCES})

if(TARGET tbb)
target_link_libraries(minsurf PUBLIC tbb)
endif()
if(TARGET TBB::tbb)
target_link_libraries(minsurf PUBLIC TBB::tbb)
endif()

target_link_libraries(minsurf PUBLIC Eigen3::Eigen)
target_link_libraries(minsurf PUBLIC fmt::fmt)
target_link_libraries(minsurf PUBLIC  CGAL::CGAL)
target_link_libraries(minsurf PUBLIC CLI11::CLI11)
target_link_libraries(minsurf PUBLIC igl::core) 
target_link_libraries(minsurf PUBLIC OpenMeshCore)
target_link_libraries(minsurf PUBLIC TBB::tbbmalloc)
target_link_libraries(minsurf PUBLIC igl::core) 
target_link_libraries(minsurf PUBLIC OpenVDB::openvdb)
target_include_directories(minsurf PUBLIC MKL_INCLUDE)
target_link_libraries(minsurf PUBLIC $<LINK_ONLY:MKL::MKL>)
target_include_directories(minsurf PUBLIC ${Spectra_SOURCE_DIR}/include)
#target_link_libraries (minsurf PUBLIC SuiteSparse::CHOLMOD)

