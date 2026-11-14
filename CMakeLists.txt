# CMakeLists.txt
cmake_minimum_required(VERSION 3.16)
project(stress_resolved CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(NOT LAMMPS_SOURCE_DIR)
  message(FATAL_ERROR "Set -DLAMMPS_SOURCE_DIR=/path/to/lammps/src")
endif()

# integer sizes, must be the LAMMPS_SIZES used to build LAMMPS or tagint/bigint/imageint
# (atom IDs, image flags) have different sizes in the plugin and in LAMMPS
set(LAMMPS_SIZES "smallbig" CACHE STRING "LAMMPS integer sizes: smallbig (default), bigbig or smallsmall")
set_property(CACHE LAMMPS_SIZES PROPERTY STRINGS smallbig bigbig smallsmall)
if(NOT LAMMPS_SIZES MATCHES "^(smallbig|bigbig|smallsmall)$")
  message(FATAL_ERROR "LAMMPS_SIZES must be smallbig, bigbig or smallsmall, not ${LAMMPS_SIZES}")
endif()
string(TOUPPER ${LAMMPS_SIZES} LAMMPS_SIZES_DEFINE)

add_library(stress_resolved MODULE stress_resolved_plugin.cpp compute_stress_resolved.cpp)
target_include_directories(stress_resolved PRIVATE ${LAMMPS_SOURCE_DIR} ${LAMMPS_SOURCE_DIR}/KSPACE)
target_compile_definitions(stress_resolved PRIVATE LAMMPS_${LAMMPS_SIZES_DEFINE})
# nodelete: keep the code loaded after LAMMPS unloads the plugin at exit, because
# LAMMPS deletes the active ewald2 instance only after unloading kspace plugins
set_target_properties(stress_resolved PROPERTIES PREFIX "" SUFFIX ".so" LINK_FLAGS "-rdynamic -Wl,-z,nodelete")
