# CMakeLists.txt
cmake_minimum_required(VERSION 3.16)
project(stress_atom_resolved CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(NOT LAMMPS_SOURCE_DIR)
  message(FATAL_ERROR "Set -DLAMMPS_SOURCE_DIR=/path/to/lammps/src")
endif()

add_library(stress_atom_resolved MODULE stress_atom_resolved_plugin.cpp compute_stress_atom_resolved.cpp)
target_include_directories(stress_atom_resolved PRIVATE ${LAMMPS_SOURCE_DIR} ${LAMMPS_SOURCE_DIR}/KSPACE)
# nodelete: keep the code loaded after LAMMPS unloads the plugin at exit, because
# LAMMPS deletes the active ewald2 instance only after unloading kspace plugins
set_target_properties(stress_atom_resolved PROPERTIES PREFIX "" SUFFIX ".so" LINK_FLAGS "-rdynamic -Wl,-z,nodelete")
