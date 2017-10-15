cmake_minimum_required(VERSION 2.8)
project(FMM)

set(FMM_MAJOR_VERSION 0)
set(FMM_MINOR_VERSION 1)
set(FMM_PATCH_VERSION 0)
set(FMM_VERSION
  ${FMM_MAJOR_VERSION}.${FMM_MINOR_VERSION}.${FMM_PATCH_VERSION})

# Select flags.
SET(CMAKE_CXX_FLAGS "-std=c++11 -march=native -fPIC -mavx")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fno-finite-math-only")
SET(CMAKE_CXX_FLAGS_DEBUG  "-Wall -Wno-unused-local-typedefs -g")

# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH
  "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/CMake/FMM)
endif()
set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()


# set up include-directories
include_directories(
  ${PROJECT_SOURCE_DIR}/console
  ${PROJECT_BINARY_DIR}/ndgridmap
  ${PROJECT_BINARY_DIR}/fmm
  ${PROJECT_BINARY_DIR}/fmm/fmdata 
  ${PROJECT_BINARY_DIR}/fm2
  ${PROJECT_BINARY_DIR}/io
  ${PROJECT_BINARY_DIR}/gradientdescent)


add_library(fmm SHARED console/console.cpp
                    ndgridmap/cell.cpp
                    fmm/fmdata/fmcell.cpp)
                    
target_link_libraries (fmm X11 pthread ) 

SET(FMM_INCLUDE_DIRECTORIES
  console
  ndgridmap
  fmm 
  fm2
  io
  gradientdescent)

INSTALL( DIRECTORY ${FMM_INCLUDE_DIRECTORIES} DESTINATION ${INSTALL_INCLUDE_DIR} )
  
install(TARGETS fmm
  # IMPORTANT: Add the foo library to the "export-set"
  EXPORT FMMTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}" COMPONENT lib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/fmm"
    COMPONENT dev)
    
add_executable (fmm_demo main.cpp) 
target_link_libraries (fmm_demo fmm X11 pthread) 
  
install(TARGETS fmm_demo
  # IMPORTANT: Add the fmm_demo executable to the "export-set"
  EXPORT FMMTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
  
# The interesting stuff goes here
# ===============================

# Add all targets to the build-tree export set
export(TARGETS fmm_demo fmm
  FILE "${PROJECT_BINARY_DIR}/FMMTargets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE FMM)

# Create the FMMConfig.cmake and FMMConfigVersion files
file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(FMMConfig.cmake.in
  "${PROJECT_BINARY_DIR}/FMMConfig.cmake" @ONLY)
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${FMM_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(FMMConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FMMConfig.cmake" @ONLY)
# ... for both
configure_file(FMMConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/FMMConfigVersion.cmake" @ONLY)

# Install the FMMConfig.cmake and FMMConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FMMConfig.cmake"
  "${PROJECT_BINARY_DIR}/FMMConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

# Install the export set for use with the install-tree
install(EXPORT FMMTargets DESTINATION
  "${INSTALL_CMAKE_DIR}" COMPONENT dev)
