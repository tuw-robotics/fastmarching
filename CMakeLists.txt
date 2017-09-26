cmake_minimum_required (VERSION 2.6) 
project (FMM) 

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")
include(FindGperftools.cmake)


# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

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

find_package(Gperftools REQUIRED)

# Self-made includes
include_directories (${FMM_SOURCE_DIR}/console)
include_directories (${FMM_SOURCE_DIR}/ndgridmap)
include_directories (${FMM_SOURCE_DIR}/fmm)
include_directories (${FMM_SOURCE_DIR}/fmm/fmdata) 
include_directories (${FMM_SOURCE_DIR}/fm2)
include_directories (${FMM_SOURCE_DIR}/io)
include_directories (${FMM_SOURCE_DIR}/gradientdescent)


# Third party internal includes
include_directories (${FMM_SOURCE_DIR}/thirdparty) 


# External dependencies
include_directories (~/boost_1_55_0)

add_library(fmmlib console/console.cpp
                    ndgridmap/cell.cpp
                    fmm/fmdata/fmcell.cpp)

add_executable (fmm main.cpp
                    console/console.cpp
                    ndgridmap/cell.cpp
                    fmm/fmdata/fmcell.cpp
                    ) 

# Linking CImg dependencies.
target_link_libraries (fmm X11 pthread -lprofiler) 
target_link_libraries (fmmlib X11 pthread ) 

# install(TARGETS
#   fmm
#   DESTINATION ${FMM_SOURCE_DIR}/build/bin
#   LIBRARY DESTINATION ${FMM_SOURCE_DIR}/build/lib
# )
# 
# install(DIRECTORY fmm/ ndgridmap/
#   DESTINATION ${FMM_SOURCE_DIR}/build/include
#   PATTERN ".svn" EXCLUDE
# )

# Add all targets to the build-tree export set
export(TARGETS fmm fmmlib
  FILE "${PROJECT_BINARY_DIR}/FMMTargets.cmake")
 
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE FMM)
 
# Create the FooBarConfig.cmake and FooBarConfigVersion files
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
 
# Install the FooBarConfig.cmake and FooBarConfigVersion.cmake
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/FMMConfig.cmake"
  "${PROJECT_BINARY_DIR}/FMMConfigVersion.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
 
# Install the export set for use with the install-tree
install(EXPORT FMMTargets DESTINATION
  "${INSTALL_CMAKE_DIR}" COMPONENT dev)
