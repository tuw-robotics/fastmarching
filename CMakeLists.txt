cmake_minimum_required (VERSION 2.6) 
project (FMM) 

# set a default build type for single-configuration
# CMake generators if no build type is set.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

set(BUILD_EXAMPLES true CACHE STRING "True (default) to build examples")
if(BUILD_EXAMPLES)
  message(STATUS "Examples are being built (NOT installed).")
endif(BUILD_EXAMPLES)

# Select flags.
set(CMAKE_CXX_FLAGS "-std=c++11")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fno-finite-math-only")
set(CMAKE_CXX_FLAGS_DEBUG  "-Wall -Wno-unused-local-typedefs -g")

# Finding Boost
find_package(Boost REQUIRED COMPONENTS system filesystem program_options)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "Boost NOT FOUND. Please install it following the instructions on the README.")
endif()

include(GNUInstallDirs)

link_directories(/usr/local/)
# Self-made includes
include_directories (
  console
  ndgridmap
  fmm
  fmm/fmdata
  fm2
  io
  gradientdescent
  thirdparty
)

# Create main library
add_library (fastmarching SHARED
  console/console.cpp
  ndgridmap/cell.cpp
  fmm/fmdata/fmcell.cpp
)

# Linking 
target_link_libraries (fastmarching
  ${Boost_LIBRARIES}
  #CImg dependencies.
  X11 
  pthread
)

# Add benchmarking capabilities
add_subdirectory(benchmark)

# Create main example
if(BUILD_EXAMPLES)
    add_executable (fmm
      main.cpp
    )
    target_link_libraries (fmm fastmarching)

    add_subdirectory(examples)
endif(BUILD_EXAMPLES)

# Install section
install (TARGETS fastmarching
         DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY ${CMAKE_SOURCE_DIR} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING PATTERN "*.hpp"
    REGEX "/build$" EXCLUDE
    REGEX ".git$" EXCLUDE
    REGEX "/examples$" EXCLUDE
    REGEX "/doc$" EXCLUDE
    REGEX "/data$" EXCLUDE
)


# Uninstall target
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
