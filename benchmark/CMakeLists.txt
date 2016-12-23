add_executable (fmmbenchmark fmm_benchmark.cpp)  
target_link_libraries (fmmbenchmark fastmarching)
install (TARGETS fmmbenchmark
         DESTINATION ${CMAKE_INSTALL_BINDIR}
)
