cmake_minimum_required(VERSION 3.18)
set(CMAKE_CXX_STANDARD 17)

# Explicitly list the test source code and headers. The Catch header-only unit
# test framework is stored in with the test source.
set(Proj_TEST_SRC
    dummy.cpp
)
set(Proj_TEST_HEADER
    dummy.hpp
)

PREPEND(Proj_TEST_SRC)

# Make an executable target that depends on the test source code we specified
# above.
add_executable(TestSurfaces ${Proj_TEST_SRC} ${Proj_TEST_HEADER})
set_property(TARGET TestSurfaces
             PROPERTY CUDA_SEPARABLE_COMPILATION ON)

# Enable testing via CTest
enable_testing()
# Add our test as runnable via CTest
add_test(NAME TestSurfaces  COMMAND TestSurfaces)

# Link our unit tests against the library we compiled
target_link_libraries(TestSurfaces MolecularSurfaces)

if(APPLE)
  # We need to add the path to the driver (libcuda.dylib) as an rpath,
  # so that the static cuda runtime can find it at runtime.
  set_property(TARGET TestSurfaces
               PROPERTY
               BUILD_RPATH ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
endif()