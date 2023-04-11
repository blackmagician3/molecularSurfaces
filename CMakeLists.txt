cmake_minimum_required(VERSION 3.18 FATAL_ERROR)
project(MolecularSurfaces LANGUAGES CXX CUDA)

# enable unit testing via "make test" once the code has been compiled.
include(CTest)

# Function to prepend the subdirectory to source files in subdirectories
FUNCTION(PREPEND var )
   SET(listVar "")
   FOREACH(f ${${var}})
       LIST(APPEND listVar "${CMAKE_CURRENT_SOURCE_DIR}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)

# After a normal build, we can specify the location of various outputs of the
# build. We put executables and static libraries outside the build directory in
# bin/ and lib/, respectively.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build/lib")

# Include source code and headers. This runs the CMakeLists.txt in each
# subdirectory. These can define their own libraries, executables, etc. as targets,
# but here we define all exportable targets in the root CMakeLists.txt.
add_subdirectory(app)
add_subdirectory(src)
add_subdirectory(include)

add_library(MolecularSurfaces STATIC
            app/main.cu
            include/example.h
            src/example.cpp
            )

target_include_directories(MolecularSurfaces PUBLIC include/)

# Request that particles be built with -std=c++11
# As this is a public compile feature anything that links to
# particles will also build with -std=c++17
target_compile_features(MolecularSurfaces PUBLIC cxx_std_17)

# Set CUDA architecture
# default value: 75 to match testing system with nvidia 2080 card
if(NOT DEFINED ${CMAKE_CUDA_ARCHITECTURES})
  set(CMAKE_CUDA_ARCHITECTURES 75 CACHE STRING "CUDA architectures" FORCE)
endif()
message(STATUS "CUDA architectures set to ${CMAKE_CUDA_ARCHITECTURES}")


# We need to explicitly state that we need all CUDA files in the
# particle library to be built with -dc as the member functions
# could be called by other libraries and executables
set_target_properties( MolecularSurfaces
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON
                       CUDA_RESOLVE_DEVICE_SYMBOLS ON
                       )

add_executable(Surfaces app/main.cu)

set_property(TARGET Surfaces
             PROPERTY CUDA_SEPARABLE_COMPILATION ON
             CUDA_RESOLVE_DEVICE_SYMBOLS ON
             )

target_link_libraries(Surfaces PRIVATE MolecularSurfaces)