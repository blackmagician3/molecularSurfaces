####################################################################
### project setup ##################################################
####################################################################
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

# Define project and executable name
set(PROJECT_NAME MolecularSurfaces)
set(EXECUTABLE_NAME Surfaces)
project(${PROJECT_NAME} LANGUAGES CXX CUDA)


# Make sure binary directory is not the same as source directory
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there."
  )
endif()


# enable unit testing via "make test" once the code has been compiled.
include(CTest)

# OpenGL
find_package(OpenGL REQUIRED)


# # Define common paths
set(DEP_DIR ${CMAKE_CURRENT_BINARY_DIR}/dep)
set(GLAD_INSTALL_DIR ${DEP_DIR}/glad)
set(GLFW_INSTALL_DIR ${DEP_DIR}/glfw)
set(GLM_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/ext/glm)

####################################################################
### external projects ##############################################
####################################################################
# Add glfw, glad as external project
include(ExternalProject)

SET(CMAKE_EXE_LINKER_FLAGS /NODEFAULTLIB:\"libcmt.lib\")

set(libGLFW glfw)
ExternalProject_Add(${libGLFW}
PREFIX          ${DEP_DIR}/${libGLFW}
GIT_REPOSITORY  https://github.com/glfw/glfw.git
GIT_TAG         3.3.4
GIT_SHALLOW     ON
UPDATE_DISCONNECTED 1 # Do not update the external project if it is already downloaded
CMAKE_ARGS      -DCMAKE_INSTALL_PREFIX:PATH=${GLFW_INSTALL_DIR}
                -DGLFW_BUILD_DOCS:BOOL=OFF
                -DGLFW_BUILD_EXAMPLES:BOOL=OFF
                -DGLFW_BUILD_TESTS:BOOL=OFF
)

set(libGLAD glad)
ExternalProject_Add(${libGLAD}
PREFIX          ${DEP_DIR}/${libGLAD}
GIT_REPOSITORY  https://github.com/Dav1dde/glad.git
GIT_TAG         origin/master
GIT_SHALLOW     ON
UPDATE_DISCONNECTED 1 # Do not update the external project if it is already downloaded
CMAKE_ARGS      -DCMAKE_INSTALL_PREFIX:PATH=${GLAD_INSTALL_DIR}
                -DGLAD_INSTALL:BOOL=ON
                -DGLAD_PROFILE:STRING="core"
                -DGLAD_ALL_EXTENSIONS:BOOL=ON
                -DUSE_MSVC_RUNTIME_LIBRARY_DLL:BOOL=OFF
)

# set(libGLM glm)
# ExternalProject_Add(${libGLM}
#     PREFIX ${DEP_DIR}/${libGLM}
#     GIT_REPOSITORY https://github.com/g-truc/glm.git
#     LOG_DOWNLOAD ON
#     BUILD_COMMAND ""
#     INSTALL_COMMAND   ""
#     UPDATE_DISCONNECTED 1 # Do not update the external project if it is already downloaded
#     CMAKE_ARGS  -DCMAKE_INSTALL_PREFIX:PATH=${GLM_INSTALL_DIR}
#                 -DGLAD_INSTALL:BOOL=ON
#                 -DGLAD_PROFILE:STRING="core"
#                 -DGLAD_ALL_EXTENSIONS:BOOL=ON
#                 -DUSE_MSVC_RUNTIME_LIBRARY_DLL:BOOL=OFF
#     )


# SETUP GLFW
add_library(GLFW_LIBRARY OBJECT IMPORTED)
set_target_properties(GLFW_LIBRARY PROPERTIES IMPORTED_OBJECTS ${GLFW_INSTALL_DIR}/lib/glfw3.lib)
set_property(TARGET GLFW_LIBRARY PROPERTY IMPORTED_LOCATION ${GLFW_INSTALL_DIR})

# SETUP GLAD
add_library(GLAD_LIBRARY OBJECT IMPORTED)
set_target_properties(GLAD_LIBRARY PROPERTIES IMPORTED_OBJECTS ${GLAD_INSTALL_DIR}/lib/glad.lib)
set_property(TARGET GLAD_LIBRARY PROPERTY IMPORTED_LOCATION ${GLAD_INSTALL_DIR})

# SETUP IMGUI
file(GLOB_RECURSE IMGUI_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui/*.cpp")
add_library(IMGUI STATIC ${IMGUI_SOURCES})
target_include_directories(IMGUI PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui"
    "${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui/backends"
    "${GLFW_INSTALL_DIR}/include"
)

####################################################################
### project files ##################################################
####################################################################

# list project files
# TODO when finishing project: list files explicitly as GLOB is not recommended
file(GLOB_RECURSE app CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/app/*.cu")
file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")
file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/ext/glm/glm/*.hpp")
file(GLOB_RECURSE headers CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.cuh")
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB_RECURSE sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")



# Create a single .lib file containing our compiled objects, and the compiled objects of all other dependencies
add_library(${PROJECT_NAME} STATIC ${app} ${headers} ${sources})


# Add all include file paths
target_include_directories(${PROJECT_NAME}
    PUBLIC ${GLAD_INSTALL_DIR}/include
    PUBLIC ${GLFW_INSTALL_DIR}/include
    PUBLIC ${GLM_INCLUDE}
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/app
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ext/imgui
)

# copy shaders
add_custom_target(copy_shaders ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/include/shaders
    ${CMAKE_CURRENT_BINARY_DIR}/dep/shaders
    COMMENT "Copying shaders into binary directory")

  # Specify the order in which libs depend on each other
add_dependencies(${PROJECT_NAME} ${libGLFW} ${libGLAD} copy_shaders)
target_link_libraries(${PROJECT_NAME} GLAD_LIBRARY GLFW_LIBRARY IMGUI)

####################################################################
### cuda configurations ############################################
####################################################################
# Request that particles be built with -std=c++17
# As this is a public compile feature anything that links to
# particles will also build with -std=c++17
# TODO: test rdc (shouldnt be required while header only)
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)
target_compile_options(${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
                        # --generate-line-info
                        --device-debug
                        --use_fast_math
                        # --relocatable-device-code=true
                        >)
target_compile_options(${PROJECT_NAME} PRIVATE $<$<CONFIG:DEBUG>:-G>)


# Set CUDA architecture
# default value: 75 to match testing system with nvidia 2080 card
# TODO variable with architecture/set by user
if(NOT DEFINED ${CMAKE_CUDA_ARCHITECTURES})
  set(CMAKE_CUDA_ARCHITECTURES 75 CACHE STRING "CUDA architectures" FORCE)
endif()
message(STATUS "CUDA architectures set to ${CMAKE_CUDA_ARCHITECTURES}")

# We need to explicitly state that we need all CUDA files in the
# particle library to be built with -dc as the member functions
# could be called by other libraries and executables
set_target_properties( ${PROJECT_NAME}
                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON
                       CUDA_RESOLVE_DEVICE_SYMBOLS ON
                       )



# add executable
add_executable(${EXECUTABLE_NAME}  app/main.cu)

target_include_directories(${EXECUTABLE_NAME}
    PUBLIC ${GLAD_INSTALL_DIR}/include
    PUBLIC ${GLFW_INSTALL_DIR}/include
    PUBLIC ${GLM_INCLUDE}
    PUBLIC ${GLAD_C_DIR}
)

set_property(TARGET ${EXECUTABLE_NAME}
             PROPERTY CUDA_SEPARABLE_COMPILATION ON
             CUDA_RESOLVE_DEVICE_SYMBOLS ON
             )
target_link_libraries(${EXECUTABLE_NAME} ${PROJECT_NAME})