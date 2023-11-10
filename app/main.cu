/////////////////////////////////////////////////////////////////////////////
//////// INCLUDES ///////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// includes, system
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes OpenGL
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "helper_gl.h"
#include "shader.h"
#include "camera.h"

// includes cuda
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

// Utility functions
#include <helper_math.h> // includes vector types

// #include <timer.h>               // timing functions

// CUDA helper functions
#include <helper_cuda.h> // helper functions for CUDA error check
// #include <helper_cuda_gl.h>      // helper functions for CUDA/GL interop

#include <application_settings.hpp>
#include <raymarching_kernel.cuh>
#include <inputWindow.hpp>

/////////////////////////////////////////////////////////////////////////////
//////// CONSTANTS //////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// TODO check if defines are necessary
#ifndef M_PI
#define M_PI 3.141592653589793
#endif
#ifndef INFINITY
#define INFINITY 1e8
#endif
bool first_debug = true;

// const float near = 0.1;
// const float far = 100.0f;

// GL stuff
// TODO: different resolutions from monitor
const unsigned int window_width = 1920;
const unsigned int window_height = 1080;

const char application_name[] = "Molecular Surfaces";

// timing
float deltaTime = 0.0f; // time between current frame and last frame
float lastFrame = 0.0f;

GLuint texture_1; // textures used as buffer between CUDA and OpenGL

/////////////////////////////////////////////////////////////////////////////
//////// FUNCTIONS //////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void createTextureBuffers()
{
  glEnable(GL_TEXTURE_2D);

  ////////////////////////////////////////////////////////////////////////////////
  // initialize texture 1
  glGenTextures(1, &texture_1);
  glBindTexture(GL_TEXTURE_2D, texture_1);

  // initialize texture with glTexImage2D
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F_ARB, window_width, window_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);

  // texture wrapping
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  // texture filtering
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

  // TODO implement second texture as buffer for read (OpenGL) and write (CUDA)

  glBindTexture(GL_TEXTURE_2D, 0);

  // register OpenGL textures with CUDA
  checkCudaErrors(cudaGraphicsGLRegisterImage(&cuda_resource_1, texture_1, GL_TEXTURE_2D, cudaGraphicsRegisterFlagsWriteDiscard));

  createResourceDesciption();
}

void cleanup(AppSettings *settings)
{
  cudaGraphicsUnregisterResource(cuda_resource_1);

  glDeleteTextures(1, &texture_1);

  delete settings;
  // cudaDeviceReset causes the driver to clean up all state. While
  // not mandatory in normal operation, it is good practice.  It is also
  // needed to ensure correct operation when the application is being
  // profiled. Calling cudaDeviceReset causes all profile data to be
  // flushed before the application exits
  cudaDeviceReset();
}
/////////////////////////////////////////////////////////////////////////////
//////// MAIN ///////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // init GL
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  InputWindow *iWindow = new InputWindow(window_width, window_height, application_name);
  if (iWindow->checkWindowCreation())
  {
    return -1;
  }
  GLFWwindow *window = iWindow->getWindowPointer();
  Camera *cam = getCameraPointer();

  // register callbacks
  glfwSetCursorPosCallback(window, mouse_callback);
  glfwSetScrollCallback(window, scroll_callback);
  glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // initializations
  int kernel_call = 0;
  AppSettings *settings = new AppSettings(window_width, window_height, texture_1); // initialize application settings
  iWindow->linkSettings(settings);

  createTextureBuffers();                         // create textures
  Shader quadTex;                                 // create Shader
  quadTex.init("fb_vertex.vs", "fb_fragment.fs"); // load shaders
  quadTex.use();                                  // use shader
  glm::mat4 model(0);                             // create model matrix
  int VAO = quadBuffers();                        // create and bind buffers for quad (used to display the molecule using a texture)

  settings->update();       // copy parameters to GPU
  settings->loadMolecule(); // load molecule data
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // start rendering mainloop
  while (!glfwWindowShouldClose(window))
  {
    //////////////////////////////////////////////////////////////////////////////////////
    // time settings
    kernel_call++;
    // std::cout << "kernel call " << kernel_call << std::endl;
    // use delta time
    float currentFrame = glfwGetTime();
    deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;
    //////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////
    // handle inputs
    processInput(window, settings, deltaTime); // handle keyboard events

    // glm::vec2 cam_mov = mouse_offset(window); // get mouse movement
    // cam->ProcessMouseMovement(cam_mov.x, cam_mov.y);
    // glm::mat4 view = cam->GetViewMatrix(); // view matrix
    // glm::mat4 proj = glm::perspective(glm::radians(cam->Zoom), (float)window_width / (float)window_height, near, far); // projection matrix
    // glm::mat4 proj = glm::mat4(1.0);
    // quadTex.setMat("view", view);
    // quadTex.setMat("projection", proj);
    //////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////
    // run CUDA kernel to generate vertex positions
    runCuda(cam, settings->getAllHostParams(), settings->getDeviceMolecule());
    //////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////
    // display molecule
    draw(window_width, window_height, VAO, &texture_1);

    //////////////////////////////////////////////////////////////////////////////////////#

    glfwSwapBuffers(window);
    glfwPollEvents();

    // TODO: implement fps counter
  }
  cleanup(settings);
  delete iWindow;
  std::cout << "task completed" << std::endl;

  glfwTerminate();
  return 0;
}
