/////////////////////////////////////////////////////////////////////////////
//////// INCLUDES ///////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// includes, system
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <windows.h>

// includes OpenGL
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "helper_gl.h"
#include "shader.h"
#include "camera.hpp"

// includes GUI
#include <nanogui/nanogui.h>

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
#include <performance.hpp>

/////////////////////////////////////////////////////////////////////////////
//////// CONSTANTS //////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// TODO check if defines are necessary

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

extern cudaGraphicsResource_t cuda_resource_1, cuda_resource_2;

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
  // nanogui::init();
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

  // setup GUI
  nanogui::Screen *screen = new nanogui::Screen;
  screen->initialize(window, true);
  // nanogui::Screen app{{(int)window_width, (int)window_height}, application_name};

  // nanogui::Window window{&app, ""};

  // printf("DEBUG: before GUI setup\n");
  screen->set_position({15, 15});
  screen->set_layout(new nanogui::GroupLayout(5, 5, 0, 0));

  nanogui::Label *l = new nanogui::Label(screen, "MODULATION", "sans-bold");
  l->set_font_size(10);
  nanogui::Slider *slider = new nanogui::Slider(screen);
  slider->set_value(0.5f);
  float modulation = 5.0f;
  slider->set_callback([&modulation](float value)
                       { modulation = value * 10.0f; });

  // Do the layout calculations based on what was added to the GUI
  screen->perform_layout();
  // printf("DEBUG: gui draw init\n");
  screen->draw_all();

  screen->set_visible(true);

  // register callbacks
  glfwSetCursorPosCallback(window, mouse_callback);
  glfwSetScrollCallback(window, scroll_callback);
  // glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

  // iWindow->setupGUI();
  // screenPtr = iWindow->getScreenPointer();
  // screenPtr->draw_all();
  // screenPtr->set_visible(true);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // initializations
  AppSettings *settings = new AppSettings(window_width, window_height, texture_1); // initialize application settings
  iWindow->linkSettings(settings);

  createTextureBuffers();                         // create textures
  Shader quadTex;                                 // create Shader
  quadTex.init("fb_vertex.vs", "fb_fragment.fs"); // load shaders
  quadTex.use();                                  // use shader
  int VAO = quadBuffers();                        // create and bind buffers for quad (used to display the molecule using a texture)

  // testing molecules
  std::string moleculePath = "C:/Users/lukab/OneDrive/Dokumente/repos/molecularSurfaces/testMolecules/3i40.pdb"; // 446 atoms
  // std::string moleculePath = "C:/Users/lukab/OneDrive/Dokumente/repos/molecularSurfaces/testMolecules/6hvb.pdb"; // 1410 atoms
  // std::string moleculePath = "C:/Users/lukab/OneDrive/Dokumente/repos/molecularSurfaces/testMolecules/5bny.pdb"; // 12460 atoms
  // std::string moleculePath = "C:/Users/lukab/OneDrive/Dokumente/repos/molecularSurfaces/testMolecules/2btv.pdb1";
  // std::string moleculePath = "";
  settings->setColorScheme(1);                              // chooses color scheme for visualization
  settings->loadMolecule(moleculePath);                     // load molecule data
  cam->intializeCameraPosition(settings->getCameraFocus()); // setup camera

  settings->setKnearest(10);
  settings->setVoxelUsage(true);
  settings->setEpsilon(0.01);

  // performance measuring (disabled by default)
  settings->changePerformanceDisplay(false); // activate performance measuring
  PerformanceCounter performance(60);        // add performance counter

  // frames per second limit (0 for unlimited)
  settings->changeFrameLimit(20);

  settings->update();
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  int frame = 0;
  // start rendering mainloop
  while (!glfwWindowShouldClose(window))
  {
    //////////////////////////////////////////////////////////////////////////////////////
    // time settings
    float currentFrame = glfwGetTime();
    deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;

    //////////////////////////////////////////////////////////////////////////////////////
    // performance measuring
    if (settings->getPerformanceDisplay())
    {
      performance.runMeasuring(deltaTime);
    }

    //////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////
    // handle inputs
    processInput(window, settings, deltaTime); // handle keyboard events

    //////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////
    // run CUDA kernel to generate vertex positions
    runCuda(cam, settings->getAllHostParams(), settings->getDeviceMolecule(), settings->getDeviceColors(), settings->getDeviceVoxelData(), settings->getDeviceVoxelCount(), frame);
    //////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////
    // display molecule

    // draw(window_width, window_height, VAO, &texture_1, &app);
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    // glBindVertexArray(VAO);
    // glDisable(GL_DEPTH_TEST);
    // glBindTexture(GL_TEXTURE_2D, texture_1);
    // glDrawArrays(GL_TRIANGLES, 0, 6);
    // glBindBuffer(GL_ARRAY_BUFFER, 0);
    // glBindVertexArray(0);
    // glBindTexture(GL_TEXTURE_2D, 0);
    // glDeleteTextures(1, &texture_1);
    screen->draw_contents();
    screen->draw_widgets();

    //////////////////////////////////////////////////////////////////////////////////////

    // limit frames per second
    // TODO: sleep function differs under ubuntu and windows
    if (settings->getFrameLimit() > 0)
    {
      float endTime = glfwGetTime();
      float deltaFrame = endTime - lastFrame;
      float sleepTime = ((1.0f / settings->getFrameLimit()) - deltaFrame) * 1000.0f;
      if (sleepTime > 0)
        Sleep(sleepTime);
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // debugging
    if (settings->getDebugMode())
    {
      int value = settings->getDebugFrame();
      value++;
      settings->setDebugFrame(value);
      settings->update();
    }
    //////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////

    // swap buffers
    glfwSwapBuffers(window);
    glfwPollEvents();

    frame++;
    if (frame > INFINITY)
      frame = 0;
  }
  cleanup(settings);
  delete iWindow;
  std::cout << "task completed" << std::endl;

  // // glfwTerminate();
  // nanogui::shutdown();
  // exit(EXIT_SUCCESS);
  return 0;
}
