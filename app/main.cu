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
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

// includes cuda
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

// Utility functions
#include <helper_math.h> // includes vector types

// CUDA helper functions
#include <helper_cuda.h> // helper functions for CUDA error check
#include <cuda_inits.cuh>

#include <application_settings.hpp>
#include <raymarching_kernel.cuh>
#include <inputWindow.hpp>
// #include <performance.hpp>

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

// extern cudaGraphicsResource_t cuda_resource_1, cuda_resource_2;
extern cudaGraphicsResource_t *cuda_resources;

// timing
float deltaTime = 0.0f; // time between current frame and last frame
float lastFrame = 0.0f;

std::chrono::system_clock::time_point start_time = std::chrono::system_clock::now();

extern GLuint *cuda_textures; // cuda_textures used as buffer between CUDA and OpenGL

/////////////////////////////////////////////////////////////////////////////
//////// FUNCTIONS //////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

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
  // glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

  // setup imgui
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGuiIO &io = ImGui::GetIO();
  (void)io;
  // io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls
  ImGui::StyleColorsDark();
  ImGui_ImplGlfw_InitForOpenGL(window, true);
  ImGui_ImplOpenGL3_Init();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // initializations
  AppSettings *settings = new AppSettings(window_width, window_height, cuda_textures[0]); // initialize application settings
  iWindow->linkSettings(settings);

  createTextureBuffers(settings->getWindowWidth(), settings->getWindowHeight()); // create textures
  Shader quadTex;                                                                // create Shader
  quadTex.init("fb_vertex.vs", "fb_fragment.fs");                                // load shaders
  quadTex.use();                                                                 // use shader
  int VAO = quadBuffers();                                                       // create and bind buffers for quad (used to display the molecule using a texture)

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
  settings->show_performance_tool = false; // activate performance measuring
  // PerformanceCounter performance(60);      // add performance counter
  settings->show_gui = true;

  settings->update();
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // start rendering mainloop
  while (!glfwWindowShouldClose(window))
  {
    // time settings
    float currentFrame = glfwGetTime();
    deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;

    // handle inputs
    glfwPollEvents();
    processInput(window, settings, deltaTime); // handle keyboard events

    // start imgui frame
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    if (settings->show_performance_tool)
    {
      setupMetricsTool(settings, start_time);
      // ImGui::ShowMetricsWindow(&settings->show_performance_tool);
    }
    if (settings->show_gui)
    {
      setupImgui(settings);
    }

    // ImGui::ShowDemoWindow();

    // run CUDA kernel to generate vertex positions
    runCuda(cam, settings->getAllHostParams(), settings->getDeviceMolecule(), settings->getDeviceColors(), settings->getDeviceVoxelData(), settings->getDeviceVoxelCount(), settings->frame);

    if (settings->frame == 10)
    {
      test_kernel<<<1, 1>>>(settings->frame);
    }

    // display molecule

    glClearColor(0, 0, 0, 1);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);

    glBindVertexArray(VAO);

    glBindTexture(GL_TEXTURE_2D, cuda_textures[0]);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // display imgui
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

    // limit frames per second
    // TODO: sleep function differs under ubuntu and windows
    if (settings->limit_framerate)
    {
      float endTime = glfwGetTime();
      float deltaFrame = endTime - lastFrame;
      float sleepTime = ((1.0f / settings->frameLimit) - deltaFrame) * 1000.0f;
      if (sleepTime > 0)
        Sleep(sleepTime);
    }

    // debugging
    if (settings->getDebugMode())
    {
      int value = settings->getDebugFrame();
      value++;
      settings->setDebugFrame(value);
      settings->update();
    }

    // swap buffers
    glfwSwapBuffers(window);

    settings->frame++;
    if (settings->frame > INFINITY)
      settings->frame = 0;
  }

  cuda_cleanup();

  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  ImGui::DestroyContext();

  delete iWindow;
  delete settings;

  cudaDeviceReset();

  glfwTerminate();
  exit(EXIT_SUCCESS);
  return 0;
}
