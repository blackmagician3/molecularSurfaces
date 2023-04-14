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

// includes cuda
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

// Utility functions
// #include <helper_functions.h>    // includes cuda.h and cuda_runtime_api.h
// #include <timer.h>               // timing functions

// CUDA helper functions
// #include <helper_cuda.h>         // helper functions for CUDA error check
// #include <helper_cuda_gl.h>      // helper functions for CUDA/GL interop


/////////////////////////////////////////////////////////////////////////////
//////// CONSTANTS //////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// TODO: different resolutions from monitor
const unsigned int window_width  = 1920;
const unsigned int window_height = 1080;


const unsigned int mesh_width    = 128;
const unsigned int mesh_height   = 128;

const char application_name[] = "Molecular Surfaces";



/////////////////////////////////////////////////////////////////////////////
//////// STEPS //////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// The steps are:
// 1. Create an empty vertex buffer object (VBO)
// 2. Register the VBO with Cuda
// 3. Map the VBO for writing from Cuda
// 4. Run Cuda kernel to modify the vertex positions
// 5. Unmap the VBO
// 6. Render the results using OpenGL
/////////////////////////////////////////////////////////////////////////////

GLuint vbo;
struct cudaGraphicsResource *cuda_vbo;
void *d_vbo_buffer = NULL;
float g_fAnim = 0.0; // TODO g_fAnim, what's that? Animation speed?

/////////////////////////////////////////////////////////////////////////////
//////// DECLARATIONS ///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

bool runTest(int argc, char **argv, char *ref_file);
void cleanup();

// GL functionality
bool initGL();
void createVBO(GLuint *vbo, struct cudaGraphicsResource **vbo_res, unsigned int vbo_res_flags);
void deleteVBO(GLuint *vbo, struct cudaGraphicsResource *vbo_res);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);

// Cuda functionality
void runCuda(struct cudaGraphicsResource **vbo_resource);
void runAutoTest(int devID, char **argv, char *ref_file);
void checkResultCuda(int argc, char **argv, const GLuint &vbo);

/////////////////////////////////////////////////////////////////////////////
//////// KERNEL /////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

__global__ void simple_vbo_kernel(float4 *pos, unsigned int width, unsigned int height, float time)
{
    unsigned int x = blockIdx.x*blockDim.x + threadIdx.x;
    unsigned int y = blockIdx.y*blockDim.y + threadIdx.y;

    // calculate uv coordinates
    float u = x / (float) width;
    float v = y / (float) height;
    u = u*2.0f - 1.0f;
    v = v*2.0f - 1.0f;

    // calculate simple sine wave pattern
    float freq = 4.0f;
    float w = sinf(u*freq + time) * cosf(v*freq + time) * 0.5f;

    // write output vertex
    pos[y*width+x] = make_float4(u, w, v, 1.0f);
}

// to launch the kernel
void launch_kernel(float4 *pos, unsigned int mesh_width, unsigned int mesh_height, float time)
{
dim3 block(8, 8, 1);
dim3 grid(mesh_width / block.x, mesh_height / block.y, 1);
simple_vbo_kernel<<< grid, block>>>(pos, mesh_width, mesh_height, time);
}

// determine capable gpu hardware
int findGraphicsGPU(char *name)
{
  //TODO: add code to determine installed graphics card(s). Requires helper functions currently not included
    int nGraphicsGPU = 1;

    return nGraphicsGPU;
}

/////////////////////////////////////////////////////////////////////////////
//////// MAIN ///////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  //char *ref_file = NULL;

  //init GL
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  GLFWwindow* window = glfwCreateWindow(window_width, window_height, application_name, NULL, NULL);
  if (window == NULL)
  {
    std::cout << "Failed to create GLFW window" << std::endl;
    glfwTerminate();
    return -1;
  }
  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
  {
    std::cout << "Failed to initialize GLAD" << std::endl;
    return -1;
  }

  // default initialization
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glDisable(GL_DEPTH_TEST);


  // create VBO
  createVBO(&vbo, &cuda_vbo, cudaGraphicsMapFlagsWriteDiscard);

  // run the cuda part
  runCuda(&cuda_vbo);

  // start rendering mainloop
  while (!glfwWindowShouldClose(window))
  {
      // run CUDA kernel to generate vertex positions
      runCuda(&cuda_vbo);

      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


      // TODO: set view matrix, camera, ...


      // render from the vbo
      glBindBuffer(GL_ARRAY_BUFFER, vbo);

      glClearColor(0.2f, 0.3f, 0.4f, 1.0f);


      glfwSwapBuffers(window);
      glfwPollEvents();


      g_fAnim += 0.01f;

      //TODO: implement fps counter

  }
  std::cout << "task completed" << std::endl;

  glfwTerminate();
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//////// Initialize GL //////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
bool initGL()
{

  return true;
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
  glViewport(0, 0, width, height);
}

/////////////////////////////////////////////////////////////////////////////
//////// Run test ///////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
bool runTest(int argc, char **argv, char *ref_file)
{


  return true;
}

/////////////////////////////////////////////////////////////////////////////
//////// Cuda Part //////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void runCuda(struct cudaGraphicsResource **vbo_resource)
{
  // map OpenGL buffer object for writing from CUDA
  float4 *dptr;
  cudaGraphicsMapResources(1, vbo_resource, 0);
  size_t num_bytes;
  cudaGraphicsResourceGetMappedPointer((void **)&dptr, &num_bytes, *vbo_resource);

  //printf("CUDA mapped VBO: May access %ld bytes\n", num_bytes);

  launch_kernel(dptr, mesh_width, mesh_height, g_fAnim);

  // unmap buffer object
  cudaGraphicsUnmapResources(1, vbo_resource, 0);
}

/////////////////////////////////////////////////////////////////////////////
//////// Create VBO /////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void createVBO(GLuint *vbo, struct cudaGraphicsResource **vbo_res, unsigned int vbo_res_flags)
{
//assert(vbo);

// create buffer object
glGenBuffers(1, vbo);
glBindBuffer(GL_ARRAY_BUFFER, *vbo);

// initialize buffer object
unsigned int size = mesh_width * mesh_height * 4 * sizeof(float);
glBufferData(GL_ARRAY_BUFFER, size, 0, GL_DYNAMIC_DRAW);

glBindBuffer(GL_ARRAY_BUFFER, 0);

// register this buffer object with CUDA
cudaGraphicsGLRegisterBuffer(vbo_res, *vbo, vbo_res_flags);

//SDK_CHECK_ERROR_GL();
}

/////////////////////////////////////////////////////////////////////////////
//////// Delete VBO /////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void deleteVBO(GLuint *vbo, struct cudaGraphicsResource *vbo_res)
{
  // unregister this buffer object with CUDA
  cudaGraphicsUnregisterResource(vbo_res);

  glBindBuffer(1, *vbo);
  glDeleteBuffers(1, vbo);

  *vbo = 0;
}

/////////////////////////////////////////////////////////////////////////////
//////// Display Callback ///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void display()
{

}


void cleanup()
{
    if (vbo)
    {
        deleteVBO(&vbo, cuda_vbo);
    }

    // cudaDeviceReset causes the driver to clean up all state. While
    // not mandatory in normal operation, it is good practice.  It is also
    // needed to ensure correct operation when the application is being
    // profiled. Calling cudaDeviceReset causes all profile data to be
    // flushed before the application exits
    cudaDeviceReset();
}