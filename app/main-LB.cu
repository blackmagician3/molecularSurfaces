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
const unsigned int window_width = 1920;
const unsigned int window_height = 1080;

const unsigned int mesh_width = 128;
const unsigned int mesh_height = 128;

const char application_name[] = "Molecular Surfaces";

// TODO set in program
const float R = 10;        // constant solvent radius
const unsigned int k = 20; // maximum number of nearby atoms

/////////////////////////////////////////////////////////////////////////////
//////// STEPS //////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// The steps are:
// 1. Create an empty vertex buffer object (VBO)
// 2. Register the VBO with Cuda
// 3. Map the VBO for writing from Cuda
// 4. Run Cuda kernel to calculate rays, points and normals
// 5. Unmap the VBO
// 6. Render the results using OpenGL
/////////////////////////////////////////////////////////////////////////////

GLuint vbo;
struct cudaGraphicsResource *cuda_vbo;
void *d_vbo_buffer = NULL;
float g_fAnim = 0.0; // animation time

/////////////////////////////////////////////////////////////////////////////
//////// DECLARATIONS ///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

bool runTest(int argc, char **argv, char *ref_file);
void cleanup();

// GL functionality
void createVBO(GLuint *vbo, struct cudaGraphicsResource **vbo_res, unsigned int vbo_res_flags);
void deleteVBO(GLuint *vbo, struct cudaGraphicsResource *vbo_res);
void framebuffer_size_callback(GLFWwindow *window, int width, int height);

// Cuda functionality
void runCuda(struct cudaGraphicsResource **vbo_resource);

/////////////////////////////////////////////////////////////////////////////
//////// KERNEL /////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Method by Parulek/Viola:
// 1. Basis: sample point p of a ray
// 2. Retrieve k clostest points to point p (in ascending order)
// 3. Generate f_SES function (??? paper wants to generate the function twice, but intersection points are unknown as of here)
// 4. Test intersection (solvent extended) surfaces of atom pairs
// 5. Test intersection (solvent extended) surfaces of atom triplets
// 6. Generate f_SES with all stored intersection points

__global__ void compute_SES(float4 ray_pos, unsigned int width, unsigned int height)
{
  unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;

  // 1 //
  unsigned int atom_count = 100;

  // 2 //
  // allocate memory for k nearby atoms
  // TODO: compare cuda Mappoc speed with new
  float *group_S_id;
  float *group_S_dist;
  cudaMallocArray(group_S_id, k * sizeof(float));
  cudaMallocArray(group_S_dist, k * sizeof(float));

  // initialize distances with 2R
  for (unsigned int i = 0; i < k; i++)
  {
    group_S_dist[i] = 2 * R;
  }

  float max_dist = 2 * R;
  float max_id = 0;
  unsigned int num_participants = 0;

  // calculate distance from point p to all atoms
  // TODO: paper uses cubic voxels (3D grid to place atoms -> only calculate distances for 3x3x3 voxels instead of all atoms)
  for (unsigned int i = 0; i < atom_count; i++)
  {

    float dist = 0f;
    // TODO implement distance function
    // dist = ray_pos*ray_pos-a_pos*a_pos
    if (dist < max_dist)
    {
      group_S_id[max_id] = i;
      group_S_dist[max_id] = dist;

      // TODO update max_dist + max_id
      // max_dist = ... return max of array
      // max_id = ...   return id corresponding to maximal array value
      num_participants++;
    }
  }
  int a_size;
  if (num_participants < k)
  {
    // TODO: resize array
    a_size = num_participants;
  }
  else
  {
    a_size = k;
  }

  // TODO sort array in ascending order of distance
  cudaFree(group_S_dist); // distances not necessary anymore after sort

  // 3 //

  // 4 // test pairwise intersection
  // TODO correct calculation
  float4 intersec_2D;
  float intersec_2D_dist = 2 * R;

  float4 intersec_3D;
  float intersec_3D_dist = 2 * R;

  for (unsigned int i = 0; i < a_size; i++)
  {
    for (unsigned int j = 0; j < i; j++)
    {
      // TODO: perform pairwise test, input: i, j // return true/false
      bool test_2D;
      if (test_2D)
      {
        // TODO calculate intersection point (2D) name: i_pos
        // TODO calculate distance to point p (2D) name: dist
        float i_dist;
        float4 i_pos;
        if (i_dist < intersec_2D_dist)
        {
          intersec_2D_dist = i_dist;
          intersec_2D = i_pos;
        }

        for (unsigned int m = 0; m < j; m++)
        {
          // TODO: perform triplet test, input: i, j, m // return true/false
          bool test_3D;
          if (test_3D)
          {
            // TODO calculate intersection point (3D)
            // TODO calculate distance to point p (3D)
            // float i_dist; // set with 3D data
            // float4 i_pos; // set with 3D data
            if (i_dist < intersec_3D_dist)
            {
              intersec_3D_dist = i_dist;
              intersec_3D = i_pos;
            }
          }
        }
      }
    }
  }
}

__global__ void simple_vbo_kernel(float4 *pos, unsigned int width, unsigned int height, float time)
{
  unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;

  // calculate uv coordinates
  float u = x / (float)width;
  float v = y / (float)height;
  u = u * 2.0f - 1.0f;
  v = v * 2.0f - 1.0f;

  // calculate simple sine wave pattern
  float freq = 4.0f;
  float w = sinf(u * freq + time) * cosf(v * freq + time) * 0.5f;

  // write output vertex
  pos[y * width + x] = make_float4(u, w, v, 1.0f);
}

// to launch the kernel
void launch_kernel(float4 *pos, unsigned int mesh_width, unsigned int mesh_height, float time)
{
  dim3 block(8, 8, 1);
  dim3 grid(mesh_width / block.x, mesh_height / block.y, 1);
  simple_vbo_kernel<<<grid, block>>>(pos, mesh_width, mesh_height, time);
}

// determine capable gpu hardware
int findGraphicsGPU(char *name)
{
  // TODO: add code to determine installed graphics card(s). Requires helper functions currently not included
  int nGraphicsGPU = 1;

  return nGraphicsGPU;
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
  GLFWwindow *window = glfwCreateWindow(window_width, window_height, application_name, NULL, NULL);
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

    // TODO: implement fps counter
  }
  std::cout << "task completed" << std::endl;

  glfwTerminate();
  return 0;
}

/////////////////////////////////////////////////////////////////////////////
//////// Initialize GL //////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void framebuffer_size_callback(GLFWwindow *window, int width, int height)
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

  // printf("CUDA mapped VBO: May access %ld bytes\n", num_bytes);

  launch_kernel(dptr, mesh_width, mesh_height, g_fAnim);

  // unmap buffer object
  cudaGraphicsUnmapResources(1, vbo_resource, 0);
}

/////////////////////////////////////////////////////////////////////////////
//////// Create VBO /////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
void createVBO(GLuint *vbo, struct cudaGraphicsResource **vbo_res, unsigned int vbo_res_flags)
{
  // assert(vbo);

  // create buffer object
  glGenBuffers(1, vbo);
  glBindBuffer(GL_ARRAY_BUFFER, *vbo);

  // initialize buffer object
  unsigned int size = mesh_width * mesh_height * 4 * sizeof(float);
  glBufferData(GL_ARRAY_BUFFER, size, 0, GL_DYNAMIC_DRAW);

  glBindBuffer(GL_ARRAY_BUFFER, 0);

  // register this buffer object with CUDA
  cudaGraphicsGLRegisterBuffer(vbo_res, *vbo, vbo_res_flags);

  // SDK_CHECK_ERROR_GL();
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