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
#include <helper_math.h> // includes vector types
#include <functions.cuh> // additional CUDA functions
// #include <timer.h>               // timing functions

// CUDA helper functions
// #include <helper_cuda.h>         // helper functions for CUDA error check
// #include <helper_cuda_gl.h>      // helper functions for CUDA/GL interop

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

// TODO: different resolutions from monitor
const unsigned int window_width = 1920;
const unsigned int window_height = 1080;

const unsigned int mesh_width = 128;
const unsigned int mesh_height = 128;

const float epsilon = 0.001;
const float epsilon2 = epsilon * epsilon;

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
// 1. Basis: sample point p of a ray / set data...
// 2. Retrieve k clostest points to point p (in ascending order)
// 3. Test intersection with single atoms (point p in single SAS)
// 4. Test intersection with atom pairs (point p in 2 SAS)
// 5. Test intersection with atom triplets

__global__ void marching_kernel(float4 *pos, unsigned int width, unsigned int height)
{

  /////////////////////////////////////////////////
  // 0 // Preperations
  /////////////////////////////////////////////////
  unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;

  // calculate uv coordinates
  float u = x / (float)width;
  float v = y / (float)height;
  u = u * 2.0f - 1.0f;
  v = v * 2.0f - 1.0f;

  // distance to next surface
  float f = epsilon + 1;

  // starting position
  // TODO calculate starting position of ray
  float4 ray_pos;

  // ray direction
  // TODO calculate and normalize ray direction
  float4 ray_dir;

  // information to be filled on surface hit
  float4 point;
  float4 normal;

  unsigned int atom_count = 10;

  // set molecule data
  // TODO implement class for atom data structure (shared float4 array as pointer to be set in class constructor?)
  __shared__ float4 molecule[atom_count];

  molecule[0] = make_float4(-0.4, 0.0, 0.0, 0.4);
  molecule[1] = make_float4(0.5, 0.0, 0.0, 0.6);
  molecule[2] = make_float4(0.7, 0.7, 0.0, 0.8);
  molecule[3] = make_float4(-0.1, 0.7, -0.2, 0.5);
  molecule[4] = make_float4(-0.3, 0.2, 0.5, 1.0);
  molecule[5] = make_float4(0.1, 0.7, 0.3, 0.2);
  molecule[6] = make_float4(0.9, -0.2, 0.3, 0.1);
  molecule[7] = make_float4(-0.1, 0.8, -0.6, 1.1);
  molecule[8] = make_float4(-0.6, -0.7, 0.1, 0.6);
  molecule[9] = make_float4(0.1, -0.2, -0.3, 0.8);

  // ray marching loop
  while (f > epsilon)
  {
    // distance to ses for current marching step
    float fses = compute_SES(ray_pos, molecule, atom_count);

    // calculate new ray position
  }
  // end:
  // update pos to transfer data to cuda
  // TODO add second buffer for normals?
  pos[y * width + x] = make_float4(u, w, v, 1.0f);
}

// calculate distance to ses based on position ray_pos
// TODO finish returns
// TODO free memory of allocated arrays
__device__ void compute_SES(float4 ray_pos, float4 *molecule_data, const int atom_count)
{
  // hold current distance to SES
  float fses;

  /////////////////////////////////////////////////
  // 1 // determine k closest points
  /////////////////////////////////////////////////
  // allocate float array for all (squared) distances to point p
  float *distances;
  cudaMallocArray(distances, atom_count * sizeof(float));

  // allocate array for ids of k atoms closest to p
  unsigned int *k_nearest_keys;
  unsigned int nearest_count = min(k, atom_count);
  cudaMallocArray(k_nearest_keys, nearest_count * sizeof(int));

  // minimum distance to solvent accessible surface over all atoms
  float min_sas;

  // calculate distance
  // TODO: paper uses cubic voxels (3D grid to place atoms -> only calculate distances for 3x3x3 voxels instead of all atoms)
  for (unsigned int i = 0; i < atom_count; i++)
  {
    // calculate distance from sample point p to atom centers
    distances[i] = distance_squared(ray_pos, molecule[i]);

    // calculate minimum distance to SAS
    float sphere_sas = distances[i] - molecule[i].w * molecule[i].w;
    if (sphere_sas < min_sas || i == 0)
    {
      min_sas = sphere_sas;
    }

    // set key
    k_nearest_keys[i] = i;
  }

  // exit calculation if position of sample point p outside of SAS
  if (min_sas > 0)
  {
    fses = min_sas;
    return fses;
  }

  // determine k atoms with smallest distance to sample point p (sorted)
  k_nearest_keys = k_smallest_by_key(distances, atom_count, nearest_count);

  /////////////////////////////////////////////////
  // 2 // intersection with one atom
  /////////////////////////////////////////////////
  // TODO case 0 - when using voxels for nearest neighbours, there might be zero neighbours
  // Case 1 - test if intersection occurs with one atom surface directly

  // group of participating atoms is ordered, k_nearest_keys[0] contains the id of the closest atom
  float4 dir = normalize(ray_pos - molecule[k_nearest_keys[0]]);
  float4 surface_hit = molecule[k_nearest_keys[0]] + molecule[k_nearest_keys[0]].w * dir;

  float surface_min_dist = INFINITY;
  // check, if point is in SAS of another atom
  bool case1 = true;
  for (unsigned int i = 0; i < nearest_count; i++)
  {
    float dist_sas = distance_squared_sas(surface_hit, molecule[k_nearest_keys[i]]);
    if (dist_sas < surface_min_dist)
    {
      surface_min_dist = dist_sas;
      // if closest atom differs from atom with minimal distance to SAS -> case2
      if (k_nearest_keys[i] != 0)
      {
        case1 = false;
      }
      else
      {
        case1 = true;
      }
    }
    surface_min_dist = fmin(dist_sas, surface_min_dist);
  }

  // case1 = true -> case 2 and 3 dont have to be calculated

  // TODO check distance to surface (f_ses). Only when f_ses < epsilon calculate point & normal
  // TODO check whether to use epsilon or epsilon2
  if (case1 && abs(surface_min_dist) < epsilon)
  {
    // surface point is surface_hit
    // TODO case1 hit -> Calculate normal and send to shader
    fses = -surface_min_dist;
    return fses;
  }

  /////////////////////////////////////////////////
  // 3 // Intersection with atom pair
  /////////////////////////////////////////////////
  // Case 2 - test if intersection with atom pair occurs (arc)
  float4 arc_hit;

  for (unsigned int i = 1; i < nearest_count; i++)
  {
    for (unsigned int j = 0; j < i; j++)
    {
      // test predicate (p has to be in extended radius of both atom centers)
      bool case2_potential = case2_predicate(ray_pos, molecule[k_nearest_keys[i]], molecule[k_nearest_keys[j]], R);

      // if predicate adhered, calculate nearest surface point on arc between atoms
      if (case2_potential)
      {
        float4 intersec = sphere_sphere_intersec(ray_pos, molecule[k_nearest_keys[i]], molecule[k_nearest_keys[j]], R);
        // TODO: store intersection point for case3
        arc_hit = arc_surface_point(ray_pos, intersec, R);

        // test distance from sample point p to surface
        // TODO check whether to use epsilon or epsilon2
        if (distance_squared(arc_hit, ray_pos) < epsilon2)
        {
          // TODO exit on first point or continue for other possible surface hits? (-> case3)
          // TODO case2 hit -> Calculate normal and send to shader
          return;
        }
      }
    }
  }

  /////////////////////////////////////////////////
  // 4 // Intersection with atom triplet
  /////////////////////////////////////////////////
  // Case 3 - test if intersection with atom triple occurs
  float4 triplet_hit;

  // idea from paper:
  //  - only iterate over all intersection points from case2
  //  - test predicate for case3
  //  - calculate surface point

  return fses;
}

// template
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
void launch_kernel(float4 *pos, unsigned int mesh_width, unsigned int mesh_height)
{
  dim3 block(8, 8, 1);
  dim3 grid(mesh_width / block.x, mesh_height / block.y, 1);

  marching_kernel<<<grid, block>>>(float4 pos, mesh_width, mesh_height);
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