/**
 * @file raymarching.cuh
 * @author Luka Blümler (luka.bluemler@uni-jena.de)
 * @brief device code for raymarching focused on the kernel setup
 *
 */

#ifndef RAYMARCHING_KERNEL_HPP
#define RAYMARCHING_KERNEL_HPP

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <helper_cuda.h> // helper functions for CUDA error check
#include <camera.h>
#include <surfaces.hpp>
#include <application_settings.hpp>
#include <simulation_config.hpp>

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// CUDA resources to write to OpenGL texture
cudaGraphicsResource_t cuda_resource_1, cuda_resource_2;     // cuda resource to write to
cudaArray *cuda_tex_array_1, *cuda_tex_array_2;              // cuda array to access texture
cudaResourceDesc cuda_resource_desc_1, cuda_resource_desc_2; // resource description for surface
/////////////////////////////////////////////////////////////////////////////////////////////////////////
const int TEST_SIZE = 10;
// simulation parameters in constant memory
__constant__ SimulationParams params;

void createResourceDesciption()
{
    // initialize resource desciption for surface
    memset(&cuda_resource_desc_1, 0, sizeof(cuda_resource_desc_1));
    cuda_resource_desc_1.resType = cudaResourceTypeArray;
}

void setParameters(SimulationParams hostParams)
{
    // copy parameters to constant memory
    checkCudaErrors(cudaMemcpyToSymbol(params, &hostParams, sizeof(SimulationParams)));
}

__device__ void setFunc(float4 *array, int array_length)
{
    float4 value = make_float4(0.0f, 0.42f, 0.42f, 10.3f);
    for (int i = 0; i < array_length; i++)
    {
        array[i] = value;
    }
}
__device__ void printFunc(float4 *array, int array_length, int x, int y)
{
    if (x == 100 && y == 100)
    {
        for (int i = 0; i < array_length; i++)
        {
            printf("array[%d]: %f, %f, %f, %f\n", i, array[i].x, array[i].y, array[i].z, array[i].w);
        }
    }
}

__device__ float4 calculateLighting(SimulationParams params, float4 colour, float4 position, float4 normal, float4 cam_pos)
{
    // colors
    float4 light_colour = make_float4(1.0f, 1.0f, 1.0f, 0.0f);
    float4 light_position = make_float4(5.0f, 0.0f, 10.0f, 0.0f);

    // ambient
    float ambientStrength = 0.1f;
    float4 ambient = ambientStrength * light_colour;

    // diffuse
    float4 norm = normalize(normal);
    float4 light_dir = normalize(light_position - position);

    float diff = max(dot(norm, light_dir), 0.0);
    float4 diffuse = diff * light_colour;

    // specular
    float specularStrength = 0.6f;
    float4 viewDir = normalize(cam_pos - position);
    float4 reflectDir = reflect(-light_dir, norm);

    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
    float4 specular = specularStrength * spec * light_colour;

    // final color
    float4 result = (ambient + diffuse + specular) * colour;

    return result;
}
__device__ float4 calculateNormal(float4 p, float4 surface_orientation, int case_id, bool debug)
{
    switch (case_id)
    {
    case 1:
        return normalize(p - surface_orientation);
        break;
    case 2:
        return normalize(surface_orientation - p);
        break;
    case 3:
        return normalize(surface_orientation - p);
        break;
    default:
        if (debug)
            printf("normal calculation missing defined case\n");
        break;
    }
}

__global__ void marching_kernel(cudaSurfaceObject_t surface, float4 *molecule, float4 screen_center, float4 cam_focus, float4 cam_right, float4 cam_up, float4 cam_position, float4 cam_front)
{
    /////////////////////////////////////////////////
    // 0 // Preperations
    /////////////////////////////////////////////////
    int x = (blockIdx.x * blockDim.x) + threadIdx.x;
    int y = (blockIdx.y * blockDim.y) + threadIdx.y;
    if (x > params.window_width || y > params.window_height)
        return;

    // calculate uv coordinates
    float u = x / (float)params.window_width;
    float v = y / (float)params.window_height;
    u = u * 2.0f - 1.0f;
    v = v * 2.0f - 1.0f;

    // determine ray
    float4 ray_origin = screen_center + u * normalize(cam_right) + v * normalize(cam_up);
    float4 ray = normalize(ray_origin - cam_position);
    float4 ray_pos;

    float depth = params.depth_min;
    bool hit = false;
    int steps_max = 250;
    float4 b_colour = make_float4(0.3f, 0.3f, 0.3f, 1.0f);
    float4 colour;
    float4 surface_orientation;
    int case_id = 0;

    // ray marching loop
    for (int i = 0; i < steps_max; i++)
    {

        // calculate new ray position
        ray_pos = cam_position + depth * ray;

        // compute distance to surface

        float f_sdf = computeSurface(ray_pos, molecule, params, &hit, &colour, &surface_orientation, &case_id, x, y);

        // for SES extend distance by solvent radius
        f_sdf += params.solvent_radius;

        depth += f_sdf;
        if (f_sdf < params.epsilon)
        {
            if (params.debug_mode && (x == (int)params.mouse_x_pos) && (y == (int)params.mouse_y_pos) && case_id == 4)
            {
                printf("case 4 with f_sdf = 0\n");
            }
            break;
        }
        if (depth > params.depth_max)
        {
            colour = b_colour;
            break;
        }
    }

    if (case_id != 0 && case_id != 4)
    {
        float4 normal = calculateNormal(ray_pos, surface_orientation, case_id, params.debug_mode);
        colour = calculateLighting(params, colour, ray_pos, normal, cam_position);
    }

    surf2Dwrite(colour, surface, x * sizeof(float4), y);
}
void runCuda(Camera *cam, SimulationParams host_params, float4 *molecule)
{
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // thread/grid setup

    dim3 block(host_params.thread_x, host_params.thread_y);
    dim3 grid(host_params.window_width / block.x, host_params.window_height / block.y);
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // camera setup
    // TODO: organise as array (cam_properties)
    float4 cam_center = cam->getCenterAsFloat4();
    float4 cam_focus = cam->getFocusAsFloat4();
    float4 cam_right = cam->getRightAsFloat4();
    float4 cam_up = cam->getUpAsFloat4();
    float4 cam_pos = cam->getPosAsFloat4();
    float4 cam_front = cam->getFrontAsFloat4();

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // map OpenGL buffer object for writing from CUDA
    // size_t num_bytes_1, num_bytes_2;
    // cudaArray *resArray_1;             //*resArray_2;
    cudaSurfaceObject_t surf_object_1; // surf_object_2;

    checkCudaErrors(cudaGraphicsMapResources(1, &cuda_resource_1, (cudaStream_t)0));                  // map resource 1 for even frames
    checkCudaErrors(cudaGraphicsSubResourceGetMappedArray(&cuda_tex_array_1, cuda_resource_1, 0, 0)); // get pointer to texture 1 for cuda
    cuda_resource_desc_1.res.array.array = cuda_tex_array_1;                                          // set cuda resource description for surface 1
    checkCudaErrors(cudaCreateSurfaceObject(&surf_object_1, &cuda_resource_desc_1));                  // create surface object 1
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // run kernel and update texture 1 via surface object
    marching_kernel<<<grid, block>>>(surf_object_1,
                                     molecule,
                                     cam_center,
                                     cam_focus,
                                     cam_right,
                                     cam_up,
                                     cam_pos,
                                     cam_front);
    cudaDeviceSynchronize();
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // unmap resources
    checkCudaErrors(cudaGraphicsUnmapResources(1, &cuda_resource_1, (cudaStream_t)0));

    ////////////////////////////////////////////////////////////////////////////////////////////////
}

#endif