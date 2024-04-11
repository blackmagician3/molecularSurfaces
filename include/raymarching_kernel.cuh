/**
 * @file raymarching.cuh
 * @author Luka Blümler (luka.bluemler@uni-jena.de)
 * @brief device code for raymarching focused on the kernel setup
 *
 */

#ifndef RAYMARCHING_KERNEL_CUH
#define RAYMARCHING_KERNEL_CUH

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <helper_cuda.h> // helper functions for CUDA error check
#include <camera.hpp>
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
const int MAX_STEP = 10000;

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

__device__ float4 convertHexToRGB(uint hexColor)
{
    // convert hex color to rgb
    float4 color;
    color.x = ((hexColor >> 16) & 0xFF) / 255.0f; // Extract the RR byte
    color.y = ((hexColor >> 8) & 0xFF) / 255.0f;  // Extract the GG byte
    color.z = ((hexColor) & 0xFF) / 255.0f;       // Extract the BB byte
    color.w = 1.0f;
    return color;
}
__device__ float4 convertAndWeight(uint hexColor, float weight)
{
    // convert hex color to rgb
    float4 color;
    color.x = weight * (((hexColor >> 16) & 0xFF) / 255.0f); // Extract the RR byte
    color.y = weight * (((hexColor >> 8) & 0xFF) / 255.0f);  // Extract the GG byte
    color.z = weight * (((hexColor) & 0xFF) / 255.0f);       // Extract the BB byte
    color.w = 1.0f;

    return color;
}

__device__ float4 blendColors(uint colorA, uint colorB, float weight1, float weight2)
{
    float4 cA = convertAndWeight(colorA, weight1);
    float4 cB = convertAndWeight(colorB, weight2);
    float4 color = cA + cB;
    color.w = 1.0f;

    return color;
}

__device__ float4 blendColors(uint colorA, uint colorB, uint colorC, float weights[])
{
    float4 cA = convertAndWeight(colorA, weights[0]);
    float4 cB = convertAndWeight(colorB, weights[1]);
    float4 cC = convertAndWeight(colorC, weights[2]);

    float4 color = cA + cB + cC;

    color.w = 1.0f;
    return color;
}

__device__ float4 getSurfaceColor(hitInfo *surfData, float4 *molecule, uint *colors, uint colorScheme, float solventRadius, float x, float y, SimulationParams params, float4 rayDir)
{
    // float4 color = make_float4(0.2f, 0.2f, 0.2f, 1.0f);
    float4 color = make_float4(1.0f, 1.0f, 1.0f, 1.0f);

    switch (surfData->collisionType)
    {
    case 0:
        color = make_float4(0.0f, 0.0f, 1.0f, 1.0f);
        break;
    case 1:
        if (colorScheme == 0)
        {
            color = make_float4(0.2f, 0.2f, 1.0f, 1.0f);
        }
        else
        {
            color = convertHexToRGB(colors[surfData->bondId1]);
        }
        break;
    case 2:
        if (colorScheme == 0)
        {
            color = make_float4(0.0f, 1.0f, 0.0f, 1.0f);
        }
        else
        {
            uint id1 = surfData->bondId1;
            uint id2 = surfData->bondId2;
            float r1 = molecule[id1].w;
            float r2 = molecule[id2].w;

            // use distance from atoms to p and normalize with length between atoms
            float atom1ToP = (length(surfData->rayPos - molecule[id1]) - (r1 - solventRadius));
            float atom2ToP = (length(surfData->rayPos - molecule[id2]) - (r2 - solventRadius));
            float sum = atom1ToP + atom2ToP;
            float weight2 = atom1ToP / sum;
            float weight1 = atom2ToP / sum;

            color = blendColors(colors[surfData->bondId1], colors[surfData->bondId2], weight1, weight2);
        }
        break;
    case 3:
        if (colorScheme == 0)
        {
            color = make_float4(0.0f, 1.0f, 0.0f, 1.0f);
        }
        else
        {
            int ids[3];
            ids[0] = surfData->bondId1;
            ids[1] = surfData->bondId2;
            ids[2] = surfData->bondId3;

            float r1 = molecule[ids[0]].w - solventRadius;
            float r2 = molecule[ids[1]].w - solventRadius;
            float r3 = molecule[ids[2]].w - solventRadius;

            float4 p1 = molecule[ids[0]] + normalize(surfData->surfaceHit - molecule[ids[0]]) * r1;
            float4 p2 = molecule[ids[1]] + normalize(surfData->surfaceHit - molecule[ids[1]]) * r2;
            float4 p3 = molecule[ids[2]] + normalize(surfData->surfaceHit - molecule[ids[2]]) * r3;

            float weights[3];

            // baricentric coordinates
            float4 u = p2 - p1;
            float4 v = p3 - p1;
            float4 w = surfData->rayPos - p1;
            float4 n = cross(u, v);

            weights[2] = dot(cross(u, w), n) / dot(n, n);
            weights[1] = dot(cross(w, v), n) / dot(n, n);
            weights[0] = 1 - weights[2] - weights[1];

            color = blendColors(colors[ids[0]], colors[ids[1]], colors[ids[2]], weights);
        }
        break;
    case 4:
        color = make_float4(1.0f, 0.0f, 1.0f, 1.0f);
        break;
    default:
        break;
    }
    return color;
}

__device__ float4 calculateLighting(SimulationParams params, float4 color, float4 position, float4 normal, float4 cam_pos)
{
    // TODO: remove point light source

    //  colors
    float4 light_color = make_float4(1.0f, 1.0f, 1.0f, 0.0f);
    // float4 light_position = make_float4(5.0f, 0.0f, 10.0f, 0.0f);
    float4 light_position = cam_pos;

    // ambient
    float ambientStrength = 0.1f;
    float4 ambient = ambientStrength * light_color;

    // diffuse
    float4 norm = normalize(normal);
    float4 light_dir = normalize(light_position - position);

    float diff = max(dot(norm, light_dir), 0.0);
    float4 diffuse = diff * light_color;

    // specular
    float specularStrength = 0.6f;
    float4 viewDir = normalize(cam_pos - position);
    float4 reflectDir = reflect(-light_dir, norm);

    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
    float4 specular = specularStrength * spec * light_color;

    // final color
    float4 result = (ambient + diffuse + specular) * color;

    return result;
}
__device__ float4 calculateNormal(hitInfo *surfacePointData, float4 *molecule, bool debug)
{
    switch (surfacePointData->collisionType)
    {
    case 0:
        return normalize(surfacePointData->rayPos - molecule[surfacePointData->bondId1]);
        break;
    case 1:
        return normalize(surfacePointData->rayPos - molecule[surfacePointData->bondId1]);
        break;
    case 2:
        return normalize(surfacePointData->surfaceHit - surfacePointData->rayPos);
        break;
    case 3:
        return normalize(surfacePointData->surfaceHit - surfacePointData->rayPos);
        break;
    default:
        if (debug)
            printf("ERROR: normal calculation missing defined case\n");
        break;
    }
}

// TODO: summarize camera values in struct
__global__ void marching_kernel(cudaSurfaceObject_t surface, float4 *molecule, uint *colors, int *voxel_data, int *voxel_count, float4 screen_center, float4 cam_focus, float4 cam_right, float4 cam_up, float4 cam_position, float4 cam_front, int frame = 0)
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
    // TODO: implement ray as struct for ease of use
    float4 ray_origin = screen_center + u * normalize(cam_right) + v * normalize(cam_up);
    float4 ray = normalize(ray_origin - cam_position);
    float4 ray_pos;

    float depth = params.depth_min;
    // TODO: max number of steps needs to be more precise
    int steps_max = 100;

    float4 b_color = make_float4(0.2f, 0.2f, 0.2f, 1.0f);
    // float4 b_color = make_float4(1.0f, 1.0f, 1.0f, 1.0f);
    float4 surfaceColor = b_color;

    // if (frame == 0 && x == 0 && y == 0)
    // {
    //     printf("Size of params: %i\n", (int)sizeof(params));
    // }

    hitInfo surfacePointData;
    surfacePointData.collisionType = 0;
    surfacePointData.isInGrid = false;
    surfacePointData.traversedGrid = false;

    surfaceColor = b_color;
    surfacePointData.collisionType = 42;

    // ray marching loop
    int step = 0;
    while (step < MAX_STEP)
    {
        // calculate new ray position
        ray_pos = cam_position + depth * ray;

        // compute distance to surface
        float f_sdf = computeSurface(ray_pos, ray, molecule, colors, voxel_data, voxel_count, params, &surfacePointData, x, y, frame, step);

        // for SES extend distance by solvent radius
        f_sdf += params.solvent_radius;

        depth += f_sdf;
        if (f_sdf < params.epsilon)
        {
            surfacePointData.rayPos = ray_pos + f_sdf * ray;
            // surfacePointData.rayPos = ray_pos;
            break;
        }
        if (depth > params.depth_max)
        {
            surfaceColor = b_color;
            surfacePointData.collisionType = 42;
            break;
        }
        if (surfacePointData.traversedGrid)
        {
            surfaceColor = b_color;
            surfacePointData.collisionType = 42;
            break;
        }
        if (step > steps_max)
        {
            surfaceColor = b_color;
            surfacePointData.collisionType = 42;
            break;
        }

        step++;
    }

    bool in_highlight = false;
    if (params.debug_mode)
    {
        float x_diff = params.mouse_x_pos - x;
        float y_diff = params.mouse_y_pos - y;
        in_highlight = ((x_diff * x_diff + y_diff * y_diff) <= (params.highlight_radius * params.highlight_radius));
        if (in_highlight)
        {
            surfaceColor = make_float4(0.89f, 0.17f, 0.74f, 1.0f);
        }
    }
    if (!params.debug_mode || !in_highlight)
    {
        surfaceColor = getSurfaceColor(&surfacePointData, molecule, colors, params.colorScheme, params.solvent_radius, x, y, params, ray);

        if (surfacePointData.collisionType < 4)
        {
            float4 normal = calculateNormal(&surfacePointData, molecule, params.debug_mode);
            surfaceColor = calculateLighting(params, surfaceColor, ray_pos, normal, cam_position);
        }
    }

    // if (surfacePointData.collisionType < 4)
    // {
    //     float4 normal = calculateNormal(&surfacePointData, molecule, params.debug_mode);
    //     surfaceColor = calculateLighting(params, surfaceColor, ray_pos, normal, cam_position);
    // }

    surf2Dwrite(surfaceColor, surface, x * sizeof(float4), y);

    // if (params.debug_mode && params.debug_frame == 0 && x == 728 && y == 423)
    // {
    //     printf("------DEBUG after write -------\n");
    //     printf("Pixel: %i, %i\n", x, y);
    //     printf("SurfaceHit: %.2f, %.2f, %.2f\n", surfacePointData.surfaceHit.x, surfacePointData.surfaceHit.y, surfacePointData.surfaceHit.z);
    //     printf("RayPos: %.2f, %.2f, %.2f\n", surfacePointData.rayPos.x, surfacePointData.rayPos.y, surfacePointData.rayPos.z);
    //     printf("RayDir: %.4f, %.4f, %.4f\n", ray.x, ray.y, ray.z);
    //     printf("Box_min: %.2f, %.2f, %.2f\n", params.box_start.x, params.box_start.y, params.box_start.z);
    //     printf("Box_max: %.2f, %.2f, %.2f\n", params.box_end.x, params.box_end.y, params.box_end.z);
    //     printf("BondID1: %i\n", surfacePointData.bondId1);
    //     printf("Colour: %.2f, %.2f, %.2f\n", surfaceColor.x, surfaceColor.y, surfaceColor.z);
    //     printf("CollisionType: %i\n", surfacePointData.collisionType);
    //     printf("-------------------------------\n");
    // }
}
void runCuda(Camera *cam, SimulationParams host_params, float4 *molecule, uint *colors, int *voxel_data, int *voxel_count, int frame = 0)
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

    if (host_params.debug_mode && host_params.debug_frame == 0)
    {
        printf("Cam Position: %.2f, %.2f, %.2f\n", cam_pos.x, cam_pos.y, cam_pos.z);
    }
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
    cudaDeviceSynchronize();
    marching_kernel<<<grid, block>>>(surf_object_1,
                                     molecule,
                                     colors,
                                     voxel_data,
                                     voxel_count,
                                     cam_center,
                                     cam_focus,
                                     cam_right,
                                     cam_up,
                                     cam_pos,
                                     cam_front,
                                     frame);
    cudaDeviceSynchronize();
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // unmap resources
    checkCudaErrors(cudaGraphicsUnmapResources(1, &cuda_resource_1, (cudaStream_t)0));

    ////////////////////////////////////////////////////////////////////////////////////////////////
}

#endif