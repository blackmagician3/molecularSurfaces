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
    // color.x = sqrt((1 - weight) * cA.x * cA.x + weight * cB.x * cB.x);
    // color.y = sqrt((1 - weight) * cA.y * cA.y + weight * cB.y * cB.y);
    // color.z = sqrt((1 - weight) * cA.z * cA.z + weight * cB.z * cB.z);
    // color.w = 1.0f;
    return color;
}
// __device__ float4 blendColors(uint colorA, uint colorB, uint colorC, float weights[])
// {
//     // float4 cA = convertHexToRGB(colorA);
//     // float4 cB = convertHexToRGB(colorB);
//     // float4 cC = convertHexToRGB(colorC);

//     float4 cA = convertAndWeight(colorA, weights[0]);
//     float4 cB = convertAndWeight(colorB, weights[1]);

//     float4 color1 = cA + cB;
//     float4 color = weights[2] * (color1) + convertAndWeight(colorC, weights[3]);
//     color.w = 1.0f;
//     // color.x = pow(weightA * pow(cA.x, 3) + weightB * pow(cB.x, 3) + weightC * pow(cC.x, 3), 1 / 3);
//     // color.y = pow(weightA * pow(cA.y, 3) + weightB * pow(cB.y, 3) + weightC * pow(cC.y, 3), 1 / 3);
//     // color.z = pow(weightA * pow(cA.z, 3) + weightB * pow(cB.z, 3) + weightC * pow(cC.z, 3), 1 / 3);
//     return color;
// }
__device__ float4 blendColors(uint colorA, uint colorB, uint colorC, float weights[])
{
    float4 cA = convertAndWeight(colorA, weights[0]);
    float4 cB = convertAndWeight(colorB, weights[1]);
    float4 cC = convertAndWeight(colorC, weights[2]);

    float4 color = cA + cB + cC;
    // float4 color = weights[0] * cA + weights[1] * cB + weights[2] * cC;

    color.w = 1.0f;
    return color;
}

__device__ float4 getSurfaceColor(hitInfo *surfData, float4 *molecule, uint *colors, uint colorScheme, float solventRadius, float x, float y, SimulationParams params, float4 rayDir)
{
    float4 color = make_float4(0.2f, 0.2f, 0.2f, 1.0f);

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
            if (params.debug_mode && (x == 768) && (y == 432))
            {
                // float4 c1 = convertHexToRGB(colors[surfData->bondId1]);
                // float4 c2 = convertHexToRGB(colors[surfData->bondId2]);
                // float4 c3 = convertHexToRGB(colors[surfData->bondId3]);

                // printf("weights: [%.3f, %.3f, %.3f], colors: [(%.2f, %.2f, %.2f), (%.2f, %.2f, %.2f), (%.2f, %.2f, %.2f)]\n",
                //        weights[0], weights[1], weights[2], c1.x, c1.y, c1.z, c2.x, c2.y, c3.z, c3.x, c3.y, c3.z);
                printf("pixel %i, %i:\n", (int)x, (int)y);
                printf("-------------------------------------------------------------------------\n");
                printf("p1.x/p1.y/p1.z/p2.x/p2.y/p2.z/p3.x/p3.y/p3.z/ray.x/ray.y/ray.z\n");
                printf("%f/%f/%f", p1.x, p1.y, p1.z);
                printf("/%f/%f/%f", p2.x, p2.y, p2.z);
                printf("/%f/%f/%f", p3.x, p3.y, p3.z);
                printf("/%f/%f/%f\n", surfData->rayPos.x, surfData->rayPos.y, surfData->rayPos.z);
                printf("-------------------------------------------------------------------------\n");
            }
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
    float4 light_position = make_float4(5.0f, 0.0f, 10.0f, 0.0f);

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
            printf("normal calculation missing defined case\n");
        break;
    }
}

__global__ void marching_kernel(cudaSurfaceObject_t surface, float4 *molecule, uint *colors, float4 screen_center, float4 cam_focus, float4 cam_right, float4 cam_up, float4 cam_position, float4 cam_front)
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
    // TODO: max number of steps needs to be more precise
    int steps_max = 250;

    float4 b_color = make_float4(0.2f, 0.2f, 0.2f, 1.0f);
    float4 surfaceColor = b_color;

    hitInfo surfacePointData;
    surfacePointData.collisionType = 0;

    // ray marching loop
    for (int i = 0; i < steps_max; i++)
    {

        // calculate new ray position
        ray_pos = cam_position + depth * ray;

        // compute distance to surface
        float f_sdf = computeSurface(ray_pos, molecule, colors, params, &surfacePointData, x, y);

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
    }

    surfaceColor = getSurfaceColor(&surfacePointData, molecule, colors, params.colorScheme, params.solvent_radius, x, y, params, ray);

    // if (params.debug_mode && (x == 704) && (y == 432))
    // {
    //     if (surfacePointData.collisionType == 3)
    //     {
    //         printf("testing pixel %i, %i:\n", x, y);
    //         printf("-----------------------------------------------------------------------------------------\n");
    //         printf("atom1.x/atom1.y/atom1.z/atom1.r/atom1.color/atom2.x/atom2.y/atom2.z/atom2.r/atom2.color/atom3.x/atom3.y/atom3.z/atom3.r/atom3.color");
    //         printf("/collision_type/ray_position.x/ray_position.y/ray_position.z/surface_hit.x/surface_hit.y/surface_hit.z");
    //         printf("/focus.x/focus.y/focus.z");
    //         printf("-----------------------------------------------------------------------------------------\n");
    //         // atom1
    //         printf("%f/%f/%f", molecule[surfacePointData.bondId1].x, molecule[surfacePointData.bondId1].y, molecule[surfacePointData.bondId1].z);
    //         printf("/%f/%i", molecule[surfacePointData.bondId1].w - params.solvent_radius, colors[surfacePointData.bondId1]);

    //         // atom2
    //         printf("/%f/%f/%f", molecule[surfacePointData.bondId2].x, molecule[surfacePointData.bondId2].y, molecule[surfacePointData.bondId2].z);
    //         printf("/%f/%i", molecule[surfacePointData.bondId2].w - params.solvent_radius, colors[surfacePointData.bondId2]);

    //         // atom3
    //         printf("/%f/%f/%f", molecule[surfacePointData.bondId3].x, molecule[surfacePointData.bondId3].y, molecule[surfacePointData.bondId3].z);
    //         printf("/%f/%i", molecule[surfacePointData.bondId3].w - params.solvent_radius, colors[surfacePointData.bondId3]);

    //         // general information
    //         printf("/%i", surfacePointData.collisionType);
    //         printf("/%f/%f/%f", surfacePointData.rayPos.x, surfacePointData.rayPos.y, surfacePointData.rayPos.z);
    //         printf("/%f/%f/%f", surfacePointData.surfaceHit.x, surfacePointData.surfaceHit.y, surfacePointData.surfaceHit.z);
    //         printf("/%f/%f/%f\n", cam_position.x, cam_position.y, cam_position.z);
    //         printf("-----------------------------------------------------------------------------------------\n");
    //     }
    //     else
    //     {
    //         printf("wrong pixel, no type 3 collision!\n");
    //     }
    //     // if (surfacePointData.collisionType == 2)
    //     // {
    //     //     printf("testing pixel %i, %i:\n", x, y);
    //     //     printf("-----------------------------------------------------------------------------------------\n");
    //     //     printf("atom1.x/atom1.y/atom1.z/atom1.r/atom1.color/atom2.x/atom2.y/atom2.z/atom2.r/atom2.color/atom3.x/atom3.y/atom3.z/atom3.r/atom3.color");
    //     //     printf("/collision_type/ray_position.x/ray_position.y/ray_position.z/surface_hit.x/surface_hit.y/surface_hit.z\n");
    //     //     printf("-----------------------------------------------------------------------------------------\n");
    //     //     // atom1
    //     //     printf("%f/%f/%f", molecule[surfacePointData.bondId1].x, molecule[surfacePointData.bondId1].y, molecule[surfacePointData.bondId1].z);
    //     //     printf("/%f/%i", molecule[surfacePointData.bondId1].w - params.solvent_radius, colors[surfacePointData.bondId1]);

    //     //     // atom2
    //     //     printf("/%f/%f/%f", molecule[surfacePointData.bondId2].x, molecule[surfacePointData.bondId2].y, molecule[surfacePointData.bondId2].z);
    //     //     printf("/%f/%i", molecule[surfacePointData.bondId2].w - params.solvent_radius, colors[surfacePointData.bondId2]);

    //     //     // general information
    //     //     printf("/%i", surfacePointData.collisionType);
    //     //     printf("/%f/%f/%f", surfacePointData.rayPos.x, surfacePointData.rayPos.y, surfacePointData.rayPos.z);
    //     //     printf("/%f/%f/%f", surfacePointData.surfaceHit.x, surfacePointData.surfaceHit.y, surfacePointData.surfaceHit.z);
    //     //     printf("/%f/%f/%f\n", cam_focus.x, cam_focus.y, cam_focus.z);
    //     //     printf("-----------------------------------------------------------------------------------------\n");
    //     // }
    //     // else
    //     // {
    //     //     printf("wrong pixel, no type 2 collision!\n");
    //     // }
    // }

    if (surfacePointData.collisionType < 4)
    {
        float4 normal = calculateNormal(&surfacePointData, molecule, params.debug_mode);
        surfaceColor = calculateLighting(params, surfaceColor, ray_pos, normal, cam_position);
    }

    surf2Dwrite(surfaceColor, surface, x * sizeof(float4), y);
}
void runCuda(Camera *cam, SimulationParams host_params, float4 *molecule, uint *colors)
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
    cudaDeviceSynchronize();
    marching_kernel<<<grid, block>>>(surf_object_1,
                                     molecule,
                                     colors,
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