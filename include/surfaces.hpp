/**
 * @file surfaces.hpp
 * @author Luka Blümler (luka.bluemler@uni-jena.de)
 * @brief device functions to calculate the molecule surface
 *
 */

#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#ifndef INFINITY
#define INFINITY 1e8
#endif

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <assert.h>

#include <helper_cuda.h> // helper functions for CUDA error check
#include <camera.hpp>
#include <helper_math.h> // includes vector types
#include <sorter.cuh>
#include <simulation_config.hpp>

typedef unsigned int uint;

// TODO: search better solution; confirm size threshold
// const GROUP_SIZE intended to allocate fast array memory in registers
// !!! has to be greater than the number of nearest atoms considered (k)
// !!! must not exceed available register size GROUP_SIZE*3*sizeof(int)< 1000
const int GROUP_SIZE = 50;

struct pairFloat4
{
    float4 values[2];
};

struct hitInfo
{
    uint bondId1;
    uint bondId2;
    uint bondId3;
    uint collisionType;
    float4 surfaceHit;
    float4 rayPos;
    bool grid_reached;
};

// struct voxel
// {
//     int start_cell;
//     int cell_count;
// };

__host__ __device__ float calcSignedDistanceSphere(float4 p, float4 atom)
{
    return length(p - atom) - atom.w;
}
/**
 * @brief returns the shortest distance of a point p to any atom surface
 *
 * @param p sample point p
 * @param atoms array of atoms
 * @param atomCount number of atoms
 * @return distance as float
 */
__host__ __device__ float calcSignedDistanceOuter(float4 p, Atom *atoms, int atomCount)
{
    float dist = calcSignedDistanceSphere(p, atoms[0].location);
    for (unsigned int i = 1; i < atomCount; ++i)
    {
        dist = min(dist, calcSignedDistanceSphere(p, atoms[i].location));
    }
    return dist;
}

__host__ __device__ float calcSignedDistanceOuter(float4 p, float4 *atoms, int atomCount)
{
    float dist = calcSignedDistanceSphere(p, atoms[0]);
    for (unsigned int i = 1; i < atomCount; ++i)
    {
        dist = min(dist, calcSignedDistanceSphere(p, atoms[i]));
    }
    return dist;
}
__host__ __device__ float calcSignedDistanceOuterWithNearestAtom(float4 p, float4 *atoms, int atomCount, hitInfo *surfacePointData)
{
    int nearId = 0;
    float dist = calcSignedDistanceSphere(p, atoms[0]);
    for (unsigned int i = 1; i < atomCount; i++)
    {
        float tmp = dist;
        dist = min(dist, calcSignedDistanceSphere(p, atoms[i]));
        if (tmp > dist)
            nearId = i;
    }
    surfacePointData->bondId1 = nearId;

    return dist;
}

/**
 * @brief determine the id of the max element of an array
 *
 * @param array array of floats
 * @param length length of the array
 * @return int
 */
__host__ __device__ int maxId(Atom array[], int length)
{
    int max_id = 0;
    float max_value = array[0].distance;

    for (int i = 1; i < length; i++)
    {
        if (array[i].distance > max_value)
        {
            max_value = array[i].distance;
            max_id = i;
        }
    }
    return max_id;
}

__device__ int castVoxelToId(int4 voxel_dim, int voxel_x, int voxel_y, int voxel_z)
{
    int id = voxel_x + (voxel_y * voxel_dim.x) + (voxel_z * voxel_dim.x * voxel_dim.y);

    return id;
}
/**
 * @brief calculate the intersection point between two spheres, that is closest to a sample point/position p
 *
 * @param p
 * @param atom1
 * @param atom2
 */
__device__ float4 intersectTwoSpheres(float4 p, float4 atom1, float4 atom2, float x, float y)
{
    // calculate distance to midpoint between atoms (2nd intercept theorem)
    float d = length(atom1 - atom2);
    float mLength = (d * d - atom2.w * atom2.w + atom1.w * atom1.w) / (2 * d);

    // calculate distance from midpoint to intersection point
    float h = sqrt(atom1.w * atom1.w - mLength * mLength);

    // calculate vector from atom1 to midpoint
    float4 b = normalize(atom2 - atom1);
    float4 m = b * mLength + atom1;

    // calculate unit vector from midpoint between atom centers to intersection point
    float4 a = normalize(p - atom1);
    float4 q = normalize(a - dot(a, b) * b);

    // calculate nearest intersection point
    float4 intersec = m + q * h;

    return intersec;
}

__device__ bool in_triangle(float4 p, float4 c1, float4 c2, float4 intersect)
{
    c1 -= p;
    c2 -= p;
    intersect -= p;

    float4 vecU = cross(c1, c2);
    float4 vecV = cross(c2, intersect);
    float4 vecW = cross(intersect, c1);
    if (dot(vecU, vecV) < 0)
        return false;
    if (dot(vecV, vecW) < 0)
        return false;
    else
        return true;
}

__device__ pairFloat4 calculate_case3(float4 p, float4 c1, float4 c2, float4 c3, float4 surface_hit, Atom atoms[], int atom_count, float epsilon)
{
    float4 base1 = c2 - c1;
    float4 base2 = c3 - c1;
    float b1 = length(base1);
    float b2 = length(base2);

    float4 i = normalize(base1);
    float4 j = normalize((base2 - dot(base2, i) * i));

    float r_squared1 = c1.w * c1.w;
    float r_squared2 = c2.w * c2.w;
    float r_squared3 = c3.w * c3.w;

    float x = (r_squared1 - r_squared2 + b1 * b1) / (2 * b1);
    float y = (r_squared1 - r_squared3 + b2 * b2 - 2 * dot(base2, i) * x) / (2 * dot(base2, j));
    float z = sqrt(r_squared1 - x * x - y * y);

    pairFloat4 hit;

    hit.values[0] = c1 + i * x + j * y + cross(i, j) * z;
    hit.values[1] = c1 + i * x + j * y - cross(i, j) * z;

    return hit;
}

__device__ float calculateDistanceToGrid(float4 ray_pos, float4 ray_dir, SimulationParams params, int frame, int x = 0, int y = 0)
{
    // TODO: optimize code: use only start and end point of box for the faces (makes first array unnecessary)
    // find the closest face of the grid
    float4 box_dimensions = (params.box_end - params.box_start) / 2;
    float4 grid_center = params.box_start + box_dimensions;

    float distance_to_face = INFINITY;
    float4 face_midpoint, face_normal;

    // calculate midpoint, normal and midpoint distance for every face of the box
    // continue with data of closest face
    for (int k = -1; k < 2; k = k + 2)
    {
        for (int i = 0; i < 3; i++)
        {
            float4 normal = makeFloat4Component((float)k * getFloat4Component(&box_dimensions, i), i);
            float4 midpoint = grid_center + normal;

            float distance = length(ray_pos - midpoint);
            if (distance < distance_to_face)
            {
                distance_to_face = distance;
                face_midpoint = midpoint;
                face_normal = normal;
            }
        }
    }

    // calculate distance to face along ray_dir
    float numerator = dot((face_midpoint - ray_pos), face_normal);
    float denominator = dot(normalize(ray_dir), face_normal);

    float traveldist = 2 * params.solvent_max; // minimum distance to be traveled in undefined cases
    if (denominator == 0)
    {
        if (numerator == 0)
            return params.epsilon; // point already lies on the face; set epsilon as distance to avoid infinite loops
        else
            return traveldist; // moving parallel to the box -> ray can't hit the molecule
    }
    float result = numerator / denominator;

    // distance to grid is measured along ray direction
    // -> negative distances implicate the box has already been traversed
    if (result < 0)
    {
        return traveldist;
    }
    return result;
}

__device__ void findNearestAtoms(float4 *molecule_data, uint number_atoms, Atom nearest_atoms[], uint nearest_count, float4 position, int x = 0, int y = 0)
{
    for (uint i = 0; i < nearest_count; i++)
    {
        float dist = calcSignedDistanceSphere(position, molecule_data[i]);

        nearest_atoms[i].location = molecule_data[i];
        nearest_atoms[i].id = i;
        nearest_atoms[i].distance = dist;
    }
    for (uint i = nearest_count; i < number_atoms; i++)
    {
        float dist = calcSignedDistanceSphere(position, molecule_data[i]);

        int max_element = findFurthestByIndex(nearest_atoms, nearest_count);
        if (nearest_atoms[max_element].distance > dist)
        {
            nearest_atoms[max_element].location = molecule_data[i];
            nearest_atoms[max_element].id = i;
            nearest_atoms[max_element].distance = dist;
        }
    }
    if (nearest_count > 1)
    {
        selectionSort(nearest_count, nearest_atoms);
    }
}
__device__ unsigned int nearestOverVoxel(float4 *molecule_data, int *voxel_data, int *voxel_count, Atom nearest_atoms[], float4 position, float4 direction, SimulationParams params, int x = 0, int y = 0)
{
    unsigned int nearest_count = 0;

    // calculate current voxel
    int4 current_voxel = castf2i(floorf((position - params.box_start) / params.voxel_size));

    // TODO: calculate nearest atoms with voxels
    // determine closest points for 3x3x3 voxels
    int4 range_min = max(current_voxel - make_int4(1, 1, 1, 0), make_int4(0, 0, 0, 0));
    int4 range_max = min(current_voxel + make_int4(1, 1, 1, 0), params.voxel_dim - make_int4(1, 1, 1, 0));

    bool k_full = false;

    for (int i = range_min.x; i <= range_max.x; i++)
    {
        for (int j = range_min.y; j <= range_max.y; j++)
        {
            for (int k = range_min.z; k <= range_max.z; k++)
            {
                // 1. calculate voxel id
                int id = castVoxelToId(params.voxel_dim, i, j, k);

                int voxel_start = id * params.atomsPerVoxel;
                int voxel_end = (id + 1) * params.atomsPerVoxel;

                for (int v = voxel_start; v < voxel_end; v++)
                {
                    int m = voxel_data[v];
                    if (m >= 0)
                    {
                        float dist = calcSignedDistanceSphere(position, molecule_data[m]);

                        if (!k_full)
                        {
                            // 3. copy all elements to nearest_atoms
                            nearest_atoms[nearest_count].location = molecule_data[m];
                            nearest_atoms[nearest_count].id = m;
                            nearest_atoms[nearest_count].distance = dist;
                            nearest_count++;

                            if (nearest_count >= params.k_nearest)
                            {
                                k_full = true;
                            }
                        }
                        else
                        {
                            int max_element = findFurthestByIndex(nearest_atoms, nearest_count);
                            if (nearest_atoms[max_element].distance > dist)
                            {
                                nearest_atoms[max_element].location = molecule_data[m];
                                nearest_atoms[max_element].id = m;
                                nearest_atoms[max_element].distance = dist;
                            }
                        }
                    }
                }
            }
        }
    }

    // swap the closest atom to be at the start of the array(only applicable with multiple atoms in neighbourhood)
    if (nearest_count > 1)
    {
        selectionSort(nearest_count, nearest_atoms);
    }
    return nearest_count;
}
__device__ float computeSurface(float4 ray_pos, float4 ray_dir, float4 *molecule, uint *colors, int *voxel_data, int *voxel_count, SimulationParams params, hitInfo *surfacePointData, int x = 0, int y = 0, int frame = 0, int step = 0)
{

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // // 0 // variables/constants
    float f_sdf = 0; // current distance to SAS

    surfacePointData->collisionType = 0;
    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // // 1 // determine signed distance to closest atom

    // if (params.debug_mode && params.debug_frame == 0 && x == params.mouse_x_pos && y == params.mouse_y_pos)
    // {
    //     printf("DEBUG|%i|%i|%i|%i|%.3f|%.3f|%.3f\n", x, y, params.debug_frame, step, ray_pos.x, ray_pos.y, ray_pos.z);
    //     // printf("DEBUG: pixel %i, %i at frame %i, step %i -> pos %.2f, %.2f, %.2f\n", x, y, params.debug_frame, step, ray_pos.x, ray_pos.y, ray_pos.z);
    // }
    if (params.use_voxel)
    {
        // position outside grid (position cannot lie on molecule surface)
        if (ray_pos < (params.box_start) || ray_pos > (params.box_end))
        {
            if (surfacePointData->grid_reached)
                return INFINITY;

            f_sdf = calculateDistanceToGrid(ray_pos, ray_dir, params, frame, x, y);
            // TODO: test if a constant travel distance (at most 1x solvent radius) would be faster than calculating the distance to the grid
            // f_sdf = params.solvent_radius;

            // if (params.debug_mode && params.debug_frame == 0 && x == 728 && y == 423)
            // {
            //     printf("DEBUG: OUTSIDE GRID (%.2f) pixel %i, %i, ray_pos: %.2f, %.2f, %.2f, debugframe %i, step: %i\n",
            //            f_sdf, x, y, ray_pos.x, ray_pos.y, ray_pos.z, params.debug_frame, step);
            // }

            return f_sdf;
        }
        surfacePointData->grid_reached = true;
    }
    if (!params.use_voxel)
    {

        f_sdf = calcSignedDistanceOuterWithNearestAtom(ray_pos, molecule, params.numAtoms, surfacePointData);
        if (f_sdf >= 0)
        {
            // if (params.debug_mode && params.debug_frame == 0 && x == 728 && y == 423)
            // {
            //     printf("DEBUG: OUTSIDE GRID (%.2f) pixel %i, %i, ray_pos: %.2f, %.2f, %.2f, debugframe %i, step: %i\n",
            //            f_sdf, x, y, ray_pos.x, ray_pos.y, ray_pos.z, params.debug_frame, step);
            // }
            return f_sdf;
        }
    }

    Atom nearest_atoms[GROUP_SIZE];
    unsigned int nearest_count = 0;

    if (params.use_voxel)
    {
        // determine neighbourhood of atoms over 3 x 3 x 3 voxels
        nearest_count = nearestOverVoxel(molecule, voxel_data, voxel_count, nearest_atoms, ray_pos, ray_dir, params, x, y);

        if (nearest_count == 0)
        {
            f_sdf = 2 * params.solvent_max;
            return f_sdf;
        }

        f_sdf = calcSignedDistanceOuter(ray_pos, nearest_atoms, nearest_count);
        if (f_sdf >= 0)
        {
            return f_sdf;
        }
    }

    if (f_sdf < 0)
    {
        // if (params.debug_mode && params.debug_frame == 0 && x == params.mouse_x_pos && y == params.mouse_y_pos)
        // if (params.debug_mode && params.debug_frame == 0 && x == 728 && y == 423)
        // {
        //     printf("DEBUG: INSIDE GRID (%.2f) pixel %i, %i, ray_pos: %.2f, %.2f, %.2f, debugframe %i, step: %i\n",
        //            f_sdf, x, y, ray_pos.x, ray_pos.y, ray_pos.z, params.debug_frame, step);
        // }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 2 // determine k closest points
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // allocate array for ids of k atoms closest to p
        // TODO: add checks in settings, that k_nearest is smaller or equal to numAtoms
        // find the k nearest atoms

        if (!params.use_voxel)
        {
            // nearest_count = min(params.k_nearest, params.numAtoms);
            nearest_count = (uint)params.k_nearest;
            findNearestAtoms(molecule, params.numAtoms, nearest_atoms, nearest_count, ray_pos);
        }

        // if (params.debug_mode && params.debug_frame == 0 && x == 693 && y == 546)
        // {
        //     printf("DEBUG: pixel %i, %i at frame %i, step %i -> pos %.2f, %.2f, %.2f\n", x, y, params.debug_frame, step, ray_pos.x, ray_pos.y, ray_pos.z);
        // }

        // if (params.debug_mode && params.debug_frame == 0 && step == 5 && x == 728 && y == 423)
        // {
        //     printf("-----------------------------------------------\n");
        //     printf("DEBUG: INSIDE GRID NEAREST ATOMS:\n");
        //     printf("rank|id|pos.x|pos.y|pos.z|radius|distance\n");
        //     for (int i = 0; i < nearest_count; i++)
        //     {
        //         printf("%i|%i|%.3f|%.3f|%.3f|%.3f|%.3f\n",
        //                i, nearest_atoms[i].id, nearest_atoms[i].location.x, nearest_atoms[i].location.y, nearest_atoms[i].location.z,
        //                nearest_atoms[i].location.w - params.solvent_radius, nearest_atoms[i].distance);
        //     }
        //     printf("-----------------------------------------------\n");
        // }
        // TODO: DEBUG appproach I:
        //       ->write Debug method to activate in range of specific raypos (custom delta)
        //       -> get all k-nearest atoms
        //       -> compare grid and non grid version
        if (params.debug_mode && params.debug_frame == 0 && x == 1062 && y == 396)
        {
            float epsilon_pos = 0.01;
            float4 target = make_float4(-14.112, -10.838, -6.265, 0);
            if ((ray_pos - target) < make_float4(epsilon_pos))
            {
                for (int i = 0; i < nearest_count; i++)
                {
                    if (params.use_voxel)
                    {
                        int4 current_voxel = castf2i(floorf((ray_pos - params.box_start) / params.voxel_size));
                        int voxel = castVoxelToId(params.voxel_dim, current_voxel.x, current_voxel.y, current_voxel.z);
                        printf("%i|%i|%i|%.2f|%.2f|%.2f|%i|%i|%.3f\n",
                               i, (int)true, step, ray_pos.x, ray_pos.y, ray_pos.z,
                               nearest_atoms[i].id, voxel, nearest_atoms[i].distance);
                    }
                    else
                    {
                        printf("%i|%i|%i|%.2f|%.2f|%.2f|%i|%i|%.3f\n",
                               i, (int)false, step, ray_pos.x, ray_pos.y, ray_pos.z,
                               nearest_atoms[i].id, 0, nearest_atoms[i].distance);
                    }
                }
            }
        }

        // TODO: DEBUG appproach II:
        //       ->write Debug method to activate in range of specific raypos (custom delta)
        //       -> calculate point for each voxel
        //       -> get all atoms for each voxel & compare with no grid version (k-nearest)
        //       -> compare grid and non grid version

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 3 // intersection with one atom

        // group of participating atoms is ordered, nearest_atoms[0] contains the id of the closest atom

        float4 dir = normalize(ray_pos - nearest_atoms[0].location);
        float4 surface_hit = nearest_atoms[0].location + nearest_atoms[0].location.w * dir;

        surfacePointData->collisionType = 1;
        surfacePointData->bondId1 = nearest_atoms[0].id;

        // if (params.debug_mode && params.debug_frame == 0 && x == params.mouse_x_pos && y == params.mouse_y_pos)
        // {
        //     printf("------DEBUG case 1 (nearest atoms)-------\n");
        //     printf("raypos = %.2f, %.2f, %.2f\n", ray_pos.x, ray_pos.y, ray_pos.z);
        //     printf("rank|id|x|y|z|radius|distance\n");
        //     for (int i = 0; i < nearest_count; i++)
        //     {
        //         printf("%i|%i|%.2f|%.2f|%.2f|%.2f|%.2f\n",
        //                i, nearest_atoms[i].id,
        //                nearest_atoms[i].location.x, nearest_atoms[i].location.y, nearest_atoms[i].location.z,
        //                nearest_atoms[i].location.w - params.solvent_radius, nearest_atoms[i].distance);
        //     }

        //     printf("---------------------------------------\n");
        // }
        if (abs(calcSignedDistanceOuter(surface_hit, nearest_atoms, params.k_nearest)) > params.epsilon)
        {
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // 4 // intersection with two atoms (arc)

            bool surface_found = false;
            if (params.use_iterative_solver)
            {
                /* code */
            }
            else
            {
                // calculate all intersection points for the group of nearest atoms
                for (unsigned int i = 0; i < params.k_nearest - 1; i++)
                {
                    for (unsigned int j = i + 1; j < params.k_nearest; j++)
                    {
                        surface_hit = intersectTwoSpheres(ray_pos, nearest_atoms[i].location, nearest_atoms[j].location, x, y);

                        // check if points lie in correct section (arc) of the surface
                        bool triangle = in_triangle(ray_pos, nearest_atoms[i].location, nearest_atoms[j].location, surface_hit);
                        // check if point lies on the surface
                        bool surface = (abs(calcSignedDistanceOuter(surface_hit, nearest_atoms, params.k_nearest)) < params.epsilon);
                        if (triangle && surface)
                        {
                            surface_found = true;
                            surfacePointData->collisionType = 2;
                            surfacePointData->bondId1 = nearest_atoms[j].id;
                            surfacePointData->bondId2 = nearest_atoms[i].id;
                            surfacePointData->surfaceHit = surface_hit;

                            break;
                        }
                    }
                    if (surface_found)
                    {
                        break;
                    }
                }
            }

            //////////////////////////////////////////////////////////////////////////////////////////////////////

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // 5 // intersection with three atoms (spherical triangle)
            if (!surface_found)
            {
                if (params.use_iterative_solver)
                {
                    /* code */
                }
                else
                {
                    bool first = true;
                    for (unsigned int i = 0; i < params.k_nearest - 2; i++)
                    {
                        for (unsigned int j = i + 1; j < params.k_nearest - 1; j++)
                        {
                            for (unsigned int k = j + 1; k < params.k_nearest; k++)
                            {
                                pairFloat4 s_pot = calculate_case3(ray_pos, nearest_atoms[i].location, nearest_atoms[j].location, nearest_atoms[k].location, surface_hit, nearest_atoms, params.k_nearest, params.epsilon);

                                // every case 3 has two potential surface points. loop over both and determine it's distance to the current position.
                                // -> bool first ensures that both distances are tested, if both points of a potential pair lie on the surface
                                // -> the 4th component of the final surface point stores the distance of the previous calculated point
                                for (unsigned int n = 0; n < 2; n++)
                                {
                                    if (abs(calcSignedDistanceOuter(s_pot.values[n], nearest_atoms, params.k_nearest)) < params.epsilon)
                                    {

                                        surface_found = true;
                                        float dist = length(ray_pos - s_pot.values[n]);
                                        if (first)
                                        {
                                            surfacePointData->bondId1 = nearest_atoms[i].id;
                                            surfacePointData->bondId2 = nearest_atoms[j].id;
                                            surfacePointData->bondId3 = nearest_atoms[k].id;

                                            surface_hit = s_pot.values[n];
                                            surface_hit.w = length(ray_pos - s_pot.values[n]);
                                            first = false;
                                        }
                                        else if (dist < surface_hit.w)
                                        {
                                            surfacePointData->bondId1 = nearest_atoms[i].id;
                                            surfacePointData->bondId2 = nearest_atoms[j].id;
                                            surfacePointData->bondId3 = nearest_atoms[k].id;

                                            surface_hit = s_pot.values[n];
                                            surface_hit.w = length(ray_pos - s_pot.values[n]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (surface_found)
                    {
                        surfacePointData->surfaceHit = surface_hit;
                        surfacePointData->collisionType = 3;
                        // surfacePointData->bondId1 = nearest_atoms[0].id;
                        // surfacePointData->bondId2 = nearest_atoms[1].id;
                        // surfacePointData->bondId3 = nearest_atoms[2].id;
                    }
                    else
                    {
                        surfacePointData->collisionType = 4;
                        if (params.debug_mode && (x == (int)params.mouse_x_pos) && (y == (int)params.mouse_y_pos))
                        {
                            printf("no viable case3 point found\n");
                        }
                    }
                }
            }
            f_sdf = -1 * length(surface_hit - ray_pos);
        }
    }
    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    return f_sdf;
}
#endif