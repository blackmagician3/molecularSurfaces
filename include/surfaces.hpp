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
const int GROUP_SIZE = 20;

struct pairFloat4
{
    float4 values[2];
};

// TODO: change wording case_id -> collisionType
struct hitInfo
{
    uint bondId1;
    uint bondId2;
    uint bondId3;
    uint collisionType;
    float4 surfaceHit;
    float4 rayPos;
};

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

/**
 * @brief takes an array of atoms and calculates their distance to the current position. The k first atoms (smallest distance) are written back to nearest_atoms
 *
 * @param molecule_data pointer to input array containing all atoms
 * @param number_atoms
 * @param nearest_atoms
 * @param nearest_ids
 * @param params.k_nearest
 * @param position
 */
__host__ __device__ void findNearestAtoms(float4 *molecule_data, uint number_atoms, Atom nearest_atoms[], uint nearest_count, float4 position, int x = 0, int y = 0)
{
    for (uint i = 0; i < nearest_count; i++)
    {
        float dist = calcSignedDistanceSphere(position, molecule_data[i]);

        nearest_atoms[i].id = i;
        nearest_atoms[i].distance = dist;
    }
    for (uint i = nearest_count; i < number_atoms; i++)
    {
        selectionSort(nearest_count, nearest_atoms);
        float dist = calcSignedDistanceSphere(position, molecule_data[i]);

        if ((dist < nearest_atoms[nearest_count - 1].distance))
        {
            nearest_atoms[nearest_count - 1].distance = dist;
            nearest_atoms[nearest_count - 1].id = i;
        }
    }
    selectionSort(nearest_count, nearest_atoms);

    for (uint i = 0; i < nearest_count; i++)
    {
        nearest_atoms[i].location = molecule_data[nearest_atoms[i].id];
    }
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

__device__ float computeSurface(float4 ray_pos, float4 *molecule, uint *colors, SimulationParams params, hitInfo *surfacePointData, int x = 0, int y = 0)
{

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // // 0 // variables/constants
    float f_sdf = 0; // current distance to SAS

    surfacePointData->collisionType = 0;
    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // // 1 // determine signed distance to closest atom

    f_sdf = calcSignedDistanceOuterWithNearestAtom(ray_pos, molecule, params.numAtoms, surfacePointData);

    // // if (x == 0 && y > 826 && y < 831)
    // // {
    // //     printf("test at y = %i-a\n", y);
    // // }

    // if (x == 0 && y == 827)
    // {
    //     for (int i = 0; i < params.numAtoms; i++)
    //     {
    //         printf("atom %i is at: %.3f, %.3f, %.3f with a radius of %.3f\n", i, molecule[i].x, molecule[i].y, molecule[i].z, molecule[i].w);
    //     }
    // }

    if (f_sdf < 0)
    {
        //     // if (x == 0 && y > 826 && y < 831)
        //     // {
        //     //     printf("test at y = %i-b\n", y);
        //     // }
        //     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        //     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //     // 2 // determine k closest points

        //     // allocate array for ids of k atoms closest to p
        // TODO: add checks in settings, that k_nearest is smaller or equal to numAtoms
        // unsigned int nearest_count = min(params.k_nearest, params.numAtoms);
        unsigned int nearest_count = (uint)params.k_nearest;
        //     // find the k nearest atoms
        Atom nearest_atoms[GROUP_SIZE];
        //     // if (x == 0 && y > 826 && y < 831)
        //     // {
        //     //     printf("test at y = %i-c\n", y);
        //     // }
        findNearestAtoms(molecule, params.numAtoms, nearest_atoms, nearest_count, ray_pos);

        //     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //     // 3 // intersection with one atom

        //     // TODO case 0 - when using voxels for nearest neighbours, there might be zero neighbours
        //     // group of participating atoms is ordered, nearest_atoms[0] contains the id of the closest atom

        float4 dir = normalize(ray_pos - nearest_atoms[0].location);
        float4 surface_hit = nearest_atoms[0].location + nearest_atoms[0].location.w * dir;

        surfacePointData->collisionType = 1;
        //     // if (params.debug_mode && (x == (int)params.mouse_x_pos) && (y == (int)params.mouse_y_pos))
        //     // {
        //     //     printf("nearest atom radius:  %.3f vs solvent radius: %.3f\n", nearest_atoms[0].location.w, params.solvent_radius);
        //     // }

        // TODO: check validity: signed distance calculation limited to k nearest atoms
        if (abs(calcSignedDistanceOuter(surface_hit, nearest_atoms, params.k_nearest)) > params.epsilon)
        {
            //         ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

            //         ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //         // 4 // intersection with two atoms (arc)

            bool surface_found = false;
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

            ////////////////////////////////////////////////////////////////////////////////////////////////////////

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // 5 // intersection with three atoms (spherical triangle)
            if (!surface_found)
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
            f_sdf = -1 * length(surface_hit - ray_pos);
        }
    }
    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    return f_sdf;
}
#endif