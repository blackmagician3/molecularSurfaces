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
#include <camera.h>
#include <helper_math.h> // includes vector types
#include <sorter.cuh>
#include <simulation_config.hpp>

// TODO: search better solution; confirm size threshold
// const GROUP_SIZE intended to allocate fast array memory in registers
// !!! has to be greater than the number of nearest atoms considered (k)
// !!! must not exceed available register size GROUP_SIZE*3*sizeof(int)< 1000
const int GROUP_SIZE = 20;

struct pairFloat4
{
    float4 values[2];
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
__host__ __device__ float calcSignedDistanceOuterWithNearestAtom(float4 p, float4 *atoms, int atomCount, float4 *surface_orientation)
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
    *surface_orientation = atoms[nearId];

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
    float max_value = array[0].id;

    for (int i = 1; i < length; i++)
    {
        if (array[i].id > max_value)
        {
            max_value = array[i].id;
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
 * @param nearest_count
 * @param position
 */
__host__ __device__ void findNearestAtoms(float4 *molecule_data, uint number_atoms, Atom nearest_atoms[], int nearest_count, float4 position, int x = 0, int y = 0)
{
    for (int i = 0; i < nearest_count; i++)
    {
        float dist = calcSignedDistanceSphere(position, molecule_data[i]);

        nearest_atoms[i].id = i;
        nearest_atoms[i].location.w = dist;
    }
    for (int i = nearest_count; i < number_atoms; i++)
    {
        int id = maxId(nearest_atoms, nearest_count);
        float dist = calcSignedDistanceSphere(position, molecule_data[i]);

        if ((dist < nearest_atoms[id].location.w))
        {
            nearest_atoms[id].location.w = dist;
            nearest_atoms[id].id = i;
        }
    }
    selectionSort(nearest_count, nearest_atoms);

    for (int i = 0; i < nearest_count; i++)
    {
        nearest_atoms[i].location = molecule_data[nearest_atoms[i].id];
    }
}

/**
 * @brief asserts, that point p lies in the extended radius of both atoms
 *
 * @param p
 * @param atom1
 * @param atom2
 */
__device__ bool test_predicate2(float4 p, float4 atom1, float4 atom2)
{
    bool predicate = (length(p - atom1) > atom1.w) && (length(p - atom2) > atom2.w);
    return predicate;
}

/**
 * @brief calculate the intersection point between two spheres, that is closest to a sample point/position p
 *
 * @param p
 * @param atom1
 * @param atom2
 */
__device__ float4 intersectTwoSpheres(float4 p, float4 atom1, float4 atom2)
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

/**
 * @brief assert, that the calculated intersection point for case 2 (sphere-sphere intersection) does...
 * a) lie within the surface spanned by the atom centers and the current position
 * b) lie on the surface
 * @param p current position/sample point
 * @param c1 first atom
 * @param c2 second atom
 * @param intersect intersection point
 * @param atoms array of atoms
 * @param epsilon
 */
__device__ bool in_triangle(float4 p, float4 c1, float4 c2, float4 intersect)
{
    float4 v1 = c1 - p;
    float4 v2 = c2 - p;
    float4 v3 = intersect - p;

    float4 vecU = cross(v1, v2);
    float4 vecV = cross(v2, v3);
    float4 vecW = cross(v3, v1);

    // bool triangle = (dot(vecU, vecV) >= 0) && (dot(vecV, vecW) >= 0);
    // return triangle;
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

__device__ float computeSurface(float4 ray_pos, float4 *molecule, SimulationParams params, bool *hit, float4 *colour, float4 *surface_orientation, int *case_id, int x = 0, int y = 0)
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 0 // variables/constants
    float f_sdf = 0;                                      // current distance to SAS
    float4 c_case0 = make_float4(0.0f, 0.0f, 1.0f, 1.0f); // color for case 0 (SAS)
    float4 c_case1 = make_float4(0.0f, 0.0f, 1.0f, 1.0f); // color for case 1
    float4 c_case2 = make_float4(0.0f, 1.0f, 0.0f, 1.0f); // color for case 2
    float4 c_case3 = make_float4(1.0f, 0.0f, 0.0f, 1.0f); // color for case 3
    float4 c_error = make_float4(1.0f, 0.0f, 1.0f, 1.0f); // color for debugging

    *case_id = 0;
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1 // determine signed distance to closest atom

    f_sdf = calcSignedDistanceOuterWithNearestAtom(ray_pos, molecule, params.numAtoms, surface_orientation);

    *colour = c_case0;
    // exit calculation if position of sample point p outside of SAS
    if (f_sdf >= 0)
        return f_sdf;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2 // determine k closest points

    // allocate array for ids of k atoms closest to p
    unsigned int nearest_count = min(params.k_nearest, params.numAtoms);

    // find the k nearest atoms
    Atom nearest_atoms[GROUP_SIZE];
    findNearestAtoms(molecule, params.numAtoms, nearest_atoms, nearest_count, ray_pos);

    // if (params.debug_mode && (x == (int)params.mouse_x_pos) && (y == (int)params.mouse_y_pos))
    // {

    //     printf("nearest atom: atom %d\n", nearest_atoms[0].id);
    //     printf("second nearest atom: atom %d\n", nearest_atoms[1].id);
    // }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 3 // intersection with one atom

    // TODO case 0 - when using voxels for nearest neighbours, there might be zero neighbours
    // group of participating atoms is ordered, nearest_atoms[0] contains the id of the closest atom

    float4 dir = normalize(ray_pos - nearest_atoms[0].location);
    float4 surface_hit = nearest_atoms[0].location + nearest_atoms[0].location.w * dir;

    *colour = c_case1;
    *case_id = 1;

    if (abs(calcSignedDistanceOuter(surface_hit, nearest_atoms, params.numAtoms)) <= params.epsilon)
    {
        return f_sdf;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 4 // intersection with two atoms (arc)

    bool surface_found = false;
    // calculate all intersection points for the group of nearest atoms
    for (unsigned int i = 0; i < nearest_count - 1; i++)
    {
        for (unsigned int j = i + 1; j < nearest_count; j++)
        {
            surface_hit = intersectTwoSpheres(ray_pos, nearest_atoms[i].location, nearest_atoms[j].location);

            // check if points lie in correct section (arc) of the surface
            bool triangle = in_triangle(ray_pos, nearest_atoms[i].location, nearest_atoms[j].location, surface_hit);
            // check if point lies on the surface
            bool surface = (abs(calcSignedDistanceOuter(surface_hit, nearest_atoms, nearest_count)) < params.epsilon);
            if (triangle && surface)
            {
                surface_found = true;
                *colour = c_case2;
                *surface_orientation = surface_hit;
                *case_id = 2;
                break;
            }
        }
        if (surface_found)
        {
            *surface_orientation = surface_hit;
            break;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 5 // intersection with three atoms (spherical triangle)
    if (!surface_found)
    {
        bool first = true;
        // float surface_dist = INFINITY;
        for (unsigned int i = 0; i < nearest_count - 2; i++)
        {
            for (unsigned int j = i + 1; j < nearest_count - 1; j++)
            {
                for (unsigned int k = j + 1; k < nearest_count; k++)
                {
                    pairFloat4 s_pot = calculate_case3(ray_pos, nearest_atoms[i].location, nearest_atoms[j].location, nearest_atoms[k].location, surface_hit, nearest_atoms, nearest_count, params.epsilon);

                    // every case 3 has two potential surface points. loop over both and determine it's distance to the current position.
                    // -> bool first ensures that both distances are tested, if both points of a potential pair lies on the surface
                    // -> the 4th component of the final surface point stores the distance of the previous calculated point to the position
                    for (unsigned int n = 0; n < 2; n++)
                    {
                        if (abs(calcSignedDistanceOuter(s_pot.values[n], nearest_atoms, nearest_count)) < params.epsilon)
                        {
                            surface_found = true;
                            float dist = length(ray_pos - s_pot.values[n]);
                            if (first)
                            {
                                surface_hit = s_pot.values[n];
                                surface_hit.w = length(ray_pos - s_pot.values[n]);
                                first = false;
                            }
                            else if (dist < surface_hit.w)
                            {
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
            *surface_orientation = surface_hit;
            *colour = c_case3;
            *case_id = 3;
            if (params.debug_mode && (x == (int)params.mouse_x_pos) && (y == (int)params.mouse_y_pos))
            {
                printf("case3 found at (%f, %f, %f) for pos (%f, %f, %f)\n", surface_hit.x, surface_hit.y, surface_hit.z, ray_pos.x, ray_pos.y, ray_pos.z);
            }
        }
        else
        {
            *colour = c_error;
            *case_id = 4;
            if (params.debug_mode && (x == (int)params.mouse_x_pos) && (y == (int)params.mouse_y_pos))
            {
                printf("no viable case3 point found\n");
            }
        }
    }

    f_sdf = -1 * length(surface_hit - ray_pos);

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    return f_sdf;
}
#endif