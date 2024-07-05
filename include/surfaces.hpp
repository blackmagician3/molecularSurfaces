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

#ifndef ZERO
#define ZERO 1e-8
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
#include <solver.hpp>

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
    bool isInGrid;
    bool traversedGrid;
};

/**
 * @brief collection of local pointers and variables
 *
 */
struct functionArgs
{
    const float4 ray_pos;       // copy of ray position
    const float4 ray_dir;       // copy of ray direction
    SimulationParams *p_params; // device parameters that can be configured from the UI
    float4 *atom_data_global;   // points to the complete array of atoms
    int *voxel_data;            // points to array of sorted atoms based on their respective voxel in the grid
    Atom *atom_data_local;      // points to a shortend array of close atoms
    hitInfo *hit_feedback;      // contains information of surface hit
    float shortest_distance;

    // debug information
    const int x, y;  // current pixel
    const int frame; // currently generated frame
    int step;        // can be used as counter for loops
    int calcs1;      // can be used as counter for loops
    int calcs2;      // can be used as counter for loops
    int calcs3;      // can be used as counter for loops
};

__host__ __device__ float calcSignedDistanceSphere(float4 p, float4 atom)
{
    return length(p - atom) - atom.w;
}

// __host__ __device__ float calcSignedDistanceOuter(float4 p, Atom *atoms, int atomCount)
// {
//     float dist = calcSignedDistanceSphere(p, atoms[0].location);
//     for (unsigned int i = 1; i < atomCount; ++i)
//     {
//         dist = min(dist, calcSignedDistanceSphere(p, atoms[i].location));
//     }
//     return dist;
// }

__host__ __device__ float localSDouter(float4 p, const functionArgs *args)
{

    float dist = calcSignedDistanceSphere(p, args->atom_data_local[0].location);
    for (unsigned int i = 1; i < args->p_params->k_nearest; ++i)
    {
        dist = min(dist, calcSignedDistanceSphere(p, args->atom_data_local[i].location));
    }
    return dist;
}

// __host__ __device__ float calcSignedDistanceOuter(float4 p, float4 *atoms, int atomCount)
// {
//     float dist = calcSignedDistanceSphere(p, atoms[0]);
//     for (unsigned int i = 1; i < atomCount; ++i)
//     {
//         dist = min(dist, calcSignedDistanceSphere(p, atoms[i]));
//     }
//     return dist;
// }
__host__ __device__ float globalSDouter(float4 p, const functionArgs *args)
{

    float dist = calcSignedDistanceSphere(p, args->atom_data_global[0]);
    for (unsigned int i = 1; i < args->p_params->numAtoms; ++i)
    {
        dist = min(dist, calcSignedDistanceSphere(p, args->atom_data_global[i]));
    }
    return dist;
}

// __host__ __device__ float calcSignedDistanceOuterWithNearestAtom(float4 p, float4 *atoms, int atomCount, hitInfo *surfacePointData)
// {
//     int nearId = 0;
//     float dist = calcSignedDistanceSphere(p, atoms[0]);
//     for (unsigned int i = 1; i < atomCount; i++)
//     {
//         float tmp = dist;
//         dist = min(dist, calcSignedDistanceSphere(p, atoms[i]));
//         if (tmp > dist)
//             nearId = i;
//     }
//     surfacePointData->bondId1 = nearId;

//     return dist;
// }

__host__ __device__ float globalSDouterWithFeedback(float4 p, const functionArgs *args)
{
    int nearId = 0;
    float dist = calcSignedDistanceSphere(p, args->atom_data_global[0]);
    for (unsigned int i = 1; i < args->p_params->numAtoms; i++)
    {
        float tmp = dist;
        dist = min(dist, calcSignedDistanceSphere(p, args->atom_data_global[i]));
        if (tmp > dist)
            nearId = i;
    }

    args->hit_feedback->bondId1 = nearId;

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
 * @brief calculate the intersection point between two spheres (analytically), that is closest to a sample point/position p
 *
 */
__device__ float4 intersectTwoSpheres1(float4 atom1, float4 atom2, functionArgs *args)
{
    float4 p = args->ray_pos;
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
 * @brief calculate the intersection point between two spheres (newton iterations), that is closest to a sample point/position p
 *
 * @param p
 * @param atom1
 * @param atom2
 * @param x // for debugging purposes
 * @param y // for debugging purposes
 */
__device__ float4 intersectTwoSpheres2a(float4 atom1, float4 atom2, functionArgs *args)
{
    float4 p = args->ray_pos;
    float delta = INFINITY;
    float4 intersec_prev = make_float4(p.x, p.y, p.z, 0);
    float4 intersec = intersec_prev;
    int solver_iter = 0; // for debugging
    int max_iter_newton = 20;

    // // derivation from paper
    // while (!near_enough)
    // {
    //     float3 v = make_float3(atom1.w - length(intersec_prev - atom1), atom2.w - length(intersec_prev - atom2), 0);
    //     if (abs(v.x) < p_params->epsilon && abs(v.y) < p_params->epsilon)
    //     {
    //         near_enough = true;
    //         break;
    //     }

    //     float4 a = -(atom1 - intersec_prev) / (length(atom1 - intersec_prev));
    //     float4 b = -(atom2 - intersec_prev) / (length(atom2 - intersec_prev));
    //     float4 c = cross(a, b);

    //     // float m_det = dot(c, c);
    //     float m_det = a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) - a.z * (b.x * c.y - b.y * c.x);
    //     if (abs(m_det) < ZERO)
    //         m_det = ZERO;
    //     float4 m1 = make_float4(b.y * c.z - b.z * c.y, a.z * c.y - a.y * c.z, a.y * b.z - a.z * b.y, 0);
    //     float4 m2 = make_float4(b.z * c.x - b.x * c.z, a.x * c.z - a.z * c.x, a.z * b.x - a.x * b.z, 0);
    //     float4 m3 = make_float4(b.x * c.y - b.y * c.x, a.y * c.x - a.x * c.y, a.x * b.y - a.y * b.x, 0); // last column of matrix not needed, as v.z = 0

    //     // float4 inter_change = make_float4(v.x * m1.x + v.y * m2.x, v.x * m1.y + v.y * m2.y, v.x * m1.z + v.y * m2.z, 0) / m_det;  // no z-components, as v.z = 0
    //     // float4 inter_change = make_float4(v.x * m1.x + v.y * m2.x + v.z * m3.x, v.x * m1.y + v.y * m2.y + v.z * m3.y, v.x * m1.z + v.y * m2.z + v.z * m3.z, 0) / m_det;

    //     float4 inter_change = make_float4(v.x * m1.x + v.y * m1.y + v.z * m1.z, v.x * m2.x + v.y * m2.y + v.z * m2.z, v.x * m3.x + v.y * m3.y + v.z * m3.z, 0) / m_det;
    //     // float4 inter_change = make_float4(v.x * m1.x + v.y * m1.y + v.z * m1.z, v.x * m2.x + v.y * m2.y + v.z * m2.z, 0, 0) / m_det;
    //     // float4 inter_change = make_float4(v.x * m1.x + v.y * m1.y + v.z * m1.z, v.x * m2.x + v.y * m2.y + v.z * m2.z, v.x * m3.x + v.y * m3.y + v.z * m3.z, 0) / m_det;

    //     intersec = intersec_prev + inter_change;

    //     delta = length(intersec - intersec_prev);

    //     // for debugging
    //     solver_iter++;

    //     if (frame == 200 && x == 982 && y == 653)
    //     {
    //         // x|y|frame|calculation|step|delta|prev.x|prev.y|prev.z|next.x|next.y|next.z|det|change.x|change.y|change.z
    //         printf("%i|%i|%i|%i|%i|%i|%.8f|%.3f|%.3f|%.3f|%.3f|%.3f|%.3f|%.5f|%.3f|%.3f|%.3f|%.3f|%.3f\n", x, y, frame, step, calculations, solver_iter - 1, delta, intersec_prev.x, intersec_prev.y, intersec_prev.z,
    //                intersec.x, intersec.y, intersec.z, m_det, inter_change.x, inter_change.y, inter_change.z, v.x, v.y);
    //         // printf("DEBUG: step: %i, delta: %.3f, x: %.3f, %.3f, %.3f, det: %.3f\n", solver_iter - 1, delta, inter_change.x, inter_change.y, inter_change.z, m_det);
    //     }
    //     intersec_prev = intersec;

    //     if (solver_iter > max_iter_newton)
    //         break;
    // }
    while (delta > args->p_params->epsilon)
    // while (!near_enough)
    {
        float3 v = make_float3(atom1.w - length(intersec_prev - atom1), atom2.w - length(intersec_prev - atom2), 0);

        float4 a = -(atom1 - intersec_prev) / length(atom1 - intersec_prev);
        float4 b = -(atom2 - intersec_prev) / length(atom2 - intersec_prev);
        float4 c = cross(a, b);

        float m_det = a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) - a.z * (b.x * c.y - b.y * c.x);
        if (abs(m_det) < ZERO)
            m_det = ZERO;
        float4 m1 = make_float4(b.y * c.z - b.z * c.y, a.z * c.y - a.y * c.z, a.y * b.z - a.z * b.y, 0);
        float4 m2 = make_float4(b.z * c.x - b.x * c.z, a.x * c.z - a.z * c.x, a.z * b.x - a.x * b.z, 0);
        float4 m3 = make_float4(b.x * c.y - b.y * c.x, a.y * c.x - a.x * c.y, a.x * b.y - a.y * b.x, 0);

        float4 inter_change = make_float4(v.x * m1.x + v.y * m1.y + v.z * m1.z, v.x * m2.x + v.y * m2.y + v.z * m2.z, v.x * m3.x + v.y * m3.y + v.z * m3.z, 0) / m_det;

        // backtracking line search
        float alpha = 1.0;    // initial step size
        float armijo_c = 0.5; // Armijo condition constant

        float4 next_point;
        float3 f_next;
        float3 f = v;

        while (true)
        {
            next_point = intersec_prev + alpha * inter_change;
            f_next = make_float3(atom1.w - length(next_point - atom1), atom2.w - length(next_point - atom2), 0);
            if (dot(f_next, f_next) <= dot(f, f) + armijo_c * alpha * dot(make_float3(inter_change.x, inter_change.y, inter_change.z), f))
                break;

            alpha *= 0.5;
        }

        intersec = next_point;
        delta = length(intersec - intersec_prev);

        // for debugging
        solver_iter++;

        // if (args->frame == 200 && args->x == 987 && args->y == 508)
        // {
        //     // x|y|frame|calculation|step|delta|prev.x|prev.y|prev.z|next.x|next.y|next.z|det|change.x|change.y|change.z
        //     printf("%i|%i|%i|%i|%i|%i|%.8f|%.3f|%.3f|%.3f|%.3f|%.3f|%.3f|%.5f|%.3f|%.3f|%.3f|%.3f|%.3f\n", args->x, args->y, args->frame, args->step, args->calcs3, solver_iter - 1, delta, intersec_prev.x, intersec_prev.y, intersec_prev.z,
        //            intersec.x, intersec.y, intersec.z, m_det, inter_change.x, inter_change.y, inter_change.z, v.x, v.y);
        // }
        intersec_prev = intersec;

        if (solver_iter > max_iter_newton)
            break;
    }

    return intersec;
}

__device__ float4 intersectTwoSpheres2b(float4 atom1, float4 atom2, functionArgs *args)
{
    float4 p = args->ray_pos;
    float3 intersec_prev = make_float3(p.x, p.y, p.z);
    float3 intersec;

    int solver_iter = 0;
    int max_iter_newton = args->p_params->solver_iter_max;

    float delta = INFINITY;

    const MatrixDim3 identity = {make_float3(1.0f, 0.0f, 0.0f), make_float3(0.0f, 1.0f, 0.0f), make_float3(0.0f, 0.0f, 1.0f)};

    const float3 atoms[2] = {make_float3(atom1), make_float3(atom2)};
    const float r[2] = {atom1.w, atom2.w};

    while (delta >= args->p_params->solver_threshold)
    {
        float3 grad = make_float3(.0f);

        MatrixDim3 hessian = {make_float3(.0f), make_float3(.0f), make_float3(.0f)};
        MatrixDim3 inverse = {make_float3(.0f), make_float3(.0f), make_float3(.0f)};

        bool quasi = false;
        bool posDef[2];
        for (int i = 0; i < 2; i++)
        {
            float3 v = intersec_prev - atoms[i];
            float l = length(v);
            posDef[i] = (l > r[i]);
            if (posDef[i])
                quasi = true;
        }

        for (int i = 0; i < 2; i++)
        {
            float3 v = intersec_prev - atoms[i];
            float l = length(v);
            grad += 2 * v * (1 - (r[i] / l));

            if (quasi)
            {
                if (posDef[i])
                {
                    hessian += (r[i] / (l * l * l)) * outer_product(v, v) + ((2 * (1 - (r[i] / l))) * identity); // Quasi-Newton
                }
                else
                {
                    hessian += (r[i] / (l * l * l)) * outer_product(v, v);
                }
            }
            else
            {
                hessian += (r[i] / (l * l * l)) * outer_product(v, v) + ((1 - (r[i] / l)) * identity); // Newton
            }
        }

        delta = length(grad);
        float det = 0;
        inverse = matrix_inverse1(hessian, true, &det);

        float3 inter_change = inverse * (-grad);

        intersec = intersec_prev + args->p_params->solver_step_size * inter_change;

        solver_iter++;
        // if (args->frame == 200 && args->x == 992 && args->y == 654)
        // {
        //     printf("%i|%i|%i|%i|%i|%i|%i|%.8f|%.8f|%.5f|%.5f|%.5f|newton 2\n", 3, args->params->solver, args->x, args->y, args->frame, args->calcs3, solver_iter,
        //            delta, length(inter_change), intersec.x, intersec.y, intersec.z);
        // }
        if (solver_iter > max_iter_newton)
            break;

        intersec_prev = intersec;
    }
    args->calcs2 = solver_iter;

    return make_float4(intersec, .0f);
}

/**
 * @brief intersection of 3 spheres using newton iteration with direct matrix inverse calculation
 *
 * @param atom1
 * @param atom2
 * @param atom3
 * @param args
 * @return __device__
 */
__device__ float4 intersectThreeSpheres2b(float4 atom1, float4 atom2, float4 atom3, functionArgs *args)
{
    float4 p = args->ray_pos;
    float delta = INFINITY;
    float4 intersec_prev = make_float4(p.x, p.y, p.z, 0);
    float4 intersec;
    int solver_iter = 0; // for debugging
    int max_iter_newton = args->p_params->solver_iter_max;
    float f_length = INFINITY;
    while (delta > args->p_params->solver_threshold)
    {
        float4 a = -(atom1 - intersec_prev) / (length(atom1 - intersec_prev));
        float4 b = -(atom2 - intersec_prev) / (length(atom2 - intersec_prev));
        float4 c = -(atom3 - intersec_prev) / (length(atom3 - intersec_prev));
        float3 v = make_float3(atom1.w - length(intersec_prev - atom1), atom2.w - length(intersec_prev - atom2), atom3.w - length(intersec_prev - atom3));
        float m_det = a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) - a.z * (b.x * c.y - b.y * c.x);
        if (abs(m_det) < ZERO)
            m_det = ZERO;
        float4 m1 = make_float4(b.y * c.z - b.z * c.y, a.z * c.y - a.y * c.z, a.y * b.z - a.z * b.y, 0);
        float4 m2 = make_float4(b.z * c.x - b.x * c.z, a.x * c.z - a.z * c.x, a.z * b.x - a.x * b.z, 0);
        float4 m3 = make_float4(b.x * c.y - b.y * c.x, a.y * c.x - a.x * c.y, a.x * b.y - a.y * b.x, 0);

        float4 inter_change = make_float4(v.x * m1.x + v.y * m1.y + v.z * m1.z, v.x * m2.x + v.y * m2.y + v.z * m2.z, v.x * m3.x + v.y * m3.y + v.z * m3.z, 0) / m_det;

        intersec = inter_change + intersec_prev;
        delta = length(intersec - intersec_prev);

        // for debugging
        solver_iter++;

        intersec_prev = intersec;

        if (solver_iter > max_iter_newton)
        {
            args->calcs1 = 42;
            break;
        }
    }
}

return intersec;
}
__device__ float4 intersectThreeSpheres2c(float4 atom1, float4 atom2, float4 atom3, functionArgs *args)
{

    float4 p = args->ray_pos;
    float delta = INFINITY;
    float4 intersec_prev = make_float4(p.x, p.y, p.z, 0);
    float4 intersec;
    int solver_iter = 0; // for debugging
    int max_iter_newton = args->p_params->solver_iter_max;
    float f_length = INFINITY;
    while (f_length > args->p_params->solver_threshold)
    // while (delta > args->p_params->epsilon)
    {
        float4 v_a1 = intersec_prev - atom1;
        float4 v_a2 = intersec_prev - atom2;
        float4 v_a3 = intersec_prev - atom3;

        float3 f = make_float3(atom1.w * atom1.w - square(v_a1), atom2.w * atom2.w - square(v_a2), atom3.w * atom3.w - square(v_a3));
        // f_length = length(f); // for debugging

        float4 a = -2 * v_a1;
        float4 b = -2 * v_a2;
        float4 c = -2 * v_a3;

        float m_det = a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) - a.z * (b.x * c.y - b.y * c.x);
        if (abs(m_det) < ZERO)
            m_det = ZERO;
        float3 m1 = make_float3(b.y * c.z - b.z * c.y, a.z * c.y - a.y * c.z, a.y * b.z - a.z * b.y);
        float3 m2 = make_float3(b.z * c.x - b.x * c.z, a.x * c.z - a.z * c.x, a.z * b.x - a.x * b.z);
        float3 m3 = make_float3(b.x * c.y - b.y * c.x, a.y * c.x - a.x * c.y, a.x * b.y - a.y * b.x);

        // Backtracking Line Search
        float alpha = 1.0;    // Initial step size
        float armijo_c = 0.5; // Armijo condition parameter
                              // Evaluate the function value at the current point

        // Compute the direction using the inverse of the Jacobian
        float4 inter_change = make_float4(f.x * m1.x + f.y * m1.y + f.z * m1.z,
                                          f.x * m2.x + f.y * m2.y + f.z * m2.z,
                                          f.x * m3.x + f.y * m3.y + f.z * m3.z, 0) /
                              m_det;

        float4 intersec_temp;
        float f_temp;
        while (true)
        {
            intersec_temp = intersec_prev + alpha * inter_change;
            float4 v1_temp = intersec_temp - atom1;
            float4 v2_temp = intersec_temp - atom2;
            float4 v3_temp = intersec_temp - atom3;

            f_temp = length(make_float3(atom1.w * atom1.w - square(v1_temp), atom2.w * atom2.w - square(v2_temp), atom3.w * atom3.w - square(v3_temp)));
            // float f_temp = length(atom1 - intersec_temp) + length(atom2 - intersec_temp) + length(atom3 - intersec_temp);
            float f_prev = length(f);
            if (f_temp <= f_prev + armijo_c * alpha * square(inter_change))
                break;

            alpha *= 0.5;
        }

        intersec = intersec_temp;

        delta = length(intersec - intersec_prev);
        f_length = f_temp;

        solver_iter++;
        // if (args->frame == 200 && args->x == 987 && args->y == 508)
        // {
        //     // x|y|frame|step|calculation|delta|prev.x|prev.y|prev.z|next.x|next.y|next.z|det|change.x|change.y|change.z
        //     printf("%i|%i|%i|%i|%i|%.8f|%.3f|%.3f|%.3f|NA|NA|NA|%.8f\n", args->x, args->y, args->frame, args->calcs3, solver_iter - 1, delta,
        //            intersec.x, intersec.y, intersec.z, m_det);
        // }
        intersec_prev = intersec;

        if (solver_iter > max_iter_newton)
        {
            args->calcs1 = 42;
            break;
        }
    }
}

return intersec;
}
__device__ float4 intersectThreeSpheres2d(float4 atom1, float4 atom2, float4 atom3, functionArgs *args)
{
    float4 p = args->ray_pos;
    float delta = INFINITY;
    float4 intersec_prev = make_float4(p.x, p.y, p.z, 0);
    float4 intersec;
    int solver_iter = 0;
    int max_iter_newton = args->p_params->solver_iter_max;
    float f_length = INFINITY;

    while (delta > args->p_params->solver_threshold)
    {
        float4 v_a1 = intersec_prev - atom1;
        float4 v_a2 = intersec_prev - atom2;
        float4 v_a3 = intersec_prev - atom3;

        float3 f = make_float3(atom1.w * atom1.w - square(v_a1), atom2.w * atom2.w - square(v_a2), atom3.w * atom3.w - square(v_a3));

        // matrix A:
        float3 A[3] = {make_float3(-2 * v_a1), make_float3(-2 * v_a2), make_float3(-2 * v_a3)};
        int P[3] = {0, 0, 0};
        lu_factorization(A, P);

        float3 inv_A[3];

        invert_matrix(A, P, inv_A);
        float4 inter_change = make_float4(dot(f, inv_A[0]), dot(f, inv_A[1]), dot(f, inv_A[2]), 0);

        // float4 inter_change = make_float4(f.x * m1.x + f.y * m1.y + f.z * m1.z,
        //                                   f.x * m2.x + f.y * m2.y + f.z * m2.z,
        //                                   f.x * m3.x + f.y * m3.y + f.z * m3.z, 0);

        intersec = intersec_prev - inter_change;

        delta = length(inter_change);
        f_length = length(f);

        solver_iter++;
        // if (args->frame == 200 && args->x == 987 && args->y == 508)
        // {
        //     // x|y|frame|step|calculation|delta|prev.x|prev.y|prev.z|next.x|next.y|next.z|det|change.x|change.y|change.z
        //     printf("%i|%i|%i|%i|%i|%.8f|%.3f|%.3f|%.3f|NA|NA|NA|%.8f\n", args->x, args->y, args->frame, args->calcs3, solver_iter - 1, delta,
        //            intersec.x, intersec.y, intersec.z, m_det);
        // }
        intersec_prev = intersec;
        if (solver_iter > max_iter_newton)
        {
            args->calcs1 = 42;
            break;
        }
    }
}

return intersec;
}

/**
 * @brief newton iteration @Oliver Kolossoski paper (using Quasi-Newtn iteration for non positive definite hessian matrix)
 *
 * @param atom1
 * @param atom2
 * @param atom3
 * @param args
 * @return __device__
 */
__device__ float4 intersectThreeSpheres2e(float4 atom1, float4 atom2, float4 atom3, functionArgs *args)
{
    float4 p = args->ray_pos;
    float3 intersec_prev = make_float3(p.x, p.y, p.z);
    float3 intersec;

    int solver_iter = 0;
    int max_iter_newton = args->p_params->solver_iter_max;

    float delta = INFINITY;

    const MatrixDim3 identity = {make_float3(1.0f, 0.0f, 0.0f), make_float3(0.0f, 1.0f, 0.0f), make_float3(0.0f, 0.0f, 1.0f)};
    const float3 atoms[3] = {make_float3(atom1), make_float3(atom2), make_float3(atom3)};
    const float r[3] = {atom1.w, atom2.w, atom3.w};

    while (delta >= args->p_params->solver_threshold)
    {
        float3 grad = make_float3(.0f);

        MatrixDim3 hessian = {make_float3(.0f), make_float3(.0f), make_float3(.0f)};
        MatrixDim3 inverse = {make_float3(.0f), make_float3(.0f), make_float3(.0f)};

        bool quasi = false;
        bool posDef[3];
        for (int i = 0; i < 3; i++)
        {
            float3 v = intersec_prev - atoms[i];
            float l = length(v);
            posDef[i] = (l > r[i]);
            if (posDef[i])
                quasi = true;
        }

        for (int i = 0; i < 3; i++)
        {
            float3 v = intersec_prev - atoms[i];
            float l = length(v);
            grad += 2 * v * (1 - (r[i] / l));

            if (quasi)
            {
                if (posDef[i])
                {
                    hessian += (r[i] / (l * l * l)) * outer_product(v, v) + ((2 * (1 - (r[i] / l))) * identity); // Quasi-Newton
                }
                else
                {
                    hessian += (r[i] / (l * l * l)) * outer_product(v, v);
                }
            }
            else
            {
                hessian += (r[i] / (l * l * l)) * outer_product(v, v) + ((1 - (r[i] / l)) * identity); // Newton
            }
        }

        delta = length(grad);
        float det = 0;
        inverse = matrix_inverse1(hessian, true, &det);

        float3 inter_change = inverse * (-grad);

        intersec = intersec_prev + args->p_params->solver_step_size * inter_change;

        solver_iter++;
        // if (args->frame == 200 && args->x == 992 && args->y == 654)
        // {
        //     printf("%i|%i|%i|%i|%i|%i|%i|%.8f|%.8f|%.5f|%.5f|%.5f|newton 2\n", 3, args->params->solver, args->x, args->y, args->frame, args->calcs3, solver_iter,
        //            delta, length(inter_change), intersec.x, intersec.y, intersec.z);
        // }
        if (solver_iter > max_iter_newton)
        {
            args->calcs1 = 42;
            break;
        }
    }

    intersec_prev = intersec;
}
args->calcs2 = solver_iter;

return make_float4(intersec, .0f);
}

__device__ float4 intersectThreeSpheres3(float4 atom1, float4 atom2, float4 atom3, functionArgs *args)
{
    float4 p = args->ray_pos;
    float delta = INFINITY;
    float4 intersec_prev = make_float4(p.x, p.y, p.z, 0);
    float4 intersec;
    int solver_iter = 0;

    int max_iter_gd = args->p_params->solver_iter_max;

    float initial_alpha = 1.0; // initial step size
    float beta = 0.5;          // reduction factor
    float armijo_c = 0.5;      // sufficient decrease parameter
    float f_length = INFINITY; // alternative exit criteria for debugging

    // while (delta > args->params->epsilon) // use delta
    while (f_length > args->p_params->solver_threshold) // use f_length as exit criteria
    {
        float4 v_a1 = intersec_prev - atom1;
        float4 v_a2 = intersec_prev - atom2;
        float4 v_a3 = intersec_prev - atom3;
        float3 r = make_float3(atom1.w * atom1.w, atom2.w * atom2.w, atom3.w * atom3.w);
        float3 f = make_float3(r.x - square(v_a1),
                               r.y - square(v_a2),
                               r.z - square(v_a3));

        // calculate jacobian
        float4 a = -2 * v_a1;
        float4 b = -2 * v_a2;
        float4 c = -2 * v_a3;

        float4 inter_change = make_float4(
            f.x * a.x + f.y * b.x + f.z * c.x,
            f.x * a.y + f.y * b.y + f.z * c.y,
            f.x * a.z + f.y * b.z + f.z * c.z,
            0);

        // Backtracking line search
        float alpha = initial_alpha;
        float4 intersec_candidate;
        while (true)
        {
            intersec_candidate = intersec_prev - alpha * inter_change;
            float f_candidate = length(make_float3(r.x - square(intersec_candidate - atom1), r.y - square(intersec_candidate - atom2), r.z - square(intersec_candidate - atom3)));

            // Armijo condition
            if (f_candidate <= length(f) - armijo_c * alpha * length(inter_change))
            {
                break;
            }
            alpha *= beta;
        }
        intersec = intersec_candidate;
        delta = length(intersec - intersec_prev);
        f_length = length(f);

        // for debugging
        solver_iter++;
        intersec_prev = intersec;

        if (solver_iter > max_iter_newton)
        {
            args->calcs1 = 42;
            break;
        }
    }
}
args->calcs2 = solver_iter;

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

__device__ pairFloat4 intersectThreeSpheres1(float4 atom1, float4 atom2, float4 atom3, functionArgs *args)
{
    float4 p = args->ray_pos;
    float4 base1 = atom2 - atom1;
    float4 base2 = atom3 - atom1;
    float b1 = length(base1);
    float b2 = length(base2);

    float4 i = normalize(base1);
    float4 j = normalize((base2 - dot(base2, i) * i));

    float r_squared1 = atom1.w * atom1.w;
    float r_squared2 = atom2.w * atom2.w;
    float r_squared3 = atom3.w * atom3.w;

    float v1 = (r_squared1 - r_squared2 + b1 * b1) / (2 * b1);
    float v2 = (r_squared1 - r_squared3 + b2 * b2 - 2 * dot(base2, i) * v1) / (2 * dot(base2, j));
    float v3 = sqrt(r_squared1 - v1 * v1 - v2 * v2);

    pairFloat4 hit;

    hit.values[0] = atom1 + i * v1 + j * v2 + cross(i, j) * v3;
    hit.values[1] = atom1 + i * v1 + j * v2 - cross(i, j) * v3;

    if (args->frame == 200 && (modf(args->x, 12) + 0) == 0)
    {
        // float4 v_a1 = p - atom1;
        // float4 v_a2 = p - atom2;
        // float4 v_a3 = p - atom3;
        // float3 f = make_float3(-square(v_a1) + atom1.w * atom1.w, -square(v_a2) + atom2.w * atom2.w, -square(v_a3) + atom3.w * atom3.w);
        // float f_length = length(f);
        // // printf("case|solver|pixel.x|pixel.y|frame|calculation|step|length(f)|delta|intersection.x|intersection.y|intersection.z|note\n");
        // printf("%i|%i|%i|%i|%i|%i|%i|%.8f|%.8f|%.5f|%.5f|%.5f|point 1\n", 3, args->p_params->solver, args->x, args->y, args->frame, args->calcs3, args->calcs2,
        //        f_length, .0f, hit.values[0].x, hit.values[0].y, hit.values[0].z);
        // printf("%i|%i|%i|%i|%i|%i|%i|%.8f|%.8f|%.5f|%.5f|%.5f|point 2\n", 3, args->p_params->solver, args->x, args->y, args->frame, args->calcs3, args->calcs2,
        //        f_length, .0f, hit.values[1].x, hit.values[1].y, hit.values[1].z);
        printf("frame|x|y|p.x|p.y|p.z|c1.x|c1.y|c1.z|c2.x|c2.y|c2.z|c3.x|c3.y|c3.z|comp.x|comp.y|comp.z\n");
        float d = INFINITY;
        int final = 0;
        for (int k = 0; k < 2; k++)
        {
            float dist = length(p - hit.values[k]);
            if (dist < d)
            {
                d = dist;
                final = k;
            }
        }

        printf("%i|%i|%i|%.8f|%.8f|%.8f|%.8f|%.8f|%.8f|%.8f|%.8f|%.8f|%.8f|%.8f|%.8f|%.8f|%.8f|%.8f\n", args->frame, args->x, args->y,
               p.x, p.y, p.z,
               atom1.x, atom1.y, atom1.z,
               atom2.x, atom2.y, atom2.z,
               atom3.x, atom3.y, atom3.z,
               hit.values[final].x, hit.values[final].y, hit.values[final].z);
    }
    args->calcs1 = 0;
}
return hit;
}

__device__ float4 intersectThreeSpheres2a(float4 atom1, float4 atom2, float4 atom3, functionArgs *args)
{
    float4 p = args->ray_pos;
    float delta = INFINITY;
    float4 intersec_prev = make_float4(p.x, p.y, p.z, 0);
    float4 intersec;
    int solver_iter = 0;
    int max_iter_newton = args->p_params->solver_iter_max;
    float f_length = INFINITY;

    while (delta > args->p_params->solver_threshold)
    {
        float4 a = -(atom1 - intersec_prev) / (length(atom1 - intersec_prev));
        float4 b = -(atom2 - intersec_prev) / (length(atom2 - intersec_prev));
        float4 c = -(atom3 - intersec_prev) / (length(atom3 - intersec_prev));
        float3 v = make_float3(atom1.w - length(intersec_prev - atom1), atom2.w - length(intersec_prev - atom2), atom3.w - length(intersec_prev - atom3));
        float m_det = a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) - a.z * (b.x * c.y - b.y * c.x);
        if (abs(m_det) < ZERO)
            m_det = ZERO;
        float4 m1 = make_float4(b.y * c.z - b.z * c.y, a.z * c.y - a.y * c.z, a.y * b.z - a.z * b.y, 0);
        float4 m2 = make_float4(b.z * c.x - b.x * c.z, a.x * c.z - a.z * c.x, a.z * b.x - a.x * b.z, 0);
        float4 m3 = make_float4(b.x * c.y - b.y * c.x, a.y * c.x - a.x * c.y, a.x * b.y - a.y * b.x, 0);

        float4 inter_change = make_float4(v.x * m1.x + v.y * m1.y + v.z * m1.z, v.x * m2.x + v.y * m2.y + v.z * m2.z, v.x * m3.x + v.y * m3.y + v.z * m3.z, 0) / m_det;

        intersec = inter_change + intersec_prev;
        delta = length(intersec - intersec_prev);

        // for debugging
        solver_iter++;

        // if (frame == 200 && x == 1000 && y == 650)
        // {
        //     // x|y|frame|calculation|delta|prev.x|prev.y|prev.z|next.x|next.y|next.z|det|change.x|change.y|change.z
        //     printf("%i|%i|%i|%i|%i|%.8f|%.3f|%.3f|%.3f|%.3f|%.3f|%.3f|%.5f|%.3f|%.3f|%.3f|%.3f|%.3f\n", x, y, frame, calculations, solver_iter - 1, delta, intersec_prev.x, intersec_prev.y, intersec_prev.z,
        //            intersec.x, intersec.y, intersec.z, m_det, inter_change.x, inter_change.y, inter_change.z, v.x, v.y);
        //     // printf("DEBUG: step: %i, delta: %.3f, x: %.3f, %.3f, %.3f, det: %.3f\n", solver_iter - 1, delta, inter_change.x, inter_change.y, inter_change.z, m_det);
        // }
        intersec_prev = intersec;

        if (solver_iter > max_iter_newton)
        {
            args->calcs1 = 42;
            break;
        }
    }

    return intersec;
}

__device__ bool testIntersectionType3(const float4 intersection_point, float4 &hit, const int i, const int j, const int k, functionArgs *args)
{
    if (abs(localSDouter(intersection_point, args)) < args->p_params->epsilon)
    {
        float dist = length(args->ray_pos - intersection_point);

        if (dist < args->shortest_distance)
        {
            args->hit_feedback->bondId1 = args->atom_data_local[i].id;
            args->hit_feedback->bondId2 = args->atom_data_local[j].id;
            args->hit_feedback->bondId3 = args->atom_data_local[k].id;

            hit = intersection_point;
            args->shortest_distance = dist;
            return true;
        }
    }
    return false;
}
__device__ bool testIntersectionType3(const pairFloat4 intersection_points, float4 &hit, const int i, const int j, const int k, functionArgs *args)
{
    bool surface_found = false;
    for (int n = 0; n < 2; n++)
    {
        if (abs(localSDouter(intersection_points.values[n], args)) < args->p_params->epsilon)
        {
            float dist = length(args->ray_pos - intersection_points.values[n]);
            if (dist < args->shortest_distance)
            {
                args->hit_feedback->bondId1 = args->atom_data_local[i].id;
                args->hit_feedback->bondId2 = args->atom_data_local[j].id;
                args->hit_feedback->bondId3 = args->atom_data_local[k].id;

                hit = intersection_points.values[n];
                args->shortest_distance = dist;
                surface_found = true;
            }
        }
    }
    return surface_found;
}
__device__ float calculateDistanceToGrid(functionArgs *args)
{
    // TODO: optimize code: use only start and end point of box for the faces (makes first array unnecessary)
    // find the closest face of the grid
    float4 box_dimensions = (args->p_params->box_end - args->p_params->box_start) / 2;
    float4 grid_center = args->p_params->box_start + box_dimensions;

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

            float distance = length(args->ray_pos - midpoint);
            if (distance < distance_to_face)
            {
                distance_to_face = distance;
                face_midpoint = midpoint;
                face_normal = normal;
            }
        }
    }

    // calculate distance to face along ray_dir
    float numerator = dot((face_midpoint - args->ray_pos), face_normal);
    float denominator = dot(normalize(args->ray_dir), face_normal);

    float traveldist = 2 * args->p_params->solvent_max; // minimum distance to be traveled in undefined cases
    if (denominator == 0)
    {
        if (numerator == 0)
            return args->p_params->epsilon; // point already lies on the face; set epsilon as distance to avoid infinite loops
        else
            return traveldist; // moving parallel to the box -> ray can't hit the molecule
    }
    float result = (numerator / denominator);

    // distance to grid is measured along ray direction
    // -> negative distances implicate the box has already been traversed
    if (result < 0)
    {
        return traveldist;
    }

    return result + args->p_params->epsilon;
}

// __device__ void findNearestAtoms(float4 *molecule_data, uint number_atoms, Atom nearest_atoms[], uint nearest_count, float4 position, int x = 0, int y = 0)
// {
//     for (uint i = 0; i < nearest_count; i++)
//     {
//         float dist = calcSignedDistanceSphere(position, molecule_data[i]);

//         nearest_atoms[i].location = molecule_data[i];
//         nearest_atoms[i].id = i;
//         nearest_atoms[i].distance = dist;
//     }
//     for (uint i = nearest_count; i < number_atoms; i++)
//     {
//         float dist = calcSignedDistanceSphere(position, molecule_data[i]);

//         int max_element = findFurthestByIndex(nearest_atoms, nearest_count);
//         if (nearest_atoms[max_element].distance > dist)
//         {
//             nearest_atoms[max_element].location = molecule_data[i];
//             nearest_atoms[max_element].id = i;
//             nearest_atoms[max_element].distance = dist;
//         }
//     }
//     if (nearest_count > 1)
//     {
//         selectionSort(nearest_count, nearest_atoms);
//     }
// }

__device__ void findNearestAtoms(functionArgs *args)
{
    for (uint i = 0; i < args->p_params->k_nearest; i++)
    {
        float dist = calcSignedDistanceSphere(args->ray_pos, args->atom_data_global[i]);

        args->atom_data_local[i].location = args->atom_data_global[i];
        args->atom_data_local[i].id = i;
        args->atom_data_local[i].distance = dist;
    }
    for (uint i = args->p_params->k_nearest; i < args->p_params->numAtoms; i++)
    {
        float dist = calcSignedDistanceSphere(args->ray_pos, args->atom_data_global[i]);

        int max_element = findFurthestByIndex(args->atom_data_local, args->p_params->k_nearest);
        if (args->atom_data_local[max_element].distance > dist)
        {
            args->atom_data_local[max_element].location = args->atom_data_global[i];
            args->atom_data_local[max_element].id = i;
            args->atom_data_local[max_element].distance = dist;
        }
    }
    if (args->p_params->k_nearest > 1)
    {
        selectionSort(args->p_params->k_nearest, args->atom_data_local);
    }
}

// __device__ unsigned int nearestOverVoxel(float4 *molecule_data, int *voxel_data, int *voxel_count, Atom nearest_atoms[], float4 position, float4 direction, SimulationParams params, int x = 0, int y = 0)
// {
//     unsigned int nearest_count = 0;

//     // calculate current voxel
//     int4 current_voxel = castf2i(floorf((position - params.box_start) / params.voxel_size));

//     // TODO: calculate nearest atoms with voxels
//     // determine closest points for 3x3x3 voxels
//     int4 range_min = max(current_voxel - make_int4(1, 1, 1, 0), make_int4(0, 0, 0, 0));
//     int4 range_max = min(current_voxel + make_int4(1, 1, 1, 0), params.voxel_dim - make_int4(1, 1, 1, 0));

//     bool k_full = false;

//     for (int i = range_min.x; i <= range_max.x; i++)
//     {
//         for (int j = range_min.y; j <= range_max.y; j++)
//         {
//             for (int k = range_min.z; k <= range_max.z; k++)
//             {
//                 // 1. calculate voxel id
//                 int id = castVoxelToId(params.voxel_dim, i, j, k);

//                 int voxel_start = id * params.atomsPerVoxel;
//                 int voxel_end = (id + 1) * params.atomsPerVoxel;

//                 for (int v = voxel_start; v < voxel_end; v++)
//                 {
//                     int m = voxel_data[v];
//                     if (m >= 0)
//                     {
//                         float dist = calcSignedDistanceSphere(position, molecule_data[m]);

//                         if (!k_full)
//                         {
//                             // 3. copy all elements to nearest_atoms
//                             nearest_atoms[nearest_count].location = molecule_data[m];
//                             nearest_atoms[nearest_count].id = m;
//                             nearest_atoms[nearest_count].distance = dist;
//                             nearest_count++;

//                             if (nearest_count >= params.k_nearest)
//                             {
//                                 k_full = true;
//                             }
//                         }
//                         else
//                         {
//                             int max_element = findFurthestByIndex(nearest_atoms, nearest_count);
//                             if (nearest_atoms[max_element].distance > dist)
//                             {
//                                 nearest_atoms[max_element].location = molecule_data[m];
//                                 nearest_atoms[max_element].id = m;
//                                 nearest_atoms[max_element].distance = dist;
//                             }
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     // swap the closest atom to be at the start of the array(only applicable with multiple atoms in neighbourhood)
//     if (nearest_count > 1)
//     {
//         selectionSort(nearest_count, nearest_atoms);
//     }
//     return nearest_count;
// }

__device__ unsigned int nearestOverVoxel(functionArgs *args)
{
    unsigned int nearest_count = 0;

    // calculate current voxel
    int4 current_voxel = castf2i(floorf((args->ray_pos - args->p_params->box_start) / args->p_params->voxel_size));

    // TODO: calculate nearest atoms with voxels
    // determine closest points for 3x3x3 voxels
    int4 range_min = max(current_voxel - make_int4(1, 1, 1, 0), make_int4(0, 0, 0, 0));
    int4 range_max = min(current_voxel + make_int4(1, 1, 1, 0), args->p_params->voxel_dim - make_int4(1, 1, 1, 0));

    bool k_full = false;

    for (int i = range_min.x; i <= range_max.x; i++)
    {
        for (int j = range_min.y; j <= range_max.y; j++)
        {
            for (int k = range_min.z; k <= range_max.z; k++)
            {
                // 1. calculate voxel id
                int id = castVoxelToId(args->p_params->voxel_dim, i, j, k);

                int voxel_start = id * args->p_params->atomsPerVoxel;
                int voxel_end = (id + 1) * args->p_params->atomsPerVoxel;

                for (int v = voxel_start; v < voxel_end; v++)
                {
                    int m = args->voxel_data[v];
                    if (m >= 0)
                    {
                        float dist = calcSignedDistanceSphere(args->ray_pos, args->atom_data_global[m]);

                        if (!k_full)
                        {
                            // 3. copy all elements to nearest_atoms
                            args->atom_data_local[nearest_count].location = args->atom_data_global[m];
                            args->atom_data_local[nearest_count].id = m;
                            args->atom_data_local[nearest_count].distance = dist;
                            nearest_count++;

                            if (nearest_count >= args->p_params->k_nearest)
                            {
                                k_full = true;
                            }
                        }
                        else
                        {
                            int max_element = findFurthestByIndex(args->atom_data_local, nearest_count);
                            if (args->atom_data_local[max_element].distance > dist)
                            {
                                args->atom_data_local[max_element].location = args->atom_data_global[m];
                                args->atom_data_local[max_element].id = m;
                                args->atom_data_local[max_element].distance = dist;
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
        selectionSort(nearest_count, args->atom_data_local);
    }
    return nearest_count;
}
__device__ float computeSurface(const float4 ray_pos, const float4 ray_dir, float4 *molecule, int *voxel_data, int *voxel_count, SimulationParams params, hitInfo *surfacePointData, int x = 0, int y = 0, int frame = 0, int step = 0)
{
    // const float4 ray_pos;     // copy of ray position
    // const float4 ray_dir;     // copy of ray direction
    // SimulationParams *params; // device parameters that can be configured from the UI
    // float4 *atom_data_global; // points to the complete array of atoms
    // int *voxel_data;          // points to array of sorted atoms based on their respective voxel in the grid
    // Atom *atom_data_local;    // points to a shortend array of close atoms
    // hitInfo *hit_feedback;    // contains information of surface hit

    // const int x, y;  // current pixel
    // const int frame; // currently generated frame
    // int step;        // can be used as counter for loops
    // int calcs1;      // can be used as counter for loops
    // int calcs2;      // can be used as counter for loops
    // int calcs3;      // can be used as counter for loops

    functionArgs args = {ray_pos, ray_dir, &params, molecule, voxel_data, nullptr, surfacePointData, INFINITY, x, y, frame, step, 0, 0, 0};

    float f_sdf = 0; // current distance to SAS

    surfacePointData->collisionType = 0;

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // // 1 // determine signed distance to closest atom
    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (params.use_voxel)
    {
        // position outside grid (position cannot lie on molecule surface)
        if (ray_pos < (params.box_start) || ray_pos > (params.box_end))
        {

            if (surfacePointData->isInGrid)
            {
                surfacePointData->isInGrid = false;
                surfacePointData->traversedGrid = true;
                return INFINITY;
            }

            if (step == 0)
            {
                float4 diag = params.box_end - params.box_start;
                float4 gridcenter = params.box_start + 0.5 * diag;
                float traveldist = length(gridcenter - ray_pos) - 0.5 * length(diag);
                return max(traveldist, 0.0f);
            }

            // f_sdf = calculateDistanceToGrid(ray_pos, ray_dir, params, frame, x, y);
            // f_sdf = params.solvent_radius;

            return 0;
        }
        surfacePointData->isInGrid = true;
    }
    if (!params.use_voxel)
    {
        f_sdf = globalSDouterWithFeedback(ray_pos, &args);
        // f_sdf = calcSignedDistanceOuterWithNearestAtom(ray_pos, molecule, params.numAtoms, surfacePointData);
        if (f_sdf >= 0)
            return f_sdf;
    }

    Atom nearest_atoms[GROUP_SIZE];
    args.atom_data_local = nearest_atoms;

    unsigned int nearest_count = 0;

    if (params.use_voxel)
    {
        // determine neighbourhood of atoms over 3 x 3 x 3 voxels
        nearest_count = nearestOverVoxel(&args);
        // nearest_count = nearestOverVoxel(molecule, voxel_data, voxel_count, nearest_atoms, ray_pos, ray_dir, params, x, y);

        if (nearest_count == 0)
        {
            f_sdf = 2 * params.solvent_max;
            return f_sdf;
        }

        f_sdf = localSDouter(ray_pos, &args);
        // f_sdf = calcSignedDistanceOuter(ray_pos, nearest_atoms, nearest_count);

        if (f_sdf >= 0)
            return f_sdf;
    }

    if (f_sdf < 0)
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 2 // determine k closest points
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // allocate array for ids of k atoms closest to p
        // find the k nearest atoms

        if (!params.use_voxel)
        {
            // nearest_count = min(params.k_nearest, params.numAtoms);
            nearest_count = (uint)params.k_nearest;
            findNearestAtoms(&args);
            // findNearestAtoms(molecule, params.numAtoms, nearest_atoms, nearest_count, ray_pos);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 3 // intersection with one atom

        // group of participating atoms is ordered, nearest_atoms[0] contains the id of the closest atom

        float4 dir = normalize(ray_pos - nearest_atoms[0].location);
        float4 surface_hit = nearest_atoms[0].location + nearest_atoms[0].location.w * dir;
        surfacePointData->surfaceHit = surface_hit;
        surfacePointData->collisionType = 1;
        surfacePointData->bondId1 = nearest_atoms[0].id;

        if (abs(localSDouter(surface_hit, &args)) > params.epsilon)
        // if (abs(calcSignedDistanceOuter(surface_hit, nearest_atoms, params.k_nearest)) > params.epsilon)
        {
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

            ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // 4 // intersection with two atoms (arc)

            bool surface_found = false;
            int calcs = 0;
            // calculate all intersection points for the group of nearest atoms
            for (unsigned int i = 0; i < params.k_nearest - 1; i++)
            {
                for (unsigned int j = i + 1; j < params.k_nearest; j++)
                {
                    // surface_hit = intersectTwoSpheres1( nearest_atoms[i].location, nearest_atoms[j].location);
                    switch (params.solver)
                    {
                    default:
                        surface_hit = intersectTwoSpheres1(nearest_atoms[i].location, nearest_atoms[j].location, &args);
                        break;
                        // case 0:
                        //     surface_hit = intersectTwoSpheres1(nearest_atoms[i].location, nearest_atoms[j].location, &args);
                        //     break;
                        // case 1:
                        //     surface_hit = intersectTwoSpheres2b(nearest_atoms[i].location, nearest_atoms[j].location, &args);
                        //     calcs++;
                        //     break;
                        // case 2:
                        //     surface_hit = intersectTwoSpheres3(nearest_atoms[i].location, nearest_atoms[j].location, &args);
                        //     calcs++;
                        //     break;

                        // default:
                        //     if (frame == 0 && x == 0 && y == 0)
                        //         printf("ERROR: invalid solver. Solver %i is not defined.\n", params.solver);
                        //     break;
                    }

                    // check if points lie in correct section (arc) of the surface
                    bool triangle = in_triangle(ray_pos, nearest_atoms[i].location, nearest_atoms[j].location, surface_hit);
                    // check if point lies on the surface

                    bool surface = (abs(localSDouter(surface_hit, &args)) < params.epsilon);
                    // bool surface = (abs(calcSignedDistanceOuter(surface_hit, nearest_atoms, params.k_nearest)) < params.epsilon);
                    if (triangle && surface)
                    {
                        surface_found = true;
                        surfacePointData->collisionType = 2;
                        surfacePointData->surfaceHit = surface_hit;
                        surfacePointData->bondId1 = nearest_atoms[i].id;
                        surfacePointData->bondId2 = nearest_atoms[j].id;
                        break;
                    }
                }
                if (surface_found)
                {
                    break;
                }
            }

            //////////////////////////////////////////////////////////////////////////////////////////////////////

            /////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // 5 // intersection with three atoms (spherical triangle)
            if (!surface_found)
            {
                for (unsigned int i = 0; i < params.k_nearest - 2; i++)
                {
                    for (unsigned int j = i + 1; j < params.k_nearest - 1; j++)
                    {
                        for (unsigned int k = j + 1; k < params.k_nearest; k++)
                        {
                            switch (params.solver)
                            {
                            case 0:
                                // default:
                                {
                                    pairFloat4 s_pot = intersectThreeSpheres1(nearest_atoms[i].location, nearest_atoms[j].location, nearest_atoms[k].location, &args);
                                    args.calcs3++; // for debugging
                                    if (testIntersectionType3(s_pot, surface_hit, i, j, k, &args))
                                        surface_found = true;
                                    break;
                                }
                            case 1:
                            {
                                // TODO: adjust intersection functions to take integer (ids) for atoms
                                float4 s_pot = intersectThreeSpheres2e(nearest_atoms[i].location, nearest_atoms[j].location, nearest_atoms[k].location, &args);
                                args.calcs3++;
                                if (testIntersectionType3(s_pot, surface_hit, i, j, k, &args))
                                    surface_found = true;
                                // surface_hit = make_float4(INFINITY);
                                // surface_found = true;
                                // surfacePointData->collisionType = 4;
                                // return INFINITY;
                                break;
                            }
                            case 2:
                            {
                                float4 s_pot = intersectThreeSpheres3(nearest_atoms[i].location, nearest_atoms[j].location, nearest_atoms[k].location, &args);
                                args.calcs3++;
                                if (testIntersectionType3(s_pot, surface_hit, i, j, k, &args))
                                    surface_found = true;
                                break;
                            }
                            default:
                            {
                                if (frame == 0 && x == 0 && y == 0)
                                    printf("ERROR: invalid solver. Solver %i is not defined.\n", params.solver);
                                break;
                            }
                            }
                        }
                    }
                }

                if (surface_found)
                {
                    surfacePointData->surfaceHit = surface_hit;
                    surfacePointData->collisionType = 3;

                    // if (frame == 199 && x == 761 && y == 634)
                    // {
                    //     // // test matrix inverse
                    //     // MatrixDim3 A;
                    //     // // A.rows[0] = make_float3(4, 7, 2);
                    //     // // A.rows[1] = make_float3(3, 5, 1);
                    //     // // A.rows[2] = make_float3(2, 4, 3);
                    //     // float3 s = make_float3(3, 9, 2);
                    //     // float3 q = make_float3(1, 2, -4);
                    //     // A = outer_product(s, q);
                    //     // MatrixDim3 expected;

                    //     // expected.rows[0] = make_float3(3, 6, -12);
                    //     // expected.rows[1] = make_float3(9, 18, -36);
                    //     // expected.rows[2] = make_float3(2, 4, -8);

                    //     // printf("DEBUG: testing inverse\n");
                    //     // for (int l = 0; l < 3; l++)
                    //     // {
                    //     //     float3 test = A.rows[l] - expected.rows[l];
                    //     //     printf("delta row %i: %.8f, %.8f,%.8f\n", l, test.x, test.y, test.z);
                    //     // }

                    //     // printf("case|solver|pixel.x|pixel.y|frame|calculation|step|length(f)|delta|intersection.x|intersection.y|intersection.z|note\n");
                    //     // printf("%i|%i|%i|%i|%i|%i|%i|%.8f|%.8f|%.5f|%.5f|%.5f|final surface hit\n", 3, params.solver, x, y, frame, calcs, args.calcs2,
                    //     //        .0f, .0f, surface_hit.x, surface_hit.y, surface_hit.z);
                    //     // printf("%i|%i|%i|%i|%i|%i|%i|%.8f|%.8f|%.5f|%.5f|%.5f|final ray position\n", 3, params.solver, x, y, frame, calcs, 0,
                    //     //        .0f, .0f, ray_pos.x, ray_pos.y, ray_pos.z);
                    //     // printf("%i|%i|%i|%i|%i|%i|%i|%.8f|%.8f|%.5f|%.5f|%.5f|final ray direction\n", 3, params.solver, x, y, frame, calcs, 0,
                    //     //        .0f, .0f, ray_dir.x, ray_dir.y, ray_dir.z);
                    // }
                }
                else
                {
                    // collision type 4 occurs when passing through the molecule without hitting anything until reaching max depth
                    surfacePointData->collisionType = 4;
                }
            }
            f_sdf = -1 * length(surface_hit - ray_pos);
        }
    }
    // if (frame == 199 && x > 1013 && x < 1047 && y > 443 && y < 461)
    // {
    //     printf("%i|%i|%i|%i|%i|%i|%i|%.8f|%.8f|%.5f|%.5f|%.5f|final surface hit\n", 3, params.solver, x, y, frame, surfacePointData->collisionType, args.calcs2,
    //            .0f, .0f, surfacePointData->surfaceHit.x, surfacePointData->surfaceHit.y, surfacePointData->surfaceHit.z);
    // }

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    return f_sdf;
}

const int TEST_NUM = 500;
const int TEST_METHODS = 5;
struct intersection_tester
{
    float4 position[TEST_NUM];
    float4 point1[TEST_NUM];
    float4 point2[TEST_NUM];
    float4 point3[TEST_NUM];
    float4 comparison[TEST_NUM];
};

__global__ void test_kernel(SimulationParams params, int frame = 0)
{
    float epsilon = 0.001f;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    intersection_tester data = {};

    functionArgs args = {make_float4(.0f), make_float4(.0f), &params, nullptr, nullptr, nullptr, nullptr, INFINITY, x, y, frame, 0, 0, 0, 0};
    float accumulated_delta = 0.0f;

    int count_hit[TEST_METHODS];          // countup if intersection point found
    int count_miss[TEST_METHODS];         // countup if no intersection point found
    int count_correct[TEST_METHODS];      // countup if delta to comparison correct
    int count_incorrect[TEST_METHODS];    // countup if delta to comparison not correct
    int count_converged[TEST_METHODS];    // countup if intersection calculation terminates before reaching the iteration threshold
    int count_iter_reached[TEST_METHODS]; // countup if calculation reaches iteration threshold
    float time[TEST_METHODS];             // for time measurement (in milliseconds)

    for (int i = 0; i < TEST_METHODS; i++)
    {
        count_hit[i] = 0;          // countup if intersection point found
        count_miss[i] = 0;         // countup if no intersection point found
        count_correct[i] = 0;      // countup if delta to comparison correct
        count_incorrect[i] = 0;    // countup if delta to comparison not correct
        count_converged[i] = 0;    // countup if intersection calculation terminates before reaching the iteration threshold
        count_iter_reached[i] = 0; // countup if calculation reaches iteration threshold
    }

    for (int j = 0; j < TEST_METHODS; j++)
    {
        cudaEventRecord(start);
        for (int i = 0; i < TEST_NUM; i++)
        {
            args.ray_pos = data.position;
            float4 intersection;

            switch (j)
            {
            case 0:
                pairFloat4 temp = intersectThreeSpheres1(data.point1, data.point2, data.point3, &args);
                float d = INFINITY;
                int finalRes = 0;
                for (int k = 0; k < 2; k++)
                {
                    float dist = length(data.comparison - temp.values[k]);

                    if (dist < d)
                    {
                        finalRes = k;
                    }
                }
                intersection = temp.values[finalRes];
                break;

            case 1:
                intersection = intersectThreeSpheres2a(data.point1, data.point2, data.point3, &args);
                break;
            case 2:
                intersection = intersectThreeSpheres2b(data.point1, data.point2, data.point3, &args);
                break;
            case 3:
                intersection = intersectThreeSpheres2c(data.point1, data.point2, data.point3, &args);
                break;
            case 4:
                intersection = intersectThreeSpheres2d(data.point1, data.point2, data.point3, &args);
                break;
            case 5:
                intersection = intersectThreeSpheres2e(data.point1, data.point2, data.point3, &args);
                break;
            case 6:
                intersection = intersectThreeSpheres3(data.point1, data.point2, data.point3, &args);
                break;
            default:
                break;
            }
            // calculate statistics
            float delta = length(data.comparison - intersection);

            (intersection = make_float4(.0f)) ? count_miss[j]++ : count_hit[j]++;
            (delta < epsilon) ? count_correct[j]++ : count_incorrect[j]++;
            (args.calcs1 == 42) ? count_iter_reached[j]++ : count_converged[j]++;
        }
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&time[j], start, stop);
    }
    printf("---------------------------------------------------------------\n");
    printf("Intersection report (epsilon = %.8f\n", epsilon);
    printf("---------------------------------------------------------------\n");
    printf("method|time|total_calculations|share hits|share misses|share correct|share incorrect|share_converged|share_terminated\n");
    printf("---------------------------------------------------------------\n");
    for (int i = 0; i < TEST_METHODS; i++)
    {
        printf("%i|%.0f|%i|%i|%i|%i|%i|%i|%i\n", i, time[i], TEST_NUM,
               count_hit[i],
               count_miss[i],
               count_correct[i],
               count_incorrect[i],
               count_converged[i],
               count_iter_reached[i]);
    }
    printf("---------------------------------------------------------------\n");
}
#endif