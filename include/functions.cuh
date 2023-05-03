#ifndef FUNCTIONS_CUH
#define FUNCTIONS_CUH

#ifndef INFINITY
#define INFINITY 1e8
#endif

#include <math.h>
#include <thrust/sort.h>

// get first 3 coordinates from float4
inline __host__ __device__ float3 get_pos(float4 a)
{
    float3 b = make_float3(a.x, a.y, a.z);
    return b;
}

inline __host__ __device__ float distance_squared(float3 a, float3 b)
{
    float dist = a.x * a.x + a.y * a.y + a.z * a.z;
    return dist;
}

// takes first 3 coordinates of two float4 vectors to calculate euclidian distance squared
inline __host__ __device__ float distance_squared(float4 a, float4 b)
{
    float dist = a.x * a.x + a.y * a.y + a.z * a.z - b.x * b.x - b.y * b.y - b.z * b.z;
    return dist;
}

__device__ float distance_from_sas(float4 vec, float4 *atoms, int atom_count, bool with_keys = false, int *keys = atoms)
{
    float dist = INFINITY;
    if (with_keys)
    {
        for (unsigned int i = 0; i < atom_count; i++)
        {
            dist = min(dist, distance_squared(vec, atoms[keys[i]]));
        }
    }
    else
    {
        for (unsigned int i = 0; i < atom_count; i++)
        {
            dist = min(dist, distance_squared(vec, atoms[i]));
        }
    }

    return dist;
}
// takes first 3 coordinates of two float4 vectors to calculate
// euclidian distance squared minus the square of the second radius
inline __host__ __device__ float distance_squared_sas(float4 a, float4 b)
{
    float dist = a.x * a.x + a.y * a.y + a.z * a.z - b.x * b.x - b.y * b.y - b.z * b.z - b.w * b.w;
    return dist;
}

// calculate k smallest elements based on value array
// return array of corresponding keys
__device__ int *k_smallest_by_key(float *values_in[], int arr_length, int k)
{
    // generate keys
    int keys_in[arr_length];
    int keys_out[k];

    for (unsigned int i = 0; i < arr_length; i++)
    {
        keys_in[i] = i;
    }

    // sort
    thrust::sort_by_key(keys_in, keys_in + arr_length, values);

    // set output
    for (unsigned int i = 0; i < k; i++)
    {
        keys_out[i] = keys_in[i];
    }

    return keys_out;
}

// test predicate for case2:
// position of sample point has to be within extended atom radius of both involved atoms
// (extended radius refers to atom radius plus solvent radius)
__device__ bool case2_predicate(float4 p, float4 atom1, float4 atom2, float solventRadius)
{
    float g1_squared = (atom1.w + solventRadius) * (atom1.w + solventRadius);
    float g2_squared = (atom2.w + solventRadius) * (atom2.w + solventRadius);
    if (distance_squared(p, atom1) > g1_squared && distance_squared(p, atom2) > g2_squared)
    {
        return true;
    }
    return false;
}

// calculates the intersection point between two spheres that is closest to the sample point p
__device__ float4 sphere_sphere_intersec(float4 p, float4 atom1, float4 atom2, float solventRadius)
{
    // extend radii
    float r1 = atom1.w + solventRadius;
    float r2 = atom2.w + solventRadius;
    // calculate distance to midpoint between atoms (2nd intercept theorem)
    float d = length(atom1, atom2);
    float m_length = (d * d - r2 * r2 + r1 * r1) / (2 * d);
    // calculate distance from midpoint to intersection point
    float h = sqrt(r1 * r1 - m_length * m_length);

    // calculate vector from atom1 to midpoint
    float4 dir = normalize(atom2 - atom1);
    float4 m = dir * m_length + atom1;

    // calculate vector perpendicular to the line between the atom centers in the direction of p
    //  - vector direction is obtained from perpendicular through p
    //  - vector start point is midpoint between atoms
    float4 q = p - m + dot((atom1 - p), dir) * dir;

    // calculate nearest intersection point
    float4 intersec = normalize(q) * h + m;
    return intersec;
}

// takes the intersection point of two spheres and calculates the nearest surface point of the arc connecting
// the two corresponding circles
__device__ float4 arc_surface_point(float4 p, float4 intersec, float solventRadius)
{
    float4 surfacePoint = normalize(p - intersec) * solventRadius + intersec;
    return surfacePoint;
}

// tests for given surface point:
// a) if the sample point p lies in the triangle spanned by the atom centers and the surface point
// b) if the surface point lies on the surface
__device__ bool test_toroid(float4 p, float4 c1, float4 c2, float3 surface_point, float4 *atoms, int atom_count, int *keys)
{
    float3 v0 = get_pos(c1) - surface_point;
    float3 v1 = get_pos(c2) - surface_point;
    float3 v2 = get_pos(p) - surface_point;

    float dot00 = dot(v0, v0);
    float dot01 = dot(v0, v1);
    float dot02 = dot(v0, v2);
    float dot11 = dot(v1, v1);
    float dot12 = dot(v1, v2);

    // Compute barycentric coordinates
    float denominator = 1 / (dot00 * dot11 - dot01 * dot01);
    float u = (dot11 * dot02 - dot01 * dot12) * denominator;
    float v = (dot00 * dot12 - dot01 * dot02) * denominator;

    // Check if point is in triangle
    bool triangle = (u >= 0) &&
                    (v >= 0) &&
                    (u + v < 1);

    if (distance_from_sas(surface_point, atoms, atom_count, true, keys) < epsilon && triangle)
    {
        return true;
    }

    return false;
}
// test, if sample point p lies on the toroid surface generated by atom a and b
__device__ float3 toroid_surface_hit(float4 p, float4 c1, float4 c2)
{
    // calculate the midpoint between both circles
    float d_squared = distance_squared(c1, c2);
    float d = (d_squared - c1.w * c1.w + c2.w * c2.w) / sqrt(d_squared);

    a = normalize(get_pos(p) - get_pos(c1));
    b = normalize(get_pos(c1) - get_pos(c2));
    float3 m = get_pos(c1) + d * b;
    float h = sqrt(c1.w * c1.w - d * d);
    float3 f = normalize(a - dot(a, b) * b);

    return true;
}
