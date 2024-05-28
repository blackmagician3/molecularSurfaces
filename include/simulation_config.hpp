#ifndef SIMULATION_CONFIG_HPP
#define SIMULATION_CONFIG_HPP

#include <vector>
#include <helper_math.h> // includes vector types

typedef unsigned int uint;
struct ThreadBlock
{
    uint threadX;
    uint threadY;
};
struct SimulationParams
{
    // window size
    uint window_width;
    uint window_height;

    // kernel config
    uint thread_x;
    uint thread_y;

    // textures as buffer between CUDA and OpenGL
    uint texture_1;
    uint texture_2;

    // constants & variables for raymarching
    uint numAtoms;
    uint k_nearest;
    bool is_SES;
    bool use_voxel;
    int solver;
    float solvent_radius;
    float solvent_max;
    float epsilon;
    float epsilon_squared;
    float depth_min;
    float depth_max;
    int4 voxel_dim; // (number voxels in x, number voxels in y, number voxels in z, 0)
    float voxel_size;
    float4 box_start; // offset for grid starting point
    float4 box_end;   // offset for grid ending point
    int atomsPerVoxel;

    // coloring
    uint colorScheme;
    // for debugging
    bool debug_mode;
    double mouse_x_pos;
    double mouse_y_pos;
    float highlight_radius;
    int debug_frame;
};

/**
 * @brief takes an array of atoms and adds a value to their 4th component (radius)
 *
 * @param arr array of atoms (float4)
 * @param arrLength number of atoms
 * @param value value to increase the radius by
 */
inline __host__ __device__ void addToRadius(float4 *arr, int arrLength, float value)
{
    for (unsigned int i = 0; i < arrLength; i++)
    {
        arr[i].w += value;
    }
}

/**
 * @brief find all divisors of a given integer
 *
 * @param number
 * @return std::vector<unsigned int>
 */
std::vector<unsigned int> findDivs(unsigned int number)
{
    std::vector<unsigned int> result;
    for (int i = 1; i <= sqrt(number); i++)
    {
        // if 'i' is factor of n
        if (number % i == 0)
        {
            // check if divisors are equal
            if (number / i == i)
                result.push_back(i);
            else
            {
                result.push_back(i);
                result.push_back(number / i);
            }
        }
    }
    return result;
}

/**
 * @brief determine threads per block for the raymarching kernel; maximizing threads per block for current resolution
 *
 * @param threadX pointer to the number of threads in x direction
 * @param threadY pointer to the number of threads in y direction
 * @param width window width for current resolution
 * @param height window height for current resolution
 */
ThreadBlock calcThreadBlock(unsigned int width, unsigned int height)
{
    // get hardware specs
    // TODO: make sure suitable hardware is used

    // get current device
    int device;
    cudaGetDevice(&device);

    // get hardware specs
    struct cudaDeviceProp props;
    cudaGetDeviceProperties(&props, device);

    unsigned int maxT = floor(props.maxThreadsPerBlock / 2);

    std::vector<unsigned int> divWidth = findDivs(width);
    std::vector<unsigned int> divHeight = findDivs(height);
    unsigned int threads = 0;
    unsigned int tx, ty;
    for (unsigned int x : divWidth)
    {
        for (unsigned int y : divHeight)
        {
            unsigned int prod = x * y;
            if ((prod <= maxT) && (prod > threads))
            {
                threads = prod;
                tx = x;
                ty = y;
            }
        }
    }
    ThreadBlock threadB;
    threadB.threadX = tx;
    threadB.threadY = ty;

    return threadB;
}

#endif