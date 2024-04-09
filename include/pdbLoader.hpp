
#ifndef PDB_LOADER_HPP
#define PDB_LOADER_HPP

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream> // for file inputs
#include <vector>
#include <typeinfo>

#include <cooperative_groups.h>
#include <cub/cub.cuh>

#include <simulation_config.hpp>
#include <atomdata.cuh>

#ifndef ELEMENTS
#define ELEMENTS 112
#endif

float castElement(char element[3])
{
    unsigned int result = (element[0] - '0') + 100 * (element[1] - '0');

    return (float)result;
}
bool readFile(std::string path, float4 **molecule, uint *atom_count)
{
    std::ifstream input_file(path);
    if (!input_file)
    {
        return false;
    }
    std::vector<float4> tempMolecule;
    unsigned int atom = 0;

    std::string delimiter = " ";
    std::string word;

    printf("DEBUG: symbols:\n");
    for (std::string line; std::getline(input_file, line);)
    {
        std::string word = line.substr(0, 6);
        if ((word == "ATOM  ") || (word == "HETATM"))
        {

            // read coordinates
            float x_coord = std::stof(line.substr(30, 8));
            float y_coord = std::stof(line.substr(38, 8));
            float z_coord = std::stof(line.substr(46, 8));
            char element_symbol[3];
            std::strcpy(element_symbol, line.substr(76, 2).c_str());
            float4 data = make_float4(x_coord, y_coord, z_coord, castElement(element_symbol));
            tempMolecule.push_back(data);

            ++atom;
        }
    }

    *atom_count = atom;
    // molecule = nullptr;
    assert(molecule != NULL);
    *molecule = (float4 *)malloc(atom * sizeof(float4));

    std::copy(tempMolecule.begin(), tempMolecule.end(), *molecule);

    return true;
}

__device__ int castVoxelToIdTemp(int4 voxel_dim, int voxel_x, int voxel_y, int voxel_z)
{
    int id = voxel_x + (voxel_y * voxel_dim.x) + (voxel_z * voxel_dim.x * voxel_dim.y);

    return id;
}

__device__ float4 testHexColor(uint hexColor)
{
    // convert hex color to rgb
    float4 color;
    color.x = ((hexColor >> 16) & 0xFF) / 255.0f; // Extract the RR byte
    color.y = ((hexColor >> 8) & 0xFF) / 255.0f;  // Extract the GG byte
    color.z = ((hexColor) & 0xFF) / 255.0f;       // Extract the BB byte
    color.w = 1.0f;
    return color;
}

unsigned int nextPow2(unsigned int x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

void getNumBlocksAndThreads(int n, int maxBlocks, int maxThreads, int &blocks, int &threads)
{
    // get device capability, to avoid block/grid size exceed the upper bound
    cudaDeviceProp prop;
    int device;
    checkCudaErrors(cudaGetDevice(&device));
    checkCudaErrors(cudaGetDeviceProperties(&prop, device));

    threads = (n < maxThreads * 2) ? nextPow2((n + 1) / 2) : maxThreads;
    blocks = (n + (threads * 2 - 1)) / (threads * 2);

    if ((float)threads * blocks > (float)prop.maxGridSize[0] * prop.maxThreadsPerBlock)
    {
        printf("ERROR: molecule size is too large!\n");
    }

    if (blocks > prop.maxGridSize[0])
    {
        printf("Grid size <%d> exceeds the device capability <%d>, block size is set as %d (original %d)\n",
               blocks, prop.maxGridSize[0], threads * 2, threads);

        blocks /= 2;
        threads *= 2;
    }
}

/**
 * @brief kernel to add colour information and radius
 * By default, the van-der-Waals radius is used ( required by the problem definition)
 *
 * @param molecule
 * @param colour
 * @param params
 * @return __global__
 */
__global__ void supplement(float4 *molecule, uint *colour, SimulationParams params)
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id < (params.numAtoms))
    {
        unsigned int symbol = (unsigned int)molecule[id].w;
        findEntry(symbol, params.colorScheme, id, molecule, colour);
    }
}

__global__ void measureGrid(float4 *molecule, float4 *min, float4 *max, int offset, SimulationParams params)
{
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * (blockDim.x * 2) + threadIdx.x;
    extern __shared__ float4 sbuffer[];

    float4 *tempMin = &sbuffer[0];
    float4 *tempMax = &sbuffer[offset];

    float4 redMin, redMax;

    if (i < params.numAtoms)
    {
        float4 adj = make_float4(molecule[i].w, molecule[i].w, molecule[i].w, 0);
        redMin = molecule[i] - adj;
        redMax = molecule[i] + adj;
        molecule[i].w += params.solvent_radius;
        // printf("%i|%.3f|%.3f|%.3f|%.3f\n", i, molecule[i].x, molecule[i].y,
        //        molecule[i].z, molecule[i].w);
    }
    else
    {
        redMin = make_float4(INFINITY, INFINITY, INFINITY, INFINITY);
        redMax = make_float4(-INFINITY, -INFINITY, -INFINITY, 0);
    }

    if (i + blockDim.x < params.numAtoms)
    {
        float4 adj = make_float4(molecule[i + blockDim.x].w, molecule[i + blockDim.x].w, molecule[i + blockDim.x].w, 0);

        redMin = fminf(redMin, molecule[i + blockDim.x] - adj);
        redMax = fmaxf(redMax, molecule[i + blockDim.x] + adj);
        molecule[i + blockDim.x].w += params.solvent_radius;
    }
    tempMin[tid] = redMin;
    tempMax[tid] = redMax;
    __syncthreads();

    // do reduction in shared mem
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s)
        {
            tempMin[tid] = redMin = fminf(redMin, tempMin[tid + s]);
            tempMax[tid] = redMax = fmaxf(redMax, tempMax[tid + s]);
        }

        __syncthreads();
    }
    // write result for this block to global mem
    if (tid == 0)
    {
        *min = redMin;
        *max = redMax;
    }
}
template <typename T>
__global__ void initArray(T *array, int n, T value)
{
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n)
        array[i] = value;
}
__global__ void countVoxels(float4 *molecule, int *voxel_count, SimulationParams params)
{
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < params.numAtoms)
    {
        int4 current_voxel = castf2i(floorf((molecule[i] - params.box_start) / params.voxel_size));
        int id = castVoxelToIdTemp(params.voxel_dim, current_voxel.x, current_voxel.y, current_voxel.z);
        atomicAdd(&(voxel_count[id]), 1);
    }
}
__global__ void countPerVoxel(int *voxel_count, int *max_count, int n)
{
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * (blockDim.x * 2) + threadIdx.x;
    extern __shared__ int sVoxelCount[];

    int maximum = 0;
    if (i < n)
        maximum = voxel_count[i];

    // if (i < n)
    // {
    //     printf("%i|%i\n", i, voxel_count[i]);
    // }

    if (i + blockDim.x < n)
        maximum = max(maximum, voxel_count[i + blockDim.x]);

    // if (i + blockDim.x < n)
    // {
    //     printf("%i|%i\n", i + blockDim.x, voxel_count[i + blockDim.x]);
    // }

    sVoxelCount[tid] = maximum;
    __syncthreads();

    // do reduction in shared mem
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s)
            sVoxelCount[tid] = maximum = max(maximum, sVoxelCount[tid + s]);

        __syncthreads();
    }
    // write result for this block to global mem
    if (tid == 0)
        *max_count = maximum;
}
__global__ void initializeGrid(float4 *molecule, int *voxel_data, int *voxel_count, SimulationParams params)
{
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < params.numAtoms)
    {

        int4 current_voxel = castf2i(floorf((molecule[i] - params.box_start) / params.voxel_size));
        int id = castVoxelToIdTemp(params.voxel_dim, current_voxel.x, current_voxel.y, current_voxel.z);
        int current_count = atomicAdd(&(voxel_count[id]), -1) - 1;
        voxel_data[id * params.atomsPerVoxel + current_count] = i;
    }
}

__global__ void printDeviceArray(int *array, int n, int startOffset = 0, int valuesOffset = 0)
{
    unsigned int start = 0;
    unsigned int end = n;
    if (startOffset > 0 || valuesOffset > 0)
    {
        (startOffset < n - 1) ? (start = startOffset) : (start = 0);
        (start + valuesOffset < n - 1) ? (end = start + valuesOffset) : (end = n);
    }

    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    printf("id|array.value\n");
    if (i == 0)
    {
        for (unsigned int l = start; l < end; l++)
            printf("%i|%i\n", l, array[l]);
    }
}
__global__ void printDeviceArray(float *array, int n, int startOffset = 0, int valuesOffset = 0)
{
    unsigned int start = 0;
    unsigned int end = n;
    if (startOffset > 0 && valuesOffset > 0)
    {
        (startOffset < n - 1) ? (start = startOffset) : (start = 0);
        (start + valuesOffset < n - 1) ? (end = start + valuesOffset) : (end = n);
    }
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    printf("id|array.value\n");
    if (i == 0)
    {
        for (unsigned int l = start; l < end; l++)
            printf("%i|%.3f\n", l, array[l]);
    }
}
__global__ void printDeviceArray(int4 *array, int n, int startOffset = 0, int valuesOffset = 0)
{
    unsigned int start = 0;
    unsigned int end = n;
    if (startOffset > 0 && valuesOffset > 0)
    {
        (startOffset < n - 1) ? (start = startOffset) : (start = 0);
        (start + valuesOffset < n - 1) ? (end = start + valuesOffset) : (end = n);
    }
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    printf("id|array.x|array.y|array.z|array.w\n");
    if (i == 0)
    {
        for (unsigned int l = start; l < end; l++)
            printf("%i|%i|%i|%i|%i\n", l, array[l].x, array[l].y, array[l].z, array[l].w);
    }
}
__global__ void printDeviceArray(float4 *array, int n, int startOffset = 0, int valuesOffset = 0)
{
    unsigned int start = 0;
    unsigned int end = n;
    if (startOffset > 0 && valuesOffset > 0)
    {
        (startOffset < n - 1) ? (start = startOffset) : (start = 0);
        (start + valuesOffset < n - 1) ? (end = start + valuesOffset) : (end = n);
    }
    unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;

    printf("id|array.x|array.y|array.z|array.w\n");
    if (i == 0)
    {
        for (unsigned int l = start; l < end; l++)
            printf("%i|%.3f|%.3f|%.3f|%.3f\n", l, array[l].x, array[l].y, array[l].z, array[l].w);
    }
}
template <typename T>
__global__ void printDeviceArray(T *array, int n, int startOffset = 0, int valuesOffset = 0)
{
    printf("INFO: the argument types are not intended for this method\n");
}

#endif