#ifndef SORTER_CUH
#define SORTER_CUH

#include <stdio.h>

typedef struct
{
    float4 location;
    int id;
    float distance;
} Atom;

__host__ __device__ void swap(Atom *xp, Atom *yp)
{
    Atom temp = *xp;
    *xp = *yp;
    *yp = temp;
}

__host__ __device__ void swap(float *xp, float *yp)
{
    float temp = *xp;
    *xp = *yp;
    *yp = temp;
}

__host__ __device__ void swap(int *xp, int *yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

__host__ __device__ void swap(unsigned int *xp, unsigned int *yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

__host__ __device__ int findFurthestByIndex(Atom array[], int size)
{
    int result = 0;
    float dist_max = array[0].distance;
    for (int i = 1; i < size; i++)
    {
        if (array[i].distance > dist_max)
        {
            dist_max = array[i].distance;
            result = i;
        }
    }
    return result;
}

// Function to perform Selection Sort
__host__ __device__ void selectionSort(int n, Atom array[])
{
    int min_idx;

    // One by one move boundary of unsorted subarray
    for (int i = 0; i < n - 1; i++)
    {
        // Find the minimum element in unsorted array
        min_idx = i;
        for (int j = i + 1; j < n; j++)
            if (array[j].distance < array[min_idx].distance)
                min_idx = j;

        // Swap the found minimum element with the first element
        swap(&array[min_idx], &array[i]);
    }
}
__host__ __device__ void sortVector(int length, float array[])
{
    if (length != 3)
        return;
    if (array[0] > array[2])
        swap(&array[0], &array[2]);
    if (array[0] > array[1])
        swap(&array[0], &array[1]);
    if (array[1] > array[2])
        swap(&array[1], &array[2]);
}

__host__ __device__ void sortParallelVector(int length, float array[], int ids[])
{
    if (length != 3)
        return;
    if (array[0] > array[2])
    {
        swap(&array[0], &array[2]);
        swap(&ids[0], &ids[2]);
    }

    if (array[0] > array[1])
    {
        swap(&array[0], &array[1]);
        swap(&ids[0], &ids[1]);
    }

    if (array[1] > array[2])
    {
        swap(&array[1], &array[2]);
        swap(&ids[1], &ids[2]);
    }
}

#endif