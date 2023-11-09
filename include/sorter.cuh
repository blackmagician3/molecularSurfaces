#ifndef SORTER_CUH
#define SORTER_CUH

#include <stdio.h>

typedef struct
{
    float4 location;
    int id;
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
            if (array[j].location.w < array[min_idx].location.w)
                min_idx = j;

        // Swap the found minimum element with the first element
        swap(&array[min_idx], &array[i]);

    }
}

#endif