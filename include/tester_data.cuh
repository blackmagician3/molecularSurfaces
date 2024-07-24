#ifndef TESTER_DATA_CUH
#define TESTER_DATA_CUH

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <surfaces.hpp>
#include <application_settings.hpp>
#include <simulation_config.hpp>

#ifndef INFINITY
#define INFINITY 1e8
#endif

#ifndef ZERO
#define ZERO 1e-8
#endif

const int TEST_METHODS = 8;

__global__ void test_kernel(SimulationParams params, float3 *pos, float4 *p1, float4 *p2, float4 *p3, float3 *comp1, float3 *comp2, int test_num, float clockRate, int frame = 0)
{

    printf("DEBUG: kernel started\n");
    float epsilon = 0.01f;
    clock_t start, stop;
    functionArgs args = {make_float4(.0f), make_float4(.0f), &params, nullptr, nullptr, nullptr, nullptr, INFINITY, true, 0, 0, frame, 0, 0, 0, 0};
    args.p_params->solver_iter_max = 100;
    float accumulated_delta = 0.0f;

    int count_hit[TEST_METHODS];              // countup if intersection point found
    int count_miss[TEST_METHODS];             // countup if no intersection point found
    int count_exist_correct[TEST_METHODS];    // countup if delta to comparison correct
    int count_exist_incorrect[TEST_METHODS];  // countup if delta to comparison not correct
    int count_absent_correct[TEST_METHODS];   // countup if intersection calculation terminates before reaching the iteration threshold
    int count_absent_incorrect[TEST_METHODS]; // countup if calculation reaches iteration threshold
    float min_delta[TEST_METHODS];            // stores the lowest delta value
    float time[TEST_METHODS];                 // for time measurement (in milliseconds)
    int count_no_comparison[TEST_METHODS];    // countup if the comparison value is 420 (set manually and unique from infinite/nan values)

    for (int i = 0; i < TEST_METHODS; i++)
    {
        count_hit[i] = 0;              // countup if intersection point found
        count_miss[i] = 0;             // countup if no intersection point found
        count_exist_correct[i] = 0;    // countup if delta to comparison correct
        count_exist_incorrect[i] = 0;  // countup if delta to comparison not correct
        count_absent_correct[i] = 0;   // countup if intersection calculation terminates before reaching the iteration threshold
        count_absent_incorrect[i] = 0; // countup if calculation reaches iteration threshold
        min_delta[i] = INFINITY;       // stores the lowest delta value
        count_no_comparison[i] = 0;    // countup if the comparison value is 420 (set manually and unique from infinite/nan values)
    }

    printf("DEBUG: in kernel > variables set\n");

    printf("DEBUG: in kernel > start testing\n");

    for (int j = 0; j < TEST_METHODS; j++)
    {
        printf("DEBUG: in kernel > start testing method %i\n", j);
        start = clock();
        for (int i = 0; i < test_num; i++)
        {
            args.ray_pos = make_float4(pos[i], .0f);
            args.calcs1 = 0;
            float3 intersection;

            switch (j)
            {
            case 0:
            {
                // if (i == 200)
                // {
                //     args.calcs2 = 42424242;
                // }

                pairFloat4 temp = intersectThreeSpheres1(p1[i], p2[i], p3[i], &args);
                float d = INFINITY;
                int finalRes = 0;
                for (int k = 0; k < 2; k++)
                {
                    float3 comp;
                    (k < 1) ? comp = comp1[i] : comp = comp2[i];
                    float dist = length(comp - make_float3(temp.values[k]));
                    if (dist < d)
                    {
                        d = dist;
                        finalRes = k;
                    }
                }

                intersection = make_float3(temp.values[finalRes]);
                break;
            }

            case 1:
                intersection = make_float3(intersectThreeSpheres2a(p1[i], p2[i], p3[i], &args));
                break;
            case 2:
                intersection = make_float3(intersectThreeSpheres2b(p1[i], p2[i], p3[i], &args));
                break;
            case 3:
                intersection = make_float3(intersectThreeSpheres2c(p1[i], p2[i], p3[i], &args));
                break;
            case 4:
                intersection = make_float3(intersectThreeSpheres2d(p1[i], p2[i], p3[i], &args));
                break;
            case 5:
                intersection = make_float3(intersectThreeSpheres2e(p1[i], p2[i], p3[i], &args));
                break;
            case 6:
                intersection = make_float3(intersectThreeSpheres2f(p1[i], p2[i], p3[i], &args));
                break;
            case 7:
                intersection = make_float3(intersectThreeSpheres3(p1[i], p2[i], p3[i], &args));
                break;
            default:
                break;
            }
            // calculate statistics
            float delta = INFINITY;
            int finalId = 0;
            for (int k = 0; k < 2; k++)
            {
                float3 comp;
                (k < 1) ? comp = (comp1[i]) : comp = (comp2[i]);
                float dist = length(comp - intersection);
                if (dist < delta)
                {
                    delta = dist;
                    finalId = k;
                }
            }
            if (delta < min_delta[j])
                min_delta[j] = delta;

            if (args.existence)
            {
                count_hit[j]++;
                (delta < epsilon) ? count_exist_correct[j]++ : count_exist_incorrect[j]++;
            }
            else
            {
                count_miss[j]++;
                (delta < epsilon) ? count_absent_correct[j]++ : count_absent_incorrect[j]++;
            }
        }
        stop = clock();

        time[j] = (float)(((long long)abs(stop - start)) / (clockRate));
    }
    printf("---------------------------------------------------------------\n");
    printf("Intersection report (epsilon = %.8f)\n", epsilon);
    printf("clockRate = %f kilohertz\n", clockRate);
    printf("---------------------------------------------------------------\n");
    printf("method|time[miliseconds]|total_calculations|share hits|share misses|lowest delta|exist & correct|exist & incorrect|absent & correct|absent & incorrect\n");
    for (int i = 0; i < TEST_METHODS; i++)
    {
        printf("%i|%.3f|%i|%i|%i|%.8f|%i|%i|%i|%i\n", i, time[i], test_num,
               count_hit[i],
               count_miss[i],
               min_delta[i],
               count_exist_correct[i],
               count_exist_incorrect[i],
               count_absent_correct[i],
               count_absent_incorrect[i]);
    }
    printf("---------------------------------------------------------------\n");
}

// if (settings->frame == 10)
//     {
//       const int TEST_NUM = 100000;
//       float3 *t_pos;
//       float4 *t_p1;
//       float4 *t_p2;
//       float4 *t_p3;
//       float3 *t_comp1;
//       float3 *t_comp2;
//       printf("DEBUG: before malloc\n");
//       t_pos = new float3[TEST_NUM];
//       t_p1 = new float4[TEST_NUM];
//       t_p2 = new float4[TEST_NUM];
//       t_p3 = new float4[TEST_NUM];
//       t_comp1 = new float3[TEST_NUM];
//       t_comp2 = new float3[TEST_NUM];

//       printf("DEBUG: after malloc\n");
//       // device pointers
//       float3 *d_pos, *d_comp1, *d_comp2;
//       float4 *d_p1, *d_p2, *d_p3;

//       // load data from csv
//       const int MAX_COLS = 21;
//       // std::ifstream file("C:/Users/lukab/OneDrive/Dokumente/repos/molecularSurfaces/Debugging/test_data.csv");
//       std::ifstream file("C:/Users/lukab/OneDrive/Dokumente/repos/molecularSurfaces/Debugging/test_data_filtered.csv");
//       if (!file.is_open())
//       {
//         printf("ERROR: test file could not be opened\n");
//       }
//       printf("DEBUG: file opened\n");
//       int row = 0;
//       std::string line;
//       while (std::getline(file, line) && row < TEST_NUM)
//       {
//         std::istringstream lineStream(line);
//         std::string cell;
//         float data_row[MAX_COLS];
//         int col = 0;

//         while (std::getline(file, line) && row < TEST_NUM)
//         {
//           std::istringstream lineStream(line);
//           std::string cell;
//           float data_row[MAX_COLS];
//           int col = 0;

//           while (std::getline(lineStream, cell, ',') && col < MAX_COLS)
//           {
//             try
//             {
//               data_row[col] = std::stof(cell);
//             }
//             catch (const std::invalid_argument &e)
//             {
//               printf("ERROR: Invalid data at row %d, col %d\n", row, col);
//               data_row[col] = 0.0f; // or handle error appropriately
//             }
//             col++;
//           }

//           if (col >= 21)
//           {
//             t_pos[row] = make_float3(data_row[0], data_row[1], data_row[2]);
//             t_p1[row] = make_float4(data_row[3], data_row[4], data_row[5], data_row[6]);
//             t_p2[row] = make_float4(data_row[7], data_row[8], data_row[9], data_row[10]);
//             t_p3[row] = make_float4(data_row[11], data_row[12], data_row[13], data_row[14]);
//             t_comp1[row] = make_float3(data_row[15], data_row[16], data_row[17]);
//             t_comp2[row] = make_float3(data_row[18], data_row[19], data_row[20]);
//             row++;
//           }
//         }
//       }

//       // close the file after read opeartion is complete
//       printf("DEBUG: before file close\n");
//       file.close();

//       printf("DEBUG: import complete\n");

//       // copy data to device
//       checkCudaErrors(cudaMalloc((void **)&(d_pos), TEST_NUM * sizeof(float3)));
//       checkCudaErrors(cudaMalloc((void **)&(d_p1), TEST_NUM * sizeof(float4)));
//       checkCudaErrors(cudaMalloc((void **)&(d_p2), TEST_NUM * sizeof(float4)));
//       checkCudaErrors(cudaMalloc((void **)&(d_p3), TEST_NUM * sizeof(float4)));
//       checkCudaErrors(cudaMalloc((void **)&(d_comp1), TEST_NUM * sizeof(float3)));
//       checkCudaErrors(cudaMalloc((void **)&(d_comp2), TEST_NUM * sizeof(float3)));

//       printf("DEBUG: after cuda malloc\n");
//       checkCudaErrors(cudaMemcpy(d_pos, t_pos, TEST_NUM * sizeof(float3), cudaMemcpyHostToDevice));
//       checkCudaErrors(cudaMemcpy(d_p1, t_p1, TEST_NUM * sizeof(float4), cudaMemcpyHostToDevice));
//       checkCudaErrors(cudaMemcpy(d_p2, t_p2, TEST_NUM * sizeof(float4), cudaMemcpyHostToDevice));
//       checkCudaErrors(cudaMemcpy(d_p3, t_p3, TEST_NUM * sizeof(float4), cudaMemcpyHostToDevice));
//       checkCudaErrors(cudaMemcpy(d_comp1, t_comp1, TEST_NUM * sizeof(float3), cudaMemcpyHostToDevice));
//       checkCudaErrors(cudaMemcpy(d_comp2, t_comp2, TEST_NUM * sizeof(float3), cudaMemcpyHostToDevice));

//       printf("DEBUG: after cuda memcpy\n");

//       printf("DEBUG: before host array deletion\n");

//       // free memory on host

//       delete[] t_pos;
//       delete[] t_p1;
//       delete[] t_p2;
//       delete[] t_p3;
//       delete[] t_comp1;
//       delete[] t_comp2;

//       int nDevices;
//       float clockRate;
//       cudaGetDeviceCount(&nDevices);
//       for (int i = 0; i < nDevices; i++)
//       {
//         cudaDeviceProp prop;
//         cudaGetDeviceProperties(&prop, i);
//         clockRate = prop.clockRate;
//       }
//       // run test
//       printf("DEBUG: kernel start\n");
//       test_kernel<<<1, 1>>>(settings->getAllHostParams(), d_pos, d_p1, d_p2, d_p3, d_comp1, d_comp2, TEST_NUM, clockRate, settings->frame);

//       // free device memory
//       cudaFree(d_pos);
//       cudaFree(d_p1);
//       cudaFree(d_p2);
//       cudaFree(d_p3);
//       cudaFree(d_comp1);
//       cudaFree(d_comp2);
//     }

#endif