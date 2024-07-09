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
struct intersection_tester
{
    float4 *position;
    float4 *point1;
    float4 *point2;
    float4 *point3;
    float4 *comparison;
};

__global__ void test_kernel(SimulationParams params, float3 *pos, float4 *p1, float4 *p2, float4 *p3, float3 *comp1, float3 *comp2, int test_num, int frame = 0)
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
        time[j] = (int)(stop - start);
    }
    printf("---------------------------------------------------------------\n");
    printf("Intersection report (epsilon = %.8f)\n", epsilon);
    printf("---------------------------------------------------------------\n");
    printf("method|time[cycles]|total_calculations|share hits|share misses|lowest delta|exist & correct|exist & incorrect|absent & correct|absent & incorrect\n");
    for (int i = 0; i < TEST_METHODS; i++)
    {
        printf("%i|%.0f|%i|%i|%i|%.8f|%i|%i|%i|%i\n", i, time[i], test_num,
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

#endif