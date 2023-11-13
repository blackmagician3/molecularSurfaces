#ifndef PERFORMANCE_HPP
#define PERFORMANCE_HPP

#include <stdio.h>

#ifndef INFINITY
#define INFINITY 1e8
#endif

/**
 * @brief class to calculate performance metrics over a given time interval.
 * Measurements will be restarted every new intervall
 *
 */
class PerformanceCounter
{
public:
    //////////////////////////////////////////////////////////////////////////////////////////////
    // constructor
    PerformanceCounter(unsigned int measurement_duration_in_sec)
    {
        measurementDuration = measurement_duration_in_sec;
        currentFPS = 0.0f;
        accumulatedTime = 0.0f;
        accumulatedFPS = 0.0f;
        currentIntervall = 0.0f;
        minFPS = INFINITY;
        maxFPS = 0.0f;
        framesPerSecond = 0.0f;
    }

    ~PerformanceCounter()
    {
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // functions

    void runMeasuring(float deltaTime)
    {
        accumulatedTime += deltaTime;

        calculateFrameRate(deltaTime);

        if (accumulatedTime >= (float)measurementDuration)
        {
            displayPerformanceMetrics();
            accumulatedTime = 0.0f;
            accumulatedFPS = 0.0f;
            minFPS = INFINITY;
            maxFPS = 0.0f;
        }
    }

    void changePerformanceIntervall(unsigned int measurement_duration_in_sec)
    {
        measurementDuration = measurement_duration_in_sec;
        accumulatedTime = 0.0f;
        accumulatedFPS = 0.0f;
        currentIntervall = 0.0f;
        currentFPS = 0.0f;
        minFPS = INFINITY;
        maxFPS = 0.0f;
        framesPerSecond = 0.0f;
    }

private:
    unsigned int measurementDuration;
    float accumulatedTime;
    float accumulatedFPS;
    float currentIntervall;
    float currentFPS;
    float framesPerSecond;
    float minFPS;
    float maxFPS;

    void calculateFrameRate(float deltaTime)
    {
        currentIntervall += deltaTime;
        ++framesPerSecond;

        if (currentIntervall >= 1.0f)
        {
            currentIntervall = 0.0f;
            currentFPS = framesPerSecond;
            accumulatedFPS += framesPerSecond;
            if (framesPerSecond < minFPS)
                minFPS = framesPerSecond;
            if (framesPerSecond > maxFPS)
                maxFPS = framesPerSecond;
            framesPerSecond = 0;
        }
    }

    void displayPerformanceMetrics()
    {
        float averageFPS = accumulatedFPS / accumulatedTime;
        printf("-------------- performance display --------------\n");
        printf("average frames per second: %.2f\n", averageFPS);
        printf("minimum frames per second: %.2f\n", minFPS);
        printf("maximum frames per second: %.2f\n", maxFPS);
        printf("over a time interval of %i seconds\n", measurementDuration);
        printf("-------------------------------------------------\n");
    }
};

#endif