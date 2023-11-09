/**
 * @file application_settings.hpp
 * @author Luka Blümler (luka.bluemler@uni-jena.de)
 * @brief class for application settings containing all constants/variables for both CPU and GPU
 */

#ifndef APPLICATION_SETTINGS_HPP
#define APPLICATION_SETTINGS_HPP

#include <raymarching_kernel.cuh>
#include <simulation_config.hpp>

typedef unsigned int uint;

class AppSettings
{
private:
    bool initialized;
    bool molecule_loaded;

    // parameters for both CPU and GPU
    SimulationParams host_params;

public:
    // molecule data
    float4 *molecule;
    float4 *molecule_device;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // functions
    AppSettings(uint width, uint height, uint texture) : initialized(true), molecule_loaded(false)
    {
        host_params.window_width = width;
        host_params.window_height = height;

        host_params.numAtoms = 10;
        host_params.k_nearest = 10;
        host_params.solvent_radius = 0.2;
        host_params.epsilon = 0.01;
        host_params.epsilon_squared = host_params.epsilon * host_params.epsilon;
        host_params.depth_min = 0.1f;
        host_params.depth_max = 100.0f;
        host_params.is_SES = true;

        host_params.texture_1 = texture;

        host_params.thread_x = 4;
        host_params.thread_y = 4;

        molecule = nullptr;
        molecule_device = nullptr;

        host_params.debug_mode = false;

        initialize();
    }
    void initialize()
    {
        // TODO: currently using registers for nearest atoms grouping per thread -> register size is shared among all threads per block -> adjust threads per block
        ThreadBlock threads = calcThreadBlock(host_params.window_width, host_params.window_height);
        // host_params.thread_x = threads.threadX;
        // host_params.thread_y = threads.threadY;

        host_params.thread_x = 4;
        host_params.thread_y = 4;
        initialized = true;
    }
    void loadMolecule()
    {
        if (molecule_loaded)
        {
            delete[] molecule;
        }

        molecule = new float4[host_params.numAtoms];

        molecule[0] = make_float4(-0.4, 0.0, 0.0, 0.4);
        molecule[1] = make_float4(0.5, 0.0, 0.0, 0.6);
        molecule[2] = make_float4(0.7, 0.7, 0.0, 0.8);
        molecule[3] = make_float4(-0.1, 0.7, -0.2, 0.5);
        molecule[4] = make_float4(-0.3, 0.2, -2.0, 1.0);
        molecule[5] = make_float4(0.1, 0.7, 0.3, 0.2);
        molecule[6] = make_float4(0.9, -0.2, 0.3, 0.1);
        molecule[7] = make_float4(-0.1, 0.8, -3, 1.1);
        molecule[8] = make_float4(-0.6, -0.7, 0.1, 0.6);
        molecule[9] = make_float4(0.1, -0.2, -0.3, 0.8);

        // extend atom radi by solvent radius
        addToRadius(molecule, host_params.numAtoms, host_params.solvent_radius);

        // copy data to GPU
        checkCudaErrors(cudaMalloc((void **)&(molecule_device), host_params.numAtoms * sizeof(float4)));
        checkCudaErrors(cudaMemcpy(molecule_device, molecule, host_params.numAtoms * sizeof(float4), cudaMemcpyHostToDevice));

        // TODO: implement loader, remove test molecule
        molecule_loaded = true;
    }
    void update()
    {
        setParameters(host_params);
    }

    void reset()
    {
        host_params.window_width = 1920;
        host_params.window_height = 1080;

        host_params.numAtoms = 10;
        host_params.solvent_radius = 0.1;
        host_params.epsilon = 0.001;
        host_params.epsilon_squared = host_params.epsilon * host_params.epsilon;

        initialize();
        if (molecule_loaded)
        {
            checkCudaErrors(cudaFree(molecule_device));
            delete[] molecule;
        }

        molecule_loaded = false;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // destructor
    ~AppSettings()
    {
        if (molecule_loaded)
        {
            cudaDeviceSynchronize();
            delete[] molecule;

            checkCudaErrors(cudaFree((float4 *)molecule_device));

            molecule = nullptr;
            molecule_device = nullptr;
        }
    }
    // copy constructor
    AppSettings(const AppSettings &obj)
    {
        initialized = true;
        host_params = obj.host_params;

        molecule_loaded = obj.molecule_loaded;

        // initialize();

        if (obj.molecule_loaded)
        {
            molecule = new float4[obj.host_params.numAtoms];
            for (int i = 0; i < obj.host_params.numAtoms; i++)
            {
                molecule[i] = obj.molecule[i];
            }
            molecule_device = obj.molecule_device;
        }
    }
    // assignment operator
    AppSettings &AppSettings::operator=(AppSettings &obj)
    {
        std::swap(host_params, obj.host_params);
        std::swap(molecule, obj.molecule);
        std::swap(molecule_device, obj.molecule_device);
        std::swap(molecule_loaded, obj.molecule_loaded);
        initialized = true;
        return *this;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // getter
    uint getAtomCount()
    {
        return host_params.numAtoms;
    }

    uint getThreadX()
    {
        return host_params.thread_x;
    }
    uint getThreadY()
    {
        return host_params.thread_y;
    }

    uint getWindowWidth()
    {
        return host_params.window_width;
    }
    uint getWindowHeight()
    {
        return host_params.window_height;
    }

    uint getTexture1()
    {
        return host_params.texture_1;
    }

    float getSolventRadius()
    {
        return host_params.solvent_radius;
    }

    float getEpsilon()
    {
        return host_params.epsilon;
    }
    float getEpsilonSquared()
    {
        return host_params.epsilon_squared;
    }
    float4 *getDeviceMolecule()
    {
        return molecule_device;
    }

    bool getDebugMode()
    {
        return host_params.debug_mode;
    }

    SimulationParams getAllHostParams()
    {
        return host_params;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // setter

    void addToSolventRadius(float value)
    {
        host_params.solvent_radius += value;
    }

    void changeDebugMode(bool value)
    {
        host_params.debug_mode = value;
    }
    void setEpsilon(float value)
    {
        host_params.epsilon = value;
        host_params.epsilon_squared = value * value;
    }
    void setMousePos(double value_x, double value_y)
    {
        host_params.mouse_x_pos = value_x;
        host_params.mouse_y_pos = value_y;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
};
#endif