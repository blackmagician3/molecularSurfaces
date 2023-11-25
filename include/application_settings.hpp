/**
 * @file application_settings.hpp
 * @author Luka Blümler (luka.bluemler@uni-jena.de)
 * @brief class for application settings containing all constants/variables for both CPU and GPU
 */

#ifndef APPLICATION_SETTINGS_HPP
#define APPLICATION_SETTINGS_HPP

#include <raymarching_kernel.cuh>
#include <simulation_config.hpp>
#include <string>       // for file/path operations
#include <fstream>      // for file inputs
#include <atomdata.hpp> // contains data on atom radii and color schemes

#ifndef ELEMENTS
#define ELEMENTS 112
#endif

typedef unsigned int uint;

class AppSettings
{
private:
    bool initialized;
    bool molecule_loaded;
    bool performanceDisplay;
    unsigned int frameLimit;
    std::string recentFilePath;

    // parameters for both CPU and GPU
    SimulationParams host_params;

public:
    // molecule data
    float4 *molecule;
    float4 *molecule_device;
    uint *colors;
    uint *colors_device;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // functions
    AppSettings(uint width, uint height, uint texture) : initialized(true), molecule_loaded(false), performanceDisplay(false)
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

        host_params.colorScheme = 0;
        host_params.debug_mode = false;

        frameLimit = 0;
        recentFilePath = "";

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
    void loadMolecule(std::string path = "")
    {
        // reset color scheme and delete privious molecule
        if (molecule_loaded)
        {
            delete[] molecule;
            delete[] colors;
        }

        // default/test molecule
        if (path.length() > 0)
        {
            recentFilePath = path;
        }
        if (path.length() == 0 && recentFilePath.length() == 0)
        {
            host_params.numAtoms = 10;
            host_params.k_nearest = 10;

            molecule = new float4[host_params.numAtoms];
            molecule[0] = make_float4(-0.4, 0.0, 0.0, 0.2);
            molecule[1] = make_float4(0.5, 0.0, 0.0, 0.3);
            molecule[2] = make_float4(0.7, 0.7, 0.0, 0.4);
            molecule[3] = make_float4(-0.1, 0.7, -0.2, 0.25);
            molecule[4] = make_float4(-0.3, 0.2, -2.0, 0.5);
            molecule[5] = make_float4(0.1, 0.7, 0.3, 0.1);
            molecule[6] = make_float4(0.9, -0.2, 0.3, 0.1);
            molecule[7] = make_float4(-0.1, 0.8, -3, 0.6);
            molecule[8] = make_float4(-0.6, -0.7, 0.1, 0.3);
            molecule[9] = make_float4(0.1, -0.2, -0.3, 0.4);

            colors = new uint[host_params.numAtoms];
            colors[0] = 0x2060ff;
            colors[1] = 0x202020;
            colors[2] = 0xee2010;
            colors[3] = 0x2060ff;
            colors[4] = 0x202020;
            colors[5] = 0xffffff;
            colors[6] = 0xffffff;
            colors[7] = 0x2060ff;
            colors[8] = 0xee2010;
            colors[9] = 0x202020;

            // extend atom radi by solvent radius
            addToRadius(molecule, host_params.numAtoms, host_params.solvent_radius);
        }
        else
        {
            // open file
            std::ifstream input_file(recentFilePath);
            if (!input_file)
            {
                std::cout << "Error reading file. Loading default molecule instead." << std::endl;
                loadMolecule();
                return;
            }

            // determine molecule size
            unsigned int atom = 0;
            unsigned int hetatm = 0;
            std::string delimiter = " ";
            std::string word;
            for (std::string tmp; std::getline(input_file, tmp);)
            {
                std::string word = tmp.substr(0, 6);
                if (word == "ATOM  ")
                {
                    ++atom;
                }
                else if (word == "HETATM")
                {
                    ++hetatm;
                }
            }
            host_params.numAtoms = atom + hetatm;
            host_params.k_nearest = host_params.numAtoms;

            molecule = new float4[host_params.numAtoms];
            colors = new uint[host_params.numAtoms];

            // read file
            std::string keyword;
            // std::string atom_name;
            std::string element_symbol;
            unsigned int n = 0;
            for (std::string line; getline(input_file, line);)
            {
                keyword = line.substr(0, 6);
                if (keyword == "ATOM  " || keyword == "HETATM")
                {
                    element_symbol = line.substr(76, 2);
                    float radius;
                    uint color;
                    findEntry(element_symbol, host_params.colorScheme, &radius, &color);

                    molecule[n] = make_float4(std::stof(line.substr(30, 8)),        // x - Coordinates
                                              std::stof(line.substr(38, 8)),        // y - Coordinates
                                              std::stof(line.substr(46, 8)),        // z - Coordinates
                                              radius + host_params.solvent_radius); // atom radius
                    colors[n] = color;
                }
                ++n;
            }
        }

        // copy data to GPU
        checkCudaErrors(cudaMalloc((void **)&(molecule_device), host_params.numAtoms * sizeof(float4)));
        checkCudaErrors(cudaMemcpy(molecule_device, molecule, host_params.numAtoms * sizeof(float4), cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMalloc((void **)&(colors_device), host_params.numAtoms * sizeof(uint)));
        checkCudaErrors(cudaMemcpy(colors_device, colors, host_params.numAtoms * sizeof(uint), cudaMemcpyHostToDevice));

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
            checkCudaErrors(cudaFree(colors_device));
            delete[] molecule;
            delete[] colors;
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
            delete[] colors;
            checkCudaErrors(cudaFree((float4 *)molecule_device));
            checkCudaErrors(cudaFree((uint *)colors_device));
            molecule = nullptr;
            molecule_device = nullptr;
            colors = nullptr;
            colors_device = nullptr;
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
            colors = new uint[obj.host_params.numAtoms];
            for (int i = 0; i < obj.host_params.numAtoms; i++)
            {
                molecule[i] = obj.molecule[i];
                colors[i] = obj.colors[i];
            }
            molecule_device = obj.molecule_device;
            colors_device = obj.colors_device;
        }
    }
    // assignment operator
    AppSettings &AppSettings::operator=(AppSettings &obj)
    {
        std::swap(host_params, obj.host_params);
        std::swap(molecule, obj.molecule);
        std::swap(molecule_device, obj.molecule_device);
        std::swap(molecule_loaded, obj.molecule_loaded);

        std::swap(colors, obj.colors);
        std::swap(colors_device, obj.colors_device);
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
    uint *getDeviceColors()
    {
        return colors_device;
    }

    bool getDebugMode()
    {
        return host_params.debug_mode;
    }

    SimulationParams getAllHostParams()
    {
        return host_params;
    }
    bool getPerformanceDisplay()
    {
        return performanceDisplay;
    }
    unsigned int getFrameLimit()
    {
        return frameLimit;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // setter

    void addToSolventRadius(float value)
    {
        host_params.solvent_radius += value;
        addToRadius(molecule, host_params.numAtoms, value);
        checkCudaErrors(cudaMemcpy(molecule_device, molecule, host_params.numAtoms * sizeof(float4), cudaMemcpyHostToDevice));
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
    void changePerformanceDisplay(bool state)
    {
        performanceDisplay = state;
    }
    void changeFrameLimit(unsigned int limit)
    {
        frameLimit = limit;
    }

    /**
     * @brief Set the color scheme for the molecule.
     * !!! Forces a reload of the molecule data from most recent file path !!!
     *
     * @param value
     */
    void setColorScheme(unsigned int value)
    {
        // schemes
        // 0 - case sensitive with default colors
        // 1 - Corey
        // 2 - Koltun
        // 3 - Jmol
        // 4 - Rasmol old (prior version 2.7.3)
        // 5 - Rasmol new (version 2.7.3 and later)
        // 6 - PubChem

        // error handling
        if (value > 6)
        {
            printf("color scheme %i is not defined. Available schemes are:\n", value);
            printf("0 - case sensitive with default colors\n");
            printf("1 - Corey\n");
            printf("2 - Koltun\n");
            printf("3 - Jmol\n");
            printf("4 - Rasmol (prior version 2.7.3)\n");
            printf("5 - Rasmol  (version 2.7.3 and later)\n");
            printf("6 - PubChem\n");
        }

        // set scheme
        host_params.colorScheme = value;

        // reload molecule
        loadMolecule();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
};
#endif