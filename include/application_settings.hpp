/**
 * @file application_settings.hpp
 * @author Luka Blümler (luka.bluemler@uni-jena.de)
 * @brief class for application settings containing all constants/variables for both CPU and GPU
 */

#ifndef APPLICATION_SETTINGS_HPP
#define APPLICATION_SETTINGS_HPP

#include <raymarching_kernel.cuh>
#include <simulation_config.hpp>
#include <string>  // for file/path operations
#include <fstream> // for file inputs
#include <vector>

#include <pdbLoader.hpp>
// #include <atomdata.hpp> // contains data on atom radii and color schemes

#ifndef ELEMENTS
#define ELEMENTS 112
#endif
#ifndef MAXTHREADS
#define MAXTHREADS 512
#endif
#ifndef MAXBLOCKS
#define MAXBLOCKS 65535
#endif

// Default camera value
const float CAMOFFSET = 0;
typedef unsigned int uint;
extern const int GROUP_SIZE;

class AppSettings
{
private:
	bool initialized;
	bool molecule_loaded;
	std::string recentFilePath;
	int voxel_mem;
	float cam_offset_init;

	// parameters for initialization
	float4 focusPoint;

	// parameters for both CPU and GPU
	SimulationParams host_params;

	int castVoxelToId(int4 voxel_dim, int4 voxel_id)
	{
		int id = (voxel_id.x + voxel_dim.x * (voxel_id.y + voxel_dim.y * voxel_id.z));
		return id;
	}

public:
	// molecule data
	float4 *molecule;
	float4 *molecule_device;
	uint *colors_device;
	int *voxel_data_device;
	int *voxel_count_device;

	// parameters for gui
	bool show_gui;
	bool show_performance_tool;
	bool limit_framerate;
	int frameLimit;
	int frame;

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// functions
	AppSettings(uint width, uint height, uint texture) : initialized(true), molecule_loaded(false)
	{
		host_params.window_width = width;
		host_params.window_height = height;

		host_params.numAtoms = 10;
		host_params.k_nearest = 15;
		host_params.solvent_radius = 1.4;
		host_params.solvent_max = 2;
		host_params.epsilon = 0.01;
		host_params.epsilon_squared = host_params.epsilon * host_params.epsilon;
		host_params.use_voxel = true;
		host_params.depth_min = 0.1f;
		host_params.depth_max = 100.0f;
		host_params.is_SES = true;
		host_params.solver = 0;
		host_params.texture_1 = texture;

		host_params.thread_x = 4;
		host_params.thread_y = 4;

		molecule = nullptr;
		molecule_device = nullptr;

		host_params.colorScheme = 0;
		host_params.debug_mode = false;

		host_params.highlight_radius = 5.0f;

		frameLimit = 60;
		recentFilePath = "";

		show_gui = true;
		show_performance_tool = true;
		limit_framerate = false;
		frame = 0;

		initialize();
	}
	void initialize()
	{
		// TODO: currently using registers for nearest atoms grouping per thread -> register size is shared among all threads per block -> adjust threads per block
		// TODO: calculated threads per block currently to high
		ThreadBlock threads = calcThreadBlock(host_params.window_width, host_params.window_height);
		host_params.thread_x = threads.threadX;
		host_params.thread_y = threads.threadY;

		// printf("INFO: kernel layout set to %i threads per block (%i in x, %i in y)\n", host_params.thread_x * host_params.thread_y, host_params.thread_x, host_params.thread_y);
		initialized = true;
	}
	void loadMolecule(std::string path = "")
	{
		// reset color scheme and delete privious molecule
		if (molecule_loaded)
		{
			delete[] molecule;
		}

		// default/test molecule
		if (path.length() > 0)
		{
			recentFilePath = path;
		}
		if (path.length() == 0 && recentFilePath.length() == 0)
		{
			printf("ERROR: no path to molecule provided.\n");
		}
		else
		{
			// Error code to check return values for CUDA calls
			cudaError_t err = cudaSuccess;

			err = cudaGetLastError();
			if (err != cudaSuccess)
			{
				printf("ERROR: Failed at very start (error code %s)!\n", cudaGetErrorString(err));
				exit(EXIT_FAILURE);
			}

			// open file
			// molecule = (float4 *)malloc(446 * sizeof(float4));
			if (readFile(recentFilePath, &molecule, &host_params.numAtoms))
			{
				// printf("INFO: molecule loaded from path, size = %i atoms\n", host_params.numAtoms);
			}
			else
			{
				printf("ERROR: file could not be read. Loading default molecule instead.\n");

				delete[] molecule;
				loadMolecule();
			}

			//  set default solvent radius
			host_params.solvent_radius = 1.4;
			////////////////////////////////////////////////////////////////////////
			// add radius and colour information
			////////////////////////////////////////////////////////////////////////
			checkCudaErrors(cudaMalloc((void **)&(molecule_device), host_params.numAtoms * sizeof(float4)));
			checkCudaErrors(cudaMemcpy(molecule_device, molecule, host_params.numAtoms * 4 * sizeof(float), cudaMemcpyHostToDevice));
			checkCudaErrors(cudaMalloc((void **)&(colors_device), host_params.numAtoms * sizeof(uint)));

			int blocks = 0;
			int threads = 0;
			getNumBlocksAndThreads(host_params.numAtoms, MAXBLOCKS, MAXTHREADS, blocks, threads);
			supplement<<<2 * blocks, threads>>>(molecule_device, colors_device, host_params);
			err = cudaGetLastError();
			if (err != cudaSuccess)
			{
				printf("ERROR: failed to add radius & colour information (error code %s)!\n", cudaGetErrorString(err));
				exit(EXIT_FAILURE);
			}
			cudaDeviceSynchronize();

			// //////////////////////////////////////////////////////////////////////
			// initialize grid
			// //////////////////////////////////////////////////////////////////////
			if (host_params.use_voxel)
			{
				int smemSize = (threads <= 32) ? 2 * threads * sizeof(float4) : threads * sizeof(float4);
				smemSize *= 2;
				int offset = ((threads <= 32) ? 2 * threads : threads);

				float4 *blockMin, *blockMax;
				checkCudaErrors(cudaMallocManaged((void **)&(blockMin), blocks * sizeof(float4)));
				checkCudaErrors(cudaMallocManaged((void **)&(blockMax), blocks * sizeof(float4)));

				measureGrid<<<blocks, threads, smemSize>>>(molecule_device, blockMin, blockMax, offset, host_params);
				err = cudaGetLastError();
				if (err != cudaSuccess)
				{
					printf("ERROR: grid measurements were not successful (error code %s)!\n", cudaGetErrorString(err));
					exit(EXIT_FAILURE);
				}
				cudaDeviceSynchronize();

				float4 bMin = blockMin[0];
				float4 bMax = blockMax[0];

				// for multiple blocks reduce the min/max ararys on the host
				for (int i = 1; i < blocks; i++)
				{
					bMin = fminf(bMin, blockMin[i]);
					bMax = fmaxf(bMax, blockMax[i]);
				}

				host_params.box_start = bMin - host_params.solvent_max;
				host_params.box_end = bMax + host_params.solvent_max;

				err = cudaGetLastError();
				if (err != cudaSuccess)
				{
					printf("ERROR: Failed to launch supplement kernel (error code %s)!\n", cudaGetErrorString(err));
					exit(EXIT_FAILURE);
				}
				checkCudaErrors(cudaFree(blockMin));
				checkCudaErrors(cudaFree(blockMax));

				// count number of atoms for each voxel
				host_params.voxel_size = 2 * fabs(host_params.box_end.w + host_params.solvent_max);
				host_params.voxel_dim = castf2i(ceilf((1 / host_params.voxel_size) * fabs(host_params.box_start - host_params.box_end)));
				voxel_mem = host_params.voxel_dim.x * host_params.voxel_dim.y * host_params.voxel_dim.z;

				checkCudaErrors(cudaMalloc((void **)&(voxel_count_device), voxel_mem * sizeof(int)));

				getNumBlocksAndThreads(voxel_mem, MAXBLOCKS, MAXTHREADS, blocks, threads);
				initArray<<<2 * blocks, threads>>>(voxel_count_device, voxel_mem, 0);
				cudaDeviceSynchronize();
				getNumBlocksAndThreads(host_params.numAtoms, MAXBLOCKS, MAXTHREADS, blocks, threads);
				countVoxels<<<2 * blocks, threads>>>(molecule_device, voxel_count_device, host_params);
				err = cudaGetLastError();
				if (err != cudaSuccess)
				{
					printf("ERROR: Failed to determine atoms per voxel (error code %s)!\n", cudaGetErrorString(err));
					exit(EXIT_FAILURE);
				}
				cudaDeviceSynchronize();

				// determine maximum number of atoms per voxel
				int maxVoxelCount;
				int *d_maxVoxelCount;
				checkCudaErrors(cudaMallocManaged(&d_maxVoxelCount, sizeof(int)));
				initArray<<<1, 1>>>(d_maxVoxelCount, 1, 0);
				smemSize = (threads <= 32) ? 2 * threads * sizeof(int) : threads * sizeof(int);

				getNumBlocksAndThreads(voxel_mem, MAXBLOCKS, MAXTHREADS, blocks, threads);
				countPerVoxel<<<blocks, threads, smemSize>>>(voxel_count_device, d_maxVoxelCount, voxel_mem);
				err = cudaGetLastError();
				if (err != cudaSuccess)
				{
					printf("ERROR: Failed to determine maximum number of atoms per voxel (error code %s)!\n", cudaGetErrorString(err));
					exit(EXIT_FAILURE);
				}
				cudaDeviceSynchronize();
				checkCudaErrors(cudaMemcpy(&maxVoxelCount, d_maxVoxelCount, sizeof(int), cudaMemcpyDeviceToHost));
				// int atomsPerDim = (int)(floorf(host_params.voxel_size / host_params.box_start.w));
				host_params.atomsPerVoxel = maxVoxelCount;
				// host_params.atomsPerVoxel = 100;
				checkCudaErrors(cudaFree(d_maxVoxelCount));

				checkCudaErrors(cudaMalloc((void **)&(voxel_data_device), host_params.atomsPerVoxel * voxel_mem * sizeof(int)));

				// initialize grid -> sort atom ids in voxel_data array based on ascending voxel number
				getNumBlocksAndThreads(host_params.atomsPerVoxel * voxel_mem, MAXBLOCKS, MAXTHREADS, blocks, threads);
				initArray<<<2 * blocks, threads>>>(voxel_data_device, host_params.atomsPerVoxel * voxel_mem, -1);
				cudaDeviceSynchronize();

				initializeGrid<<<2 * blocks, threads>>>(molecule_device, voxel_data_device, voxel_count_device, host_params);
				err = cudaGetLastError();
				if (err != cudaSuccess)
				{
					printf("ERROR: Failed to initialize grid (error code %s)!\n", cudaGetErrorString(err));
					exit(EXIT_FAILURE);
				}
			}

			// setup initial camera position to cover the whole molecule
			// TODO: synchronize Zoom Default Value with Camera class
			float Zoom = 10;
			focusPoint = 0.5 * (host_params.box_start + host_params.box_end);
			float zOffset = max(abs(host_params.box_end.x - host_params.box_start.x) * Zoom, abs(host_params.box_end.y - host_params.box_start.y) * Zoom);
			zOffset = max(zOffset, 0.5 * abs(host_params.box_end.z - host_params.box_start.z));
			cam_offset_init = zOffset + CAMOFFSET;
			focusPoint.w = cam_offset_init;
			host_params.depth_max = 4 * zOffset + CAMOFFSET;
		}

		update();
		molecule_loaded = true;
	}
	void update()
	{
		// test thresholds
		if (host_params.k_nearest > GROUP_SIZE)
		{
			host_params.k_nearest = GROUP_SIZE;
		}

		// synchronize CPU & GPU parameters
		setParameters(host_params);
	}

	void reset_view()
	{
		host_params.window_width = 1920;
		host_params.window_height = 1080;
		host_params.k_nearest = 10;
		host_params.solvent_radius = 1.4;
		host_params.solvent_max = 2;
		host_params.epsilon = 0.01;
		host_params.epsilon_squared = host_params.epsilon * host_params.epsilon;
		host_params.use_voxel = true;
		host_params.depth_min = 0.1f;
		host_params.is_SES = true;
		host_params.solver = 0;
		host_params.debug_mode = false;
		host_params.highlight_radius = 5.0f;
		frameLimit = 60;

		show_gui = true;
		show_performance_tool = true;
		limit_framerate = false;
		frame = 0;

		if (molecule_loaded)
		{
			focusPoint = 0.5 * (host_params.box_start + host_params.box_end);
			focusPoint.w = cam_offset_init;
		}

		update();
	}
	void reload_molecule()
	{
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// destructor
	~AppSettings()
	{
		if (molecule_loaded)
		{
			checkCudaErrors(cudaFree((float4 *)molecule_device));
			checkCudaErrors(cudaFree((uint *)colors_device));

			molecule_device = nullptr;
			colors_device = nullptr;

			if (host_params.use_voxel)
			{
				checkCudaErrors(cudaFree((int *)voxel_data_device));
				checkCudaErrors(cudaFree((int *)voxel_count_device));
				voxel_data_device = nullptr;
				voxel_count_device = nullptr;
			}

			delete[] molecule;
			molecule = nullptr;
		}
	}
	// copy constructor
	AppSettings(const AppSettings &obj)
	{
		initialized = true;
		host_params = obj.host_params;
		molecule_loaded = obj.molecule_loaded;

		if (obj.molecule_loaded)
		{
			for (int i = 0; i < host_params.numAtoms; i++)
			{
				molecule[i] = obj.molecule[i];
			}

			molecule_device = obj.molecule_device;

			colors_device = obj.colors_device;
			voxel_data_device = obj.voxel_data_device;
			voxel_count_device = obj.voxel_count_device;
		}
	}
	// assignment operator
	AppSettings &AppSettings::operator=(AppSettings &obj)
	{
		std::swap(host_params, obj.host_params);
		std::swap(molecule, obj.molecule);
		std::swap(molecule_device, obj.molecule_device);
		std::swap(molecule_loaded, obj.molecule_loaded);

		std::swap(colors_device, obj.colors_device);

		std::swap(voxel_data_device, obj.voxel_data_device);
		std::swap(voxel_count_device, obj.voxel_count_device);
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

	float getMaxSolventRadius()
	{
		return host_params.solvent_max;
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

	int *getDeviceVoxelData()
	{
		return voxel_data_device;
	}
	int *getDeviceVoxelCount()
	{
		return voxel_count_device;
	}

	bool getDebugMode()
	{
		return host_params.debug_mode;
	}

	SimulationParams getAllHostParams()
	{
		return host_params;
	}
	SimulationParams *getParamsPointer()
	{
		return &host_params;
	}

	float4 getCameraFocus()
	{
		return focusPoint;
	}

	int getDebugFrame()
	{
		return host_params.debug_frame;
	}
	std::string getRecentFilePath()
	{
		return recentFilePath;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setter

	void setWindowWidth(uint value)
	{
		host_params.window_width = value;
	}
	void setWindowHeight(uint value)
	{
		host_params.window_height = value;
	}
	void addToSolventRadius(float value)
	{
		host_params.solvent_radius += value;
		if (host_params.solvent_radius > host_params.solvent_max)
		{
			value = host_params.solvent_radius - host_params.solvent_max;
			host_params.solvent_radius = host_params.solvent_max;
		}
		if (host_params.solvent_radius < 0)
		{
			value = -host_params.solvent_radius;
			host_params.solvent_radius = host_params.solvent_max;
		}
		addToRadius(molecule, host_params.numAtoms, value);
		checkCudaErrors(cudaMemcpy(molecule_device, &molecule, host_params.numAtoms * sizeof(float4), cudaMemcpyHostToDevice));
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
	void setKnearest(unsigned int value)
	{
		host_params.k_nearest = value;
	}

	void setVoxelUsage(bool value)
	{
		host_params.use_voxel = value;
	}

	void setDebugFrame(int value)
	{
		host_params.debug_frame = value;
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

		// TODO: on color change, figure out how to reload molecule
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////
};
#endif