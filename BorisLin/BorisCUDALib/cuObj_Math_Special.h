#pragma once

#include "alloc_cpy.h"
#include "cuFuncs_Math.h"
#include "cuObject.h"

#include <vector>

//This mirrors the CPU version in Obj_Math_Special.h
//Just a wrapper to store analogous data on the GPU so we can use it in TEquation_CUDA

class ManagedFuncs_Special_CUDA {

	//smallest value to calculate in values : store it at values[0]
	cuBReal start;

	//number of points per unit to store in values
	int resolution;

	cuBReal* values;
	size_t size;

public:

	///////////////////////////////////////////////////////////////
	
	__host__ void construct_cu_obj(void)
	{
		set_gpu_value(start, (cuBReal)0.0);
		set_gpu_value(resolution, (int)1);
		
		set_gpu_value(size, (size_t)1);
		gpu_alloc_managed(values, 1);
		gpu_set_managed(values, (cuBReal)0.0, (size_t)1);
	}
	
	__host__ void destruct_cu_obj(void)
	{
		gpu_free_managed(values);
	}

	///////////////////////////////////////////////////////////////

	//set gpu data from cpu data
	__host__ bool set_data(std::vector<double>& data_cpu, double start_cpu, int resolution_cpu)
	{
		if (data_cpu.size() != get_gpu_value(size)) {

			cudaError_t error = gpu_alloc_managed(values, data_cpu.size());
			if (error != cudaSuccess) return false;
		}

		cudaError_t error = cpu_to_gpu_managed(values, data_cpu.data(), data_cpu.size());
		if (error != cudaSuccess) return false;

		set_gpu_value(start, (cuBReal)start_cpu);
		set_gpu_value(resolution, (int)resolution_cpu); 
		set_gpu_value(size, data_cpu.size());

		return true;
	}

	///////////////////////////////////////////////////////////////

	__device__ cuBReal evaluate(cuBReal x)
	{
		cuBReal findex = (x - start) * resolution;
		int index = (int)cu_floor_epsilon(findex);

		if (index + 1 < (int)size && index >= 0) {

			return values[index] * (cuBReal(index + 1) - findex) + values[index + 1] * (findex - cuBReal(index));
		}
		else {

			if (index >= 0) return values[size - 1];
			else return values[0];
		}
	}
};
