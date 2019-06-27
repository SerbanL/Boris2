#pragma once

#include "cuVEC.h"

template <typename VType>
class cuTransfer {

	//array with transfer information from mesh_in : full cell indexes for transfer : i - mesh index, j - mesh cell index, k - supermesh (this mesh) cell index
	cuPair<cuINT3, cuReal>* transfer_in_info;
	size_t transfer_in_info_size;

	//array with transfer information to mesh_out : full cell indexes for transfer : i - mesh index, j - mesh cell index, k - supermesh (this mesh) cell index
	cuPair<cuINT3, cuReal>* transfer_out_info;
	size_t transfer_out_info_size;

	//array of meshes for transfer in and out from / to
	cuVEC<VType>* mesh_in;
	cuVEC<VType>* mesh_out;

private:
	
	//set and allocate memory for transfer_in_info array
	__host__ bool set_transfer_in_info_size(size_t size);

	//set and allocate memory for transfer_out_info array
	__host__ bool set_transfer_out_info_size(size_t size);
	
public:

	__host__ void construct_cu_obj(void)
	{
		nullgpuptr(transfer_in_info);
		nullgpuptr(transfer_out_info);
		nullgpuptr(mesh_in);
		nullgpuptr(mesh_out);

		set_transfer_in_info_size(0);
		set_transfer_out_info_size(0);
	}

	__host__ void destruct_cu_obj(void)
	{
		gpu_free_managed(transfer_in_info);
		gpu_free_managed(transfer_out_info);
		gpu_free_managed(mesh_in);
		gpu_free_managed(mesh_out);
	}

	//reset to zero all transfer info data
	__host__ void clear_transfer_data(void);

	//--------------------------------------------MESH TRANSFER COPY

	template <typename cpuVEC>
	__host__ bool copy_transfer_info(cu_arr<cuVEC<VType>>& mesh_in_arr, cu_arr<cuVEC<VType>>& mesh_out_arr, cpuVEC& cpuVEC);

	//--------------------------------------------MESH TRANSFER

	//do the actual transfer of values to and from this mesh using these - pass in the size of quantity (= get_gpu_value(n).dim()), and transfer_info size to speed up call
	void transfer_in(size_t size_transfer, VType*& sMesh_quantity);

	//transfer to output meshes. setOutput = true means set mesh values, not add. Pass in size_transfer (transfer_info_size) and number of output meshes to speed up call
	void transfer_out(size_t size_transfer, VType*& sMesh_quantity, int mesh_out_num);
};

//--------------------------------------------HELPERS

//reset to zero all transfer info data
template <typename VType>
__host__ void cuTransfer<VType>::clear_transfer_data(void)
{
	set_transfer_in_info_size(0);
	set_transfer_out_info_size(0);

	gpu_free_managed(transfer_in_info);
	gpu_free_managed(transfer_out_info);
	gpu_free_managed(mesh_in);
	gpu_free_managed(mesh_out);
}

//set and allocate memory for transfer_in_info array
template <typename VType>
__host__ bool cuTransfer<VType>::set_transfer_in_info_size(size_t size)
{
	//size_t local_value = size;
	set_gpu_value(transfer_in_info_size, size);

	if (size == 0) {

		gpu_free_managed(transfer_in_info);
		return true;
	}

	//new size value set : must also adjust memory allocation
	cudaError_t error = gpu_alloc_managed(transfer_in_info, size);
	if (error != cudaSuccess) {

		gpu_free_managed(transfer_in_info);
		set_gpu_value(transfer_in_info_size, (size_t)0);

		return false;
	}

	return true;
}

//set and allocate memory for transfer_out_info array
template <typename VType>
__host__ bool cuTransfer<VType>::set_transfer_out_info_size(size_t size)
{
	//size_t local_value = size;
	set_gpu_value(transfer_out_info_size, size);

	if (size == 0) {

		gpu_free_managed(transfer_out_info);
		return true;
	}

	//new size value set : must also adjust memory allocation
	cudaError_t error = gpu_alloc_managed(transfer_out_info, size);
	if (error != cudaSuccess) {

		gpu_free_managed(transfer_out_info);
		set_gpu_value(transfer_out_info_size, (size_t)0);

		return false;
	}

	return true;
}

//--------------------------------------------MESH TRANSFER COPY

//copy precalculated transfer info from cpu memory
template <typename VType>
template <typename cpuVEC>
__host__ bool cuTransfer<VType>::copy_transfer_info(cu_arr<cuVEC<VType>>& mesh_in_arr, cu_arr<cuVEC<VType>>& mesh_out_arr, cpuVEC& cpuVEC)
{
	clear_transfer_data();

	//------

	//copy cuVEC pointers from mesh_in_arr to the managed mesh_in array for later usage
	size_t size_in = mesh_in_arr.size();

	if (size_in > 0) {

		gpu_alloc_managed(mesh_in, size_in);
		gpu_to_gpu_managed1st(mesh_in, (cuVEC<VType>*)mesh_in_arr, size_in);
	}

	//copy cuVEC pointers from mesh_out_arr to the managed mesh_out array for later usage
	size_t size_out = mesh_out_arr.size();

	if (size_out > 0) {

		gpu_alloc_managed(mesh_out, size_out);
		gpu_to_gpu_managed1st(mesh_out, (cuVEC<VType>*)mesh_out_arr, size_out);
	}

	//------

	//copy transfer_in_info_cpu to transfer_in_info : convert it first
	std::vector<std::pair<cuINT3, cuReal>> transfer_in_info_cpu;

	//set size
	if (!malloc_vector(transfer_in_info_cpu, cpuVEC.size_transfer_in())) return false;

	std::vector<std::pair<INT3, double>> cpuVEC_transfer_in_info = cpuVEC.get_flattened_transfer_in_info();
	if (cpuVEC_transfer_in_info.size() != cpuVEC.size_transfer_in()) return false;

	//now copy and convert
#pragma omp parallel for
	for (int idx = 0; idx < cpuVEC.size_transfer_in(); idx++) {

		transfer_in_info_cpu[idx].first = cpuVEC_transfer_in_info[idx].first;
		transfer_in_info_cpu[idx].second = cpuVEC_transfer_in_info[idx].second;
	}

	//allocate gpu memory for transfer_out_info
	if (!set_transfer_in_info_size(cpuVEC.size_transfer_in())) return false;

	//copy to transfer_in_info
	cpu_to_gpu_managed(transfer_in_info, transfer_in_info_cpu.data(), transfer_in_info_cpu.size());

	//------
	
	//copy transfer_out_info_cpu to transfer_out_info : convert it first
	std::vector<std::pair<cuINT3, cuReal>> transfer_out_info_cpu;

	//set size
	if (!malloc_vector(transfer_out_info_cpu, cpuVEC.size_transfer_out())) return false;

	std::vector<std::pair<INT3, double>> cpuVEC_transfer_out_info = cpuVEC.get_flattened_transfer_out_info();
	if (cpuVEC_transfer_out_info.size() != cpuVEC.size_transfer_out()) return false;

	//now copy and convert
#pragma omp parallel for
	for (int idx = 0; idx < cpuVEC.size_transfer_out(); idx++) {

		transfer_out_info_cpu[idx].first = cpuVEC_transfer_out_info[idx].first;
		transfer_out_info_cpu[idx].second = cpuVEC_transfer_out_info[idx].second;
	}

	//allocate gpu memory for transfer_out_info
	if (!set_transfer_out_info_size(cpuVEC.size_transfer_out())) return false;

	//copy to transfer_in_info
	cpu_to_gpu_managed(transfer_out_info, transfer_out_info_cpu.data(), transfer_out_info_cpu.size());

	return true;
}