#pragma once

#include "cuVEC.h"

template <typename VType>
class cuTransfer {

	//array with transfer information from mesh_in : full cell indexes for transfer : i - mesh index, j - mesh cell index, k - supermesh (this mesh) cell index
	cuPair<cuINT3, cuBReal>* transfer_in_info;
	size_t transfer_in_info_size;

	//array with transfer information to mesh_out : full cell indexes for transfer : i - mesh index, j - mesh cell index, k - supermesh (this mesh) cell index
	cuPair<cuINT3, cuBReal>* transfer_out_info;
	size_t transfer_out_info_size;

	//array of meshes for transfer in and out from / to
	//mesh_in and mesh_out VECs have exactly the same rectangle and cellsize for each index, but may differ in value stored (e.g. magnetization and effective field) - they could also be exactly the same VEC
	//mesh_in2 can be used if we require multiple inputs, e.g. averaging inputs or multiplying inputs
	//mesh_out2 can be used if require duplicating outputs
	//For both mesh_in2 and mesh_out2, the input averaging and output duplicating is done if the respective VECs are not empty
	//Thus when using these modes, the secondary VECs should either be empty or have exactly same size as the primary VECs.
	//In any case, if using these modes the vectors below have to have exactly the same dimensions
	cuVEC<VType>* mesh_in; 
	cuVEC<VType>* mesh_in2;
	cuVEC<cuBReal>* mesh_in2_real;
	cuVEC<VType>* mesh_out;
	cuVEC<VType>* mesh_out2;

private:
	
	//auxiliary : does the actual transfer info copy after in and out cuVEC have been set
	template <typename cpuTransfer>
	__host__ bool copy_transfer_info(cpuTransfer& vec_transfer);

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
		
		nullgpuptr(mesh_in2);
		nullgpuptr(mesh_out2);

		nullgpuptr(mesh_in2_real);

		set_transfer_in_info_size(0);
		set_transfer_out_info_size(0);
	}

	__host__ void destruct_cu_obj(void)
	{
		gpu_free_managed(transfer_in_info);
		gpu_free_managed(transfer_out_info);
		
		gpu_free_managed(mesh_in);
		gpu_free_managed(mesh_out);
		
		gpu_free_managed(mesh_in2);
		gpu_free_managed(mesh_out2);

		gpu_free_managed(mesh_in2_real);
	}

	//reset to zero all transfer info data
	__host__ void clear_transfer_data(void);

	//--------------------------------------------MESH TRANSFER COPY

	//SINGLE INPUT, SINGLE OUTPUT

	template <typename cpuTransfer>
	__host__ bool copy_transfer_info(cu_arr<cuVEC<VType>>& mesh_in_arr, cu_arr<cuVEC<VType>>& mesh_out_arr, cpuTransfer& vec_transfer);

	//MULTIPLE INPUTS, SINGLE OUTPUT

	template <typename cpuTransfer>
	__host__ bool copy_transfer_info_averagedinputs(
		cu_arr<cuVEC<VType>>& mesh_in_arr1, cu_arr<cuVEC<VType>>& mesh_in_arr2, 
		cu_arr<cuVEC<VType>>& mesh_out_arr, 
		cpuTransfer& vec_transfer);

	template <typename cpuTransfer>
	__host__ bool copy_transfer_info_multipliedinputs(
		cu_arr<cuVEC<VType>>& mesh_in_arr1, cu_arr<cuVEC<cuBReal>>& mesh_in_arr2_real,
		cu_arr<cuVEC<VType>>& mesh_out_arr,
		cpuTransfer& vec_transfer);

	//MULTIPLE INPUT, MULTIPLE OUTPUT

	template <typename cpuTransfer>
	__host__ bool copy_transfer_info_averagedinputs_duplicatedoutputs(
		cu_arr<cuVEC<VType>>& mesh_in_arr1, cu_arr<cuVEC<VType>>& mesh_in_arr2, 
		cu_arr<cuVEC<VType>>& mesh_out_arr1, cu_arr<cuVEC<VType>>& mesh_out_arr2, 
		cpuTransfer& vec_transfer);

	//--------------------------------------------MESH TRANSFER

	//SINGLE INPUT, SINGLE OUTPUT

	//do the actual transfer of values to and from this mesh using these - pass in the size of quantity (= get_gpu_value(n).dim()), and transfer_info size to speed up call
	void transfer_in(size_t size_transfer, VType*& sMesh_quantity);

	//transfer to output meshes. setOutput = true means set mesh values, not add. Pass in size_transfer (transfer_info_size) and number of output meshes to speed up call
	void transfer_out(size_t size_transfer, VType*& sMesh_quantity, int mesh_out_num);

	//AVERAGED INPUTS

	//do the actual transfer of values to and from this mesh using these - pass in the size of quantity (= get_gpu_value(n).dim()), and transfer_info size to speed up call
	void transfer_in_averaged(size_t size_transfer, VType*& sMesh_quantity);

	//MULTIPLIED INPUTS

	//do the actual transfer of values to and from this mesh using these - pass in the size of quantity (= get_gpu_value(n).dim()), and transfer_info size to speed up call
	void transfer_in_multiplied(size_t size_transfer, VType*& sMesh_quantity);

	//DUPLICATED OUTPUT

	//transfer to output meshes. setOutput = true means set mesh values, not add. Pass in size_transfer (transfer_info_size) and number of output meshes to speed up call
	void transfer_out_duplicated(size_t size_transfer, VType*& sMesh_quantity, int mesh_out_num);
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
	
	gpu_free_managed(mesh_in2);
	gpu_free_managed(mesh_out2);

	gpu_free_managed(mesh_in2_real);
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

//auxiliary : does the actual transfer info copy after in and out cuVEC have been set
template <typename VType>
template <typename cpuTransfer>
__host__ bool cuTransfer<VType>::copy_transfer_info(cpuTransfer& vec_transfer)
{
	//copy transfer_in_info_cpu to transfer_in_info : convert it first
	std::vector<std::pair<cuINT3, cuBReal>> transfer_in_info_cpu;

	//set size
	if (!malloc_vector(transfer_in_info_cpu, vec_transfer.size_transfer_in())) return false;

	std::vector<std::pair<INT3, double>> cpuVEC_transfer_in_info = vec_transfer.get_flattened_transfer_in_info();
	if (cpuVEC_transfer_in_info.size() != vec_transfer.size_transfer_in()) return false;

	//now copy and convert
#pragma omp parallel for
	for (int idx = 0; idx < vec_transfer.size_transfer_in(); idx++) {

		transfer_in_info_cpu[idx].first = cpuVEC_transfer_in_info[idx].first;
		transfer_in_info_cpu[idx].second = cpuVEC_transfer_in_info[idx].second;
	}

	//allocate gpu memory for transfer_out_info
	if (!set_transfer_in_info_size(vec_transfer.size_transfer_in())) return false;

	//copy to transfer_in_info
	cpu_to_gpu_managed(transfer_in_info, transfer_in_info_cpu.data(), transfer_in_info_cpu.size());

	//------

	//copy transfer_out_info_cpu to transfer_out_info : convert it first
	std::vector<std::pair<cuINT3, cuBReal>> transfer_out_info_cpu;

	//set size
	if (!malloc_vector(transfer_out_info_cpu, vec_transfer.size_transfer_out())) return false;

	std::vector<std::pair<INT3, double>> cpuVEC_transfer_out_info = vec_transfer.get_flattened_transfer_out_info();
	if (cpuVEC_transfer_out_info.size() != vec_transfer.size_transfer_out()) return false;

	//now copy and convert
#pragma omp parallel for
	for (int idx = 0; idx < vec_transfer.size_transfer_out(); idx++) {

		transfer_out_info_cpu[idx].first = cpuVEC_transfer_out_info[idx].first;
		transfer_out_info_cpu[idx].second = cpuVEC_transfer_out_info[idx].second;
	}

	//allocate gpu memory for transfer_out_info
	if (!set_transfer_out_info_size(vec_transfer.size_transfer_out())) return false;

	//copy to transfer_in_info
	cpu_to_gpu_managed(transfer_out_info, transfer_out_info_cpu.data(), transfer_out_info_cpu.size());

	return true;
}


//SINGLE INPUT, SINGLE OUTPUT

//copy precalculated transfer info from cpu memory
template <typename VType>
template <typename cpuTransfer>
__host__ bool cuTransfer<VType>::copy_transfer_info(cu_arr<cuVEC<VType>>& mesh_in_arr, cu_arr<cuVEC<VType>>& mesh_out_arr, cpuTransfer& vec_transfer)
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

	return copy_transfer_info(vec_transfer);
}

//MULTIPLE INPUTS, SINGLE OUTPUT

//copy precalculated transfer info from cpu memory
template <typename VType>
template <typename cpuTransfer>
__host__ bool cuTransfer<VType>::copy_transfer_info_averagedinputs(
	cu_arr<cuVEC<VType>>& mesh_in_arr1, cu_arr<cuVEC<VType>>& mesh_in_arr2, 
	cu_arr<cuVEC<VType>>& mesh_out_arr,
	cpuTransfer& vec_transfer)
{
	clear_transfer_data();

	//------

	//copy cuVEC pointers from mesh_in_arr to the managed mesh_in array for later usage
	size_t size_in = mesh_in_arr1.size();
	size_t size_in2 = mesh_in_arr2.size();

	if (size_in != size_in2) return false;

	if (size_in > 0) {

		gpu_alloc_managed(mesh_in, size_in);
		gpu_to_gpu_managed1st(mesh_in, (cuVEC<VType>*)mesh_in_arr1, size_in);

		gpu_alloc_managed(mesh_in2, size_in);
		gpu_to_gpu_managed1st(mesh_in2, (cuVEC<VType>*)mesh_in_arr2, size_in);
	}

	//copy cuVEC pointers from mesh_out_arr to the managed mesh_out array for later usage
	size_t size_out = mesh_out_arr.size();

	if (size_out > 0) {

		gpu_alloc_managed(mesh_out, size_out);
		gpu_to_gpu_managed1st(mesh_out, (cuVEC<VType>*)mesh_out_arr, size_out);
	}

	//------

	return copy_transfer_info(vec_transfer);
}

//copy precalculated transfer info from cpu memory
template <typename VType>
template <typename cpuTransfer>
__host__ bool cuTransfer<VType>::copy_transfer_info_multipliedinputs(
	cu_arr<cuVEC<VType>>& mesh_in_arr1, cu_arr<cuVEC<cuBReal>>& mesh_in_arr2_real,
	cu_arr<cuVEC<VType>>& mesh_out_arr,
	cpuTransfer& vec_transfer)
{
	clear_transfer_data();

	//------

	//copy cuVEC pointers from mesh_in_arr to the managed mesh_in array for later usage
	size_t size_in = mesh_in_arr1.size();
	size_t size_in2 = mesh_in_arr2_real.size();

	if (size_in != size_in2) return false;

	if (size_in > 0) {

		gpu_alloc_managed(mesh_in, size_in);
		gpu_to_gpu_managed1st(mesh_in, (cuVEC<VType>*)mesh_in_arr1, size_in);

		gpu_alloc_managed(mesh_in2_real, size_in);
		gpu_to_gpu_managed1st(mesh_in2_real, (cuVEC<cuBReal>*)mesh_in_arr2_real, size_in);
	}

	//copy cuVEC pointers from mesh_out_arr to the managed mesh_out array for later usage
	size_t size_out = mesh_out_arr.size();

	if (size_out > 0) {

		gpu_alloc_managed(mesh_out, size_out);
		gpu_to_gpu_managed1st(mesh_out, (cuVEC<VType>*)mesh_out_arr, size_out);
	}

	//------

	return copy_transfer_info(vec_transfer);
}

//MULTIPLE INPUTS, MULTIPLE OUTPUTS

//copy precalculated transfer info from cpu memory
template <typename VType>
template <typename cpuTransfer>
__host__ bool cuTransfer<VType>::copy_transfer_info_averagedinputs_duplicatedoutputs(
	cu_arr<cuVEC<VType>>& mesh_in_arr1, cu_arr<cuVEC<VType>>& mesh_in_arr2, 
	cu_arr<cuVEC<VType>>& mesh_out_arr1, cu_arr<cuVEC<VType>>& mesh_out_arr2, 
	cpuTransfer& vec_transfer)
{
	clear_transfer_data();

	//------

	//copy cuVEC pointers from mesh_in_arr to the managed mesh_in array for later usage
	size_t size_in = mesh_in_arr1.size();
	size_t size_in2 = mesh_in_arr2.size();

	size_t size_out = mesh_out_arr1.size();
	size_t size_out2 = mesh_out_arr2.size();

	if (size_in != size_in2 || size_out != size_out2) return false;

	if (size_in > 0) {

		gpu_alloc_managed(mesh_in, size_in);
		gpu_to_gpu_managed1st(mesh_in, (cuVEC<VType>*)mesh_in_arr1, size_in);

		gpu_alloc_managed(mesh_in2, size_in);
		gpu_to_gpu_managed1st(mesh_in2, (cuVEC<VType>*)mesh_in_arr2, size_in);
	}

	//copy cuVEC pointers from mesh_out_arr to the managed mesh_out array for later usage

	if (size_out > 0) {

		gpu_alloc_managed(mesh_out, size_out);
		gpu_to_gpu_managed1st(mesh_out, (cuVEC<VType>*)mesh_out_arr1, size_out);

		gpu_alloc_managed(mesh_out2, size_out);
		gpu_to_gpu_managed1st(mesh_out2, (cuVEC<VType>*)mesh_out_arr2, size_out);
	}

	//------

	return copy_transfer_info(vec_transfer);
}