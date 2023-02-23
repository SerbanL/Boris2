#pragma once

#include "cuArray.h"

//------------------------------------------- RESIZING cu_obj managed

//------------------------------------------- RESIZING

template <typename VType>
__host__ bool cu_arr<VType>::resize(size_t size_)
{
	if (size_ == size()) return true;

	if (!size_) clear();
	else {

		cudaError_t error = gpu_alloc(cu_array, size_);
		if (error == cudaSuccess) {

			arr_size = size_;

			if (!pcu_array) gpu_alloc(pcu_array);
			cpu_to_gpu(pcu_array, &cu_array);

			return true;
		}
		else {

			clear();
			return false;
		}
	}

	return true;
}

template <typename VType>
__host__ void cu_arr<VType>::clear(void)
{
	gpu_free(cu_array);
	cu_array = nullptr;

	gpu_free(pcu_array);
	pcu_array = nullptr;

	arr_size = 0;
}

//------------------------------------------- STORE NEW ENTRIES : cuArray_sizing.h

//new_entry is a pointer in cpu memory to an object in gpu memory
template <typename VType>
__host__ void cu_arr<VType>::push_back(VType*& new_entry)
{
	//allocate new memory size in a temporary array
	VType* new_array = nullptr;
	cudaError_t error = gpu_alloc(new_array, arr_size + 1);

	if (error != cudaSuccess) {

		gpu_free(new_array);
		return;
	}

	//copy data currently in array to temporary array (if any)
	if (arr_size > 0) {

		gpu_to_gpu(new_array, cu_array, arr_size);
	}

	//add new entry to end of temporary array
	gpu_to_gpu(new_array + arr_size, new_entry);

	//swap pointers so array now points to newly constructed memory
	gpu_swap(cu_array, new_array);

	//free old memory
	gpu_free(new_array);

	//set new size
	arr_size++;

	if (!pcu_array) gpu_alloc(pcu_array);
	cpu_to_gpu(pcu_array, &cu_array);
}