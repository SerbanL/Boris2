#pragma once

#include "cuArray.h"

//------------------------------------------- GPU <-> CPU TRANSFER : cuArray_transfer.h

//copy values from a std::vector into gpu memory. Doesn't set size, but copies up to currently allocated size.
template <typename VType>
template <typename Type>
__host__ void cu_arr<VType>::copy_from_vector(std::vector<Type>& cpuvec)
{
	size_t size_transfer = (cpuvec.size() < arr_size ? cpuvec.size() : arr_size);

	if (size_transfer) cpu_to_gpu(cu_array, cpuvec.data(), (cpuvec.size() < arr_size ? cpuvec.size() : arr_size));
}

//copy values to a std::vector into cpu memory. Doesn't set size, but copies up to currently allocated size, starting at given offset in cpuvec
template <typename VType>
template <typename Type>
__host__ void cu_arr<VType>::copy_to_vector(std::vector<Type>& cpuvec, size_t offset)
{
	size_t size_transfer = (cpuvec.size() - offset < arr_size ? cpuvec.size() - offset : arr_size);

	if (size_transfer) gpu_to_cpu(cpuvec.data() + offset, cu_array, size_transfer);
}