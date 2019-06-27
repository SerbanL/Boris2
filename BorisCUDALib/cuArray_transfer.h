#pragma once

#include "cuArray.h"

//------------------------------------------- GPU <-> CPU TRANSFER : cuArray_transfer.h

//copy values from a std::vector into gpu memory. Doesn't set size, but copies up to currently allocated size.
template <typename VType>
template <typename Type>
__host__ void cu_arr<VType>::copy_from_cpuvector(std::vector<Type>& cpuvec)
{
	cpu_to_gpu(cu_array, cpuvec.data(), (cpuvec.size() < size_cpu ? cpuvec.size() : size_cpu));
}

//copy values to a std::vector into cpu memory. Doesn't set size, but copies up to currently allocated size.
template <typename VType>
template <typename Type>
__host__ void cu_arr<VType>::copy_to_cpuvector(std::vector<Type>& cpuvec)
{
	gpu_to_cpu(cpuvec.data(), cu_array, (cpuvec.size() < size_cpu ? cpuvec.size() : size_cpu));
}