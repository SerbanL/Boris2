#pragma once

#include "cuArray.h"


//------------------------------------------- SET VALUE : cuArray_aux.h

template <typename VType>
__host__ void cu_arr<VType>::set(VType value)
{
	gpu_set(cu_array, value, size_cpu);
}

//------------------------------------------- GET SIZE : cuArray_aux.h

template <typename VType>
__host__ size_t cu_arr<VType>::size(void)
{
	return size_cpu;
}