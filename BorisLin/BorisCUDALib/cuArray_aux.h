#pragma once

#include "cuArray.h"


//------------------------------------------- SET VALUE : cuArray_aux.h

template <typename VType>
__host__ void cu_arr<VType>::set(VType value)
{
	gpu_set(cu_array, value, arr_size);
}

//set single value at given index
template <typename VType>
__host__ void cu_arr<VType>::setvalue(int index, VType value)
{
	cpu_to_gpu(cu_array + index, &value);
}

//------------------------------------------- GET SIZE : cuArray_aux.h

template <typename VType>
__host__ size_t cu_arr<VType>::size(void)
{
	return arr_size;
}