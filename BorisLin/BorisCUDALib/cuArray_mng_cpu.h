#pragma once

#include "cuArray.h"

//------------------------------------------- CONSTRUCTORS

//void constructor
template <typename VType>
__host__ cu_arr<VType>::cu_arr(void)
{
	cu_array = nullptr;
	pcu_array = nullptr;

	arr_size = 0;
}

//size constructor
template <typename VType>
__host__ cu_arr<VType>::cu_arr(size_t size_)
{
	cu_array = nullptr;
	pcu_array = nullptr;

	arr_size = 0;

	resize(size_);
}

//------------------------------------------- DESTRUCTOR

//destructor
template <typename VType>
__host__ cu_arr<VType>::~cu_arr()
{
	gpu_free(cu_array);
	gpu_free(pcu_array);
}

	