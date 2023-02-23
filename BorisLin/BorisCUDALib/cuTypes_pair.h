#pragma once

#include <cuda_runtime.h>

#include "alloc_cpy.h"

#include <utility>

////////////////////////////////////////////////////////////////////////////////////////////////// cuPair, similar to a std::pair but for cuda use
//
//

template <typename FType, typename SType> 
struct cuPair {

	//----------------------------- DATA

	FType first;
	SType second;

	//----------------------------- cu_obj MANAGED CONSTRUCTORS / DESTRUCTOR

	__host__ void construct_cu_obj(void)
	{
		set_gpu_value(first, FType());
		set_gpu_value(second, FType());
	}

	__host__ void construct_cu_obj(FType first_, SType second_)
	{
		set_gpu_value(first, first_);
		set_gpu_value(second, second_);
	}

	__host__ void construct_cu_obj(const cuPair& copyThis)
	{
		assign_cu_obj(copyThis);
	}

	__host__ void assign_cu_obj(const cuPair& copyThis)
	{
		gpu_to_gpu(first, copyThis.first);
		gpu_to_gpu(second, copyThis.second);
	}

	__host__ void destruct_cu_obj(void)
	{
	}

	//----------------------------- CONVERSION TO/FROM NON-CUDA VERSION

	template <typename FType_, typename SType_>
	__host__ operator std::pair<FType_, SType_>() const
	{
		return std::pair<FType_, SType_>((FType_)first, (SType_)second);
	}

	template <typename FType_, typename SType_>
	__host__ cuPair<FType, SType>& operator=(const std::pair<FType_, SType_> &rhs)
	{
		first = (FType)rhs.first; 
		second = (SType)rhs.second;
		return *this;
	}

	template <typename FType_, typename SType_>
	__host__ cuPair(const std::pair<FType_, SType_> &rhs)
	{
		first = (FType)rhs.first;
		second = (SType)rhs.second;
	}
};