#pragma once

#include "cuVEC.h"
#include "launchers.h"

//------------------------------------------------------------------- ADDITION

template <typename VType>
__global__ void add_kernel(cuVEC<VType>& add_to_this, cuVEC<VType>& add_this)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < add_this.linear_size() && idx < add_to_this.linear_size()) {

		add_to_this[idx] = add_to_this[idx] + add_this[idx];
	}
}

template void cuVEC<float>::add_values(size_t size, cu_obj<cuVEC<float>>& add_this);
template void cuVEC<double>::add_values(size_t size, cu_obj<cuVEC<double>>& add_this);

template void cuVEC<cuFLT3>::add_values(size_t size, cu_obj<cuVEC<cuFLT3>>& add_this);
template void cuVEC<cuDBL3>::add_values(size_t size, cu_obj<cuVEC<cuDBL3>>& add_this);

template void cuVEC<cuFLT4>::add_values(size_t size, cu_obj<cuVEC<cuFLT4>>& add_this);
template void cuVEC<cuDBL4>::add_values(size_t size, cu_obj<cuVEC<cuDBL4>>& add_this);

template void cuVEC<cuReIm3>::add_values(size_t size, cu_obj<cuVEC<cuReIm3>>& add_this);

//add to this vec the values in add_this : must have same size : size
template <typename VType>
__host__ void cuVEC<VType>::add_values(size_t size, cu_obj<cuVEC<VType>>& add_this)
{
	add_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, (cuVEC<VType>&)add_this);
}

template <typename VType>
__host__ void cuVEC<VType>::add_values(size_t size, cuVEC<VType>& add_this)
{
	add_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, add_this);
}

//------------------------------------------------------------------- SUBTRACTION

template <typename VType>
__global__ void sub_kernel(cuVEC<VType>& sub_from_this, cuVEC<VType>& sub_this)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < sub_this.linear_size() && idx < sub_from_this.linear_size()) {

		sub_from_this[idx] = sub_from_this[idx] + sub_this[idx];
	}
}

template void cuVEC<float>::sub_values(size_t size, cu_obj<cuVEC<float>>& sub_this);
template void cuVEC<double>::sub_values(size_t size, cu_obj<cuVEC<double>>& sub_this);

template void cuVEC<cuFLT3>::sub_values(size_t size, cu_obj<cuVEC<cuFLT3>>& sub_this);
template void cuVEC<cuDBL3>::sub_values(size_t size, cu_obj<cuVEC<cuDBL3>>& sub_this);

template void cuVEC<cuFLT4>::sub_values(size_t size, cu_obj<cuVEC<cuFLT4>>& sub_this);
template void cuVEC<cuDBL4>::sub_values(size_t size, cu_obj<cuVEC<cuDBL4>>& sub_this);

template void cuVEC<cuReIm3>::sub_values(size_t size, cu_obj<cuVEC<cuReIm3>>& sub_this);

//subtract from this vec the values in sub_this : must have same size : size
template <typename VType>
__host__ void cuVEC<VType>::sub_values(size_t size, cu_obj<cuVEC<VType>>& sub_this)
{
	sub_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (*this, (cuVEC<VType>&)sub_this);
}

template <typename VType>
__host__ void cuVEC<VType>::sub_values(size_t size, cuVEC<VType>& sub_this)
{
	sub_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, sub_this);
}