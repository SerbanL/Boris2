#pragma once

#include "cuVEC_VC.h"
#include "cuFuncs_Math.h"
#include "launchers.h"
#include "Reduction.cuh"

//------------------------------------------------------------------- SCALE VALUES

template <typename VType>
__global__ void scale_values_kernel(cuSZ3& n, int*& ngbrFlags, VType*& quantity, cuBReal constant)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < n.dim()) {

		if (ngbrFlags[idx] & NF_NOTEMPTY) {

			quantity[idx] *= constant;
		}
	}
}

template void cuVEC_VC<float>::scale_values(size_t size, cuBReal constant);
template void cuVEC_VC<double>::scale_values(size_t size, cuBReal constant);

template void cuVEC_VC<cuFLT3>::scale_values(size_t size, cuBReal constant);
template void cuVEC_VC<cuDBL3>::scale_values(size_t size, cuBReal constant);

//scale all stored values by the given constant
template <typename VType>
void cuVEC_VC<VType>::scale_values(size_t size, cuBReal constant)
{
	scale_values_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, ngbrFlags, cuVEC<VType>::quantity, constant);
}