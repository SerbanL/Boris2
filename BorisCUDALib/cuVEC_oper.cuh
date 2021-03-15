#pragma once

#include "cuVEC.h"
#include "cuFuncs_Math.h"
#include "launchers.h"
#include "Reduction.cuh"

//------------------------------------------------------------------- SETBOX

template <typename VType>
__global__ void setbox_kernel(cuSZ3& n, cuBox box, VType value, VType*& quantity)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	if (idx < n.dim()) {

		if (box.Contains(ijk)) quantity[idx] = value;
	}
}

template void cuVEC<char>::setbox(cuBox box, char value);
template void cuVEC<int>::setbox(cuBox box, int value);
template void cuVEC<unsigned>::setbox(cuBox box, unsigned value);
template void cuVEC<long>::setbox(cuBox box, long value);
template void cuVEC<size_t>::setbox(cuBox box, size_t value);
template void cuVEC<float>::setbox(cuBox box, float value);
template void cuVEC<double>::setbox(cuBox box, double value);

template void cuVEC<cuINT3>::setbox(cuBox box, cuINT3 value);
template void cuVEC<cuSZ3>::setbox(cuBox box, cuSZ3 value);
template void cuVEC<cuFLT3>::setbox(cuBox box, cuFLT3 value);
template void cuVEC<cuDBL3>::setbox(cuBox box, cuDBL3 value);

template void cuVEC<cuINT4>::setbox(cuBox box, cuINT4 value);
template void cuVEC<cuSZ4>::setbox(cuBox box, cuSZ4 value);
template void cuVEC<cuFLT4>::setbox(cuBox box, cuFLT4 value);
template void cuVEC<cuDBL4>::setbox(cuBox box, cuDBL4 value);

template void cuVEC<cuINT33>::setbox(cuBox box, cuINT33 value);
template void cuVEC<cuFLT33>::setbox(cuBox box, cuFLT33 value);
template void cuVEC<cuDBL33>::setbox(cuBox box, cuDBL33 value);

template void cuVEC<cuReIm>::setbox(cuBox box, cuReIm value);
template void cuVEC<cuReIm3>::setbox(cuBox box, cuReIm3 value);

template <typename VType>
__host__ void cuVEC<VType>::setbox(cuBox box, VType value)
{
	setbox_kernel <<< (get_gpu_value(n).dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, box, value, quantity);
}

//------------------------------------------------------------------- SET

template <typename VType>
__global__ void set_kernel(size_t size, VType value, VType*& quantity)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		quantity[idx] = value;
	}
}

template void cuVEC<char>::set(size_t size, char value);
template void cuVEC<int>::set(size_t size, int value);
template void cuVEC<unsigned>::set(size_t size, unsigned value);
template void cuVEC<long>::set(size_t size, long value);
template void cuVEC<size_t>::set(size_t size, size_t value);
template void cuVEC<float>::set(size_t size, float value);
template void cuVEC<double>::set(size_t size, double value);

template void cuVEC<cuINT3>::set(size_t size, cuINT3 value);
template void cuVEC<cuSZ3>::set(size_t size, cuSZ3 value);
template void cuVEC<cuFLT3>::set(size_t size, cuFLT3 value);
template void cuVEC<cuDBL3>::set(size_t size, cuDBL3 value);

template void cuVEC<cuINT4>::set(size_t size, cuINT4 value);
template void cuVEC<cuSZ4>::set(size_t size, cuSZ4 value);
template void cuVEC<cuFLT4>::set(size_t size, cuFLT4 value);
template void cuVEC<cuDBL4>::set(size_t size, cuDBL4 value);

template void cuVEC<cuINT33>::set(size_t size, cuINT33 value);
template void cuVEC<cuFLT33>::set(size_t size, cuFLT33 value);
template void cuVEC<cuDBL33>::set(size_t size, cuDBL33 value);

template void cuVEC<cuReIm>::set(size_t size, cuReIm value);
template void cuVEC<cuReIm3>::set(size_t size, cuReIm3 value);

template <typename VType>
__host__ void cuVEC<VType>::set(size_t size, VType value)
{
	set_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (size, value, quantity);
}

//------------------------------------------------------------------- RENORMALIZE

template <typename VType, typename PType>
__global__ void renormalize_kernel(cuSZ3& n, VType*& quantity, PType new_norm)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < n.dim()) {

		PType curr_norm = cu_GetMagnitude(quantity[idx]);

		if (cuIsNZ(curr_norm)) quantity[idx] *= new_norm / curr_norm;
	}
}

template void cuVEC<float>::renormalize(size_t arr_size, float new_norm);
template void cuVEC<double>::renormalize(size_t arr_size, double new_norm);

template void cuVEC<cuFLT3>::renormalize(size_t arr_size, float new_norm);
template void cuVEC<cuDBL3>::renormalize(size_t arr_size, double new_norm);

template void cuVEC<cuFLT4>::renormalize(size_t arr_size, float new_norm);
template void cuVEC<cuDBL4>::renormalize(size_t arr_size, double new_norm);

template <typename VType>
template <typename PType>
__host__ void cuVEC<VType>::renormalize(size_t arr_size, PType new_norm)
{
	renormalize_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, quantity, new_norm);
}