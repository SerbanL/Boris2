#pragma once

#include <cuda_runtime.h>

#include "launchers.h"
#include "cuTypes.h"

//------------------------------------------------------------------- SET

template <typename DType>
__global__ void set_cuda_array(size_t size, DType* cu_dest_ptr, DType value)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cu_dest_ptr[idx] = value;
	}
}

template void set_cuda_array_launcher<bool>(size_t size, bool* cu_dest_ptr, bool value);
template void set_cuda_array_launcher<char>(size_t size, char* cu_dest_ptr, char value);
template void set_cuda_array_launcher<int>(size_t size, int* cu_dest_ptr, int value);
template void set_cuda_array_launcher<unsigned>(size_t size, unsigned* cu_dest_ptr, unsigned value);
template void set_cuda_array_launcher<long>(size_t size, long* cu_dest_ptr, long value);
template void set_cuda_array_launcher<size_t>(size_t size, size_t* cu_dest_ptr, size_t value);
template void set_cuda_array_launcher<float>(size_t size, float* cu_dest_ptr, float value);
template void set_cuda_array_launcher<double>(size_t size, double* cu_dest_ptr, double value);

template void set_cuda_array_launcher<cuBComplex>(size_t size, cuBComplex* cu_dest_ptr, cuBComplex value);

template void set_cuda_array_launcher<cuINT2>(size_t size, cuINT2* cu_dest_ptr, cuINT2 value);
template void set_cuda_array_launcher<cuFLT2>(size_t size, cuFLT2* cu_dest_ptr, cuFLT2 value);
template void set_cuda_array_launcher<cuDBL2>(size_t size, cuDBL2* cu_dest_ptr, cuDBL2 value);

template void set_cuda_array_launcher<cuINT3>(size_t size, cuINT3* cu_dest_ptr, cuINT3 value);
template void set_cuda_array_launcher<cuSZ3>(size_t size, cuSZ3* cu_dest_ptr, cuSZ3 value);
template void set_cuda_array_launcher<cuFLT3>(size_t size, cuFLT3* cu_dest_ptr, cuFLT3 value);
template void set_cuda_array_launcher<cuDBL3>(size_t size, cuDBL3* cu_dest_ptr, cuDBL3 value);

template void set_cuda_array_launcher<cuINT4>(size_t size, cuINT4* cu_dest_ptr, cuINT4 value);
template void set_cuda_array_launcher<cuSZ4>(size_t size, cuSZ4* cu_dest_ptr, cuSZ4 value);
template void set_cuda_array_launcher<cuFLT4>(size_t size, cuFLT4* cu_dest_ptr, cuFLT4 value);
template void set_cuda_array_launcher<cuDBL4>(size_t size, cuDBL4* cu_dest_ptr, cuDBL4 value);

template void set_cuda_array_launcher<cuINT33>(size_t size, cuINT33* cu_dest_ptr, cuINT33 value);
template void set_cuda_array_launcher<cuFLT33>(size_t size, cuFLT33* cu_dest_ptr, cuFLT33 value);
template void set_cuda_array_launcher<cuDBL33>(size_t size, cuDBL33* cu_dest_ptr, cuDBL33 value);

template void set_cuda_array_launcher<cuReIm>(size_t size, cuReIm* cu_dest_ptr, cuReIm value);
template void set_cuda_array_launcher<cuReIm3>(size_t size, cuReIm3* cu_dest_ptr, cuReIm3 value);

template <typename DType>
void set_cuda_array_launcher(size_t size, DType* cu_dest_ptr, DType value)
{
	set_cuda_array <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (size, cu_dest_ptr, value);
}

//------------------------------------------------------------------- MULTIPLY

template <typename DType>
__global__ void mul_cuda_array(size_t size, DType* cu_dest_ptr, DType value)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cu_dest_ptr[idx] = cu_dest_ptr[idx] * value;
	}
}

template void mul_cuda_array_launcher<bool>(size_t size, bool* cu_dest_ptr, bool value);
template void mul_cuda_array_launcher<char>(size_t size, char* cu_dest_ptr, char value);
template void mul_cuda_array_launcher<int>(size_t size, int* cu_dest_ptr, int value);
template void mul_cuda_array_launcher<unsigned>(size_t size, unsigned* cu_dest_ptr, unsigned value);
template void mul_cuda_array_launcher<long>(size_t size, long* cu_dest_ptr, long value);
template void mul_cuda_array_launcher<size_t>(size_t size, size_t* cu_dest_ptr, size_t value);
template void mul_cuda_array_launcher<float>(size_t size, float* cu_dest_ptr, float value);
template void mul_cuda_array_launcher<double>(size_t size, double* cu_dest_ptr, double value);

template <typename DType>
void mul_cuda_array_launcher(size_t size, DType* cu_dest_ptr, DType value)
{
	mul_cuda_array << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (size, cu_dest_ptr, value);
}

//------------------------------------------------------------------- AND

template <typename DType>
__global__ void and_cuda_array(size_t size, DType* cu_dest_ptr, DType value)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cu_dest_ptr[idx] &= value;
	}
}

template void and_cuda_array_launcher<bool>(size_t size, bool* cu_dest_ptr, bool value);
template void and_cuda_array_launcher<char>(size_t size, char* cu_dest_ptr, char value);
template void and_cuda_array_launcher<int>(size_t size, int* cu_dest_ptr, int value);
template void and_cuda_array_launcher<unsigned>(size_t size, unsigned* cu_dest_ptr, unsigned value);
template void and_cuda_array_launcher<long>(size_t size, long* cu_dest_ptr, long value);
template void and_cuda_array_launcher<size_t>(size_t size, size_t* cu_dest_ptr, size_t value);

template <typename DType>
void and_cuda_array_launcher(size_t size, DType* cu_dest_ptr, DType value)
{
	and_cuda_array << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (size, cu_dest_ptr, value);
}

//------------------------------------------------------------------- OR

template <typename DType>
__global__ void or_cuda_array(size_t size, DType* cu_dest_ptr, DType value)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		cu_dest_ptr[idx] |= value;
	}
}

template void or_cuda_array_launcher<bool>(size_t size, bool* cu_dest_ptr, bool value);
template void or_cuda_array_launcher<char>(size_t size, char* cu_dest_ptr, char value);
template void or_cuda_array_launcher<int>(size_t size, int* cu_dest_ptr, int value);
template void or_cuda_array_launcher<unsigned>(size_t size, unsigned* cu_dest_ptr, unsigned value);
template void or_cuda_array_launcher<long>(size_t size, long* cu_dest_ptr, long value);
template void or_cuda_array_launcher<size_t>(size_t size, size_t* cu_dest_ptr, size_t value);

template <typename DType>
void or_cuda_array_launcher(size_t size, DType* cu_dest_ptr, DType value)
{
	or_cuda_array << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (size, cu_dest_ptr, value);
}