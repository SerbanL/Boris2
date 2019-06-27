#pragma once

#include <cuda_runtime.h>

#include "cuTypes.h"

//----------------------------------------- Reduction helpers

//----------------------------------------- Maximum

__device__ inline void atomicMax(float* result, float value)
{
	//atomicMax is not available for floating point, so this is implemented using atomicCAS.
	//This method is recommended by nvidia, see : http://docs.nvidia.com/cuda/cuda-c-programming-guide/#axzz44PKOq6jJ

#if __CUDA_ARCH__ >= 200 
	float old = *result, assumed;

	do {

		assumed = old;

		old = atomicCAS((unsigned*)result, __float_as_int(assumed), __float_as_int(fmaxf(value, assumed)));

	} while (old != assumed);
#endif
	return;
}

template <typename Type>
__device__ void atomicMax(cuVAL2<Type>* result, const cuVAL2<Type>& value)
{
	atomicMax(&(result->x), value.x);
	atomicMax(&(result->y), value.y);
}

template <typename Type>
__device__ void atomicMax(cuVAL3<Type>* result, const cuVAL3<Type>& value)
{
	atomicMax(&(result->x), value.x);
	atomicMax(&(result->y), value.y);
	atomicMax(&(result->z), value.z);
}

//----------------------------------------- Minimum

__device__ inline void atomicMin(float *result, float value)
{
#if __CUDA_ARCH__ >= 200 
	float old = *result, assumed;

	do {

		assumed = old;

		old = atomicCAS((unsigned*)result, __float_as_int(assumed), __float_as_int(fminf(value, assumed)));

	} while (old != assumed);
#endif
	return;
}

template <typename Type>
__device__ void atomicMin(cuVAL2<Type>* result, const cuVAL2<Type>& value)
{
	atomicMin(&(result->x), value.x);
	atomicMin(&(result->y), value.y);
}

template <typename Type>
__device__ void atomicMin(cuVAL3<Type>* result, const cuVAL3<Type>& value)
{
	atomicMin(&(result->x), value.x);
	atomicMin(&(result->y), value.y);
	atomicMin(&(result->z), value.z);
}

//----------------------------------------- Add

template <typename Type>
__device__ void atomicAdd(cuVAL2<Type>* result, const cuVAL2<Type>& value)
{
	atomicAdd(&(result->x), value.x);
	atomicAdd(&(result->y), value.y);
}

template <typename Type>
__device__ void atomicAdd(cuVAL3<Type>* result, const cuVAL3<Type>& value)
{
	atomicAdd(&(result->x), value.x);
	atomicAdd(&(result->y), value.y);
	atomicAdd(&(result->z), value.z);
}

#if __CUDA_ARCH__ < 600
template <typename Type_ = double>
__device__ void atomicAdd(Type_* address, Type_ val)
{
#if __CUDA_ARCH__ >= 200
	unsigned long long int* address_as_ull =
		(unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;

	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed,
			__double_as_longlong(val +
				__longlong_as_double(assumed)));

		// Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
	} while (assumed != old);

	return;
#endif
}
#endif