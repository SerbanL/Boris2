#pragma once

#include <cuda_runtime.h>

#include "atomics.cuh"
#include "cuTypes.h"

//----------------------------------------- REDUCTION : Add

//Sum elements in large_array of size array_size and store the sum in result - element_idx is the index in large_array for this thread
//The include_in_reduction flag is set by the caller after performing a check (if needed)
//Important : result must be reset to zero before calling this

template <typename FType>
__device__ void reduction_sum(int element_idx, size_t array_size, FType* large_array, FType& result, bool include_in_reduction = true)
{
	//memory shared between threads in a block
	__shared__ FType shared_memory[CUDATHREADS];

	//linearized thread id for accessing shared memory
	unsigned thread_idx = threadIdx.x;

	//store values in shared_memory so we can perform reduction on it
	if (element_idx < array_size && include_in_reduction) {

		shared_memory[thread_idx] = large_array[element_idx];
	}
	else shared_memory[thread_idx] = FType();

	//reduce values in this block using sequential addressing reduction (reduce pairs of lower half and top half elements iteratively until we are left with just the one element stored in shared_memory[0])
	
	__syncthreads();

	for (unsigned s = CUDATHREADS / 2; s > 0; s >>= 1) {

		if (thread_idx < s) {

			//summing reduction
			shared_memory[thread_idx] += shared_memory[thread_idx + s];
		}

		__syncthreads();
	}

	if (thread_idx == 0) {

		//use atomic operation to reduce all shared_memory[0] elements from different blocks to a single value. 
		atomicAdd(&result, shared_memory[0]);
	}
}

//as above, but custom block size of 1024 threads
template <typename FType>
__device__ void reduction_sum_blocksize1024(int element_idx, size_t array_size, FType* large_array, FType& result, bool include_in_reduction = true)
{
	//memory shared between threads in a block
	__shared__ FType shared_memory[1024];

	//linearized thread id for accessing shared memory
	unsigned thread_idx = threadIdx.x;

	//store values in shared_memory so we can perform reduction on it
	if (element_idx < array_size && include_in_reduction) {

		shared_memory[thread_idx] = large_array[element_idx];
	}
	else shared_memory[thread_idx] = FType();

	//reduce values in this block using sequential addressing reduction (reduce pairs of lower half and top half elements iteratively until we are left with just the one element stored in shared_memory[0])

	__syncthreads();

	for (unsigned s = blockDim.x / 2; s > 0; s >>= 1) {

		if (thread_idx < s) {

			//summing reduction
			shared_memory[thread_idx] += shared_memory[thread_idx + s];
		}

		__syncthreads();
	}

	if (thread_idx == 0) {

		//use atomic operation to reduce all shared_memory[0] elements from different blocks to a single value. 
		atomicAdd(&result, shared_memory[0]);
	}
}

//----------------------------------------- REDUCTION : Average

//Sum elements in large_array of size array_size and store the sum in result; also sum number of points that contribute - element_idx is the index in large_array for this thread
//Divide result by points_count to obtain the average value; Important : result and points_count must be reset to zero before calling this.

template <typename FType>
__device__ void reduction_avg(int element_idx, size_t array_size, FType* large_array, FType& result, size_t& points_count, bool include_in_reduction = true)
{
	//memory shared between threads in a block
	__shared__ FType shared_memory[CUDATHREADS];
	__shared__ size_t shared_memory_count[CUDATHREADS];

	//linearized thread id for accessing shared memory
	unsigned thread_idx = threadIdx.x;

	//store values in shared_memory so we can perform reduction on it
	if (element_idx < array_size && include_in_reduction) {

		shared_memory[thread_idx] = large_array[element_idx];
		shared_memory_count[thread_idx] = 1;
	}
	else {

		shared_memory[thread_idx] = FType();
		shared_memory_count[thread_idx] = 0;
	}

	//reduce values in this block using sequential addressing reduction (reduce pairs of lower half and top half elements iteratively until we are left with just the one element stored in shared_memory[0])

	__syncthreads();

	for (unsigned s = CUDATHREADS / 2; s > 0; s >>= 1) {

		if (thread_idx < s) {

			//summing reduction
			shared_memory[thread_idx] += shared_memory[thread_idx + s];
			shared_memory_count[thread_idx] += shared_memory_count[thread_idx + s];
		}

		__syncthreads();
	}

	if (thread_idx == 0) {

		//use atomic operation to reduce all shared_memory[0] elements from different blocks to a single value. 
		atomicAdd(&result, shared_memory[0]);
		atomicAdd(&points_count, shared_memory_count[0]);
	}
}

//----------------------------------------- REDUCTION : Delta

//Find maximum change from old_value as a magnitude (each thread will have a different old_value) in large_array of size array_size and store it in result - element_idx is the index in large_array for this thread
//The include_in_reduction flag is set by the caller after performing a check (if needed)
//Important : result must be reset to zero before calling this

template <typename FType>
__device__ void reduction_delta(int element_idx, size_t array_size, FType* large_array, FType old_value, cuBReal& result, bool include_in_reduction = true)
{
	//memory shared between threads in a block
	__shared__ FType shared_memory[CUDATHREADS];

	//linearized thread id for accessing shared memory
	unsigned thread_idx = threadIdx.x;

	//store values in shared_memory so we can perform reduction on it
	if (element_idx < array_size && include_in_reduction) {

		shared_memory[thread_idx] = large_array[element_idx] - old_value;
	}
	else shared_memory[thread_idx] = FType();

	//reduce values in this block using sequential addressing reduction (reduce pairs of lower half and top half elements iteratively until we are left with just the one element stored in shared_memory[0])

	__syncthreads();

	for (unsigned s = CUDATHREADS / 2; s > 0; s >>= 1) {

		if (thread_idx < s) {

			//summing reduction
			shared_memory[thread_idx] = ( cu_GetMagnitude(shared_memory[thread_idx]) > cu_GetMagnitude(shared_memory[thread_idx + s]) ? shared_memory[thread_idx] : shared_memory[thread_idx + s] );
		}

		__syncthreads();
	}

	if (thread_idx == 0) {

		//use atomic operation to reduce all shared_memory[0] elements from different blocks to a single value. 
		atomicMax(&result, (cuBReal)cu_GetMagnitude(shared_memory[0]));
	}
}

//----------------------------------------- REDUCTION : Maximum value

//Find max value in large_array of size array_size and store it in result - element_idx is the index in large_array for this thread
//Important : result must be reset to large_array[0] before calling this (not zero)

template <typename FType>
__device__ void reduction_max(int element_idx, size_t array_size, FType* large_array, FType& result, bool include_in_reduction = true)
{
	//memory shared between threads in a block
	__shared__ FType shared_memory[CUDATHREADS];

	//linearized thread id for accessing shared memory
	unsigned thread_idx = threadIdx.x;

	//store values in shared_memory so we can perform reduction on it
	if (element_idx < array_size && include_in_reduction) {

		shared_memory[thread_idx] = large_array[element_idx];
	}
	else shared_memory[thread_idx] = result;

	//reduce values in this block using sequential addressing reduction (reduce pairs of lower half and top half elements iteratively until we are left with just the one element stored in shared_memory[0])

	__syncthreads();

	for (unsigned s = CUDATHREADS / 2; s > 0; s >>= 1) {

		if (thread_idx < s) {

			//reduction : maximum
			shared_memory[thread_idx] = (shared_memory[thread_idx] > shared_memory[thread_idx + s] ? shared_memory[thread_idx] : shared_memory[thread_idx + s]);
		}

		__syncthreads();
	}

	if (thread_idx == 0) {

		//use atomic operation to reduce all shared_memory[0] elements from different blocks to a single value.
		atomicMax(&result, shared_memory[0]);
	}
}

//----------------------------------------- REDUCTION : Minimum value

//Find min value in large_array of size array_size and store it in result - element_idx is the index in large_array for this thread
//Important : result must be reset to large_array[0] before calling this (not zero)

template <typename FType>
__device__ void reduction_min(int element_idx, size_t array_size, FType* large_array, FType& result, bool include_in_reduction = true)
{
	//memory shared between threads in a block
	__shared__ FType shared_memory[CUDATHREADS];

	//linearized thread id for accessing shared memory
	unsigned thread_idx = threadIdx.x;

	//store values in shared_memory so we can perform reduction on it
	if (element_idx < array_size && include_in_reduction) {

		shared_memory[thread_idx] = large_array[element_idx];
	}
	else shared_memory[thread_idx] = result;

	//reduce values in this block using sequential addressing reduction (reduce pairs of lower half and top half elements iteratively until we are left with just the one element stored in shared_memory[0])

	__syncthreads();

	for (unsigned s = CUDATHREADS / 2; s > 0; s >>= 1) {

		if (thread_idx < s) {

			//reduction : minimum
			shared_memory[thread_idx] = (shared_memory[thread_idx] < shared_memory[thread_idx + s] ? shared_memory[thread_idx] : shared_memory[thread_idx + s]);
		}

		__syncthreads();
	}

	if (thread_idx == 0) {

		//use atomic operation to reduce all shared_memory[0] elements from different blocks to a single value. 
		atomicMin(&result, shared_memory[0]);
	}
}

//----------------------------------------- REDUCTION : Minimum and Maximum values

//Find min and max values in large_array of size array_size and store it in result - element_idx is the index in large_array for this thread
//Important : minimum and maximum output values must be reset to large_array[0] before calling this (not zero)

template <typename FType>
__device__ void reduction_minmax(int element_idx, size_t array_size, FType* large_array, FType& minimum, FType& maximum, bool include_in_reduction = true)
{
	//memory shared between threads in a block
	__shared__ FType shared_memory_min[CUDATHREADS];
	__shared__ FType shared_memory_max[CUDATHREADS];

	//linearized thread id for accessing shared memory
	unsigned thread_idx = threadIdx.x;

	//store values in shared_memory so we can perform reduction on it
	if (element_idx < array_size && include_in_reduction) {

		shared_memory_min[thread_idx] = large_array[element_idx];
		shared_memory_max[thread_idx] = large_array[element_idx];
	}
	else {

		shared_memory_min[thread_idx] = minimum;
		shared_memory_max[thread_idx] = maximum;
	}

	//reduce values in this block using sequential addressing reduction (reduce pairs of lower half and top half elements iteratively until we are left with just the one element stored in shared_memory[0])

	__syncthreads();

	for (unsigned s = CUDATHREADS / 2; s > 0; s >>= 1) {

		if (thread_idx < s) {

			//reduction : minimum
			shared_memory_min[thread_idx] = (shared_memory_min[thread_idx] < shared_memory_min[thread_idx + s] ? shared_memory_min[thread_idx] : shared_memory_min[thread_idx + s]);

			//reduction : maximum
			shared_memory_max[thread_idx] = (shared_memory_max[thread_idx] > shared_memory_max[thread_idx + s] ? shared_memory_max[thread_idx] : shared_memory_max[thread_idx + s]);
		}

		__syncthreads();
	}

	if (thread_idx == 0) {

		//use atomic operation to reduce all shared_memory[0] elements from different blocks to a single value.
		atomicMin(&minimum, shared_memory_min[0]);
		atomicMax(&maximum, shared_memory_max[0]);
	}
}
