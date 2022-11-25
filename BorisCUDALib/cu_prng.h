#pragma once

#include "alloc_cpy.h"
#include "cuBRealComplex.h"

//Very simple pseudo-random number generator based on LCG (linear congruential generator) with textbook values (modulus 2^32 used).

//linear congruential generator values obtained from "Numerical Recipes" - W.H. Press et al., Cambridge University Press, 3rd edition (2007) :
//Xn+1 = (a * Xn + c) mod m, where:
//a = 1664525
//c = 1013904223
//m = 2^32

template <typename Dummy = void>
class cuBorisRand {

//if kernel size is below this threshold then don't use atomic operations, i.e. set divisor equal to 1 by default
#define PRNG_KERNELSIZETHERSHOLD	1048576
//default divisor to use : divide the kernel size by this to determine prn_size, within which access is done using atomic operations
#define PRNG_DIVISOR	128

	//array of random numbers - access to these is done using atomic operations with the index obtained from the kernel total thread index
	unsigned* prn;

	//the prn array size
	size_t prn_size;

	//if prn_size is set to be same as the kernel size this PRNG will be used in, then atomic operations not needed
	bool use_atomics;

private:

	//access prn element at address and store next LCG value using an atomic operation
	__device__ unsigned atomicLCG(unsigned* address);

	//set starting seed values - prn array must have memory allocated
	__host__ void seed_array(unsigned seed);

public:

	//----------------------------------------- cu_obj "managed" constructor / destructor

	__host__ void construct_cu_obj(void)
	{
		clear();
	}

	//make a prng with starting seed and memory size prn_size = ker_size / divisor
	//in general recommended memory size is the kernel size (ker_size) the prng will be used in divided by 128 (PRNG_DIVISOR) - divisor.
	//at small kernel dimensions, using divisor > 1 will result in slow performance, so if kernel size is below PRNG_KERNELSIZETHERSHOLD divisor should be equal to 1, otherwise use PRNG_DIVISOR
	//default value of divisor = 0 means this pattern will be used, but user can set divisor to any value required
	//divisor = 1 is a special mode, where the prn arrray size is the same as the kernel size. This is useful if the same sequence of random numbers must be generated every time, and in this case there's no need for atomic operations.
	//(if divisor > 1, then even though same sequence of random numbers is generated for given seed, it's not guaranteed they'll be the same in every cell every time, due to access order to atomic operation)
	__host__ void construct_cu_obj(unsigned seed, size_t ker_size, size_t divisor = 0)
	{
		initialize(seed, ker_size, divisor);
	}

	__host__ void destruct_cu_obj(void)
	{
		gpu_free_managed(prn);
	}

	//----------------------------------------- resize / reseed

	//make a prng with starting seed and memory size : in general recommended mem_size is the kernel size the prng will be used in divided by 128.
	//Larger memory size might not result in speed increase, lower memory size starts to slow the prng due to atomic operations used.
	__host__ cudaError_t initialize(unsigned seed, size_t ker_size, size_t divisor = 0)
	{
		cudaError_t error = cudaSuccess;

		//implement default divisor value
		if (divisor == 0) {

			if (ker_size < PRNG_KERNELSIZETHERSHOLD) divisor = 1;
			else divisor = PRNG_DIVISOR;
		}

		//for divisor = 1 then no need to use atomic operations
		if (divisor > 1) set_gpu_value(use_atomics, true);
		else set_gpu_value(use_atomics, false);

		size_t mem_size = ker_size / divisor;
		if (!mem_size) mem_size = 1;

		if (get_gpu_value(prn_size) != mem_size) {

			error = gpu_alloc_managed(prn, mem_size);
			if (error == cudaSuccess) {

				error = set_gpu_value(prn_size, mem_size);
			}
			else {

				gpu_alloc_managed(prn, 1);
				set_gpu_value(prn_size, (size_t)1);
			}
		}

		seed_array(seed);

		return error;
	}

	//resize array to just 1 element with a default seed of 0
	__host__ void clear(void)
	{
		gpu_alloc_managed(prn);
		set_gpu_value(prn_size, (size_t)1);
	}

	//----------------------------------------- prn generation

	//unsigned integer value out : 0 to 2^32 - 1
	__device__ unsigned randi(void);
	
	//floating point value out in interval [0, 1]
	__device__ float rand(void);

	//Box-Muller transform to generate Gaussian distribution from uniform distribution
	__device__ float rand_gauss(float mean, float std);
};