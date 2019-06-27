#pragma once

#include "alloc_cpy.h"
#include "cuBRealComplex.h"

//Very simple pseudo-random number generator based on LCG (linear congruential generator) with textbook values (modulus 2^32 used).

//linear congruential generator values obtained from "Numerical Recipes" - W.H. Press et al., Cambridge University Press, 3rd edition (2007) :
//Xn+1 = (a * Xn + c) mod m, where:
//a = 1664525
//c = 1013904223
//m = 2^32

class cuBorisRand {

	//array of random numbers - access to these is done using atomic operations with the index obtained from the kernel total thread index
	unsigned* prn;

	//the prn array size
	size_t prn_size;

private:

	//set starting seed values - prn array must have memory allocated
	__host__ void seed_array(unsigned seed)
	{
		std::vector<unsigned> prn_cpu(get_gpu_value(prn_size));

		for (int idx = 0; idx < prn_cpu.size(); idx++) {

			prn_cpu[idx] = seed * (idx + 1);
		}

		cpu_to_gpu_managed(prn, prn_cpu.data(), prn_cpu.size());
	}

	//access prn element at address and store next LCG value using an atomic operation
	__device__ unsigned atomicLCG(unsigned* address);

public:

	//----------------------------------------- cu_obj "managed" constructor / destructor

	__host__ void construct_cu_obj(void)
	{
		clear();
	}

	//make a prng with starting seed and memory size : in general recommended mem_size is the kernel size the prng will be used in divided by 128.
	//Larger memory size might not result in speed increase, lower memory size starts to slow the prng due to atomic operations used.
	__host__ void construct_cu_obj(unsigned seed, size_t mem_size)
	{
		initialize(seed, mem_size);
	}

	__host__ void destruct_cu_obj(void)
	{
		gpu_free_managed(prn);
	}

	//----------------------------------------- resize / reseed

	//make a prng with starting seed and memory size : in general recommended mem_size is the kernel size the prng will be used in divided by 128.
	//Larger memory size might not result in speed increase, lower memory size starts to slow the prng due to atomic operations used.
	__host__ cudaError_t initialize(unsigned seed, size_t mem_size)
	{
		if (!mem_size) mem_size = 1;

		cudaError_t error = gpu_alloc_managed(prn, mem_size);
		if (error == cudaSuccess) {

			error = set_gpu_value(prn_size, mem_size);
		}
		else {

			gpu_alloc_managed(prn, 1);
			set_gpu_value(prn_size, (size_t)1);
		}

		seed_array(seed);

		return error;
	}

	//resize array to just 1 element with a default seed of 0
	__host__ void clear(void)
	{
		gpu_alloc_managed(prn);
		set_gpu_value(prn_size, (size_t)1);

		seed_array(0);
	}

	//----------------------------------------- prn generation

	//unsigned integer value out : 0 to 2^32 - 1
	__device__ unsigned randi(void);
	
	//floating point value out in interval [0, 1]
	__device__ float rand(void);
};