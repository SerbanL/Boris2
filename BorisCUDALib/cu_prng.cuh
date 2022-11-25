#pragma once

#include "cu_prng.h"

//-----------------------------------------

__global__ inline void seed_array_kernel(unsigned*& prn, size_t& prn_size, unsigned seed)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < prn_size) {

		prn[idx] = seed * (idx + 1);
	}
}

template void cuBorisRand<void>::seed_array(unsigned seed);

template <typename Dummy>
__host__ void cuBorisRand<Dummy>::seed_array(unsigned seed)
{
	seed_array_kernel <<< (get_gpu_value(prn_size) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (prn, prn_size, seed);
}

//-----------------------------------------

template unsigned cuBorisRand<void>::atomicLCG(unsigned* address);

template <typename Dummy>
__device__ unsigned cuBorisRand<Dummy>::atomicLCG(unsigned* address)
{
	if (!use_atomics) {

		//special case where atomic operation not required
		*address = (unsigned)1664525 * *address + (unsigned)1013904223;
		return *address;
	}

	unsigned old = *address, assumed;

	do {

		assumed = old;

		//LCG equation used to generate next number in sequence : the modulo operation is free since unsigned is 32 bits wide
		old = atomicCAS(address, assumed, (unsigned)1664525 * old + (unsigned)1013904223);

	} while (old != assumed);
	
	return old;
}

//-----------------------------------------

template unsigned cuBorisRand<void>::randi(void);

//unsigned integer value out : 0 to 2^32 - 1
template <typename Dummy>
__device__ unsigned cuBorisRand<Dummy>::randi(void)
{
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) % prn_size;

	return atomicLCG(prn + idx);
}

template float cuBorisRand<void>::rand(void);

//floating point value out in interval [0, 1]
template <typename Dummy>
__device__ float cuBorisRand<Dummy>::rand(void)
{
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) % prn_size;

	return (float)atomicLCG(prn + idx) / (unsigned)4294967295;
}

template float cuBorisRand<void>::rand_gauss(float mean, float std);

//Box-Muller transform to generate Gaussian distribution from uniform distribution
template <typename Dummy>
__device__ float cuBorisRand<Dummy>::rand_gauss(float mean, float std)
{
	//Not exactly the usual Box-Muller transform : generate a value every time, still Gaussian distribution.
	//This saves having to allocate extra spaces
	float u1, u2;
	do {
		u1 = this->rand();
		u2 = this->rand();
	} while (u1 <= 1e-10);

	float z0;
	z0 = sqrt(-2.0 * log(u1)) * cos((float)TWO_PI * u2);

	return z0 * std + mean;
}
