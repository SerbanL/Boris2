#pragma once

#include "cu_prng.h"

//-----------------------------------------

__device__ inline unsigned cuBorisRand::atomicLCG(unsigned* address)
{
	unsigned old = *address, assumed;

	do {

		assumed = old;

		//LCG equation used to generate next number in sequence : the modulo operation is free since unsigned is 32 bits wide
		old = atomicCAS(address, assumed, (unsigned)1664525 * old + (unsigned)1013904223);

	} while (old != assumed);
	
	return old;
}

//-----------------------------------------

//unsigned integer value out : 0 to 2^32 - 1
__device__ inline unsigned cuBorisRand::randi(void)
{
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) % prn_size;

	return atomicLCG(prn + idx);
}

//floating point value out in interval [0, 1]
__device__ inline float cuBorisRand::rand(void)
{
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) % prn_size;

	return (float)atomicLCG(prn + idx) / (unsigned)4294967295;
}

//Box-Muller transform to generate Gaussian distribution from uniform distribution
__device__ inline float cuBorisRand::rand_gauss(float mean, float std)
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
