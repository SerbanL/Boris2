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
