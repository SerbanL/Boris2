#include "ModulesCUDA.h"

#if COMPILECUDA == 1

__global__ void ZeroEnergy_kernel(cuReal& energy)
{
	if (threadIdx.x == 0) energy = 0.0;
}

void ModulesCUDA::ZeroEnergy(void)
{
	ZeroEnergy_kernel <<< 1, CUDATHREADS >>> (energy);
}

#endif