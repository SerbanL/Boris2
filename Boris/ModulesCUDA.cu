#include "ModulesCUDA.h"

#if COMPILECUDA == 1

__global__ void ZeroEnergy_kernel(cuBReal& energy, size_t& points_count)
{
	if (threadIdx.x == 0) energy = 0.0;
	if (threadIdx.x == 1) points_count = 0;
}

void ModulesCUDA::ZeroEnergy(void)
{
	ZeroEnergy_kernel <<< 1, CUDATHREADS >>> (energy, points_count);
}

#endif