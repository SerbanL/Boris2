#include "SDemagCUDA_Demag.h"

#if COMPILECUDA == 1

#ifdef MODULE_SDEMAG

#include "BorisCUDALib.cuh"

//-------------------Getters

__global__ void Add_Energy_Kernel(cuBReal& energy, cuBReal& total_energy)
{
	if (threadIdx.x == 0) total_energy += energy;
}

//add energy in this module to a running total
void SDemagCUDA_Demag::Add_Energy(cu_obj<cuBReal>& total_energy)
{
	Add_Energy_Kernel <<< (1 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (energy, total_energy);
}

#endif

#endif