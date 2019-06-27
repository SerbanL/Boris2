#include "SDemagCUDA_Demag.h"

#if COMPILECUDA == 1

#ifdef MODULE_SDEMAG

#include "BorisCUDALib.cuh"

#include "Mesh_FerromagneticCUDA.h"

//-------------------Getters

__global__ void Add_Energy_Kernel(cuReal& energy, cuReal& total_energy)
{
	if (threadIdx.x == 0) total_energy += energy;
}

//add energy in this module to a running total
void SDemagCUDA_Demag::Add_Energy(cu_obj<cuReal>& total_energy)
{
	Add_Energy_Kernel <<< (1 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (energy, total_energy);
}

#endif

#endif