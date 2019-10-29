#include "DiffEqAFMCUDA.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"

//-----------------------------------------

__global__ void RestoreMagnetisation_kernel(cuVEC_VC<cuReal3>& M, cuVEC<cuReal3>& sM1, cuVEC_VC<cuReal3>& M2, cuVEC<cuReal3>& sM1_2)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < M.linear_size()) {

		M[idx] = sM1[idx];
		M2[idx] = sM1_2[idx];
	}
}

//Restore magnetisation after a failed step for adaptive time-step methods
void DifferentialEquationAFMCUDA::RestoreMagnetisation(void)
{
	RestoreMagnetisation_kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->M, sM1, pMeshCUDA->M2, sM1_2);
}

//-----------------------------------------

#endif