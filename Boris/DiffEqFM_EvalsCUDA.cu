#include "DiffEqFMCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_FERROMAGNETIC

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//-----------------------------------------

__global__ void Restoremagnetization_FM_kernel(cuVEC_VC<cuReal3>& M, cuVEC<cuReal3>& sM1)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < M.linear_size()) {

		M[idx] = sM1[idx];
	}
}

//Restore magnetization after a failed step for adaptive time-step methods
void DifferentialEquationFMCUDA::Restoremagnetization(void)
{
	Restoremagnetization_FM_kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->M, sM1);
}

//-----------------------------------------

#endif
#endif