#include "Atom_DiffEqCubicCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_ATOM_CUBIC

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

//-----------------------------------------

__global__ void RestoreMoments_Cubic_kernel(cuVEC_VC<cuReal3>& M1, cuVEC<cuReal3>& sM1)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < M1.linear_size()) {

		M1[idx] = sM1[idx];
	}
}

//Restore magnetization after a failed step for adaptive time-step methods
void Atom_DifferentialEquationCubicCUDA::RestoreMoments(void)
{
	RestoreMoments_Cubic_kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->M1, sM1);
}

//-----------------------------------------

#endif
#endif