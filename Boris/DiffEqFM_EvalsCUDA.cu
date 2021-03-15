#include "DiffEqFMCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_FERROMAGNETIC

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

#include "MeshParamsControlCUDA.h"

//-----------------------------------------

__global__ void RestoreMagnetization_FM_kernel(cuVEC_VC<cuReal3>& M, cuVEC<cuReal3>& sM1)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < M.linear_size()) {

		M[idx] = sM1[idx];
	}
}

//Restore magnetization after a failed step for adaptive time-step methods
void DifferentialEquationFMCUDA::RestoreMagnetization(void)
{
	RestoreMagnetization_FM_kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->M, sM1);
}

//-----------------------------------------

__global__ void RenormalizeMagnetization_FM_kernel(ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;


	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			cuBReal Ms = *cuMesh.pMs;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms);

			(*cuMesh.pM)[idx].renormalize(Ms);
		}
	}
}

//Restore magnetization after a failed step for adaptive time-step methods
void DifferentialEquationFMCUDA::RenormalizeMagnetization(void)
{
	RenormalizeMagnetization_FM_kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh);
}

//-----------------------------------------

#endif
#endif