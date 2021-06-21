#include "DiffEqAFMCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

#include "BorisCUDALib.cuh"

#include "MeshParamsControlCUDA.h"

//-----------------------------------------

__global__ void Restoremagnetization_AFM_kernel(cuVEC_VC<cuReal3>& M, cuVEC<cuReal3>& sM1, cuVEC_VC<cuReal3>& M2, cuVEC<cuReal3>& sM1_2)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < M.linear_size()) {

		M[idx] = sM1[idx];
		M2[idx] = sM1_2[idx];
	}
}

//Restore magnetization after a failed step for adaptive time-step methods
void DifferentialEquationAFMCUDA::RestoreMagnetization(void)
{
	Restoremagnetization_AFM_kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->M, sM1, pMeshCUDA->M2, sM1_2);
}

//-----------------------------------------

__global__ void RenormalizeMagnetization_AFM_kernel(ManagedMeshCUDA& cuMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < cuMesh.pM->linear_size()) {

		if (cuMesh.pM->is_not_empty(idx)) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM);

			if (Ms_AFM.i) (*cuMesh.pM)[idx].renormalize(Ms_AFM.i);
			if (Ms_AFM.j) (*cuMesh.pM2)[idx].renormalize(Ms_AFM.j);
		}
	}
}

//Restore magnetization after a failed step for adaptive time-step methods
void DifferentialEquationAFMCUDA::RenormalizeMagnetization(void)
{
	RenormalizeMagnetization_AFM_kernel <<< (pMeshCUDA->M()->linear_size_cpu() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh);
}

//-----------------------------------------

#endif
#endif