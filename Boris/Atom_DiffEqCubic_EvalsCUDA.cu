#include "Atom_DiffEqCubicCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_ATOM_CUBIC

//defines evaluation methods kernel launchers

#include "BorisCUDALib.cuh"

#include "Atom_MeshParamsControlCUDA.h"

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

//Save current moments in sM VECs (e.g. useful to reset dM / dt calculation)
__global__ void SaveMoments_Cubic_kernel(cuVEC_VC<cuReal3>& M1, cuVEC<cuReal3>& sM1)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < M1.linear_size()) {

		sM1[idx] = M1[idx];
	}
}

//Restore magnetization after a failed step for adaptive time-step methods
void Atom_DifferentialEquationCubicCUDA::SaveMoments(void)
{
	SaveMoments_Cubic_kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->M1, sM1);
}

//-----------------------------------------

__global__ void RenormalizeMoments_Cubic_kernel(ManagedAtom_MeshCUDA& cuaMesh)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < cuaMesh.pM1->linear_size()) {

		if (cuaMesh.pM1->is_not_empty(idx)) {

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s);

			if (mu_s) (*cuaMesh.pM1)[idx].renormalize(mu_s);
		}
	}
}

//Restore magnetization after a failed step for adaptive time-step methods
void Atom_DifferentialEquationCubicCUDA::RenormalizeMoments(void)
{
	RenormalizeMoments_Cubic_kernel <<< (paMeshCUDA->M1()->linear_size_cpu() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh);
}

//-----------------------------------------

#endif
#endif