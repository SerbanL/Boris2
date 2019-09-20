#include "ZeemanCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_ZEEMAN

#include "BorisCUDALib.cuh"

#include "Mesh_FerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

__global__ void ZeemanCUDA_UpdateField(ManagedMeshCUDA& cuMesh, cuReal3& Ha, cuBReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuBReal cHA = *cuMesh.pcHA;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHA, cHA);

		Heff[idx] += (cHA * Ha);

		if (do_reduction) {

			int non_empty_cells = M.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MU0 * M[idx] * (cHA * Ha) / non_empty_cells;
		}
	}

	if(do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void ZeemanCUDA::UpdateField(void)
{
	if (pMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();

		ZeemanCUDA_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, Ha, energy, true);
	}
	else ZeemanCUDA_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, Ha, energy, false);
}

#endif

#endif