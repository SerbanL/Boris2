#include "Demag_NCUDA.h"
#include "Mesh_FerromagneticCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_DEMAG_N

#include "BorisCUDALib.cuh"

__global__ void Demag_NCUDA_UpdateField(ManagedMeshCUDA& cuMesh, cuReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	MatPCUDA<cuReal2, cuReal>& Nxy = *cuMesh.pNxy;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Heff_value = cuReal3();

		if (M.is_not_empty(idx)) {

			Heff_value = cuReal3(-cuReal2(Nxy).x * M[idx].x, -cuReal2(Nxy).y * M[idx].y, -(1 - cuReal2(Nxy).x - cuReal2(Nxy).y) * M[idx].z);

			if (do_reduction) {

				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ = -(cuReal)MU0 * M[idx] * Heff_value / (2 * non_empty_cells);
			}
		}

		Heff[idx] += Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void Demag_NCUDA::UpdateField(void)
{
	if (pMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();

		Demag_NCUDA_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, true);
	}
	else {
		
		Demag_NCUDA_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, false);
	}
}

#endif

#endif