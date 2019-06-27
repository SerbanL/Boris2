#include "RoughnessCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_ROUGHNESS

#include "BorisCUDALib.cuh"

#include "Mesh_FerromagneticCUDA.h"

__global__ void RoughnessCUDA_UpdateField_Kernel(ManagedMeshCUDA& cuMesh, cuVEC<cuReal3>& Fmul_rough, cuVEC<cuReal3>& Fomul_rough, cuReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Hrough = cuReal3();

		if (M.is_not_empty(idx)) {

			Hrough = cuReal33(
				cuReal3(Fmul_rough[idx].x, Fomul_rough[idx].x, Fomul_rough[idx].y),
				cuReal3(Fomul_rough[idx].x, Fmul_rough[idx].y, Fomul_rough[idx].z),
				cuReal3(Fomul_rough[idx].y, Fomul_rough[idx].z, Fmul_rough[idx].z)) * M[idx];

			if (do_reduction) {

				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ = -(cuReal)MU0 * M[idx] * Hrough / (2 * non_empty_cells);
			}
		}

		Heff[idx] += Hrough;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

void RoughnessCUDA::UpdateField(void)
{
	if (pMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();

		RoughnessCUDA_UpdateField_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, Fmul_rough, Fomul_rough, energy, true);
	}
	else {

		RoughnessCUDA_UpdateField_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, Fmul_rough, Fomul_rough, energy, false);
	}
}

#endif

#endif