#include "SOTFieldCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_SOTFIELD

#include "BorisCUDALib.cuh"

#include "Mesh_FerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

__global__ void SOTFieldCUDA_UpdateField(ManagedMeshCUDA& cuMesh)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC_VC<cuReal3>& Jc = *cuMesh.pJc;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Heff.linear_size()) {

		if (M.is_not_empty(idx)) {

			int idx_Jc = Jc.position_to_cellidx(M.cellidx_to_position(idx));

			cuReal Ms = *cuMesh.pMs;
			cuReal SHA = *cuMesh.pSHA;
			cuReal flSOT = *cuMesh.pflSOT;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pSHA, SHA, *cuMesh.pflSOT, flSOT);

			cuReal a_const = -((cuReal)HBAR_E / (cuReal)MU0) * SHA / (2 * Ms * Ms * (M.rect.e.z - M.rect.s.z));

			cuReal3 p_vec = cuReal3(0, 0, 1) ^ Jc[idx_Jc];

			Heff[idx] += a_const * ((M[idx] ^ p_vec) + flSOT * Ms * p_vec);
		}
	}
}

//----------------------- UpdateField LAUNCHER

void SOTFieldCUDA::UpdateField(void)
{
	if (!pMeshCUDA->EComputation_Enabled()) return;

	SOTFieldCUDA_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh);
}

#endif

#endif