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
	cuVEC_VC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Heff.linear_size()) {

		if (M.is_not_empty(idx)) {

			cuBReal grel = *cuMesh.pgrel;
			cuBReal Ms = *cuMesh.pMs;
			cuBReal SHA = *cuMesh.pSHA;
			cuBReal flSOT = *cuMesh.pflSOT;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pgrel, grel, *cuMesh.pMs, Ms, *cuMesh.pSHA, SHA, *cuMesh.pflSOT, flSOT);

			if (cuIsNZ(grel)) {

				cuBReal a_const = -(SHA * MUB_E / (GAMMA * grel)) / (Ms * Ms * (M.rect.e.z - M.rect.s.z));

				int idx_E = E.position_to_cellidx(M.cellidx_to_position(idx));
				cuReal3 p_vec = cuReal3(0, 0, 1) ^ (elC[idx_E] * E[idx_E]);

				Heff[idx] += a_const * ((M[idx] ^ p_vec) + flSOT * Ms * p_vec);
			}
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