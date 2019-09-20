#include "ExchangeCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_EXCHANGE

#include "Mesh_FerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

#include "BorisCUDALib.cuh"

__global__ void ExchangeCUDA_UpdateField(ManagedMeshCUDA& cuMesh, cuBReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Hexch = cuReal3();

		if (M.is_not_empty(idx)) {

			cuBReal Ms = *cuMesh.pMs;
			cuBReal A = *cuMesh.pA;
			
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pA, A);

			Hexch = 2 * A * M.delsq_neu(idx) / ((cuBReal)MU0 * Ms * Ms);
			
			if (do_reduction) {

				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ = -(cuBReal)MU0 * M[idx] * Hexch / (2 * non_empty_cells);
			}
		}

		Heff[idx] += Hexch;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void Exch_6ngbr_NeuCUDA::UpdateField(void)
{
	if (pMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();

		ExchangeCUDA_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, true);
	}
	else {

		ExchangeCUDA_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, false);
	}
	
	if (pMeshCUDA->GetMeshExchangeCoupling()) CalculateExchangeCoupling(energy);
}

#endif

#endif