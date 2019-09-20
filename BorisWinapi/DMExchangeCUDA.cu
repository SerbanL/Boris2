#include "DMExchangeCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_DMEXCHANGE

#include "BorisCUDALib.cuh"

#include "Mesh_FerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

__global__ void DMExchangeCUDA_UpdateField(ManagedMeshCUDA& cuMesh, cuBReal& energy, bool do_reduction)
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
			cuBReal D = *cuMesh.pD;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pA, A, *cuMesh.pD, D);

			if (M.is_interior(idx)) {

				//interior point : can use cheaper neu versions

				//direct exchange contribution
				Hexch = 2 * A * M.delsq_neu(idx) / ((cuBReal)MU0 * Ms * Ms);

				//Dzyaloshinskii-Moriya exchange contribution

				//Hdm, ex = -2D / (mu0*Ms) * curl m
				Hexch += -2 * D * M.curl_neu(idx) / ((cuBReal)MU0 * Ms * Ms);
			}
			else {

				//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
				cuReal3 bnd_dm_dx = (D / (2 * A)) * cuReal3(0, -M[idx].z, M[idx].y);
				cuReal3 bnd_dm_dy = (D / (2 * A)) * cuReal3(M[idx].z, 0, -M[idx].x);
				cuReal3 bnd_dm_dz = (D / (2 * A)) * cuReal3(-M[idx].y, M[idx].x, 0);
				cuReal33 bnd_nneu = cuReal33(bnd_dm_dx, bnd_dm_dy, bnd_dm_dz);

				//direct exchange contribution
				Hexch = 2 * A * M.delsq_nneu(idx, bnd_nneu) / ((cuBReal)MU0 * Ms * Ms);

				//Dzyaloshinskii-Moriya exchange contribution

				//Hdm, ex = -2D / (mu0*Ms) * curl m
				Hexch += -2 * D * M.curl_nneu(idx, bnd_nneu) / ((cuBReal)MU0 * Ms * Ms);
			}

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

void DMExchangeCUDA::UpdateField(void)
{
	if (pMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();

		DMExchangeCUDA_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, true);
	}
	else {

		DMExchangeCUDA_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, false);
	}

	if (pMeshCUDA->GetMeshExchangeCoupling()) CalculateExchangeCoupling(energy);
}

#endif

#endif