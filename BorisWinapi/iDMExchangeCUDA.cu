#include "iDMExchangeCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_IDMEXCHANGE

#include "BorisCUDALib.cuh"

#include "Mesh_FerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

__global__ void iDMExchangeCUDA_UpdateField(ManagedMeshCUDA& cuMesh, cuReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Hexch = cuReal3();

		if (M.is_not_empty(idx)) {

			cuReal Ms = *cuMesh.pMs;
			cuReal A = *cuMesh.pA;
			cuReal D = *cuMesh.pD;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pA, A, *cuMesh.pD, D);

			if (M.is_plane_interior(idx)) {

				//interior point : can use cheaper neu versions

				//direct exchange contribution
				Hexch = 2 * A * M.delsq_neu(idx) / ((cuReal)MU0 * Ms * Ms);

				//Dzyaloshinskii-Moriya interfacial exchange contribution

				//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
				cuReal33 Mdiff = M.grad_neu(idx);

				//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
				Hexch += -2 * D * cuReal3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y) / ((cuReal)MU0 * Ms * Ms);
			}
			else {

				//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
				cuReal3 bnd_dm_dx = (D / (2 * A * Ms)) * cuReal3(M[idx].z, 0, -M[idx].x);
				cuReal3 bnd_dm_dy = (D / (2 * A * Ms)) * cuReal3(0, M[idx].z, -M[idx].y);
				cuReal33 bnd_nneu = cuReal33(bnd_dm_dx, bnd_dm_dy, cuReal3());

				//direct exchange contribution
				Hexch = 2 * A * M.delsq_nneu(idx, bnd_nneu) / ((cuReal)MU0 * Ms * Ms);

				//Dzyaloshinskii-Moriya interfacial exchange contribution

				//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
				cuReal33 Mdiff = M.grad_nneu(idx, bnd_nneu);

				//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
				Hexch += -2 * D * cuReal3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y) / ((cuReal)MU0 * Ms * Ms);
			}

			if (do_reduction) {

				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ = -(cuReal)MU0 * M[idx] * Hexch / (2 * non_empty_cells);
			}
		}

		Heff[idx] += Hexch;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void iDMExchangeCUDA::UpdateField(void)
{
	if (pMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();

		iDMExchangeCUDA_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, true);
	}
	else {

		iDMExchangeCUDA_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, false);
	}

	if (pMeshCUDA->GetMeshExchangeCoupling()) CalculateExchangeCoupling(energy);
}

#endif

#endif