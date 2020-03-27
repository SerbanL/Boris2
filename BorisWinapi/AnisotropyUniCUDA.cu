#include "AnisotropyCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_ANIUNI

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"
#include "MeshDefs.h"

__global__ void Anisotropy_UniaxialCUDA_FM_UpdateField(ManagedMeshCUDA& cuMesh, cuBReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Heff_value = cuReal3();

		if (M.is_not_empty(idx)) {

			cuBReal Ms = *cuMesh.pMs;
			cuBReal K1 = *cuMesh.pK1;
			cuBReal K2 = *cuMesh.pK2;
			cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;

			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pK1, K1, *cuMesh.pK2, K2, *cuMesh.pmcanis_ea1, mcanis_ea1);

			//calculate m.ea dot product
			cuBReal dotprod = (M[idx] * mcanis_ea1) / Ms;

			//update effective field with the anisotropy field
			Heff_value = (2 / ((cuBReal)MU0 * Ms)) * dotprod * (K1 + 2 * K2 * (1 - dotprod * dotprod)) * mcanis_ea1;

			if (do_reduction) {

				//update energy (E/V) = K1 * sin^2(theta) + K2 * sin^4(theta) = K1 * [ 1 - dotprod*dotprod ] + K2 * [1 - dotprod * dotprod]^2
				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ = ((K1 + K2 * (1 - dotprod * dotprod)) * (1 - dotprod * dotprod)) / non_empty_cells;
			}
		}

		Heff[idx] += Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

__global__ void Anisotropy_UniaxialCUDA_AFM_UpdateField(ManagedMeshCUDA& cuMesh, cuBReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Heff_value = cuReal3();
		cuReal3 Heff2_value = cuReal3();

		if (M.is_not_empty(idx)) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuReal2 K1_AFM = *cuMesh.pK1_AFM;
			cuReal2 K2_AFM = *cuMesh.pK2_AFM;
			cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;

			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pK1_AFM, K1_AFM, *cuMesh.pK2_AFM, K2_AFM, *cuMesh.pmcanis_ea1, mcanis_ea1);

			//calculate m.ea dot product
			cuBReal dotprod = (M[idx] * mcanis_ea1) / Ms_AFM.i;
			cuBReal dotprod2 = (M2[idx] * mcanis_ea1) / Ms_AFM.j;

			//update effective field with the anisotropy field
			Heff_value = (2 / ((cuBReal)MU0 * Ms_AFM.i)) * dotprod * (K1_AFM.i + 2 * K2_AFM.i * (1 - dotprod * dotprod)) * mcanis_ea1;
			Heff2_value = (2 / ((cuBReal)MU0 * Ms_AFM.j)) * dotprod2 * (K1_AFM.j + 2 * K2_AFM.j * (1 - dotprod2 * dotprod2)) * mcanis_ea1;

			if (do_reduction) {

				//update energy (E/V) = K1 * sin^2(theta) + K2 * sin^4(theta) = K1 * [ 1 - dotprod*dotprod ] + K2 * [1 - dotprod * dotprod]^2
				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ = ((K1_AFM.i + K2_AFM.i * (1 - dotprod * dotprod)) * (1 - dotprod * dotprod) + (K1_AFM.j + K2_AFM.j * (1 - dotprod2 * dotprod2)) * (1 - dotprod2 * dotprod2)) / (2 * non_empty_cells);
			}
		}

		Heff[idx] += Heff_value;
		Heff2[idx] += Heff2_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void Anisotropy_UniaxialCUDA::UpdateField(void)
{
	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		//anti-ferromagnetic mesh

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			Anisotropy_UniaxialCUDA_AFM_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, true);
		}
		else {

			Anisotropy_UniaxialCUDA_AFM_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, false);
		}
	}
	else {

		//ferromagnetic mesh

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			Anisotropy_UniaxialCUDA_FM_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, true);
		}
		else {

			Anisotropy_UniaxialCUDA_FM_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, false);
		}
	}
}

#endif

#endif