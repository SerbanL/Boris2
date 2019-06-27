#include "AnisotropyCubiCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_ANICUBI

#include "Mesh_FerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

#include "BorisCUDALib.cuh"

__global__ void Anisotropy_CubicCUDA_UpdateField(ManagedMeshCUDA& cuMesh, cuReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Heff_value = cuReal3();

		if (M.is_not_empty(idx)) {

			cuReal Ms = *cuMesh.pMs;
			cuReal K1 = *cuMesh.pK1;
			cuReal K2 = *cuMesh.pK2;
			cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
			cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;

			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pK1, K1, *cuMesh.pK2, K2, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2);

			//vector product of ea1 and ea2 : the third orthogonal axis
			cuReal3 mcanis_ea3 = mcanis_ea1 ^ mcanis_ea2;

			//calculate m.ea1, m.ea2 and m.ea3 dot products
			cuReal d1 = (M[idx] * mcanis_ea1) / Ms;
			cuReal d2 = (M[idx] * mcanis_ea2) / Ms;
			cuReal d3 = (M[idx] * mcanis_ea3) / Ms;

			//terms for K1 contribution
			cuReal a1 = d1 * (d2*d2 + d3 * d3);
			cuReal a2 = d2 * (d1*d1 + d3 * d3);
			cuReal a3 = d3 * (d1*d1 + d2 * d2);

			//terms for K2 contribution
			cuReal d123 = d1 * d2*d3;

			cuReal b1 = d123 * d2*d3;
			cuReal b2 = d123 * d1*d3;
			cuReal b3 = d123 * d1*d2;

			//update effective field with the anisotropy field
			Heff_value = cuReal3(
				(-2 * K1 / ((cuReal)MU0*Ms)) * (mcanis_ea1.i * a1 + mcanis_ea2.i * a2 + mcanis_ea3.i * a3)
				+ (-2 * K2 / ((cuReal)MU0*Ms)) * (mcanis_ea1.i * b1 + mcanis_ea2.i * b2 + mcanis_ea3.i * b3),

				(-2 * K1 / ((cuReal)MU0*Ms)) * (mcanis_ea1.j * a1 + mcanis_ea2.j * a2 + mcanis_ea3.j * a3)
				+ (-2 * K2 / ((cuReal)MU0*Ms)) * (mcanis_ea1.j * b1 + mcanis_ea2.j * b2 + mcanis_ea3.j * b3),

				(-2 * K1 / ((cuReal)MU0*Ms)) * (mcanis_ea1.k * a1 + mcanis_ea2.k * a2 + mcanis_ea3.k * a3)
				+ (-2 * K2 / ((cuReal)MU0*Ms)) * (mcanis_ea1.k * b1 + mcanis_ea2.k * b2 + mcanis_ea3.k * b3)
			);

			if (do_reduction) {

				//update energy (E/V)		
				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ = (K1 * (d1*d1*d2*d2 + d1 * d1*d3*d3 + d2 * d2*d3*d3) + K2 * d123*d123) / non_empty_cells;
			}
		}

		Heff[idx] += Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void Anisotropy_CubicCUDA::UpdateField(void)
{
	if (pMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();

		Anisotropy_CubicCUDA_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, true);
	}
	else {

		Anisotropy_CubicCUDA_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, false);
	}
}

#endif

#endif