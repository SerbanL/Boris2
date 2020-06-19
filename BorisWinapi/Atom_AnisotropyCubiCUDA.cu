#include "Atom_AnisotropyCubiCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_ANICUBI) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "Atom_MeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"
#include "MeshDefs.h"

__global__ void Atom_Anisotropy_CubiCUDA_Cubic_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, cuBReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff1.linear_size()) {

		cuReal3 Heff_value = cuReal3();

		if (M1.is_not_empty(idx)) {

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal K = *cuaMesh.pK;
			cuReal3 mcanis_ea1 = *cuaMesh.pmcanis_ea1;
			cuReal3 mcanis_ea2 = *cuaMesh.pmcanis_ea2;
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pK, K, *cuaMesh.pmcanis_ea1, mcanis_ea1, *cuaMesh.pmcanis_ea2, mcanis_ea2);

			//vector product of ea1 and ea2 : the third orthogonal axis
			cuReal3 mcanis_ea3 = mcanis_ea1 ^ mcanis_ea2;

			//calculate m.ea1, m.ea2 and m.ea3 dot products
			cuBReal d1 = (M1[idx] * mcanis_ea1) / mu_s;
			cuBReal d2 = (M1[idx] * mcanis_ea2) / mu_s;
			cuBReal d3 = (M1[idx] * mcanis_ea3) / mu_s;

			//terms for K contribution
			cuBReal a1 = d1 * (d2*d2 + d3*d3);
			cuBReal a2 = d2 * (d1*d1 + d3*d3);
			cuBReal a3 = d3 * (d1*d1 + d2*d2);

			//update effective field with the anisotropy field
			cuReal3 Heff_value = cuReal3(
				(-2 * K / (MUB_MU0*mu_s)) * (mcanis_ea1.i * a1 + mcanis_ea2.i * a2 + mcanis_ea3.i * a3),
				(-2 * K / (MUB_MU0*mu_s)) * (mcanis_ea1.j * a1 + mcanis_ea2.j * a2 + mcanis_ea3.j * a3),
				(-2 * K / (MUB_MU0*mu_s)) * (mcanis_ea1.k * a1 + mcanis_ea2.k * a2 + mcanis_ea3.k * a3)
			);

			if (do_reduction) {

				//update energy density
				cuBReal non_empty_volume = M1.get_nonempty_cells() * M1.h.dim();
				if (non_empty_volume) energy_ = K * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) / non_empty_volume;
			}
		}

		Heff1[idx] += Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void Atom_Anisotropy_CubiCUDA::UpdateField(void)
{
	if (paMeshCUDA->GetMeshType() == MESH_ATOM_CUBIC) {

		//atomistic simple-cubic mesh

		if (paMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			Atom_Anisotropy_CubiCUDA_Cubic_UpdateField <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, energy, true);
		}
		else {

			Atom_Anisotropy_CubiCUDA_Cubic_UpdateField <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, energy, false);
		}
	}
}

//-------------------Energy density methods

__global__ void Atom_Anisotropy_CubiCUDA_Cubic_GetEnergy(ManagedAtom_MeshCUDA& cuaMesh, cuBReal& energy, size_t& points_count, cuRect avRect)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M1.linear_size()) {

		if (M1.is_not_empty(idx) && avRect.contains(M1.cellidx_to_position(idx))) {

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal K = *cuaMesh.pK;
			cuReal3 mcanis_ea1 = *cuaMesh.pmcanis_ea1;
			cuReal3 mcanis_ea2 = *cuaMesh.pmcanis_ea2;
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pK, K, *cuaMesh.pmcanis_ea1, mcanis_ea1, *cuaMesh.pmcanis_ea2, mcanis_ea2);

			//vector product of ea1 and ea2 : the third orthogonal axis
			cuReal3 mcanis_ea3 = mcanis_ea1 ^ mcanis_ea2;

			//calculate m.ea1, m.ea2 and m.ea3 dot products
			cuBReal d1 = (M1[idx] * mcanis_ea1) / mu_s;
			cuBReal d2 = (M1[idx] * mcanis_ea2) / mu_s;
			cuBReal d3 = (M1[idx] * mcanis_ea3) / mu_s;

			//update energy density
			cuBReal non_empty_volume = M1.get_nonempty_cells() * M1.h.dim();
			if (non_empty_volume) energy_ = K * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) / non_empty_volume;
			include_in_reduction = true;
		}
	}

	reduction_avg(0, 1, &energy_, energy, points_count, include_in_reduction);
}

cuBReal Atom_Anisotropy_CubiCUDA::GetEnergyDensity(cuRect avRect)
{
	ZeroEnergy();

	if (paMeshCUDA->GetMeshType() == MESH_ATOM_CUBIC) {

		//atomistic simple-cubic mesh

		Atom_Anisotropy_CubiCUDA_Cubic_GetEnergy <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, energy, points_count, avRect);
	}

	size_t points_count_cpu = points_count.to_cpu();

	if (points_count_cpu) return energy.to_cpu() / points_count_cpu;
	else return 0.0;
}

#endif

#endif