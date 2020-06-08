#include "Atom_AnisotropyCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_ANIUNI) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "Atom_MeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"
#include "MeshDefs.h"

__global__ void Atom_Anisotropy_UniaxialCUDA_Cubic_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, cuBReal& energy, bool do_reduction)
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
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pK, K, *cuaMesh.pmcanis_ea1, mcanis_ea1);

			//calculate m.ea dot product
			cuBReal dotprod = (M1[idx] * mcanis_ea1) / mu_s;

			//update effective field with the anisotropy field
			Heff_value = (2 * K / ((cuBReal)MUB_MU0*mu_s)) * dotprod * mcanis_ea1;

			if (do_reduction) {

				//energy E = -K*dotprod*dotprod 
				//update energy density
				int non_empty_volume = M1.get_nonempty_cells() * M1.h.dim();
				if (non_empty_volume) energy_ = -K * dotprod*dotprod / non_empty_volume;
			}
		}

		Heff1[idx] += Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void Atom_Anisotropy_UniaxialCUDA::UpdateField(void)
{
	if (paMeshCUDA->GetMeshType() == MESH_ATOM_CUBIC) {

		//atomistic simple-cubic mesh

		if (paMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			Atom_Anisotropy_UniaxialCUDA_Cubic_UpdateField <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (paMeshCUDA->cuaMesh, energy, true);
		}
		else {

			Atom_Anisotropy_UniaxialCUDA_Cubic_UpdateField <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (paMeshCUDA->cuaMesh, energy, false);
		}
	}
}

//-------------------Energy density methods

__global__ void Atom_Anisotropy_UniaxialCUDA_Cubic_GetEnergy(ManagedAtom_MeshCUDA& cuaMesh, cuBReal& energy, size_t& points_count, cuRect avRect)
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
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pK, K, *cuaMesh.pmcanis_ea1, mcanis_ea1);

			//calculate m.ea dot product
			cuBReal dotprod = (M1[idx] * mcanis_ea1) / mu_s;

			//energy E = -K*dotprod*dotprod 
			//update energy density
			int non_empty_volume = M1.get_nonempty_cells() * M1.h.dim();
			if (non_empty_volume) energy_ = -K * dotprod*dotprod / non_empty_volume;
			include_in_reduction = true;
		}
	}

	reduction_avg(0, 1, &energy_, energy, points_count, include_in_reduction);
}

cuBReal Atom_Anisotropy_UniaxialCUDA::GetEnergyDensity(cuRect avRect)
{
	ZeroEnergy();

	if (paMeshCUDA->GetMeshType() == MESH_ATOM_CUBIC) {

		//atomistic simple-cubic mesh

		Atom_Anisotropy_UniaxialCUDA_Cubic_GetEnergy <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, energy, points_count, avRect);
	}

	size_t points_count_cpu = points_count.to_cpu();

	if (points_count_cpu) return energy.to_cpu() / points_count_cpu;
	else return 0.0;
}

#endif

#endif