#include "Atom_ExchangeCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_EXCHANGE) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "Atom_MeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"
#include "MeshDefs.h"

//////////////////////////////////////////////////////////////////////// UPDATE FIELD

__global__ void Atom_ExchangeCUDA_Cubic_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, cuBReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff1.linear_size()) {

		cuReal3 Hexch = cuReal3();

		if (M1.is_not_empty(idx)) {

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal J = *cuaMesh.pJ;
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pJ, J);

			//update effective field with the Heisenberg exchange field
			Hexch = (J / (MUB_MU0*mu_s)) * M1.ngbr_dirsum(idx);

			if (do_reduction) {

				//energy E = -mu_s * Bex
				//update energy density
				cuBReal non_empty_volume = M1.get_nonempty_cells() * M1.h.dim();
				if (non_empty_volume) energy_ = -(cuBReal)MUB_MU0 * M1[idx] * Hexch / (2*non_empty_volume);
			}
		}

		Heff1[idx] += Hexch;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void Atom_ExchangeCUDA::UpdateField(void)
{
	if (paMeshCUDA->GetMeshType() == MESH_ATOM_CUBIC) {

		//atomistic simple cubic mesh

		if (paMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			Atom_ExchangeCUDA_Cubic_UpdateField <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, energy, true);
		}
		else {

			Atom_ExchangeCUDA_Cubic_UpdateField <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, energy, false);
		}
	}
}

//////////////////////////////////////////////////////////////////////// ENERGY DENSITY DATA METHODS

__global__ void Atom_ExchangeCUDA_Cubic_GetEnergy(ManagedAtom_MeshCUDA& cuaMesh, cuBReal& energy, size_t& points_count, cuRect avRect)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M1.linear_size()) {

		cuReal3 Hexch = cuReal3();

		if (M1.is_not_empty(idx) && avRect.contains(M1.cellidx_to_position(idx))) {

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal J = *cuaMesh.pJ;
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pJ, J);

			//update effective field with the Heisenberg exchange field
			Hexch = (J / (MUB_MU0*mu_s)) * M1.ngbr_dirsum(idx);

			//energy at this point
			energy_ = -(cuBReal)MUB_MU0 * M1[idx] * Hexch / 2;
			include_in_reduction = true;
		}
	}

	reduction_avg(0, 1, &energy_, energy, points_count, include_in_reduction);
}

__global__ void Atom_ExchangeCUDA_Cubic_GetEnergy_Max(ManagedAtom_MeshCUDA& cuaMesh, cuBReal& energy, cuRect rectangle)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M1.linear_size()) {

		cuReal3 Hexch = cuReal3();

		if (M1.is_not_empty(idx) && rectangle.contains(M1.cellidx_to_position(idx))) {

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal J = *cuaMesh.pJ;
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pJ, J);

			//update effective field with the Heisenberg exchange field
			Hexch = (J / (MUB_MU0*mu_s)) * M1.ngbr_dirsum(idx);

			//energy modulus at this point
			energy_ = fabs((cuBReal)MUB_MU0 * M1[idx] * Hexch / 2);
			include_in_reduction = true;
		}
	}

	reduction_max(0, 1, &energy_, energy, include_in_reduction);
}

cuBReal Atom_ExchangeCUDA::GetEnergyDensity(cuRect avRect)
{
	ZeroEnergy();

	if (paMeshCUDA->GetMeshType() == MESH_ATOM_CUBIC) {

		//atomistic simple-cubic mesh

		Atom_ExchangeCUDA_Cubic_GetEnergy <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, energy, points_count, avRect);
	}

	size_t points_count_cpu = points_count.to_cpu();

	//return energy density
	if (points_count_cpu) return energy.to_cpu() / (points_count_cpu * paMeshCUDA->h.dim());
	else return 0.0;
}

cuBReal Atom_ExchangeCUDA::GetEnergy_Max(cuRect rectangle)
{
	ZeroEnergy();

	if (paMeshCUDA->GetMeshType() == MESH_ATOM_CUBIC) {

		//atomistic simple-cubic mesh

		Atom_ExchangeCUDA_Cubic_GetEnergy_Max <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, energy, rectangle);
	}

	//return maximum energy density modulus
	return energy.to_cpu() / paMeshCUDA->h.dim();
}

//////////////////////////////////////////////////////////////////////// ENERGY DENSITY DISPLAY METHODS

__global__ void Atom_ExchangeCUDA_Cubic_Compute_Exchange(ManagedAtom_MeshCUDA& cuaMesh, cuVEC<cuBReal>& exchange_displayVEC)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < M1.linear_size()) {

		cuReal3 Hexch = cuReal3();

		if (M1.is_not_empty(idx)) {

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal J = *cuaMesh.pJ;
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pJ, J);

			//update effective field with the Heisenberg exchange field
			Hexch = (J / (MUB_MU0*mu_s)) * M1.ngbr_dirsum(idx);
		}

		//energy density at this point
		exchange_displayVEC[idx] = -(cuBReal)MUB_MU0 * M1[idx] * Hexch / (2 * M1.h.dim());
	}
}

void Atom_ExchangeCUDA::Compute_ExchangeCUDA(void)
{
	exchange_displayVEC()->resize(paMeshCUDA->h, paMeshCUDA->meshRect);

	if (paMeshCUDA->GetMeshType() == MESH_ATOM_CUBIC) {

		//atomistic simple-cubic mesh

		Atom_ExchangeCUDA_Cubic_Compute_Exchange <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (paMeshCUDA->cuaMesh, exchange_displayVEC);
	}
}

#endif

#endif