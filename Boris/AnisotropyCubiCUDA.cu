#include "AnisotropyCubiCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_ANICUBI

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"
#include "MeshDefs.h"

__global__ void Anisotropy_CubicCUDA_FM_UpdateField(ManagedMeshCUDA& cuMesh, cuBReal& energy, bool do_reduction)
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
			cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;

			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pK1, K1, *cuMesh.pK2, K2, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2);

			//vector product of ea1 and ea2 : the third orthogonal axis
			cuReal3 mcanis_ea3 = mcanis_ea1 ^ mcanis_ea2;

			//calculate m.ea1, m.ea2 and m.ea3 dot products
			cuBReal d1 = (M[idx] * mcanis_ea1) / Ms;
			cuBReal d2 = (M[idx] * mcanis_ea2) / Ms;
			cuBReal d3 = (M[idx] * mcanis_ea3) / Ms;

			//terms for K1 contribution
			cuBReal a1 = d1 * (d2*d2 + d3 * d3);
			cuBReal a2 = d2 * (d1*d1 + d3 * d3);
			cuBReal a3 = d3 * (d1*d1 + d2 * d2);

			//terms for K2 contribution
			cuBReal d123 = d1 * d2*d3;

			cuBReal b1 = d123 * d2*d3;
			cuBReal b2 = d123 * d1*d3;
			cuBReal b3 = d123 * d1*d2;

			//update effective field with the anisotropy field
			Heff_value = cuReal3(
				(-2 * K1 / ((cuBReal)MU0*Ms)) * (mcanis_ea1.i * a1 + mcanis_ea2.i * a2 + mcanis_ea3.i * a3)
				+ (-2 * K2 / ((cuBReal)MU0*Ms)) * (mcanis_ea1.i * b1 + mcanis_ea2.i * b2 + mcanis_ea3.i * b3),

				(-2 * K1 / ((cuBReal)MU0*Ms)) * (mcanis_ea1.j * a1 + mcanis_ea2.j * a2 + mcanis_ea3.j * a3)
				+ (-2 * K2 / ((cuBReal)MU0*Ms)) * (mcanis_ea1.j * b1 + mcanis_ea2.j * b2 + mcanis_ea3.j * b3),

				(-2 * K1 / ((cuBReal)MU0*Ms)) * (mcanis_ea1.k * a1 + mcanis_ea2.k * a2 + mcanis_ea3.k * a3)
				+ (-2 * K2 / ((cuBReal)MU0*Ms)) * (mcanis_ea1.k * b1 + mcanis_ea2.k * b2 + mcanis_ea3.k * b3)
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

__global__ void Anisotropy_CubicCUDA_AFM_UpdateField(ManagedMeshCUDA& cuMesh, cuBReal& energy, bool do_reduction)
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
			cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;

			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pK1_AFM, K1_AFM, *cuMesh.pK2_AFM, K2_AFM, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2);

			//vector product of ea1 and ea2 : the third orthogonal axis
			cuReal3 mcanis_ea3 = mcanis_ea1 ^ mcanis_ea2;

			//calculate m.ea1, m.ea2 and m.ea3 dot products
			cuBReal d1 = (M[idx] * mcanis_ea1) / Ms_AFM.i;
			cuBReal d2 = (M[idx] * mcanis_ea2) / Ms_AFM.i;
			cuBReal d3 = (M[idx] * mcanis_ea3) / Ms_AFM.i;

			//terms for K1 contribution
			cuBReal a1 = d1 * (d2*d2 + d3 * d3);
			cuBReal a2 = d2 * (d1*d1 + d3 * d3);
			cuBReal a3 = d3 * (d1*d1 + d2 * d2);

			//terms for K2 contribution
			cuBReal d123 = d1*d2*d3;

			cuBReal b1 = d123 * d2*d3;
			cuBReal b2 = d123 * d1*d3;
			cuBReal b3 = d123 * d1*d2;

			//update effective field with the anisotropy field
			Heff_value = cuReal3(
				(-2 * K1_AFM.i / ((cuBReal)MU0*Ms_AFM.i)) * (mcanis_ea1.i * a1 + mcanis_ea2.i * a2 + mcanis_ea3.i * a3)
				+ (-2 * K2_AFM.i / ((cuBReal)MU0*Ms_AFM.i)) * (mcanis_ea1.i * b1 + mcanis_ea2.i * b2 + mcanis_ea3.i * b3),

				(-2 * K1_AFM.i / ((cuBReal)MU0*Ms_AFM.i)) * (mcanis_ea1.j * a1 + mcanis_ea2.j * a2 + mcanis_ea3.j * a3)
				+ (-2 * K2_AFM.i / ((cuBReal)MU0*Ms_AFM.i)) * (mcanis_ea1.j * b1 + mcanis_ea2.j * b2 + mcanis_ea3.j * b3),

				(-2 * K1_AFM.i / ((cuBReal)MU0*Ms_AFM.i)) * (mcanis_ea1.k * a1 + mcanis_ea2.k * a2 + mcanis_ea3.k * a3)
				+ (-2 * K2_AFM.i / ((cuBReal)MU0*Ms_AFM.i)) * (mcanis_ea1.k * b1 + mcanis_ea2.k * b2 + mcanis_ea3.k * b3)
			);

			//same thing for sub-lattice B

			cuBReal d1B = (M2[idx] * mcanis_ea1) / Ms_AFM.j;
			cuBReal d2B = (M2[idx] * mcanis_ea2) / Ms_AFM.j;
			cuBReal d3B = (M2[idx] * mcanis_ea3) / Ms_AFM.j;

			cuBReal a1B = d1B * (d2B*d2B + d3B*d3B);
			cuBReal a2B = d2B * (d1B*d1B + d3B*d3B);
			cuBReal a3B = d3B * (d1B*d1B + d2B*d2B);

			cuBReal d123B = d1B*d2B*d3B;

			cuBReal b1B = d123B*d2B*d3B;
			cuBReal b2B = d123B*d1B*d3B;
			cuBReal b3B = d123B*d1B*d2B;

			Heff2_value = cuReal3(
				(-2 * K1_AFM.j / ((cuBReal)MU0*Ms_AFM.j)) * (mcanis_ea1.i * a1B + mcanis_ea2.i * a2B + mcanis_ea3.i * a3B)
				+ (-2 * K2_AFM.j / ((cuBReal)MU0*Ms_AFM.j)) * (mcanis_ea1.i * b1B + mcanis_ea2.i * b2B + mcanis_ea3.i * b3B),

				(-2 * K1_AFM.j / ((cuBReal)MU0*Ms_AFM.j)) * (mcanis_ea1.j * a1B + mcanis_ea2.j * a2B + mcanis_ea3.j * a3B)
				+ (-2 * K2_AFM.j / ((cuBReal)MU0*Ms_AFM.j)) * (mcanis_ea1.j * b1B + mcanis_ea2.j * b2B + mcanis_ea3.j * b3B),

				(-2 * K1_AFM.j / ((cuBReal)MU0*Ms_AFM.j)) * (mcanis_ea1.k * a1B + mcanis_ea2.k * a2B + mcanis_ea3.k * a3B)
				+ (-2 * K2_AFM.j / ((cuBReal)MU0*Ms_AFM.j)) * (mcanis_ea1.k * b1B + mcanis_ea2.k * b2B + mcanis_ea3.k * b3B)
			);

			if (do_reduction) {

				//update energy (E/V)		
				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ = ((K1_AFM.i * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2_AFM.i * d123*d123) + (K1_AFM.j * (d1B*d1B*d2B*d2B + d1B*d1B*d3B*d3B + d2B*d2B*d3B*d3B) + K2_AFM.j * d123B*d123B)) / (2*non_empty_cells);
			}
		}

		Heff[idx] += Heff_value;
		Heff2[idx] += Heff2_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void Anisotropy_CubicCUDA::UpdateField(void)
{
	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		//anti-ferromagnetic mesh

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			Anisotropy_CubicCUDA_AFM_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, true);
		}
		else {

			Anisotropy_CubicCUDA_AFM_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, false);
		}
	}
	else {

		//ferromagnetic mesh

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			Anisotropy_CubicCUDA_FM_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, true);
		}
		else {

			Anisotropy_CubicCUDA_FM_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, false);
		}
	}
}

//-------------------Energy density methods

__global__ void Anisotropy_CubicCUDA_FM_GetEnergy(ManagedMeshCUDA& cuMesh, cuBReal& energy, size_t& points_count, cuRect avRect)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && avRect.contains(M.cellidx_to_position(idx))) {

			cuBReal Ms = *cuMesh.pMs;
			cuBReal K1 = *cuMesh.pK1;
			cuBReal K2 = *cuMesh.pK2;
			cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
			cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;

			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pK1, K1, *cuMesh.pK2, K2, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2);

			//vector product of ea1 and ea2 : the third orthogonal axis
			cuReal3 mcanis_ea3 = mcanis_ea1 ^ mcanis_ea2;

			//calculate m.ea1, m.ea2 and m.ea3 dot products
			cuBReal d1 = (M[idx] * mcanis_ea1) / Ms;
			cuBReal d2 = (M[idx] * mcanis_ea2) / Ms;
			cuBReal d3 = (M[idx] * mcanis_ea3) / Ms;

			//terms for K2 contribution
			cuBReal d123 = d1*d2*d3;

			energy_ = (K1 * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2 * d123*d123);
			include_in_reduction = true;
		}
	}

	reduction_avg(0, 1, &energy_, energy, points_count, include_in_reduction);
}

__global__ void Anisotropy_CubicCUDA_AFM_GetEnergy(ManagedMeshCUDA& cuMesh, cuBReal& energy, size_t& points_count, cuRect avRect)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx) && avRect.contains(M.cellidx_to_position(idx))) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuReal2 K1_AFM = *cuMesh.pK1_AFM;
			cuReal2 K2_AFM = *cuMesh.pK2_AFM;
			cuReal3 mcanis_ea1 = *cuMesh.pmcanis_ea1;
			cuReal3 mcanis_ea2 = *cuMesh.pmcanis_ea2;

			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pK1_AFM, K1_AFM, *cuMesh.pK2_AFM, K2_AFM, *cuMesh.pmcanis_ea1, mcanis_ea1, *cuMesh.pmcanis_ea2, mcanis_ea2);

			//vector product of ea1 and ea2 : the third orthogonal axis
			cuReal3 mcanis_ea3 = mcanis_ea1 ^ mcanis_ea2;

			//calculate m.ea1, m.ea2 and m.ea3 dot products
			cuBReal d1 = (M[idx] * mcanis_ea1) / Ms_AFM.i;
			cuBReal d2 = (M[idx] * mcanis_ea2) / Ms_AFM.i;
			cuBReal d3 = (M[idx] * mcanis_ea3) / Ms_AFM.i;

			//terms for K2 contribution
			cuBReal d123 = d1*d2*d3;

			//same thing for sub-lattice B

			cuBReal d1B = (M2[idx] * mcanis_ea1) / Ms_AFM.j;
			cuBReal d2B = (M2[idx] * mcanis_ea2) / Ms_AFM.j;
			cuBReal d3B = (M2[idx] * mcanis_ea3) / Ms_AFM.j;

			cuBReal d123B = d1B*d2B*d3B;

			energy_ = ((K1_AFM.i * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2_AFM.i * d123*d123) + (K1_AFM.j * (d1B*d1B*d2B*d2B + d1B * d1B*d3B*d3B + d2B * d2B*d3B*d3B) + K2_AFM.j * d123B*d123B)) / 2;
			include_in_reduction = true;
		}
	}

	reduction_avg(0, 1, &energy_, energy, points_count, include_in_reduction);
}

cuBReal Anisotropy_CubicCUDA::GetEnergyDensity(cuRect avRect)
{
	ZeroEnergy();

	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		//anti-ferromagnetic mesh

		Anisotropy_CubicCUDA_AFM_GetEnergy <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, energy, points_count, avRect);
	}
	else {

		//ferromagnetic mesh

		Anisotropy_CubicCUDA_FM_GetEnergy <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, energy, points_count, avRect);
	}

	size_t points_count_cpu = points_count.to_cpu();

	if (points_count_cpu) return energy.to_cpu() / points_count_cpu;
	else return 0.0;
}

#endif

#endif