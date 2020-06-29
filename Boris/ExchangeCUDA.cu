#include "ExchangeCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_EXCHANGE

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"
#include "MeshDefs.h"

//////////////////////////////////////////////////////////////////////// UPDATE FIELD

__global__ void ExchangeCUDA_FM_UpdateField(ManagedMeshCUDA& cuMesh, cuBReal& energy, bool do_reduction)
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

__global__ void ExchangeCUDA_AFM_UpdateField(ManagedMeshCUDA& cuMesh, cuBReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuReal3 Hexch = cuReal3();
		cuReal3 Hexch2 = cuReal3();

		if (M.is_not_empty(idx)) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuReal2 A_AFM = *cuMesh.pA_AFM;
			cuReal2 Ah = *cuMesh.pAh;
			cuReal2 Anh = *cuMesh.pAnh;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pA_AFM, A_AFM, *cuMesh.pAh, Ah, *cuMesh.pAnh, Anh);

			cuReal3 delsq_M_A = M.delsq_neu(idx);
			cuReal3 delsq_M_B = M2.delsq_neu(idx);

			cuReal2 Mmag = cuReal2(M[idx].norm(), M2[idx].norm());

			Hexch = (2 * A_AFM.i / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.i)) * delsq_M_A + (-4 * Ah.i * (M[idx] ^ (M[idx] ^ M2[idx])) / (Mmag.i*Mmag.i) + Anh.i * delsq_M_B) / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.j);
			Hexch2 = (2 * A_AFM.j / ((cuBReal)MU0*Ms_AFM.j*Ms_AFM.j)) * delsq_M_B + (-4 * Ah.j * (M2[idx] ^ (M2[idx] ^ M[idx])) / (Mmag.j*Mmag.j) + Anh.j * delsq_M_A) / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.j);

			if (do_reduction) {

				int non_empty_cells = M.get_nonempty_cells();
				if (non_empty_cells) energy_ = -(cuBReal)MU0 * (M[idx] * Hexch  + M2[idx] * Hexch2) / (4 * non_empty_cells);
			}
		}

		Heff[idx] += Hexch;
		Heff2[idx] += Hexch2;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void Exch_6ngbr_NeuCUDA::UpdateField(void)
{
	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		//anti-ferromagnetic mesh

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			ExchangeCUDA_AFM_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, true);
		}
		else {

			ExchangeCUDA_AFM_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, false);
		}
	}
	else {

		//ferromagnetic mesh

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			ExchangeCUDA_FM_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, true);
		}
		else {

			ExchangeCUDA_FM_UpdateField << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (pMeshCUDA->cuMesh, energy, false);
		}
	}
	
	if (pMeshCUDA->GetMeshExchangeCoupling()) CalculateExchangeCoupling(energy);
}

//////////////////////////////////////////////////////////////////////// ENERGY DENSITY DATA METHODS

__global__ void ExchangeCUDA_FM_GetEnergy(ManagedMeshCUDA& cuMesh, cuBReal& energy, size_t& points_count, cuRect avRect)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M.linear_size()) {

		cuReal3 Hexch = cuReal3();

		if (M.is_not_empty(idx) && avRect.contains(M.cellidx_to_position(idx))) {

			cuBReal Ms = *cuMesh.pMs;
			cuBReal A = *cuMesh.pA;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pA, A);

			Hexch = 2 * A * M.delsq_neu(idx) / ((cuBReal)MU0 * Ms * Ms);

			energy_ = -(cuBReal)MU0 * M[idx] * Hexch / 2;
			include_in_reduction = true;
		}
	}

	reduction_avg(0, 1, &energy_, energy, points_count, include_in_reduction);
}

__global__ void ExchangeCUDA_AFM_GetEnergy(ManagedMeshCUDA& cuMesh, cuBReal& energy, size_t& points_count, cuRect avRect)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M.linear_size()) {

		cuReal3 Hexch = cuReal3();
		cuReal3 Hexch2 = cuReal3();

		if (M.is_not_empty(idx) && avRect.contains(M.cellidx_to_position(idx))) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuReal2 A_AFM = *cuMesh.pA_AFM;
			cuReal2 Ah = *cuMesh.pAh;
			cuReal2 Anh = *cuMesh.pAnh;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pA_AFM, A_AFM, *cuMesh.pAh, Ah, *cuMesh.pAnh, Anh);

			cuReal3 delsq_M_A = M.delsq_neu(idx);
			cuReal3 delsq_M_B = M2.delsq_neu(idx);

			cuReal2 Mmag = cuReal2(M[idx].norm(), M2[idx].norm());

			Hexch = (2 * A_AFM.i / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.i)) * delsq_M_A + (-4 * Ah.i * (M[idx] ^ (M[idx] ^ M2[idx])) / (Mmag.i*Mmag.i) + Anh.i * delsq_M_B) / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.j);
			Hexch2 = (2 * A_AFM.j / ((cuBReal)MU0*Ms_AFM.j*Ms_AFM.j)) * delsq_M_B + (-4 * Ah.j * (M2[idx] ^ (M2[idx] ^ M[idx])) / (Mmag.j*Mmag.j) + Anh.j * delsq_M_A) / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.j);

			energy_ = -(cuBReal)MU0 * (M[idx] * Hexch + M2[idx] * Hexch2) / 4;
			include_in_reduction = true;
		}
	}

	reduction_avg(0, 1, &energy_, energy, points_count, include_in_reduction);
}

__global__ void ExchangeCUDA_FM_GetEnergy_Max(ManagedMeshCUDA& cuMesh, cuBReal& energy, cuRect rectangle)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M.linear_size()) {

		cuReal3 Hexch = cuReal3();

		if (M.is_not_empty(idx) && rectangle.contains(M.cellidx_to_position(idx))) {

			cuBReal Ms = *cuMesh.pMs;
			cuBReal A = *cuMesh.pA;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pA, A);

			Hexch = 2 * A * M.delsq_neu(idx) / ((cuBReal)MU0 * Ms * Ms);

			energy_ = fabs((cuBReal)MU0 * M[idx] * Hexch / 2);
			include_in_reduction = true;
		}
	}

	reduction_max(0, 1, &energy_, energy, include_in_reduction);
}

__global__ void ExchangeCUDA_AFM_GetEnergy_Max(ManagedMeshCUDA& cuMesh, cuBReal& energy, cuRect rectangle)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	bool include_in_reduction = false;

	if (idx < M.linear_size()) {

		cuReal3 Hexch = cuReal3();
		cuReal3 Hexch2 = cuReal3();

		if (M.is_not_empty(idx) && rectangle.contains(M.cellidx_to_position(idx))) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuReal2 A_AFM = *cuMesh.pA_AFM;
			cuReal2 Ah = *cuMesh.pAh;
			cuReal2 Anh = *cuMesh.pAnh;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pA_AFM, A_AFM, *cuMesh.pAh, Ah, *cuMesh.pAnh, Anh);

			cuReal3 delsq_M_A = M.delsq_neu(idx);
			cuReal3 delsq_M_B = M2.delsq_neu(idx);

			cuReal2 Mmag = cuReal2(M[idx].norm(), M2[idx].norm());

			Hexch = (2 * A_AFM.i / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.i)) * delsq_M_A + (-4 * Ah.i * (M[idx] ^ (M[idx] ^ M2[idx])) / (Mmag.i*Mmag.i) + Anh.i * delsq_M_B) / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.j);
			Hexch2 = (2 * A_AFM.j / ((cuBReal)MU0*Ms_AFM.j*Ms_AFM.j)) * delsq_M_B + (-4 * Ah.j * (M2[idx] ^ (M2[idx] ^ M[idx])) / (Mmag.j*Mmag.j) + Anh.j * delsq_M_A) / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.j);

			energy_ = fabs((cuBReal)MU0 * (M[idx] * Hexch + M2[idx] * Hexch2) / 4);
			include_in_reduction = true;
		}
	}

	reduction_max(0, 1, &energy_, energy, include_in_reduction);
}

cuBReal Exch_6ngbr_NeuCUDA::GetEnergyDensity(cuRect avRect)
{
	ZeroEnergy();

	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		//anti-ferromagnetic mesh

		ExchangeCUDA_AFM_GetEnergy <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, energy, points_count, avRect);
	}
	else {

		//ferromagnetic mesh

		ExchangeCUDA_FM_GetEnergy <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, energy, points_count, avRect);
	}

	size_t points_count_cpu = points_count.to_cpu();

	if (points_count_cpu) return energy.to_cpu() / points_count_cpu;
	else return 0.0;
}

cuBReal Exch_6ngbr_NeuCUDA::GetEnergy_Max(cuRect rectangle)
{
	ZeroEnergy();

	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		//anti-ferromagnetic mesh

		ExchangeCUDA_AFM_GetEnergy_Max <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, energy, rectangle);
	}
	else {

		//ferromagnetic mesh

		ExchangeCUDA_FM_GetEnergy_Max <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, energy, rectangle);
	}

	return energy.to_cpu();
}

//////////////////////////////////////////////////////////////////////// ENERGY DENSITY DISPLAY METHODS

__global__ void ExchangeCUDA_FM_Compute_Exchange(ManagedMeshCUDA& cuMesh, cuVEC<cuBReal>& exchange_displayVEC)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < M.linear_size()) {

		cuReal3 Hexch = cuReal3();

		if (M.is_not_empty(idx)) {

			cuBReal Ms = *cuMesh.pMs;
			cuBReal A = *cuMesh.pA;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs, Ms, *cuMesh.pA, A);

			Hexch = 2 * A * M.delsq_neu(idx) / ((cuBReal)MU0 * Ms * Ms);
		}

		exchange_displayVEC[idx] = -(cuBReal)MU0 * (M[idx] * Hexch) / 2;
	}
}

__global__ void ExchangeCUDA_AFM_Compute_Exchange(ManagedMeshCUDA& cuMesh, cuVEC<cuBReal>& exchange_displayVEC)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < M.linear_size()) {

		cuReal3 Hexch = cuReal3();
		cuReal3 Hexch2 = cuReal3();

		if (M.is_not_empty(idx)) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuReal2 A_AFM = *cuMesh.pA_AFM;
			cuReal2 Ah = *cuMesh.pAh;
			cuReal2 Anh = *cuMesh.pAnh;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pA_AFM, A_AFM, *cuMesh.pAh, Ah, *cuMesh.pAnh, Anh);

			cuReal3 delsq_M_A = M.delsq_neu(idx);
			cuReal3 delsq_M_B = M2.delsq_neu(idx);

			cuReal2 Mmag = cuReal2(M[idx].norm(), M2[idx].norm());

			Hexch = (2 * A_AFM.i / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.i)) * delsq_M_A + (-4 * Ah.i * (M[idx] ^ (M[idx] ^ M2[idx])) / (Mmag.i*Mmag.i) + Anh.i * delsq_M_B) / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.j);
			Hexch2 = (2 * A_AFM.j / ((cuBReal)MU0*Ms_AFM.j*Ms_AFM.j)) * delsq_M_B + (-4 * Ah.j * (M2[idx] ^ (M2[idx] ^ M[idx])) / (Mmag.j*Mmag.j) + Anh.j * delsq_M_A) / ((cuBReal)MU0*Ms_AFM.i*Ms_AFM.j);
		}

		exchange_displayVEC[idx] = -(cuBReal)MU0 * (M[idx] * Hexch + M2[idx] * Hexch2) / 4;
	}
}

void Exch_6ngbr_NeuCUDA::Compute_ExchangeCUDA(void)
{
	exchange_displayVEC()->resize(pMeshCUDA->h, pMeshCUDA->meshRect);

	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		//anti-ferromagnetic mesh
		
		ExchangeCUDA_AFM_Compute_Exchange <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, exchange_displayVEC);
	}
	else {

		//ferromagnetic mesh

		ExchangeCUDA_FM_Compute_Exchange <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, exchange_displayVEC);
	}
}

#endif

#endif