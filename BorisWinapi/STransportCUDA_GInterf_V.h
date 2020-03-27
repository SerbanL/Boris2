#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "TransportCUDA_Poisson_Spin_V.h"

#include "MeshParamsControlCUDA.h"

//Use to set cmbnd V values using continuous flux only method - A_func and B_func define the the flux across the interface as flux = A + B * dV, where dV = V_pri - V_sec
class STransportCUDA_GInterf_V_Funcs {

public:

	__host__ void construct_cu_obj(void) {}
	__host__ void destruct_cu_obj(void) {}

	__device__ cuBReal A_func(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_V_Funcs& trans_sec, TransportCUDA_Spin_V_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		cuVEC_VC<cuReal3>& S_sec = *(trans_sec.pcuMesh->pS);
		cuVEC_VC<cuReal3>& S_pri = *(trans_pri.pcuMesh->pS);

		cuVEC_VC<cuBReal>& elC_sec = *(trans_sec.pcuMesh->pelC);
		cuVEC_VC<cuBReal>& elC_pri = *(trans_pri.pcuMesh->pelC);

		cuVEC_VC<cuReal3>& M_sec = *(trans_sec.pcuMesh->pM);
		cuVEC_VC<cuReal3>& M_pri = *(trans_pri.pcuMesh->pM);

		cuBReal De_sec = *(trans_sec.pcuMesh->pDe);
		cuBReal De_pri = *(trans_pri.pcuMesh->pDe);

		cuBReal Ms_sec = *(trans_sec.pcuMesh->pMs);
		cuBReal Ms_pri = *(trans_pri.pcuMesh->pMs);

		//S values at the interface obtained using interpolation based on current values
		cuReal3 S_pri_val = 1.5 * S_pri[cell1_idx] - 0.5 * S_pri[cell2_idx];
		cuReal3 S_sec_val = 1.5 * S_sec.weighted_average(relpos_m1, stencil) - 0.5 * S_sec.weighted_average(relpos_m1 + shift, stencil);

		//conductivity values at the interface obtained using interpolation
		cuBReal elC_pri_val = 1.5 * elC_pri[cell1_idx] - 0.5 * elC_pri[cell2_idx];
		cuBReal elC_sec_val = 1.5 * elC_sec.weighted_average(relpos_m1, stencil) - 0.5 * elC_sec.weighted_average(relpos_m1 + shift, stencil);

		//Get G values from top contacting mesh
		cuReal2 Gi;
		if (primary_top) {

			Gi = *(trans_pri.pcuMesh->pGi);
			trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pGi), Gi);
		}
		else {

			Gi = *(trans_sec.pcuMesh->pGi);
			trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pGi), Gi);
		}

		if (trans_pri.stsolve == STSOLVE_FERROMAGNETIC) {

			//F is on the primary side

			trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pMs), Ms_pri, *(trans_pri.pcuMesh->pDe), De_pri);
			trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pDe), De_sec);

			//M value at the interface obtained using interpolation
			int idx1_M = M_pri.position_to_cellidx(S_pri.cellidx_to_position(cell1_idx));
			int idx2_M = M_pri.position_to_cellidx(S_pri.cellidx_to_position(cell2_idx));
			cuReal3 m = (1.5 * M_pri[idx1_M] - 0.5 * M_pri[idx2_M]) / Ms_pri;

			//dVs is the difference between the primary and secondary Vs values at the interface
			cuReal3 dVs = (De_pri * S_pri_val / elC_pri_val - De_sec * S_sec_val / elC_sec_val) / (cuBReal)MUB_E;

			return (Gi.i - Gi.j) * (dVs * m);
		}
		else {

			//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

			trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pDe), De_pri);
			trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pMs), Ms_sec, *(trans_sec.pcuMesh->pDe), De_sec);

			//M value at the interface obtained using interpolation
			cuReal3 m = (1.5 * M_sec.weighted_average(relpos_m1, stencil) - 0.5 * M_sec.weighted_average(relpos_m1 + shift, stencil)) / Ms_sec;

			//dVs is the difference between the primary and secondary Vs values at the interface
			cuReal3 dVs = (De_pri * S_pri_val / elC_pri_val - De_sec * S_sec_val / elC_sec_val) / (cuBReal)MUB_E;

			return (Gi.i - Gi.j) * (dVs * m);
		}		
	}

	__device__ cuBReal B_func(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_V_Funcs& trans_sec, TransportCUDA_Spin_V_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Get G values from top contacting mesh
		cuReal2 Gi;
		if (primary_top) {

			Gi = *(trans_pri.pcuMesh->pGi);
			trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pGi), Gi);
		}
		else {

			Gi = *(trans_sec.pcuMesh->pGi);
			trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pGi), Gi);
		}

		return -(Gi.i + Gi.j);
	}
};

#endif

#endif
