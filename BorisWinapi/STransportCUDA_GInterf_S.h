#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "MaterialParameterCUDA.h"
#include "MeshParamsControlCUDA.h"

#include "TransportCUDA_Poisson_Spin_S.h"

#include "ManagedDiffEqCUDA.h"

//Use to set cmbnd S values using discontinuous method - A_func_sec, B_func_sec and A_func_pri, B_func_pri define the fluxes either side of the interface as flux = A + B * dVs, where dVs = Vs_pri - Vs_sec
//This is for an NF interface : N secondary, F primary
class STransportCUDA_GInterf_S_NF_Funcs {

public:

	__host__ void construct_cu_obj(void) {}
	__host__ void destruct_cu_obj(void) {}
	
	//A-Function for N
	__device__ cuReal3 A_func_sec(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_S_Funcs& trans_sec, TransportCUDA_Spin_S_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Secondary is an N here

		cuVEC_VC<cuReal>& V_sec = *(trans_sec.pcuMesh->pV);
		cuVEC_VC<cuReal>& V_pri = *(trans_pri.pcuMesh->pV);

		cuVEC_VC<cuReal3>& M_pri = *(trans_pri.pcuMesh->pM);

		cuReal Ms_pri = *(trans_pri.pcuMesh->pMs);
		cuReal pump_eff_pri = *(trans_pri.pcuMesh->ppump_eff);

		//Get G values from top contacting mesh
		cuReal2 Gi, Gmix;
		if (primary_top) {

			Gi = *(trans_pri.pcuMesh->pGi);
			Gmix = *(trans_pri.pcuMesh->pGmix);
			trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pGi), Gi, *(trans_pri.pcuMesh->pGmix), Gmix);
		}
		else {

			Gi = *(trans_sec.pcuMesh->pGi);
			Gmix = *(trans_sec.pcuMesh->pGmix);
			trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pGi), Gi, *(trans_sec.pcuMesh->pGmix), Gmix);
		}

		cuReal3 m;

		//on the N side must include spin pumping current : subtract it from the A value.
		cuReal3 Js_pump = cuReal3();

		if (M_pri.linear_size()) {

			//F is on the primary side

			trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pMs), Ms_pri, *(trans_pri.pcuMesh->ppump_eff), pump_eff_pri);

			//M value at the interface obtained using interpolation
			int idx1_M = M_pri.position_to_cellidx(V_pri.cellidx_to_position(cell1_idx));
			int idx2_M = M_pri.position_to_cellidx(V_pri.cellidx_to_position(cell2_idx));
			m = (1.5 * M_pri[idx1_M] - 0.5 * M_pri[idx2_M]) / Ms_pri;
			
			if (cuIsNZ((cuReal)pump_eff_pri) && trans_pri.pcuDiffEq) {

				cuReal3 dmdt = (1.5 * (trans_pri.pcuDiffEq)->dMdt(idx1_M) - 0.5 * (trans_pri.pcuDiffEq)->dMdt(idx2_M)) / Ms_pri;

				Js_pump = pump_eff_pri * (cuReal)MUB_E * (cuReal)HBAR_E * (Gmix.i * (m ^ dmdt) + Gmix.j * m);
			}
		}

		//V values at the interface obtained using interpolation
		cuReal V_pri_val = 1.5 * V_pri[cell1_idx] - 0.5 * V_pri[cell2_idx];
		cuReal V_sec_val = 1.5 * V_sec.weighted_average(relpos_m1, stencil) - 0.5 * V_sec.weighted_average(relpos_m1 + shift, stencil);
		cuReal dV = V_pri_val - V_sec_val;

		return (cuReal)MUB_E * (Gi.i - Gi.j) * dV * m - Js_pump;
	}

	//A-Function for F
	__device__ cuReal3 A_func_pri(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_S_Funcs& trans_sec, TransportCUDA_Spin_S_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Primary is an F here

		cuVEC_VC<cuReal>& V_sec = *(trans_sec.pcuMesh->pV);
		cuVEC_VC<cuReal>& V_pri = *(trans_pri.pcuMesh->pV);
		
		cuVEC_VC<cuReal3>& M_pri = *(trans_pri.pcuMesh->pM);
		
		cuReal Ms_pri = *(trans_pri.pcuMesh->pMs);

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

		cuReal3 m;

		if (M_pri.linear_size()) {

			//F is on the primary side

			trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pMs), Ms_pri);

			//M value at the interface obtained using interpolation
			int idx1_M = M_pri.position_to_cellidx(V_pri.cellidx_to_position(cell1_idx));
			int idx2_M = M_pri.position_to_cellidx(V_pri.cellidx_to_position(cell2_idx));
			m = (1.5 * M_pri[idx1_M] - 0.5 * M_pri[idx2_M]) / Ms_pri;
		}

		//V values at the interface obtained using interpolation
		cuReal V_pri_val = 1.5 * V_pri[cell1_idx] - 0.5 * V_pri[cell2_idx];
		cuReal V_sec_val = 1.5 * V_sec.weighted_average(relpos_m1, stencil) - 0.5 * V_sec.weighted_average(relpos_m1 + shift, stencil);

		cuReal dV = V_pri_val - V_sec_val;

		return MUB_E * (Gi.i - Gi.j) * dV * m;
	}

	//B-Function for N
	__device__ cuReal33 B_func_sec(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_S_Funcs& trans_sec, TransportCUDA_Spin_S_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Secondary is an N here

		cuVEC_VC<cuReal3>& S_pri = *(trans_pri.pcuMesh->pS);

		cuVEC_VC<cuReal3>& M_pri = *(trans_pri.pcuMesh->pM);
		
		cuReal Ms_pri = *(trans_pri.pcuMesh->pMs);

		//Get G values from top contacting mesh
		cuReal2 Gi, Gmix;
		if (primary_top) {

			Gi = *(trans_pri.pcuMesh->pGi);
			Gmix = *(trans_pri.pcuMesh->pGmix);
			trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pGi), Gi, *(trans_pri.pcuMesh->pGmix), Gmix);
		}
		else {

			Gi = *(trans_sec.pcuMesh->pGi);
			Gmix = *(trans_sec.pcuMesh->pGmix);
			trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pGi), Gi, *(trans_sec.pcuMesh->pGmix), Gmix);
		}

		cuReal3 m;

		if (M_pri.linear_size()) {

			//F is on the primary side

			trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pMs), Ms_pri);

			//M value at the interface obtained using interpolation
			int idx1_M = M_pri.position_to_cellidx(S_pri.cellidx_to_position(cell1_idx));
			int idx2_M = M_pri.position_to_cellidx(S_pri.cellidx_to_position(cell2_idx));
			m = (1.5 * M_pri[idx1_M] - 0.5 * M_pri[idx2_M]) / Ms_pri;
		}

		cuReal gR = 2 * (cuReal)MUB_E * Gmix.i;
		cuReal gI = 2 * (cuReal)MUB_E * Gmix.j;

		cuReal33 eps_m = cu_epsilon3(m);

		return (-(cuReal)MUB_E * (Gi.i + Gi.j) * (m | m)) + (eps_m * (gR * eps_m + gI * cu_ident<cuReal33>()));
	}

	//B-Function for F
	__device__ cuReal33 B_func_pri(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_S_Funcs& trans_sec, TransportCUDA_Spin_S_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Primary is an F here

		cuVEC_VC<cuReal3>& S_pri = *(trans_pri.pcuMesh->pS);

		cuVEC_VC<cuReal3>& M_pri = *(trans_pri.pcuMesh->pM);

		cuReal Ms_pri = *(trans_pri.pcuMesh->pMs);

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

		cuReal3 m;

		if (M_pri.linear_size()) {

			//F is on the primary side

			trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pMs), Ms_pri);

			//M value at the interface obtained using interpolation
			int idx1_M = M_pri.position_to_cellidx(S_pri.cellidx_to_position(cell1_idx));
			int idx2_M = M_pri.position_to_cellidx(S_pri.cellidx_to_position(cell2_idx));
			m = (1.5 * M_pri[idx1_M] - 0.5 * M_pri[idx2_M]) / Ms_pri;
		}

		return -(cuReal)MUB_E * (Gi.i + Gi.j) * (m | m);
	}
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Use to set cmbnd S values using discontinuous method - A_func_sec, B_func_sec and A_func_pri, B_func_pri define the fluxes either side of the interface as flux = A + B * dVs, where dVs = Vs_pri - Vs_sec
//This is for an FN interface : F secondary, N primary
class STransportCUDA_GInterf_S_FN_Funcs {

public:

	__host__ void construct_cu_obj(void) {}
	__host__ void destruct_cu_obj(void) {}

	//A-Function for N
	__device__ cuReal3 A_func_pri(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_S_Funcs& trans_sec, TransportCUDA_Spin_S_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Primary is an N here

		cuVEC_VC<cuReal>& V_sec = *(trans_sec.pcuMesh->pV);
		cuVEC_VC<cuReal>& V_pri = *(trans_pri.pcuMesh->pV);

		cuVEC_VC<cuReal3>& M_sec = *(trans_sec.pcuMesh->pM);

		cuReal Ms_sec = *(trans_sec.pcuMesh->pMs);
		cuReal pump_eff_sec = *(trans_sec.pcuMesh->ppump_eff);

		//Get G values from top contacting mesh
		cuReal2 Gi, Gmix;
		if (primary_top) {

			Gi = *(trans_pri.pcuMesh->pGi);
			Gmix = *(trans_pri.pcuMesh->pGmix);
			trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pGi), Gi, *(trans_pri.pcuMesh->pGmix), Gmix);
		}
		else {

			Gi = *(trans_sec.pcuMesh->pGi);
			Gmix = *(trans_sec.pcuMesh->pGmix);
			trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pGi), Gi, *(trans_sec.pcuMesh->pGmix), Gmix);
		}

		cuReal3 m;

		//on the N side must include spin pumping current : subtract it from the A value.
		cuReal3 Js_pump = cuReal3();

		if (M_sec.linear_size()) {

			//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

			trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pMs), Ms_sec, *(trans_sec.pcuMesh->ppump_eff), pump_eff_sec);

			//M value at the interface obtained using interpolation
			m = (1.5 * M_sec.weighted_average(relpos_m1, stencil) - 0.5 * M_sec.weighted_average(relpos_m1 + shift, stencil)) / Ms_sec;

			if (cuIsNZ((cuReal)pump_eff_sec) && trans_sec.pcuDiffEq) {

				int idx1_M = M_sec.position_to_cellidx(relpos_m1);
				int idx2_M = M_sec.position_to_cellidx(relpos_m1 + shift);

				cuReal3 dmdt = (1.5 * (trans_sec.pcuDiffEq)->dMdt(idx1_M) - 0.5 * (trans_sec.pcuDiffEq)->dMdt(idx2_M)) / Ms_sec;

				Js_pump = pump_eff_sec * (cuReal)MUB_E * (cuReal)HBAR_E * (Gmix.i * (m ^ dmdt) + Gmix.j * m);
			}
		}

		//V values at the interface obtained using interpolation
		cuReal V_pri_val = 1.5 * V_pri[cell1_idx] - 0.5 * V_pri[cell2_idx];
		cuReal V_sec_val = 1.5 * V_sec.weighted_average(relpos_m1, stencil) - 0.5 * V_sec.weighted_average(relpos_m1 + shift, stencil);
		cuReal dV = V_pri_val - V_sec_val;

		return (cuReal)MUB_E * (Gi.i - Gi.j) * dV * m - Js_pump;
	}

	//A-Function for F
	__device__ cuReal3 A_func_sec(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_S_Funcs& trans_sec, TransportCUDA_Spin_S_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Secondary is an F here

		cuVEC_VC<cuReal>& V_sec = *(trans_sec.pcuMesh->pV);
		cuVEC_VC<cuReal>& V_pri = *(trans_pri.pcuMesh->pV);

		cuVEC_VC<cuReal3>& M_sec = *(trans_sec.pcuMesh->pM);

		cuReal Ms_sec = *(trans_sec.pcuMesh->pMs);

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

		cuReal3 m;

		if (M_sec.linear_size()) {

			//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

			trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pMs), Ms_sec);

			//M value at the interface obtained using interpolation
			m = (1.5 * M_sec.weighted_average(relpos_m1, stencil) - 0.5 * M_sec.weighted_average(relpos_m1 + shift, stencil)) / Ms_sec;
		}

		//V values at the interface obtained using interpolation
		cuReal V_pri_val = 1.5 * V_pri[cell1_idx] - 0.5 * V_pri[cell2_idx];
		cuReal V_sec_val = 1.5 * V_sec.weighted_average(relpos_m1, stencil) - 0.5 * V_sec.weighted_average(relpos_m1 + shift, stencil);

		cuReal dV = V_pri_val - V_sec_val;

		return MUB_E * (Gi.i - Gi.j) * dV * m;
	}

	//B-Function for N
	__device__ cuReal33 B_func_pri(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_S_Funcs& trans_sec, TransportCUDA_Spin_S_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Primary is an N here

		cuVEC_VC<cuReal3>& M_sec = *(trans_sec.pcuMesh->pM);

		cuReal Ms_sec = *(trans_sec.pcuMesh->pMs);

		//Get G values from top contacting mesh
		cuReal2 Gi, Gmix;
		if (primary_top) {

			Gi = *(trans_pri.pcuMesh->pGi);
			Gmix = *(trans_pri.pcuMesh->pGmix);
			trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pGi), Gi, *(trans_pri.pcuMesh->pGmix), Gmix);
		}
		else {

			Gi = *(trans_sec.pcuMesh->pGi);
			Gmix = *(trans_sec.pcuMesh->pGmix);
			trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pGi), Gi, *(trans_sec.pcuMesh->pGmix), Gmix);
		}

		cuReal3 m;

		if (M_sec.linear_size()) {

			//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

			trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pMs), Ms_sec);

			//M value at the interface obtained using interpolation
			m = (1.5 * M_sec.weighted_average(relpos_m1, stencil) - 0.5 * M_sec.weighted_average(relpos_m1 + shift, stencil)) / Ms_sec;
		}

		cuReal gR = 2 * (cuReal)MUB_E * Gmix.i;
		cuReal gI = 2 * (cuReal)MUB_E * Gmix.j;

		cuReal33 eps_m = cu_epsilon3(m);

		return (-(cuReal)MUB_E * (Gi.i + Gi.j) * (m | m)) + (eps_m * (gR * eps_m + gI * cu_ident<cuReal33>()));
	}

	//B-Function for F
	__device__ cuReal33 B_func_sec(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_S_Funcs& trans_sec, TransportCUDA_Spin_S_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Secondary is an F here

		cuVEC_VC<cuReal3>& M_sec = *(trans_sec.pcuMesh->pM);

		cuReal Ms_sec = *(trans_sec.pcuMesh->pMs);

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

		cuReal3 m;

		if (M_sec.linear_size()) {

			//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

			trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pMs), Ms_sec);

			//M value at the interface obtained using interpolation
			m = (1.5 * M_sec.weighted_average(relpos_m1, stencil) - 0.5 * M_sec.weighted_average(relpos_m1 + shift, stencil)) / Ms_sec;
		}

		return -(cuReal)MUB_E * (Gi.i + Gi.j) * (m | m);
	}
};

#endif

#endif