#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "TransportCUDA_Poisson_Spin_S.h"

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

		cuVEC_VC<cuBReal>& V_sec = (trans_sec.pcuMesh ? *(trans_sec.pcuMesh->pV) : *(trans_sec.pcuaMesh->pV));
		cuVEC_VC<cuBReal>& V_pri = (trans_pri.pcuMesh ? *(trans_pri.pcuMesh->pV) : *(trans_pri.pcuaMesh->pV));

		cuVEC_VC<cuReal3>& M_pri = (trans_pri.pcuMesh ? *(trans_pri.pcuMesh->pM) : *(trans_pri.pcuaMesh->pM1));

		cuBReal pump_eff_pri = (trans_pri.pcuMesh ? *(trans_pri.pcuMesh->ppump_eff) : *(trans_pri.pcuaMesh->ppump_eff));

		//Get G values from top contacting mesh
		cuReal2 Gi, Gmix;
		if (primary_top) {

			if (trans_pri.pcuMesh) {

				Gi = *(trans_pri.pcuMesh->pGi);
				Gmix = *(trans_pri.pcuMesh->pGmix);
				trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pGi), Gi, *(trans_pri.pcuMesh->pGmix), Gmix);
			}
			else {

				Gi = *(trans_pri.pcuaMesh->pGi);
				Gmix = *(trans_pri.pcuaMesh->pGmix);
				trans_pri.pcuaMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuaMesh->pGi), Gi, *(trans_pri.pcuaMesh->pGmix), Gmix);
			}
		}
		else {

			if (trans_sec.pcuMesh) {

				Gi = *(trans_sec.pcuMesh->pGi);
				Gmix = *(trans_sec.pcuMesh->pGmix);
				trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pGi), Gi, *(trans_sec.pcuMesh->pGmix), Gmix);
			}
			else {

				Gi = *(trans_sec.pcuaMesh->pGi);
				Gmix = *(trans_sec.pcuaMesh->pGmix);
				trans_sec.pcuaMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuaMesh->pGi), Gi, *(trans_sec.pcuaMesh->pGmix), Gmix);
			}
		}

		cuReal3 m;

		//on the N side must include spin pumping current : subtract it from the A value.
		cuReal3 Js_pump = cuReal3();

		if (trans_pri.stsolve == STSOLVE_FERROMAGNETIC || trans_pri.stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

			//F is on the primary side

			if (trans_pri.pcuMesh) trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->ppump_eff), pump_eff_pri);
			else trans_pri.pcuaMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuaMesh->ppump_eff), pump_eff_pri);

			//M value at the interface obtained using interpolation
			int idx1_M = M_pri.position_to_cellidx(V_pri.cellidx_to_position(cell1_idx));
			int idx2_M = M_pri.position_to_cellidx(V_pri.cellidx_to_position(cell2_idx));
			m = 1.5 * cu_normalize(M_pri[idx1_M]) - 0.5 * cu_normalize(M_pri[idx2_M]);
			
			if (cuIsNZ((cuBReal)pump_eff_pri) && trans_pri.pdM_dt) {

				cuReal3 dmdt = 1.5 * cu_normalize((*(trans_pri.pdM_dt))[idx1_M], M_pri[idx1_M]) - 0.5 * cu_normalize((*(trans_pri.pdM_dt))[idx2_M], M_pri[idx2_M]);

				Js_pump = pump_eff_pri * (cuBReal)MUB_E * (cuBReal)HBAR_E * (Gmix.i * (m ^ dmdt) + Gmix.j * m);
			}
		}

		//V values at the interface obtained using interpolation
		cuBReal V_pri_val = 1.5 * V_pri[cell1_idx] - 0.5 * V_pri[cell2_idx];
		cuBReal V_sec_val = 1.5 * V_sec.weighted_average(relpos_m1, stencil) - 0.5 * V_sec.weighted_average(relpos_m1 + shift, stencil);
		cuBReal dV = V_pri_val - V_sec_val;

		return (cuBReal)MUB_E * (Gi.i - Gi.j) * dV * m - Js_pump;
	}

	//A-Function for F
	__device__ cuReal3 A_func_pri(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_S_Funcs& trans_sec, TransportCUDA_Spin_S_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Primary is an F here

		cuVEC_VC<cuBReal>& V_sec = (trans_sec.pcuMesh ? *(trans_sec.pcuMesh->pV) : *(trans_sec.pcuaMesh->pV));
		cuVEC_VC<cuBReal>& V_pri = (trans_pri.pcuMesh ? *(trans_pri.pcuMesh->pV) : *(trans_pri.pcuaMesh->pV));

		cuVEC_VC<cuReal3>& M_pri = (trans_pri.pcuMesh ? *(trans_pri.pcuMesh->pM) : *(trans_pri.pcuaMesh->pM1));

		//Get G values from top contacting mesh
		cuReal2 Gi;
		if (primary_top) {

			if (trans_pri.pcuMesh) {

				Gi = *(trans_pri.pcuMesh->pGi);
				trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pGi), Gi);
			}
			else {

				Gi = *(trans_pri.pcuaMesh->pGi);
				trans_pri.pcuaMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuaMesh->pGi), Gi);
			}
		}
		else {

			if (trans_sec.pcuMesh) {

				Gi = *(trans_sec.pcuMesh->pGi);
				trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pGi), Gi);
			}
			else {

				Gi = *(trans_sec.pcuaMesh->pGi);
				trans_sec.pcuaMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuaMesh->pGi), Gi);
			}
		}

		cuReal3 m;

		if (trans_pri.stsolve == STSOLVE_FERROMAGNETIC || trans_pri.stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

			//F is on the primary side

			//M value at the interface obtained using interpolation
			int idx1_M = M_pri.position_to_cellidx(V_pri.cellidx_to_position(cell1_idx));
			int idx2_M = M_pri.position_to_cellidx(V_pri.cellidx_to_position(cell2_idx));
			m = 1.5 * cu_normalize(M_pri[idx1_M]) - 0.5 * cu_normalize(M_pri[idx2_M]);
		}

		//V values at the interface obtained using interpolation
		cuBReal V_pri_val = 1.5 * V_pri[cell1_idx] - 0.5 * V_pri[cell2_idx];
		cuBReal V_sec_val = 1.5 * V_sec.weighted_average(relpos_m1, stencil) - 0.5 * V_sec.weighted_average(relpos_m1 + shift, stencil);

		cuBReal dV = V_pri_val - V_sec_val;

		return MUB_E * (Gi.i - Gi.j) * dV * m;
	}

	//B-Function for N
	__device__ cuReal33 B_func_sec(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_S_Funcs& trans_sec, TransportCUDA_Spin_S_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Secondary is an N here

		cuVEC_VC<cuReal3>& S_pri = (trans_pri.pcuMesh ? *(trans_pri.pcuMesh->pS) : *(trans_pri.pcuaMesh->pS));

		cuVEC_VC<cuReal3>& M_pri = (trans_pri.pcuMesh ? *(trans_pri.pcuMesh->pM) : *(trans_pri.pcuaMesh->pM1));

		//Get G values from top contacting mesh
		cuReal2 Gi, Gmix;
		if (primary_top) {

			if (trans_pri.pcuMesh) {

				Gi = *(trans_pri.pcuMesh->pGi);
				Gmix = *(trans_pri.pcuMesh->pGmix);
				trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pGi), Gi, *(trans_pri.pcuMesh->pGmix), Gmix);
			}
			else {

				Gi = *(trans_pri.pcuaMesh->pGi);
				Gmix = *(trans_pri.pcuaMesh->pGmix);
				trans_pri.pcuaMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuaMesh->pGi), Gi, *(trans_pri.pcuaMesh->pGmix), Gmix);
			}
		}
		else {

			if (trans_sec.pcuMesh) {

				Gi = *(trans_sec.pcuMesh->pGi);
				Gmix = *(trans_sec.pcuMesh->pGmix);
				trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pGi), Gi, *(trans_sec.pcuMesh->pGmix), Gmix);
			}
			else {

				Gi = *(trans_sec.pcuaMesh->pGi);
				Gmix = *(trans_sec.pcuaMesh->pGmix);
				trans_sec.pcuaMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuaMesh->pGi), Gi, *(trans_sec.pcuaMesh->pGmix), Gmix);
			}
		}

		cuReal3 m;

		if (trans_pri.stsolve == STSOLVE_FERROMAGNETIC || trans_pri.stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

			//F is on the primary side

			//M value at the interface obtained using interpolation
			int idx1_M = M_pri.position_to_cellidx(S_pri.cellidx_to_position(cell1_idx));
			int idx2_M = M_pri.position_to_cellidx(S_pri.cellidx_to_position(cell2_idx));
			m = 1.5 * cu_normalize(M_pri[idx1_M]) - 0.5 * cu_normalize(M_pri[idx2_M]);
		}

		cuBReal gR = 2 * (cuBReal)MUB_E * Gmix.i;
		cuBReal gI = 2 * (cuBReal)MUB_E * Gmix.j;

		cuReal33 eps_m = cu_epsilon3(m);

		return (-(cuBReal)MUB_E * (Gi.i + Gi.j) * (m | m)) + (eps_m * (gR * eps_m + gI * cu_ident<cuReal33>()));
	}

	//B-Function for F
	__device__ cuReal33 B_func_pri(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_S_Funcs& trans_sec, TransportCUDA_Spin_S_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Primary is an F here

		cuVEC_VC<cuReal3>& S_pri = (trans_pri.pcuMesh ? *(trans_pri.pcuMesh->pS) : *(trans_pri.pcuaMesh->pS));

		cuVEC_VC<cuReal3>& M_pri = (trans_pri.pcuMesh ? *(trans_pri.pcuMesh->pM) : *(trans_pri.pcuaMesh->pM1));

		//Get G values from top contacting mesh
		cuReal2 Gi;
		if (primary_top) {

			if (trans_pri.pcuMesh) {

				Gi = *(trans_pri.pcuMesh->pGi);
				trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pGi), Gi);
			}
			else {

				Gi = *(trans_pri.pcuaMesh->pGi);
				trans_pri.pcuaMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuaMesh->pGi), Gi);
			}
		}
		else {

			if (trans_sec.pcuMesh) {

				Gi = *(trans_sec.pcuMesh->pGi);
				trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pGi), Gi);
			}
			else {

				Gi = *(trans_sec.pcuaMesh->pGi);
				trans_sec.pcuaMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuaMesh->pGi), Gi);
			}
		}

		cuReal3 m;

		if (trans_pri.stsolve == STSOLVE_FERROMAGNETIC || trans_pri.stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

			//F is on the primary side

			//M value at the interface obtained using interpolation
			int idx1_M = M_pri.position_to_cellidx(S_pri.cellidx_to_position(cell1_idx));
			int idx2_M = M_pri.position_to_cellidx(S_pri.cellidx_to_position(cell2_idx));
			m = 1.5 * cu_normalize(M_pri[idx1_M]) - 0.5 * cu_normalize(M_pri[idx2_M]);
		}

		return -(cuBReal)MUB_E * (Gi.i + Gi.j) * (m | m);
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

		cuVEC_VC<cuBReal>& V_sec = (trans_sec.pcuMesh ? *(trans_sec.pcuMesh->pV) : *(trans_sec.pcuaMesh->pV));
		cuVEC_VC<cuBReal>& V_pri = (trans_pri.pcuMesh ? *(trans_pri.pcuMesh->pV) : *(trans_pri.pcuaMesh->pV));

		cuVEC_VC<cuReal3>& M_sec = (trans_sec.pcuMesh ? *(trans_sec.pcuMesh->pM) : *(trans_sec.pcuaMesh->pM1));

		cuBReal pump_eff_sec = (trans_sec.pcuMesh ? *(trans_sec.pcuMesh->ppump_eff) : *(trans_sec.pcuaMesh->ppump_eff));

		//Get G values from top contacting mesh
		cuReal2 Gi, Gmix;
		if (primary_top) {

			if (trans_pri.pcuMesh) {

				Gi = *(trans_pri.pcuMesh->pGi);
				Gmix = *(trans_pri.pcuMesh->pGmix);
				trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pGi), Gi, *(trans_pri.pcuMesh->pGmix), Gmix);
			}
			else {

				Gi = *(trans_pri.pcuaMesh->pGi);
				Gmix = *(trans_pri.pcuaMesh->pGmix);
				trans_pri.pcuaMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuaMesh->pGi), Gi, *(trans_pri.pcuaMesh->pGmix), Gmix);
			}
		}
		else {

			if (trans_sec.pcuMesh) {

				Gi = *(trans_sec.pcuMesh->pGi);
				Gmix = *(trans_sec.pcuMesh->pGmix);
				trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pGi), Gi, *(trans_sec.pcuMesh->pGmix), Gmix);
			}
			else {

				Gi = *(trans_sec.pcuaMesh->pGi);
				Gmix = *(trans_sec.pcuaMesh->pGmix);
				trans_sec.pcuaMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuaMesh->pGi), Gi, *(trans_sec.pcuaMesh->pGmix), Gmix);
			}
		}

		cuReal3 m;

		//on the N side must include spin pumping current : subtract it from the A value.
		cuReal3 Js_pump = cuReal3();
		
		if (trans_sec.stsolve == STSOLVE_FERROMAGNETIC || trans_sec.stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

			//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

			if (trans_sec.pcuMesh) trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->ppump_eff), pump_eff_sec);
			else trans_sec.pcuaMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuaMesh->ppump_eff), pump_eff_sec);

			//M value at the interface obtained using interpolation
			m = 1.5 * cu_normalize(M_sec.weighted_average(relpos_m1, stencil)) - 0.5 * cu_normalize(M_sec.weighted_average(relpos_m1 + shift, stencil));

			if (cuIsNZ((cuBReal)pump_eff_sec) && trans_sec.pdM_dt) {

				int idx1_M = M_sec.position_to_cellidx(relpos_m1);
				int idx2_M = M_sec.position_to_cellidx(relpos_m1 + shift);

				cuReal3 dmdt = 1.5 * cu_normalize((*(trans_sec.pdM_dt))[idx1_M], M_sec[idx1_M]) - 0.5 * cu_normalize((*(trans_sec.pdM_dt))[idx2_M], M_sec[idx2_M]);

				Js_pump = pump_eff_sec * (cuBReal)MUB_E * (cuBReal)HBAR_E * (Gmix.i * (m ^ dmdt) + Gmix.j * m);
			}
		}

		//V values at the interface obtained using interpolation
		cuBReal V_pri_val = 1.5 * V_pri[cell1_idx] - 0.5 * V_pri[cell2_idx];
		cuBReal V_sec_val = 1.5 * V_sec.weighted_average(relpos_m1, stencil) - 0.5 * V_sec.weighted_average(relpos_m1 + shift, stencil);
		cuBReal dV = V_pri_val - V_sec_val;

		return (cuBReal)MUB_E * (Gi.i - Gi.j) * dV * m - Js_pump;
	}

	//A-Function for F
	__device__ cuReal3 A_func_sec(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_S_Funcs& trans_sec, TransportCUDA_Spin_S_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Secondary is an F here

		cuVEC_VC<cuBReal>& V_sec = (trans_sec.pcuMesh ? *(trans_sec.pcuMesh->pV) : *(trans_sec.pcuaMesh->pV));
		cuVEC_VC<cuBReal>& V_pri = (trans_pri.pcuMesh ? *(trans_pri.pcuMesh->pV) : *(trans_pri.pcuaMesh->pV));

		cuVEC_VC<cuReal3>& M_sec = (trans_sec.pcuMesh ? *(trans_sec.pcuMesh->pM) : *(trans_sec.pcuaMesh->pM1));

		//Get G values from top contacting mesh
		cuReal2 Gi;
		if (primary_top) {

			if (trans_pri.pcuMesh) {

				Gi = *(trans_pri.pcuMesh->pGi);
				trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pGi), Gi);
			}
			else {

				Gi = *(trans_pri.pcuaMesh->pGi);
				trans_pri.pcuaMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuaMesh->pGi), Gi);
			}
		}
		else {

			if (trans_sec.pcuMesh) {

				Gi = *(trans_sec.pcuMesh->pGi);
				trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pGi), Gi);
			}
			else {

				Gi = *(trans_sec.pcuaMesh->pGi);
				trans_sec.pcuaMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuaMesh->pGi), Gi);
			}
		}

		cuReal3 m;

		if (trans_sec.stsolve == STSOLVE_FERROMAGNETIC || trans_sec.stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

			//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

			//M value at the interface obtained using interpolation
			m = 1.5 * cu_normalize(M_sec.weighted_average(relpos_m1, stencil)) - 0.5 * cu_normalize(M_sec.weighted_average(relpos_m1 + shift, stencil));
		}

		//V values at the interface obtained using interpolation
		cuBReal V_pri_val = 1.5 * V_pri[cell1_idx] - 0.5 * V_pri[cell2_idx];
		cuBReal V_sec_val = 1.5 * V_sec.weighted_average(relpos_m1, stencil) - 0.5 * V_sec.weighted_average(relpos_m1 + shift, stencil);

		cuBReal dV = V_pri_val - V_sec_val;

		return MUB_E * (Gi.i - Gi.j) * dV * m;
	}

	//B-Function for N
	__device__ cuReal33 B_func_pri(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_S_Funcs& trans_sec, TransportCUDA_Spin_S_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Primary is an N here

		cuVEC_VC<cuReal3>& M_sec = (trans_sec.pcuMesh ? *(trans_sec.pcuMesh->pM) : *(trans_sec.pcuaMesh->pM1));

		//Get G values from top contacting mesh
		cuReal2 Gi, Gmix;
		if (primary_top) {

			if (trans_pri.pcuMesh) {

				Gi = *(trans_pri.pcuMesh->pGi);
				Gmix = *(trans_pri.pcuMesh->pGmix);
				trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pGi), Gi, *(trans_pri.pcuMesh->pGmix), Gmix);
			}
			else {

				Gi = *(trans_pri.pcuaMesh->pGi);
				Gmix = *(trans_pri.pcuaMesh->pGmix);
				trans_pri.pcuaMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuaMesh->pGi), Gi, *(trans_pri.pcuaMesh->pGmix), Gmix);
			}
		}
		else {

			if (trans_sec.pcuMesh) {

				Gi = *(trans_sec.pcuMesh->pGi);
				Gmix = *(trans_sec.pcuMesh->pGmix);
				trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pGi), Gi, *(trans_sec.pcuMesh->pGmix), Gmix);
			}
			else {

				Gi = *(trans_sec.pcuaMesh->pGi);
				Gmix = *(trans_sec.pcuaMesh->pGmix);
				trans_sec.pcuaMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuaMesh->pGi), Gi, *(trans_sec.pcuaMesh->pGmix), Gmix);
			}
		}

		cuReal3 m;

		if (trans_sec.stsolve == STSOLVE_FERROMAGNETIC || trans_sec.stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

			//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

			//M value at the interface obtained using interpolation
			m = 1.5 * cu_normalize(M_sec.weighted_average(relpos_m1, stencil)) - 0.5 * cu_normalize(M_sec.weighted_average(relpos_m1 + shift, stencil));
		}

		cuBReal gR = 2 * (cuBReal)MUB_E * Gmix.i;
		cuBReal gI = 2 * (cuBReal)MUB_E * Gmix.j;

		cuReal33 eps_m = cu_epsilon3(m);

		return (-(cuBReal)MUB_E * (Gi.i + Gi.j) * (m | m)) + (eps_m * (gR * eps_m + gI * cu_ident<cuReal33>()));
	}

	//B-Function for F
	__device__ cuReal33 B_func_sec(int cell1_idx, int cell2_idx, cuReal3 relpos_m1, cuReal3 shift, cuReal3 stencil, TransportCUDA_Spin_S_Funcs& trans_sec, TransportCUDA_Spin_S_Funcs& trans_pri) const
	{
		//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
		bool primary_top = !(shift >= cuReal3(0));

		//Secondary is an F here

		cuVEC_VC<cuReal3>& M_sec = (trans_sec.pcuMesh ? *(trans_sec.pcuMesh->pM) : *(trans_sec.pcuaMesh->pM1));

		//Get G values from top contacting mesh
		cuReal2 Gi;
		if (primary_top) {

			if (trans_pri.pcuMesh) {

				Gi = *(trans_pri.pcuMesh->pGi);
				trans_pri.pcuMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuMesh->pGi), Gi);
			}
			else {

				Gi = *(trans_pri.pcuaMesh->pGi);
				trans_pri.pcuaMesh->update_parameters_ecoarse(cell1_idx, *(trans_pri.pcuaMesh->pGi), Gi);
			}
		}
		else {

			if (trans_sec.pcuMesh) {

				Gi = *(trans_sec.pcuMesh->pGi);
				trans_sec.pcuMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuMesh->pGi), Gi);
			}
			else {

				Gi = *(trans_sec.pcuaMesh->pGi);
				trans_sec.pcuaMesh->update_parameters_atposition(relpos_m1, *(trans_sec.pcuaMesh->pGi), Gi);
			}
		}

		cuReal3 m;

		if (trans_sec.stsolve == STSOLVE_FERROMAGNETIC || trans_sec.stsolve == STSOLVE_FERROMAGNETIC_ATOM) {

			//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

			//M value at the interface obtained using interpolation
			m = 1.5 * cu_normalize(M_sec.weighted_average(relpos_m1, stencil)) - 0.5 * cu_normalize(M_sec.weighted_average(relpos_m1 + shift, stencil));
		}

		return -(cuBReal)MUB_E * (Gi.i + Gi.j) * (m | m);
	}
};

#endif

#endif