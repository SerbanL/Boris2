#include "stdafx.h"
#include "STransport.h"

#ifdef MODULE_TRANSPORT

#include "Mesh.h"
#include "MeshParamsControl.h"

//functions to specify boundary conditions for interface conductance approach for charge : Jc_N = Jc_F = A + B * dV
//VERIFIED - CORRECT
double STransport::Afunc_V(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, Transport& trans_sec, Transport& trans_pri) const
{
	//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
	bool primary_top = !(shift >= DBL3(0));

	//primary and secondary meshes
	Mesh* pMesh_pri = trans_pri.pMesh;
	Mesh* pMesh_sec = trans_sec.pMesh;

	//S values at the interface obtained using interpolation based on current values
	DBL3 S_pri = 1.5 * pMesh_pri->S[cell1_idx] - 0.5 * pMesh_pri->S[cell2_idx];
	DBL3 S_sec = 1.5 * pMesh_sec->S.weighted_average(relpos_m1, stencil) - 0.5 * pMesh_sec->S.weighted_average(relpos_m1 + shift, stencil);

	//conductivity values at the interface obtained using interpolation
	double elC_pri = 1.5 * pMesh_pri->elC[cell1_idx] - 0.5 * pMesh_pri->elC[cell2_idx];
	double elC_sec = 1.5 * pMesh_sec->elC.weighted_average(relpos_m1, stencil) - 0.5 * pMesh_sec->elC.weighted_average(relpos_m1 + shift, stencil);

	//Get G values from top contacting mesh
	DBL2 Gi;
	if (primary_top) {

		Gi = pMesh_pri->Gi;
		pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Gi, Gi);
	}
	else {

		Gi = pMesh_sec->Gi;
		pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gi, Gi);
	}

	if (trans_pri.stsolve == STSOLVE_FERROMAGNETIC) {

		//F is on the primary side

		double Ms = pMesh_pri->Ms;
		double De_pri = pMesh_pri->De;
		double De_sec = pMesh_sec->De;
		pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Ms, Ms, pMesh_pri->De, De_pri);
		pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->De, De_sec);

		//dVs is the difference between the primary and secondary Vs values at the interface
		DBL3 dVs = (De_pri * S_pri / elC_pri - De_sec * S_sec / elC_sec) / MUB_E;
		
		//M value at the interface obtained using interpolation
		int idx1_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell1_idx));
		int idx2_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell2_idx));
		DBL3 m = (1.5 * pMesh_pri->M[idx1_M] - 0.5 * pMesh_pri->M[idx2_M]) / Ms;

		return (Gi.i - Gi.j) * (dVs * m);
	}
	else {

		//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

		double Ms = pMesh_sec->Ms;
		double De_pri = pMesh_pri->De;
		double De_sec = pMesh_sec->De;
		pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->De, De_pri);
		pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Ms, Ms, pMesh_sec->De, De_sec);

		//dVs is the difference between the primary and secondary Vs values at the interface
		DBL3 dVs = (De_pri * S_pri / elC_pri - De_sec * S_sec / elC_sec) / MUB_E;

		//M value at the interface obtained using interpolation
		DBL3 m = (1.5 * pMesh_sec->M.weighted_average(relpos_m1, stencil) - 0.5 * pMesh_sec->M.weighted_average(relpos_m1 + shift, stencil)) / Ms;

		return (Gi.i - Gi.j) * (dVs * m);
	}
}

//functions to specify boundary conditions for interface conductance approach for charge : Jc_N = Jc_F = A + B * dV
//VERIFIED - CORRECT
double STransport::Bfunc_V(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, Transport& trans_sec, Transport& trans_pri) const
{
	//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
	bool primary_top = !(shift >= DBL3(0));

	//Get G values from top contacting mesh
	DBL2 Gi;
	if (primary_top) {

		Mesh* pMesh_pri = trans_pri.pMesh;

		Gi = pMesh_pri->Gi;
		pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Gi, Gi);
	}
	else {

		Mesh* pMesh_sec = trans_sec.pMesh;
		
		Gi = pMesh_sec->Gi;
		pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gi, Gi);
	}

	return -(Gi.i + Gi.j);
}

//functions to specify boundary conditions for interface conductance approach for spin : Js_N = A_N + B_N * dVs, Js_F = A_F + B_F * dVs#
//VERIFIED - CORRECT
DBL3 STransport::Afunc_N_S(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, Transport& trans_sec, Transport& trans_pri) const
{
	//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
	bool primary_top = !(shift >= DBL3(0));

	//primary and secondary meshes
	Mesh* pMesh_pri = trans_pri.pMesh;
	Mesh* pMesh_sec = trans_sec.pMesh;
	
	//Get G values from top contacting mesh
	DBL2 Gi, Gmix;
	if (primary_top) {

		Gi = pMesh_pri->Gi;
		Gmix = pMesh_pri->Gmix;
		pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Gi, Gi, pMesh_pri->Gmix, Gmix);
	}
	else {

		Gi = pMesh_sec->Gi;
		Gmix = pMesh_sec->Gmix;
		pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gi, Gi, pMesh_sec->Gmix, Gmix);
	}


	DBL3 m;
	//on the N side must include spin pumping current : subtract it from the A value.
	DBL3 Js_pump = DBL3();

	if (trans_pri.stsolve == STSOLVE_FERROMAGNETIC) {

		//F is on the primary side

		double Ms = pMesh_pri->Ms;
		double pump_eff = pMesh_pri->pump_eff;
		pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Ms, Ms, pMesh_pri->pump_eff, pump_eff);

		//M value at the interface obtained using interpolation
		int idx1_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell1_idx));
		int idx2_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell2_idx));
		m = (1.5 * pMesh_pri->M[idx1_M] - 0.5 * pMesh_pri->M[idx2_M]) / Ms;

		if (IsNZ((double)pump_eff)) {

			DBL3 dmdt = (1.5 * pMesh_pri->dMdt(idx1_M) - 0.5 * pMesh_pri->dMdt(idx2_M)) / Ms;

			Js_pump = pump_eff * MUB_E * HBAR_E * (Gmix.i * (m ^ dmdt) + Gmix.j * m);
		}
	}
	else {

		//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

		double Ms = pMesh_sec->Ms;
		double pump_eff = pMesh_sec->pump_eff;
		pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Ms, Ms, pMesh_sec->pump_eff, pump_eff);

		//M value at the interface obtained using interpolation
		m = (1.5 * pMesh_sec->M.weighted_average(relpos_m1, stencil) - 0.5 * pMesh_sec->M.weighted_average(relpos_m1 + shift, stencil)) / Ms;

		if (IsNZ((double)pump_eff)) {

			int idx1_M = pMesh_sec->M.position_to_cellidx(relpos_m1);
			int idx2_M = pMesh_sec->M.position_to_cellidx(relpos_m1 + shift);

			DBL3 dmdt = (1.5 * pMesh_sec->dMdt(idx1_M) - 0.5 * pMesh_sec->dMdt(idx2_M)) / Ms;

			Js_pump = pump_eff * MUB_E * HBAR_E * (Gmix.i * (m ^ dmdt) + Gmix.j * m);
		}
	}

	//V values at the interface obtained using interpolation
	double V_pri = 1.5 * pMesh_pri->V[cell1_idx] - 0.5 * pMesh_pri->V[cell2_idx];
	double V_sec = 1.5 * pMesh_sec->V.weighted_average(relpos_m1, stencil) - 0.5 * pMesh_sec->V.weighted_average(relpos_m1 + shift, stencil);
	double dV = V_pri - V_sec;

	return MUB_E * (Gi.i - Gi.j) * dV * m - Js_pump;
}

//functions to specify boundary conditions for interface conductance approach for spin : Js_N = A_N + B_N * dVs, Js_F = A_F + B_F * dVs
//VERIFIED - CORRECT
DBL3 STransport::Afunc_F_S(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, Transport& trans_sec, Transport& trans_pri) const
{
	//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
	bool primary_top = !(shift >= DBL3(0));

	//primary and secondary meshes
	Mesh* pMesh_pri = trans_pri.pMesh;
	Mesh* pMesh_sec = trans_sec.pMesh;

	//Get G values from top contacting mesh
	DBL2 Gi;
	if (primary_top) {

		Gi = pMesh_pri->Gi;
		pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Gi, Gi);
	}
	else {

		Gi = pMesh_sec->Gi;
		pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gi, Gi);
	}

	DBL3 m;

	if (trans_pri.stsolve == STSOLVE_FERROMAGNETIC) {
		
		//F is on the primary side

		double Ms = pMesh_pri->Ms;
		pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Ms, Ms);

		//M value at the interface obtained using interpolation
		int idx1_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell1_idx));
		int idx2_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell2_idx));
		m = (1.5 * pMesh_pri->M[idx1_M] - 0.5 * pMesh_pri->M[idx2_M]) / Ms;
	}
	else {

		//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)
		
		double Ms = pMesh_sec->Ms;
		pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Ms, Ms);

		//M value at the interface obtained using interpolation
		m = (1.5 * pMesh_sec->M.weighted_average(relpos_m1, stencil) - 0.5 * pMesh_sec->M.weighted_average(relpos_m1 + shift, stencil)) / Ms;
	}

	//V values at the interface obtained using interpolation
	double V_pri = 1.5 * pMesh_pri->V[cell1_idx] - 0.5 * pMesh_pri->V[cell2_idx];
	double V_sec = 1.5 * pMesh_sec->V.weighted_average(relpos_m1, stencil) - 0.5 * pMesh_sec->V.weighted_average(relpos_m1 + shift, stencil);

	double dV = V_pri - V_sec;

	return MUB_E * (Gi.i - Gi.j) * dV * m;
}

//functions to specify boundary conditions for interface conductance approach for spin : Js_N = A_N + B_N * dVs, Js_F = A_F + B_F * dVs
//VERIFIED - CORRECT
DBL33 STransport::Bfunc_N_S(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, Transport& trans_sec, Transport& trans_pri) const
{
	//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
	bool primary_top = !(shift >= DBL3(0));

	//primary and secondary meshes
	Mesh* pMesh_pri = trans_pri.pMesh;
	Mesh* pMesh_sec = trans_sec.pMesh;

	//Get G values from top contacting mesh
	DBL2 Gi, Gmix;
	if (primary_top) {
		
		Gi = pMesh_pri->Gi;
		Gmix = pMesh_pri->Gmix;
		pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Gi, Gi, pMesh_pri->Gmix, Gmix);
	}
	else {

		Gi = pMesh_sec->Gi;
		Gmix = pMesh_sec->Gmix;
		pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gi, Gi, pMesh_sec->Gmix, Gmix);
	}

	DBL3 m;

	if (trans_pri.stsolve == STSOLVE_FERROMAGNETIC) {

		//F is on the primary side

		double Ms = pMesh_pri->Ms;
		pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Ms, Ms);

		//M value at the interface obtained using interpolation
		int idx1_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell1_idx));
		int idx2_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell2_idx));
		m = (1.5 * pMesh_pri->M[idx1_M] - 0.5 * pMesh_pri->M[idx2_M]) / Ms;
	}
	else {

		//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

		double Ms = pMesh_sec->Ms;
		pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Ms, Ms);

		//M value at the interface obtained using interpolation
		m = (1.5 * pMesh_sec->M.weighted_average(relpos_m1, stencil) - 0.5 * pMesh_sec->M.weighted_average(relpos_m1 + shift, stencil)) / pMesh_sec->Ms;
	}

	double gR = 2 * MUB_E * Gmix.i;
	double gI = 2 * MUB_E * Gmix.j;

	DBL33 eps_m = epsilon3(m);

	return (-MUB_E * (Gi.i + Gi.j) * (m | m)) + (eps_m * (gR * eps_m + gI * ident<DBL33>()));
}

//functions to specify boundary conditions for interface conductance approach for spin : Js_N = A_N + B_N * dVs, Js_F = A_F + B_F * dVs
//VERIFIED - CORRECT
DBL33 STransport::Bfunc_F_S(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, Transport& trans_sec, Transport& trans_pri) const
{
	//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
	bool primary_top = !(shift >= DBL3(0));

	//primary and secondary meshes
	Mesh* pMesh_pri = trans_pri.pMesh;
	Mesh* pMesh_sec = trans_sec.pMesh;

	//Get G values from top contacting mesh
	DBL2 Gi;
	if (primary_top) {

		Gi = pMesh_pri->Gi;
		pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Gi, Gi);
	}
	else {

		
		Gi = pMesh_sec->Gi;
		pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gi, Gi);
	}

	DBL3 m;

	if (trans_pri.stsolve == STSOLVE_FERROMAGNETIC) {

		//F is on the primary side

		double Ms = pMesh_pri->Ms;
		pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Ms, Ms);

		//M value at the interface obtained using interpolation
		int idx1_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell1_idx));
		int idx2_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell2_idx));
		m = (1.5 * pMesh_pri->M[idx1_M] - 0.5 * pMesh_pri->M[idx2_M]) / Ms;
	}
	else {

		//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

		double Ms = pMesh_sec->Ms;
		pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Ms, Ms);

		//M value at the interface obtained using interpolation
		m = (1.5 * pMesh_sec->M.weighted_average(relpos_m1, stencil) - 0.5 * pMesh_sec->M.weighted_average(relpos_m1 + shift, stencil)) / Ms;
	}

	return -MUB_E * (Gi.i + Gi.j) * (m | m);
}

#endif