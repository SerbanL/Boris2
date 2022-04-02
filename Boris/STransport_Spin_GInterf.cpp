#include "stdafx.h"
#include "STransport.h"

#ifdef MODULE_COMPILATION_TRANSPORT

#include "Mesh.h"
#include "MeshParamsControl.h"
#include "Atom_Mesh.h"
#include "Atom_MeshParamsControl.h"

//functions to specify boundary conditions for interface conductance approach for charge : Jc_N = Jc_F = A + B * dV
//VERIFIED - CORRECT
double STransport::Afunc_V(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, TransportBase& trans_sec, TransportBase& trans_pri) const
{
	//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
	bool primary_top = !(shift >= DBL3(0));

	//primary and secondary meshes
	Mesh *pMesh_pri = nullptr, *pMesh_sec = nullptr;
	Atom_Mesh *paMesh_pri = nullptr, *paMesh_sec = nullptr;

	if (!trans_pri.pMeshBase->is_atomistic()) pMesh_pri = dynamic_cast<Mesh*>(trans_pri.pMeshBase);
	else paMesh_pri = dynamic_cast<Atom_Mesh*>(trans_pri.pMeshBase);

	if (!trans_sec.pMeshBase->is_atomistic()) pMesh_sec = dynamic_cast<Mesh*>(trans_sec.pMeshBase);
	else paMesh_sec = dynamic_cast<Atom_Mesh*>(trans_sec.pMeshBase);

	//S values at the interface obtained using interpolation based on current values
	DBL3 S_pri = 1.5 * trans_pri.pMeshBase->S[cell1_idx] - 0.5 * trans_pri.pMeshBase->S[cell2_idx];
	DBL3 S_sec = 1.5 * trans_sec.pMeshBase->S.weighted_average(relpos_m1, stencil) - 0.5 * trans_sec.pMeshBase->S.weighted_average(relpos_m1 + shift, stencil);

	//conductivity values at the interface obtained using interpolation
	double elC_pri = 1.5 * trans_pri.pMeshBase->elC[cell1_idx] - 0.5 * trans_pri.pMeshBase->elC[cell2_idx];
	double elC_sec = 1.5 * trans_sec.pMeshBase->elC.weighted_average(relpos_m1, stencil) - 0.5 * trans_sec.pMeshBase->elC.weighted_average(relpos_m1 + shift, stencil);

	//Get G values from top contacting mesh
	DBL2 Gi;
	if (primary_top) {

		if (pMesh_pri) {

			Gi = pMesh_pri->Gi;
			pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Gi, Gi);
		}
		else {

			Gi = paMesh_pri->Gi;
			paMesh_pri->update_parameters_ecoarse(cell1_idx, paMesh_pri->Gi, Gi);
		}
	}
	else {

		if (pMesh_sec) {

			Gi = pMesh_sec->Gi;
			pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gi, Gi);
		}
		else {

			Gi = paMesh_sec->Gi;
			paMesh_sec->update_parameters_atposition(relpos_m1, paMesh_sec->Gi, Gi);
		}
	}

	if (trans_pri.Get_STSolveType() == STSOLVE_FERROMAGNETIC || trans_pri.Get_STSolveType() == STSOLVE_FERROMAGNETIC_ATOM) {

		//F is on the primary side
		double De_pri, De_sec;
		DBL3 m;

		if (pMesh_pri) {

			De_pri = pMesh_pri->De;
			pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->De, De_pri);

			//m value at the interface obtained using interpolation
			int idx1_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell1_idx));
			int idx2_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell2_idx));
			
			m = 1.5 * normalize(pMesh_pri->M[idx1_M]) - 0.5 * normalize(pMesh_pri->M[idx2_M]);
		}
		else {

			double mu_s = paMesh_pri->mu_s;
			De_pri = paMesh_pri->De;
			paMesh_pri->update_parameters_ecoarse(cell1_idx, paMesh_pri->mu_s, mu_s, paMesh_pri->De, De_pri);

			//m value at the interface obtained using interpolation
			int idx1_M = paMesh_pri->M1.position_to_cellidx(paMesh_pri->S.cellidx_to_position(cell1_idx));
			int idx2_M = paMesh_pri->M1.position_to_cellidx(paMesh_pri->S.cellidx_to_position(cell2_idx));
			m = (1.5 * paMesh_pri->M1[idx1_M] - 0.5 * paMesh_pri->M1[idx2_M]) / mu_s;
		}

		if (pMesh_sec) {

			De_sec = pMesh_sec->De;
			pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->De, De_sec);

			if (trans_sec.Get_STSolveType() == STSOLVE_TUNNELING) {

				//for tunnel barriers, set De to 1.0 if not metallic pinhole (for pinhole elecCond > 0)
				double elecCond = pMesh_sec->elecCond;
				pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->elecCond, elecCond);
				if (IsZ(elecCond)) De_sec = 1.0;
			}
		}
		else {

			De_sec = paMesh_sec->De;
			paMesh_sec->update_parameters_atposition(relpos_m1, paMesh_sec->De, De_sec);
		}

		//dVs is the difference between the primary and secondary Vs values at the interface
		DBL3 dVs = (De_pri * S_pri / elC_pri - De_sec * S_sec / elC_sec) / MUB_E;

		return (Gi.i - Gi.j) * (dVs * m);
	}
	else {

		//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

		double De_pri, De_sec;
		DBL3 m;

		if (pMesh_pri) {

			De_pri = pMesh_pri->De;
			pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->De, De_pri);

			if (trans_pri.Get_STSolveType() == STSOLVE_TUNNELING) {

				//for tunnel barriers, set De to 1.0 if not metallic pinhole (for pinhole elecCond > 0)
				double elecCond = pMesh_pri->elecCond;
				pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->elecCond, elecCond);
				if (IsZ(elecCond)) De_pri = 1.0;
			}
		}
		else {

			De_pri = paMesh_pri->De;
			paMesh_pri->update_parameters_ecoarse(cell1_idx, paMesh_pri->De, De_pri);
		}

		if (pMesh_sec) {

			De_sec = pMesh_sec->De;
			pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->De, De_sec);

			//m value at the interface obtained using interpolation
			m = 1.5 * normalize(pMesh_sec->M.weighted_average(relpos_m1, stencil)) - 0.5 * normalize(pMesh_sec->M.weighted_average(relpos_m1 + shift, stencil));
		}
		else {

			double mu_s = paMesh_sec->mu_s;
			De_sec = paMesh_sec->De;
			paMesh_sec->update_parameters_atposition(relpos_m1, paMesh_sec->mu_s, mu_s, paMesh_sec->De, De_sec);

			//m value at the interface obtained using interpolation
			m = (1.5 * paMesh_sec->M1.weighted_average(relpos_m1, stencil) - 0.5 * paMesh_sec->M1.weighted_average(relpos_m1 + shift, stencil)) / mu_s;
		}

		//dVs is the difference between the primary and secondary Vs values at the interface
		DBL3 dVs = (De_pri * S_pri / elC_pri - De_sec * S_sec / elC_sec) / MUB_E;

		return (Gi.i - Gi.j) * (dVs * m);
	}
}

//functions to specify boundary conditions for interface conductance approach for charge : Jc_N = Jc_F = A + B * dV
//VERIFIED - CORRECT
double STransport::Bfunc_V(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, TransportBase& trans_sec, TransportBase& trans_pri) const
{
	//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
	bool primary_top = !(shift >= DBL3(0));

	//primary and secondary meshes
	Mesh *pMesh_pri = nullptr, *pMesh_sec = nullptr;
	Atom_Mesh *paMesh_pri = nullptr, *paMesh_sec = nullptr;

	if (!trans_pri.pMeshBase->is_atomistic()) pMesh_pri = dynamic_cast<Mesh*>(trans_pri.pMeshBase);
	else paMesh_pri = dynamic_cast<Atom_Mesh*>(trans_pri.pMeshBase);

	if (!trans_sec.pMeshBase->is_atomistic()) pMesh_sec = dynamic_cast<Mesh*>(trans_sec.pMeshBase);
	else paMesh_sec = dynamic_cast<Atom_Mesh*>(trans_sec.pMeshBase);

	//Get G values from top contacting mesh
	DBL2 Gi;
	if (primary_top) {

		if (pMesh_pri) {

			Gi = pMesh_pri->Gi;
			pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Gi, Gi);
		}
		else {

			Gi = paMesh_pri->Gi;
			paMesh_pri->update_parameters_ecoarse(cell1_idx, paMesh_pri->Gi, Gi);
		}
	}
	else {

		if (pMesh_sec) {

			Gi = pMesh_sec->Gi;
			pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gi, Gi);
		}
		else {

			Gi = paMesh_sec->Gi;
			paMesh_sec->update_parameters_atposition(relpos_m1, paMesh_sec->Gi, Gi);
		}
	}

	return -(Gi.i + Gi.j);
}

//functions to specify boundary conditions for interface conductance approach for spin : Js_N = A_N + B_N * dVs, Js_F = A_F + B_F * dVs
//VERIFIED - CORRECT
DBL3 STransport::Afunc_N_S(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, TransportBase& trans_sec, TransportBase& trans_pri) const
{
	//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
	bool primary_top = !(shift >= DBL3(0));

	//primary and secondary meshes
	Mesh *pMesh_pri = nullptr, *pMesh_sec = nullptr;
	Atom_Mesh *paMesh_pri = nullptr, *paMesh_sec = nullptr;

	if (!trans_pri.pMeshBase->is_atomistic()) pMesh_pri = dynamic_cast<Mesh*>(trans_pri.pMeshBase);
	else paMesh_pri = dynamic_cast<Atom_Mesh*>(trans_pri.pMeshBase);

	if (!trans_sec.pMeshBase->is_atomistic()) pMesh_sec = dynamic_cast<Mesh*>(trans_sec.pMeshBase);
	else paMesh_sec = dynamic_cast<Atom_Mesh*>(trans_sec.pMeshBase);

	//Get G values from top contacting mesh
	DBL2 Gi, Gmix;
	if (primary_top) {

		if (pMesh_pri) {

			Gi = pMesh_pri->Gi;
			Gmix = pMesh_pri->Gmix;
			pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Gi, Gi, pMesh_pri->Gmix, Gmix);
		}
		else {

			Gi = paMesh_pri->Gi;
			Gmix = paMesh_pri->Gmix;
			paMesh_pri->update_parameters_ecoarse(cell1_idx, paMesh_pri->Gi, Gi, paMesh_pri->Gmix, Gmix);
		}
	}
	else {

		if (pMesh_sec) {

			Gi = pMesh_sec->Gi;
			Gmix = pMesh_sec->Gmix;
			pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gi, Gi, pMesh_sec->Gmix, Gmix);
		}
		else {

			Gi = paMesh_sec->Gi;
			Gmix = paMesh_sec->Gmix;
			paMesh_sec->update_parameters_atposition(relpos_m1, paMesh_sec->Gi, Gi, paMesh_sec->Gmix, Gmix);
		}
	}

	DBL3 m;
	//on the N side must include spin pumping current : subtract it from the A value.
	DBL3 Js_pump = DBL3();

	if (trans_pri.Get_STSolveType() == STSOLVE_FERROMAGNETIC || trans_pri.Get_STSolveType() == STSOLVE_FERROMAGNETIC_ATOM) {

		//F is on the primary side
		if (pMesh_pri) {

			double pump_eff = pMesh_pri->pump_eff;
			pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->pump_eff, pump_eff);

			//M value at the interface obtained using interpolation
			int idx1_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell1_idx));
			int idx2_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell2_idx));
			m = 1.5 * normalize(pMesh_pri->M[idx1_M]) - 0.5 * normalize(pMesh_pri->M[idx2_M]);

			if (IsNZ((double)pump_eff)) {

				DBL3 dmdt = 1.5 * normalize(pMesh_pri->dMdt(idx1_M), pMesh_pri->M[idx1_M]) - 0.5 * normalize(pMesh_pri->dMdt(idx2_M), pMesh_pri->M[idx2_M]);

				Js_pump = pump_eff * MUB_E * HBAR_E * (Gmix.i * (m ^ dmdt) + Gmix.j * m);
			}
		}
		else {

			double mu_s = paMesh_pri->mu_s;
			double pump_eff = paMesh_pri->pump_eff;
			paMesh_pri->update_parameters_ecoarse(cell1_idx, paMesh_pri->mu_s, mu_s, paMesh_pri->pump_eff, pump_eff);

			//M value at the interface obtained using interpolation
			int idx1_M = paMesh_pri->M1.position_to_cellidx(paMesh_pri->S.cellidx_to_position(cell1_idx));
			int idx2_M = paMesh_pri->M1.position_to_cellidx(paMesh_pri->S.cellidx_to_position(cell2_idx));
			m = (1.5 * paMesh_pri->M1[idx1_M] - 0.5 * paMesh_pri->M1[idx2_M]) / mu_s;

			if (IsNZ((double)pump_eff)) {

				DBL3 dmdt = (1.5 * paMesh_pri->dMdt(idx1_M) - 0.5 * paMesh_pri->dMdt(idx2_M)) / mu_s;

				Js_pump = pump_eff * MUB_E * HBAR_E * (Gmix.i * (m ^ dmdt) + Gmix.j * m);
			}
		}
	}
	else {

		//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)
		if (pMesh_sec) {

			double pump_eff = pMesh_sec->pump_eff;
			pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->pump_eff, pump_eff);

			//M value at the interface obtained using interpolation
			m = 1.5 * normalize(pMesh_sec->M.weighted_average(relpos_m1, stencil)) - 0.5 * normalize(pMesh_sec->M.weighted_average(relpos_m1 + shift, stencil));

			if (IsNZ((double)pump_eff)) {

				int idx1_M = pMesh_sec->M.position_to_cellidx(relpos_m1);
				int idx2_M = pMesh_sec->M.position_to_cellidx(relpos_m1 + shift);

				DBL3 dmdt = 1.5 * normalize(pMesh_sec->dMdt(idx1_M), pMesh_sec->M[idx1_M]) - 0.5 * normalize(pMesh_sec->dMdt(idx2_M), pMesh_sec->M[idx2_M]);

				Js_pump = pump_eff * MUB_E * HBAR_E * (Gmix.i * (m ^ dmdt) + Gmix.j * m);
			}
		}
		else {

			double mu_s = paMesh_sec->mu_s;
			double pump_eff = paMesh_sec->pump_eff;
			paMesh_sec->update_parameters_atposition(relpos_m1, paMesh_sec->mu_s, mu_s, paMesh_sec->pump_eff, pump_eff);

			//M value at the interface obtained using interpolation
			m = (1.5 * paMesh_sec->M1.weighted_average(relpos_m1, stencil) - 0.5 * paMesh_sec->M1.weighted_average(relpos_m1 + shift, stencil)) / mu_s;

			if (IsNZ((double)pump_eff)) {

				int idx1_M = paMesh_sec->M1.position_to_cellidx(relpos_m1);
				int idx2_M = paMesh_sec->M1.position_to_cellidx(relpos_m1 + shift);

				DBL3 dmdt = (1.5 * paMesh_sec->dMdt(idx1_M) - 0.5 * paMesh_sec->dMdt(idx2_M)) / mu_s;

				Js_pump = pump_eff * MUB_E * HBAR_E * (Gmix.i * (m ^ dmdt) + Gmix.j * m);
			}
		}
	}

	//V values at the interface obtained using interpolation
	double V_pri = 1.5 * trans_pri.pMeshBase->V[cell1_idx] - 0.5 * trans_pri.pMeshBase->V[cell2_idx];
	double V_sec = 1.5 * trans_sec.pMeshBase->V.weighted_average(relpos_m1, stencil) - 0.5 * trans_sec.pMeshBase->V.weighted_average(relpos_m1 + shift, stencil);
	double dV = V_pri - V_sec;

	return MUB_E * (Gi.i - Gi.j) * dV * m - Js_pump;
}

//functions to specify boundary conditions for interface conductance approach for spin : Js_N = A_N + B_N * dVs, Js_F = A_F + B_F * dVs
//VERIFIED - CORRECT
DBL3 STransport::Afunc_F_S(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, TransportBase& trans_sec, TransportBase& trans_pri) const
{
	//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
	bool primary_top = !(shift >= DBL3(0));

	//primary and secondary meshes
	Mesh *pMesh_pri = nullptr, *pMesh_sec = nullptr;
	Atom_Mesh *paMesh_pri = nullptr, *paMesh_sec = nullptr;

	if (!trans_pri.pMeshBase->is_atomistic()) pMesh_pri = dynamic_cast<Mesh*>(trans_pri.pMeshBase);
	else paMesh_pri = dynamic_cast<Atom_Mesh*>(trans_pri.pMeshBase);

	if (!trans_sec.pMeshBase->is_atomistic()) pMesh_sec = dynamic_cast<Mesh*>(trans_sec.pMeshBase);
	else paMesh_sec = dynamic_cast<Atom_Mesh*>(trans_sec.pMeshBase);

	//Get G values from top contacting mesh
	DBL2 Gi;
	if (primary_top) {

		if (pMesh_pri) {

			Gi = pMesh_pri->Gi;
			pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Gi, Gi);
		}
		else {

			Gi = paMesh_pri->Gi;
			paMesh_pri->update_parameters_ecoarse(cell1_idx, paMesh_pri->Gi, Gi);
		}
	}
	else {

		if (pMesh_sec) {

			Gi = pMesh_sec->Gi;
			pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gi, Gi);
		}
		else {

			Gi = paMesh_sec->Gi;
			paMesh_sec->update_parameters_atposition(relpos_m1, paMesh_sec->Gi, Gi);
		}
	}

	DBL3 m;

	if (trans_pri.Get_STSolveType() == STSOLVE_FERROMAGNETIC || trans_pri.Get_STSolveType() == STSOLVE_FERROMAGNETIC_ATOM) {

		//F is on the primary side
		if (pMesh_pri) {

			//M value at the interface obtained using interpolation
			int idx1_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell1_idx));
			int idx2_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell2_idx));
			m = 1.5 * normalize(pMesh_pri->M[idx1_M]) - 0.5 * normalize(pMesh_pri->M[idx2_M]);
		}
		else {

			double mu_s = paMesh_pri->mu_s;
			paMesh_pri->update_parameters_ecoarse(cell1_idx, paMesh_pri->mu_s, mu_s);

			//M value at the interface obtained using interpolation
			int idx1_M = paMesh_pri->M1.position_to_cellidx(paMesh_pri->S.cellidx_to_position(cell1_idx));
			int idx2_M = paMesh_pri->M1.position_to_cellidx(paMesh_pri->S.cellidx_to_position(cell2_idx));
			m = (1.5 * paMesh_pri->M1[idx1_M] - 0.5 * paMesh_pri->M1[idx2_M]) / mu_s;
		}
	}
	else {

		//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

		if (pMesh_sec) {

			//M value at the interface obtained using interpolation
			m = 1.5 * normalize(pMesh_sec->M.weighted_average(relpos_m1, stencil)) - 0.5 * normalize(pMesh_sec->M.weighted_average(relpos_m1 + shift, stencil));
		}
		else {

			double mu_s = paMesh_sec->mu_s;
			paMesh_sec->update_parameters_atposition(relpos_m1, paMesh_sec->mu_s, mu_s);

			//M value at the interface obtained using interpolation
			m = (1.5 * paMesh_sec->M1.weighted_average(relpos_m1, stencil) - 0.5 * paMesh_sec->M1.weighted_average(relpos_m1 + shift, stencil)) / mu_s;
		}
	}

	//V values at the interface obtained using interpolation
	double V_pri = 1.5 * trans_pri.pMeshBase->V[cell1_idx] - 0.5 * trans_pri.pMeshBase->V[cell2_idx];
	double V_sec = 1.5 * trans_sec.pMeshBase->V.weighted_average(relpos_m1, stencil) - 0.5 * trans_sec.pMeshBase->V.weighted_average(relpos_m1 + shift, stencil);

	double dV = V_pri - V_sec;

	return MUB_E * (Gi.i - Gi.j) * dV * m;
}

//functions to specify boundary conditions for interface conductance approach for spin : Js_N = A_N + B_N * dVs, Js_F = A_F + B_F * dVs
//VERIFIED - CORRECT
DBL33 STransport::Bfunc_N_S(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, TransportBase& trans_sec, TransportBase& trans_pri) const
{
	//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
	bool primary_top = !(shift >= DBL3(0));

	//primary and secondary meshes
	Mesh *pMesh_pri = nullptr, *pMesh_sec = nullptr;
	Atom_Mesh *paMesh_pri = nullptr, *paMesh_sec = nullptr;

	if (!trans_pri.pMeshBase->is_atomistic()) pMesh_pri = dynamic_cast<Mesh*>(trans_pri.pMeshBase);
	else paMesh_pri = dynamic_cast<Atom_Mesh*>(trans_pri.pMeshBase);

	if (!trans_sec.pMeshBase->is_atomistic()) pMesh_sec = dynamic_cast<Mesh*>(trans_sec.pMeshBase);
	else paMesh_sec = dynamic_cast<Atom_Mesh*>(trans_sec.pMeshBase);

	//Get G values from top contacting mesh
	DBL2 Gi, Gmix;
	if (primary_top) {

		if (pMesh_pri) {

			Gi = pMesh_pri->Gi;
			Gmix = pMesh_pri->Gmix;
			pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Gi, Gi, pMesh_pri->Gmix, Gmix);
		}
		else {

			Gi = paMesh_pri->Gi;
			Gmix = paMesh_pri->Gmix;
			paMesh_pri->update_parameters_ecoarse(cell1_idx, paMesh_pri->Gi, Gi, paMesh_pri->Gmix, Gmix);
		}
	}
	else {

		if (pMesh_sec) {

			Gi = pMesh_sec->Gi;
			Gmix = pMesh_sec->Gmix;
			pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gi, Gi, pMesh_sec->Gmix, Gmix);
		}
		else {

			Gi = paMesh_sec->Gi;
			Gmix = paMesh_sec->Gmix;
			paMesh_sec->update_parameters_atposition(relpos_m1, paMesh_sec->Gi, Gi, paMesh_sec->Gmix, Gmix);
		}
	}

	DBL3 m;

	if (trans_pri.Get_STSolveType() == STSOLVE_FERROMAGNETIC || trans_pri.Get_STSolveType() == STSOLVE_FERROMAGNETIC_ATOM) {

		//F is on the primary side

		if (pMesh_pri) {

			//M value at the interface obtained using interpolation
			int idx1_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell1_idx));
			int idx2_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell2_idx));
			m = 1.5 * normalize(pMesh_pri->M[idx1_M]) - 0.5 * normalize(pMesh_pri->M[idx2_M]);
		}
		else {

			double mu_s = paMesh_pri->mu_s;
			paMesh_pri->update_parameters_ecoarse(cell1_idx, paMesh_pri->mu_s, mu_s);

			//M value at the interface obtained using interpolation
			int idx1_M = paMesh_pri->M1.position_to_cellidx(paMesh_pri->S.cellidx_to_position(cell1_idx));
			int idx2_M = paMesh_pri->M1.position_to_cellidx(paMesh_pri->S.cellidx_to_position(cell2_idx));
			m = (1.5 * paMesh_pri->M1[idx1_M] - 0.5 * paMesh_pri->M1[idx2_M]) / mu_s;
		}
	}
	else {

		//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

		if (pMesh_sec) {

			//M value at the interface obtained using interpolation
			m = 1.5 * normalize(pMesh_sec->M.weighted_average(relpos_m1, stencil)) - 0.5 * normalize(pMesh_sec->M.weighted_average(relpos_m1 + shift, stencil));
		}
		else {

			double mu_s = paMesh_sec->mu_s;
			paMesh_sec->update_parameters_atposition(relpos_m1, paMesh_sec->mu_s, mu_s);

			//M value at the interface obtained using interpolation
			m = (1.5 * paMesh_sec->M1.weighted_average(relpos_m1, stencil) - 0.5 * paMesh_sec->M1.weighted_average(relpos_m1 + shift, stencil)) / mu_s;
		}
	}

	double gR = 2 * MUB_E * Gmix.i;
	double gI = 2 * MUB_E * Gmix.j;

	DBL33 eps_m = epsilon3(m);

	return (-MUB_E * (Gi.i + Gi.j) * (m | m)) + (eps_m * (gR * eps_m + gI * ident<DBL33>()));
}

//functions to specify boundary conditions for interface conductance approach for spin : Js_N = A_N + B_N * dVs, Js_F = A_F + B_F * dVs
//VERIFIED - CORRECT
DBL33 STransport::Bfunc_F_S(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, TransportBase& trans_sec, TransportBase& trans_pri) const
{
	//the shift value comes from CMBNDInfo and is used to find positions in secondary mesh. If negative then secondary mesh is on -ve side (bottom).
	bool primary_top = !(shift >= DBL3(0));

	//primary and secondary meshes
	Mesh *pMesh_pri = nullptr, *pMesh_sec = nullptr;
	Atom_Mesh *paMesh_pri = nullptr, *paMesh_sec = nullptr;

	if (!trans_pri.pMeshBase->is_atomistic()) pMesh_pri = dynamic_cast<Mesh*>(trans_pri.pMeshBase);
	else paMesh_pri = dynamic_cast<Atom_Mesh*>(trans_pri.pMeshBase);

	if (!trans_sec.pMeshBase->is_atomistic()) pMesh_sec = dynamic_cast<Mesh*>(trans_sec.pMeshBase);
	else paMesh_sec = dynamic_cast<Atom_Mesh*>(trans_sec.pMeshBase);

	//Get G values from top contacting mesh
	DBL2 Gi;
	if (primary_top) {

		if (pMesh_pri) {

			Gi = pMesh_pri->Gi;
			pMesh_pri->update_parameters_ecoarse(cell1_idx, pMesh_pri->Gi, Gi);
		}
		else {

			Gi = paMesh_pri->Gi;
			paMesh_pri->update_parameters_ecoarse(cell1_idx, paMesh_pri->Gi, Gi);
		}
	}
	else {

		if (pMesh_sec) {

			Gi = pMesh_sec->Gi;
			pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gi, Gi);
		}
		else {

			Gi = paMesh_sec->Gi;
			paMesh_sec->update_parameters_atposition(relpos_m1, paMesh_sec->Gi, Gi);
		}
	}

	DBL3 m;

	if (trans_pri.Get_STSolveType() == STSOLVE_FERROMAGNETIC || trans_pri.Get_STSolveType() == STSOLVE_FERROMAGNETIC_ATOM) {

		//F is on the primary side
		if (pMesh_pri) {

			//M value at the interface obtained using interpolation
			int idx1_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell1_idx));
			int idx2_M = pMesh_pri->M.position_to_cellidx(pMesh_pri->S.cellidx_to_position(cell2_idx));
			m = 1.5 * normalize(pMesh_pri->M[idx1_M]) - 0.5 * normalize(pMesh_pri->M[idx2_M]);
		}
		else {

			double mu_s = paMesh_pri->mu_s;
			paMesh_pri->update_parameters_ecoarse(cell1_idx, paMesh_pri->mu_s, mu_s);

			//M value at the interface obtained using interpolation
			int idx1_M = paMesh_pri->M1.position_to_cellidx(paMesh_pri->S.cellidx_to_position(cell1_idx));
			int idx2_M = paMesh_pri->M1.position_to_cellidx(paMesh_pri->S.cellidx_to_position(cell2_idx));
			m = (1.5 * paMesh_pri->M1[idx1_M] - 0.5 * paMesh_pri->M1[idx2_M]) / mu_s;
		}
	}
	else {

		//F is on the secondary side (this function must only be used for NF interfaces, so exactly one of the contacting meshes will be magnetic)

		if (pMesh_sec) {

			//M value at the interface obtained using interpolation
			m = 1.5 * normalize(pMesh_sec->M.weighted_average(relpos_m1, stencil)) - 0.5 * normalize(pMesh_sec->M.weighted_average(relpos_m1 + shift, stencil));
		}
		else {

			double mu_s = paMesh_sec->mu_s;
			paMesh_sec->update_parameters_atposition(relpos_m1, paMesh_sec->mu_s, mu_s);

			//M value at the interface obtained using interpolation
			m = (1.5 * paMesh_sec->M1.weighted_average(relpos_m1, stencil) - 0.5 * paMesh_sec->M1.weighted_average(relpos_m1 + shift, stencil)) / mu_s;
		}
	}

	return -MUB_E * (Gi.i + Gi.j) * (m | m);
}

#endif