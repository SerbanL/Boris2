#include "stdafx.h"
#include "Transport.h"

#ifdef MODULE_COMPILATION_TRANSPORT

#include "Mesh.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "TransportCUDA.h"
#endif

//-------------------Display Calculation Methods

//prepare displayVEC ready for calculation of display quantity
bool Transport::PrepareDisplayVEC(DBL3 cellsize)
{
	if (pSMesh->SolveSpinCurrent() && pMesh->EComputation_Enabled()) {

		//make sure memory is allocated to the correct size
		displayVEC.assign(cellsize, pMesh->meshRect, DBL3(0.0));

		return true;
	}
	else displayVEC.clear();

	return false;
}

//return x, y, or z component of spin current (component = 0, 1, or 2)
//VERIFIED - CORRECT
VEC<DBL3>& Transport::GetSpinCurrent(int component)
{
	//make sure memory is allocated to the correct size
	if (!PrepareDisplayVEC(pMesh->h_e)) return displayVEC;

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		cu_obj<cuVEC<cuReal3>>& cudisplayVEC = GetSpinCurrentCUDA(component);
		cudisplayVEC()->copy_to_cpuvec(displayVEC);

		return displayVEC;
	}
#endif

	//compute spin current and store result in displayVEC depending on required component

	bool cpump_enabled = IsNZ(pMesh->cpump_eff.get0());
	bool the_enabled = IsNZ(pMesh->the_eff.get0());

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->S.linear_size(); idx++) {

		DBL33 Js = DBL33();

		if (pMesh->S.is_not_empty(idx)) {

			if (stsolve == STSOLVE_FERROMAGNETIC) {

				//magnetic mesh terms

				double Ms = pMesh->Ms;
				double P = pMesh->P;
				double De = pMesh->De;
				pMesh->update_parameters_ecoarse(idx, pMesh->Ms, Ms, pMesh->P, P, pMesh->De, De);

				//1. drift
				int idx_M = pMesh->M.position_to_cellidx(pMesh->S.cellidx_to_position(idx));

				DBL3 Mval = pMesh->M[idx_M];
				DBL33 grad_S = pMesh->S.grad_neu(idx);

				Js = (pMesh->E[idx] | Mval) * (P * pMesh->elC[idx] / Ms) * (-MUB_E);

				//2. diffusion with homogeneous Neumann boundary condition
				Js -= grad_S * De;

				//3. charge pumping
				//4. topological Hall effect

				if (component != 2 && (cpump_enabled || the_enabled)) {

					DBL33 grad_m = pMesh->M.grad_neu(idx_M) / Ms;

					//topological Hall effect contribution
					if (the_enabled) {

						double n_density = pMesh->n_density;
						pMesh->update_parameters_ecoarse(idx, pMesh->n_density, n_density);

						DBL3 B = (grad_m.x ^ grad_m.y);
						Js += pMesh->the_eff.get0() * (HBAR_E * MUB_E * pMesh->elC[idx] * pMesh->elC[idx] / (ECHARGE * n_density)) * DBL33(-pMesh->E[idx].y * B, pMesh->E[idx].x * B, DBL3());
					}

					//charge pumping contribution
					if (cpump_enabled) {

						//value a1
						DBL3 dm_dt = dM_dt[idx_M] / Ms;
						Js += pMesh->cpump_eff.get0() * (HBAR_E * MUB_E * pMesh->elC[idx] / 2) * DBL33(dm_dt ^ grad_m.x, dm_dt ^ grad_m.y, DBL3());
					}
				}
			}
			else if (stsolve != STSOLVE_NONE) {

				//non-magnetic mesh terms

				double De = pMesh->De;
				double SHA = pMesh->SHA;
				pMesh->update_parameters_ecoarse(idx, pMesh->De, De, pMesh->SHA, SHA);

				//1. SHE contribution
				Js = epsilon3(pMesh->E[idx]) * SHA * pMesh->elC[idx] * MUB_E;

				//2. diffusion with non-homogeneous Neumann boundary condition
				Js -= pMesh->S.grad_nneu(idx, epsilon3(pMesh->E[idx]) * (SHA * pMesh->elC[idx] * MUB_E / De)) * De;
			}
		}

		switch (component) {

		case 0:
			displayVEC[idx] = Js.x;
			break;
		case 1:
			displayVEC[idx] = Js.y;
			break;
		case 2:
			displayVEC[idx] = Js.z;
			break;
		}
	}

	return displayVEC;
}

DBL3 Transport::GetAverageSpinCurrent(int component, Rect rectangle, std::vector<MeshShape> shapes)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		//update displayVEC with spin current in TransportCUDA
		GetSpinCurrentCUDA(component);

		//average spin current in displayVEC in TransportCUDA
		if (pSMesh->SolveSpinCurrent()) return cuReal3(pTransportCUDA->displayVEC()->average_nonempty(pMesh->n_e.dim(), rectangle));
		else return cuReal3(0.0);
	}
#endif

	//update displayVEC with spin current
	GetSpinCurrent(component);

	//average spin current in displayVEC
	if (!shapes.size()) return displayVEC.average_nonempty_omp(rectangle);
	else return displayVEC.shape_getaverage(shapes);
}

//Get average bulk spin torque in given rectangle - calculate it first
DBL3 Transport::GetAverageSpinTorque(Rect rectangle, std::vector<MeshShape> shapes)
{
	GetSpinTorque();
	if (displayVEC.linear_size()) {

		if (!shapes.size()) return displayVEC.average_nonempty_omp(rectangle);
		else return displayVEC.shape_getaverage(shapes);
	}

	return DBL3();
}

//Get average interfacial spin torque in given rectangle - must have been calculated already in displayVEC through the supermesh transport module
DBL3 Transport::GetAverageInterfacialSpinTorque(Rect rectangle, std::vector<MeshShape> shapes)
{ 
	if (displayVEC.linear_size()) {

		if (!shapes.size()) return displayVEC.average_nonempty_omp(rectangle);
		else return displayVEC.shape_getaverage(shapes);
	}

	return DBL3();
}

//return spin torque computed from spin accumulation
//VERIFIED - CORRECT
VEC<DBL3>& Transport::GetSpinTorque(void)
{
	//make sure memory is allocated to the correct size
	if (!PrepareDisplayVEC(pMesh->h)) return displayVEC;

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		cu_obj<cuVEC<cuReal3>>& cudisplayVEC = GetSpinTorqueCUDA();
		cudisplayVEC()->copy_to_cpuvec(displayVEC);

		return displayVEC;
	}
#endif

	if (stsolve != STSOLVE_FERROMAGNETIC) return displayVEC;

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->M.linear_size(); idx++) {

		if (pMesh->M.is_empty(idx)) {

			displayVEC[idx] = DBL3();
			continue;
		}

		double Ms = pMesh->Ms;
		double De = pMesh->De;
		double ts_eff = pMesh->ts_eff;
		double l_ex = pMesh->l_ex;
		double l_ph = pMesh->l_ph;
		pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms, pMesh->De, De, pMesh->ts_eff, ts_eff, pMesh->l_ex, l_ex, pMesh->l_ph, l_ph);

		//average S in the magnetization cell with index idx
		DBL3 Sav = pMesh->S.weighted_average(pMesh->M.cellidx_to_position(idx), pMesh->h);
		DBL3 M = pMesh->M[idx];

		displayVEC[idx] = ts_eff * ((Sav ^ M) * De / (Ms * l_ex * l_ex) + (M ^ (Sav ^ M)) * De / (Ms * Ms * l_ph * l_ph));
	}

	return displayVEC;
}

#if COMPILECUDA == 1
cu_obj<cuVEC<cuReal3>>& Transport::GetSpinCurrentCUDA(int component)
{
	return pTransportCUDA->GetSpinCurrent(component);
}

cu_obj<cuVEC<cuReal3>>& Transport::GetSpinTorqueCUDA(void)
{
	return pTransportCUDA->GetSpinTorque();
}
#endif

//Calculate the interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set), accumulating result in displayVEC
//VERIFIED - CORRECT
void Transport::CalculateDisplaySAInterfaceTorque(Transport* ptrans_sec, CMBNDInfo& contact)
{
	//the top contacting mesh sets G values
	bool isGInterface_Enabled = ((contact.IsPrimaryTop() && pMesh->GInterface_Enabled()) || (!contact.IsPrimaryTop() && ptrans_sec->GInterface_Enabled()));

	if (isGInterface_Enabled && stsolve == STSOLVE_FERROMAGNETIC && ptrans_sec->Get_STSolveType() == STSOLVE_NORMALMETAL) {

		//interface conductance method with F being the primary mesh (N-F contact): calculate and set spin torque

		//convert the cells box from S mesh to M mesh
		INT3 mbox_start = pMesh->M.cellidx_from_position(pMesh->S.cellidx_to_position(contact.cells_box.s) + pMesh->meshRect.s);
		INT3 mbox_end = pMesh->M.cellidx_from_position(pMesh->S.cellidx_to_position(contact.cells_box.e - INT3(1)) + pMesh->meshRect.s) + INT3(1);

		if ((mbox_end.i - mbox_start.i) == 0) mbox_end.i = mbox_start.i + 1;
		if ((mbox_end.j - mbox_start.j) == 0) mbox_end.j = mbox_start.j + 1;
		if ((mbox_end.k - mbox_start.k) == 0) mbox_end.k = mbox_start.k + 1;

		INT3 box_sizes = mbox_end - mbox_start;

		//the cellsize perpendicular to the contact (in the M mesh)
		double dh = (DBL3(contact.cell_shift) & pMesh->h).norm();

		Mesh* pMesh_sec = ptrans_sec->pMesh;

		//primary cells in this contact
#pragma omp parallel for
		for (int box_idx = 0; box_idx < box_sizes.dim(); box_idx++) {

			int i = (box_idx % box_sizes.x) + mbox_start.i;
			int j = ((box_idx / box_sizes.x) % box_sizes.y) + mbox_start.j;
			int k = (box_idx / (box_sizes.x * box_sizes.y)) + mbox_start.k;

			//index of magnetic cell 1
			int mcell1_idx = i + j * pMesh->n.x + k * pMesh->n.x*pMesh->n.y;

			if (pMesh->M.is_empty(mcell1_idx)) continue;

			double Ms = pMesh->Ms;
			double tsi_eff = pMesh->tsi_eff;
			pMesh->update_parameters_mcoarse(mcell1_idx, pMesh->Ms, Ms, pMesh->tsi_eff, tsi_eff);

			//position at interface relative to primary mesh
			DBL3 mhshift_primary = contact.hshift_primary.normalized() & pMesh->h;
			DBL3 relpos_interf = ((DBL3(i, j, k) + DBL3(0.5)) & pMesh->h) + mhshift_primary / 2;

			DBL3 relpos_1 = relpos_interf - contact.hshift_primary / 2;

			DBL3 relpos_m1 = pMesh->meshRect.s - pMesh_sec->meshRect.s + relpos_interf + contact.hshift_secondary / 2;

			DBL3 stencil = pMesh->h - mod(mhshift_primary) + mod(contact.hshift_secondary);

			//S values
			DBL3 S_1 = pMesh->S.weighted_average(relpos_1, stencil);
			DBL3 S_2 = pMesh->S.weighted_average(relpos_1 - contact.hshift_primary, stencil);
			DBL3 S_m1 = pMesh_sec->S.weighted_average(relpos_m1, stencil);
			DBL3 S_m2 = pMesh_sec->S.weighted_average(relpos_m1 + contact.hshift_secondary, stencil);

			//c values
			double c_m1 = ptrans_sec->cfunc_sec(relpos_m1, stencil);
			double c_m2 = ptrans_sec->cfunc_sec(relpos_m1 + contact.hshift_secondary, stencil);
			double c_1 = cfunc_sec(relpos_1, stencil);
			double c_2 = cfunc_sec(relpos_1 - contact.hshift_primary, stencil);

			//Calculate S drop at the interface
			DBL3 Vs_F = 1.5 * c_1 * S_1 - 0.5 * c_2 * S_2;
			DBL3 Vs_N = 1.5 * c_m1 * S_m1 - 0.5 * c_m2 * S_m2;
			DBL3 dVs = Vs_F - Vs_N;

			//Get G values from top contacting mesh
			DBL2 Gmix;
			if (contact.IsPrimaryTop()) {

				Gmix = pMesh->Gmix;
				pMesh->update_parameters_mcoarse(mcell1_idx, pMesh->Gmix, Gmix);
			}
			else {

				Gmix = pMesh_sec->Gmix;
				pMesh_sec->update_parameters_atposition(relpos_m1, pMesh_sec->Gmix, Gmix);
			}

			double gI = (2.0 * GMUB_2E / dh) * Gmix.j / Ms;
			double gR = (2.0 * GMUB_2E / dh) * Gmix.i / Ms;

			displayVEC[mcell1_idx] += tsi_eff * (gI * (pMesh->M[mcell1_idx] ^ dVs) + gR * (pMesh->M[mcell1_idx] ^ (pMesh->M[mcell1_idx] ^ dVs)) / Ms);
		}
	}
}

#endif