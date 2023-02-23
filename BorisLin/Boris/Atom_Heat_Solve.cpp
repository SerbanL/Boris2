#include "stdafx.h"
#include "Atom_Heat.h"

#if defined(MODULE_COMPILATION_HEAT) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

double Atom_Heat::UpdateField(void)
{
	//no contribution to total energy density
	return 0.0;
}

//-------------------Calculation Methods

//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// 1-TEMPERATURE MODEL ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

//1-temperature model
void Atom_Heat::IterateHeatEquation_1TM(double dT)
{
	//FTCS:

	/////////////////////////////////////////
	// Fixed Q set (which could be zero)
	/////////////////////////////////////////

	if (!Q_equation.is_set()) {

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
#pragma omp parallel for
		for (int idx = 0; idx < paMesh->Temp.linear_size(); idx++) {

			if (!paMesh->Temp.is_not_empty(idx) || !paMesh->Temp.is_not_cmbnd(idx)) continue;

			double density = paMesh->density;
			double shc = paMesh->shc;
			double thermCond = paMesh->thermCond;
			paMesh->update_parameters_tcoarse(idx, paMesh->density, density, paMesh->shc, shc, paMesh->thermCond, thermCond);

			double cro = density * shc;
			double K = thermCond;

			//heat equation with Robin boundaries (based on Newton's law of cooling)
			heatEq_RHS[idx] = paMesh->Temp.delsq_robin(idx, K) * K / cro;

			//add Joule heating if set
			if (paMesh->E.linear_size()) {

				double joule_eff = paMesh->joule_eff;
				paMesh->update_parameters_tcoarse(idx, paMesh->joule_eff, joule_eff);

				if (IsNZ(joule_eff)) {

					DBL3 position = paMesh->Temp.cellidx_to_position(idx);

					double elC_value = paMesh->elC.weighted_average(position, paMesh->Temp.h);
					DBL3 E_value = paMesh->E.weighted_average(position, paMesh->Temp.h);

					//add Joule heating source term
					heatEq_RHS[idx] += joule_eff * (elC_value * E_value * E_value) / cro;
				}
			}

			//add heat source contribution if set
			if (IsNZ(paMesh->Q.get0())) {

				double Q = paMesh->Q;
				paMesh->update_parameters_tcoarse(idx, paMesh->Q, Q);

				heatEq_RHS[idx] += Q / cro;
			}
		}
	}

	/////////////////////////////////////////
	// Q set using text equation
	/////////////////////////////////////////

	else {

		double time = pSMesh->GetStageTime();

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
		for (int j = 0; j < paMesh->n_t.y; j++) {
			for (int k = 0; k < paMesh->n_t.z; k++) {
				for (int i = 0; i < paMesh->n_t.x; i++) {

					int idx = i + j * paMesh->n_t.x + k * paMesh->n_t.x*paMesh->n_t.y;

					if (!paMesh->Temp.is_not_empty(idx) || !paMesh->Temp.is_not_cmbnd(idx)) continue;

					double density = paMesh->density;
					double shc = paMesh->shc;
					double thermCond = paMesh->thermCond;
					paMesh->update_parameters_tcoarse(idx, paMesh->density, density, paMesh->shc, shc, paMesh->thermCond, thermCond);

					double cro = density * shc;
					double K = thermCond;

					//heat equation with Robin boundaries (based on Newton's law of cooling)
					heatEq_RHS[idx] = paMesh->Temp.delsq_robin(idx, K) * K / cro;

					//add Joule heating if set
					if (paMesh->E.linear_size()) {

						double joule_eff = paMesh->joule_eff;
						paMesh->update_parameters_tcoarse(idx, paMesh->joule_eff, joule_eff);

						if (IsNZ(joule_eff)) {

							DBL3 position = paMesh->Temp.cellidx_to_position(idx);

							double elC_value = paMesh->elC.weighted_average(position, paMesh->Temp.h);
							DBL3 E_value = paMesh->E.weighted_average(position, paMesh->Temp.h);

							//add Joule heating source term
							heatEq_RHS[idx] += joule_eff * (elC_value * E_value * E_value) / cro;
						}
					}

					//add heat source contribution
					DBL3 relpos = DBL3(i + 0.5, j + 0.5, k + 0.5) & paMesh->h_t;
					double Q = Q_equation.evaluate(relpos.x, relpos.y, relpos.z, time);

					heatEq_RHS[idx] += Q / cro;
				}
			}
		}
	}

	//2. Now use forward time to advance by dT:
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->Temp.linear_size(); idx++) {

		if (!paMesh->Temp.is_not_empty(idx) || !paMesh->Temp.is_not_cmbnd(idx)) continue;

		paMesh->Temp[idx] += dT * heatEq_RHS[idx];
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// 2-TEMPERATURE MODEL ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

//2-temperature model : itinerant electrons <-> lattice
void Atom_Heat::IterateHeatEquation_2TM(double dT)
{
	//FTCS:

	/////////////////////////////////////////
	// Fixed Q set (which could be zero)
	/////////////////////////////////////////

	if (!Q_equation.is_set()) {

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
#pragma omp parallel for
		for (int idx = 0; idx < paMesh->Temp.linear_size(); idx++) {

			if (!paMesh->Temp.is_not_empty(idx)) continue;

			double density = paMesh->density;
			double shc = paMesh->shc;
			double shc_e = paMesh->shc_e;
			double G_el = paMesh->G_e;
			double thermCond = paMesh->thermCond;
			paMesh->update_parameters_tcoarse(idx, paMesh->density, density, paMesh->shc, shc, paMesh->shc_e, shc_e, paMesh->G_e, G_el, paMesh->thermCond, thermCond);

			double cro_e = density * shc_e;
			double K = thermCond;

			//1. Itinerant Electrons Temperature

			if (paMesh->Temp.is_not_cmbnd(idx)) {

				//heat equation with Robin boundaries (based on Newton's law of cooling) and coupling to lattice
				heatEq_RHS[idx] = (paMesh->Temp.delsq_robin(idx, K) * K - G_el * (paMesh->Temp[idx] - paMesh->Temp_l[idx])) / cro_e;

				//add Joule heating if set
				if (paMesh->E.linear_size()) {

					double joule_eff = paMesh->joule_eff;
					paMesh->update_parameters_tcoarse(idx, paMesh->joule_eff, joule_eff);

					if (IsNZ(joule_eff)) {

						DBL3 position = paMesh->Temp.cellidx_to_position(idx);

						double elC_value = paMesh->elC.weighted_average(position, paMesh->Temp.h);
						DBL3 E_value = paMesh->E.weighted_average(position, paMesh->Temp.h);

						//add Joule heating source term
						heatEq_RHS[idx] += joule_eff * (elC_value * E_value * E_value) / cro_e;
					}
				}

				//add heat source contribution if set
				if (IsNZ(paMesh->Q.get0())) {

					double Q = paMesh->Q;
					paMesh->update_parameters_tcoarse(idx, paMesh->Q, Q);

					heatEq_RHS[idx] += Q / cro_e;
				}
			}

			//2. Lattice Temperature

			//lattice specific heat capacity + electron specific heat capacity gives the total specific heat capacity
			double cro_l = density * (shc - shc_e);

			paMesh->Temp_l[idx] += dT * G_el * (paMesh->Temp[idx] - paMesh->Temp_l[idx]) / cro_l;
		}
	}

	/////////////////////////////////////////
	// Q set using text equation
	/////////////////////////////////////////

	else {

		double time = pSMesh->GetStageTime();

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
		for (int j = 0; j < paMesh->n_t.y; j++) {
			for (int k = 0; k < paMesh->n_t.z; k++) {
				for (int i = 0; i < paMesh->n_t.x; i++) {

					int idx = i + j * paMesh->n_t.x + k * paMesh->n_t.x*paMesh->n_t.y;

					if (!paMesh->Temp.is_not_empty(idx)) continue;

					double density = paMesh->density;
					double shc = paMesh->shc;
					double shc_e = paMesh->shc_e;
					double G_el = paMesh->G_e;
					double thermCond = paMesh->thermCond;
					paMesh->update_parameters_tcoarse(idx, paMesh->density, density, paMesh->shc, shc, paMesh->shc_e, shc_e, paMesh->G_e, G_el, paMesh->thermCond, thermCond);

					double cro_e = density * shc_e;
					double K = thermCond;

					//1. Itinerant Electrons Temperature

					if (paMesh->Temp.is_not_cmbnd(idx)) {

						//heat equation with Robin boundaries (based on Newton's law of cooling) and coupling to lattice
						heatEq_RHS[idx] = (paMesh->Temp.delsq_robin(idx, K) * K - G_el * (paMesh->Temp[idx] - paMesh->Temp_l[idx])) / cro_e;

						//add Joule heating if set
						if (paMesh->E.linear_size()) {

							double joule_eff = paMesh->joule_eff;
							paMesh->update_parameters_tcoarse(idx, paMesh->joule_eff, joule_eff);

							if (IsNZ(joule_eff)) {

								DBL3 position = paMesh->Temp.cellidx_to_position(idx);

								double elC_value = paMesh->elC.weighted_average(position, paMesh->Temp.h);
								DBL3 E_value = paMesh->E.weighted_average(position, paMesh->Temp.h);

								//add Joule heating source term
								heatEq_RHS[idx] += joule_eff * (elC_value * E_value * E_value) / cro_e;
							}
						}

						//add heat source contribution
						DBL3 relpos = DBL3(i + 0.5, j + 0.5, k + 0.5) & paMesh->h_t;
						double Q = Q_equation.evaluate(relpos.x, relpos.y, relpos.z, time);

						heatEq_RHS[idx] += Q / cro_e;
					}

					//2. Lattice Temperature

					//lattice specific heat capacity + electron specific heat capacity gives the total specific heat capacity
					double cro_l = density * (shc - shc_e);

					paMesh->Temp_l[idx] += dT * G_el * (paMesh->Temp[idx] - paMesh->Temp_l[idx]) / cro_l;
				}
			}
		}
	}

	//2. Now use forward time to advance by dT:
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->Temp.linear_size(); idx++) {

		if (!paMesh->Temp.is_not_empty(idx) || !paMesh->Temp.is_not_cmbnd(idx)) continue;

		paMesh->Temp[idx] += dT * heatEq_RHS[idx];
	}
}

//-------------------CMBND computation methods

//CMBND values set based on continuity of temperature and heat flux

//heat flux, f(T) = -K * grad T = a + b * grad T -> a = 0, b = -K
double Atom_Heat::afunc_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	return 0.0;
}

double Atom_Heat::afunc_pri(int cell1_idx, int cell2_idx, DBL3 shift) const
{
	return 0.0;
}

double Atom_Heat::bfunc_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	double thermCond = paMesh->thermCond;
	paMesh->update_parameters_atposition(relpos_m1, paMesh->thermCond, thermCond);

	return -1.0 * thermCond;
}

double Atom_Heat::bfunc_pri(int cell1_idx, int cell2_idx) const
{
	double thermCond = paMesh->thermCond;
	paMesh->update_parameters_tcoarse(cell1_idx, paMesh->thermCond, thermCond);

	return -1.0 * thermCond;
}

//second order differential of T at cells either side of the boundary; delsq T = -Jc^2 / K * elC - Q / K - many-temperature model coupling terms / K
double Atom_Heat::diff2_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const
{
	double thermCond = paMesh->thermCond;

	if (paMesh->E.linear_size() || IsNZ(paMesh->Q.get0())) {

		paMesh->update_parameters_atposition(relpos_m1, paMesh->thermCond, thermCond);
	}
	else return 0.0;

	double value = 0.0;

	//Joule heating
	if (paMesh->E.linear_size()) {

		double joule_eff = paMesh->joule_eff;
		paMesh->update_parameters_atposition(relpos_m1, paMesh->joule_eff, joule_eff);

		if (IsNZ(joule_eff)) {

			double elCval = paMesh->elC.weighted_average(relpos_m1, stencil);
			DBL3 Eval = paMesh->E.weighted_average(relpos_m1, stencil);
			value = -joule_eff * (elCval * Eval * Eval) / thermCond;
		}
	}

	if (!Q_equation.is_set()) {

		//heat source contribution if set
		if (IsNZ(paMesh->Q.get0())) {

			double Q = paMesh->Q;
			paMesh->update_parameters_atposition(relpos_m1, paMesh->Q, Q);

			value -= Q / thermCond;
		}
	}
	else {

		double Q = Q_equation.evaluate(relpos_m1.x, relpos_m1.y, relpos_m1.z, pSMesh->GetStageTime());

		value -= Q / thermCond;
	}

	if (paMesh->Temp_l.linear_size()) {

		double G_el = paMesh->G_e;
		paMesh->update_parameters_atposition(relpos_m1, paMesh->G_e, G_el);

		value += G_el * (paMesh->Temp.weighted_average(relpos_m1, stencil) - paMesh->Temp_l.weighted_average(relpos_m1, stencil)) / thermCond;
	}

	return value;
}

double Atom_Heat::diff2_pri(int cell1_idx, DBL3 shift) const
{
	double thermCond = paMesh->thermCond;

	if (paMesh->E.linear_size() || IsNZ(paMesh->Q.get0())) {

		paMesh->update_parameters_tcoarse(cell1_idx, paMesh->thermCond, thermCond);
	}
	else return 0.0;

	double value = 0.0;

	//Joule heating
	if (paMesh->E.linear_size()) {

		double joule_eff = paMesh->joule_eff;
		paMesh->update_parameters_tcoarse(cell1_idx, paMesh->joule_eff, joule_eff);

		if (IsNZ(joule_eff)) {

			int idx1_E = paMesh->E.position_to_cellidx(paMesh->Temp.cellidx_to_position(cell1_idx));
			value = -joule_eff * (paMesh->elC[idx1_E] * paMesh->E[idx1_E] * paMesh->E[idx1_E]) / thermCond;
		}
	}

	if (!Q_equation.is_set()) {

		//heat source contribution if set
		if (IsNZ(paMesh->Q.get0())) {

			double Q = paMesh->Q;
			paMesh->update_parameters_tcoarse(cell1_idx, paMesh->Q, Q);

			value -= Q / thermCond;
		}
	}
	else {

		DBL3 relpos = paMesh->Temp.cellidx_to_position(cell1_idx);
		double Q = Q_equation.evaluate(relpos.x, relpos.y, relpos.z, pSMesh->GetStageTime());

		value -= Q / thermCond;
	}

	if (paMesh->Temp_l.linear_size()) {

		double G_el = paMesh->G_e;
		paMesh->update_parameters_tcoarse(cell1_idx, paMesh->G_e, G_el);

		value += G_el * (paMesh->Temp[cell1_idx] - paMesh->Temp_l[cell1_idx]) / thermCond;
	}

	return value;
}

#endif