#include "stdafx.h"
#include "Heat.h"

#ifdef MODULE_COMPILATION_HEAT

#include "Mesh.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

double Heat::UpdateField(void)
{
	//no contribution to total energy density
	return 0.0;
}

//-------------------Calculation Methods

//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// 1-TEMPERATURE MODEL ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

//1-temperature model
void Heat::IterateHeatEquation_1TM(double dT)
{
	//FTCS:

	/////////////////////////////////////////
	// Fixed Q set (which could be zero)
	/////////////////////////////////////////

	if (!Q_equation.is_set()) {

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
#pragma omp parallel for
		for (int idx = 0; idx < pMesh->Temp.linear_size(); idx++) {

			if (!pMesh->Temp.is_not_empty(idx) || !pMesh->Temp.is_not_cmbnd(idx)) continue;

			double density = pMesh->density;
			double shc = pMesh->shc;
			double thermCond = pMesh->thermCond;
			pMesh->update_parameters_tcoarse(idx, pMesh->density, density, pMesh->shc, shc, pMesh->thermCond, thermCond);

			double cro = density * shc;
			double K = thermCond;

			//heat equation with Robin boundaries (based on Newton's law of cooling)
			heatEq_RHS[idx] = pMesh->Temp.delsq_robin(idx, K) * K / cro;

			//add Joule heating if set
			if (pMesh->E.linear_size()) {

				double joule_eff = pMesh->joule_eff;
				pMesh->update_parameters_tcoarse(idx, pMesh->joule_eff, joule_eff);

				if (IsNZ(joule_eff)) {

					DBL3 position = pMesh->Temp.cellidx_to_position(idx);

					double elC_value = pMesh->elC.weighted_average(position, pMesh->Temp.h);
					DBL3 E_value = pMesh->E.weighted_average(position, pMesh->Temp.h);

					//add Joule heating source term
					heatEq_RHS[idx] += joule_eff * (elC_value * E_value * E_value) / cro;
				}
			}

			//add heat source contribution if set
			if (IsNZ(pMesh->Q.get0())) {

				double Q = pMesh->Q;
				pMesh->update_parameters_tcoarse(idx, pMesh->Q, Q);

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
		for (int j = 0; j < pMesh->n_t.y; j++) {
			for (int k = 0; k < pMesh->n_t.z; k++) {
				for (int i = 0; i < pMesh->n_t.x; i++) {

					int idx = i + j * pMesh->n_t.x + k * pMesh->n_t.x*pMesh->n_t.y;

					if (!pMesh->Temp.is_not_empty(idx) || !pMesh->Temp.is_not_cmbnd(idx)) continue;

					double density = pMesh->density;
					double shc = pMesh->shc;
					double thermCond = pMesh->thermCond;
					pMesh->update_parameters_tcoarse(idx, pMesh->density, density, pMesh->shc, shc, pMesh->thermCond, thermCond);

					double cro = density * shc;
					double K = thermCond;

					//heat equation with Robin boundaries (based on Newton's law of cooling)
					heatEq_RHS[idx] = pMesh->Temp.delsq_robin(idx, K) * K / cro;

					//add Joule heating if set
					if (pMesh->E.linear_size()) {

						double joule_eff = pMesh->joule_eff;
						pMesh->update_parameters_tcoarse(idx, pMesh->joule_eff, joule_eff);

						if (IsNZ(joule_eff)) {

							DBL3 position = pMesh->Temp.cellidx_to_position(idx);

							double elC_value = pMesh->elC.weighted_average(position, pMesh->Temp.h);
							DBL3 E_value = pMesh->E.weighted_average(position, pMesh->Temp.h);

							//add Joule heating source term
							heatEq_RHS[idx] += joule_eff * (elC_value * E_value * E_value) / cro;
						}
					}

					//add heat source contribution
					DBL3 relpos = DBL3(i + 0.5, j + 0.5, k + 0.5) & pMesh->h_t;
					double Q = Q_equation.evaluate(relpos.x, relpos.y, relpos.z, time);

					heatEq_RHS[idx] += Q / cro;
				}
			}
		}
	}

	//2. Now use forward time to advance by dT:
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->Temp.linear_size(); idx++) {

		if (!pMesh->Temp.is_not_empty(idx) || !pMesh->Temp.is_not_cmbnd(idx)) continue;

		pMesh->Temp[idx] += dT * heatEq_RHS[idx];
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// 2-TEMPERATURE MODEL ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

//2-temperature model : itinerant electrons <-> lattice
void Heat::IterateHeatEquation_2TM(double dT)
{
	//FTCS:

	/////////////////////////////////////////
	// Fixed Q set (which could be zero)
	/////////////////////////////////////////

	if (!Q_equation.is_set()) {

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
#pragma omp parallel for
		for (int idx = 0; idx < pMesh->Temp.linear_size(); idx++) {

			if (!pMesh->Temp.is_not_empty(idx)) continue;

			double density = pMesh->density;
			double shc = pMesh->shc;
			double shc_e = pMesh->shc_e;
			double G_el = pMesh->G_e;
			double thermCond = pMesh->thermCond;
			pMesh->update_parameters_tcoarse(idx, pMesh->density, density, pMesh->shc, shc, pMesh->shc_e, shc_e, pMesh->G_e, G_el, pMesh->thermCond, thermCond);

			double cro_e = density * shc_e;
			double K = thermCond;

			//1. Itinerant Electrons Temperature

			if (pMesh->Temp.is_not_cmbnd(idx)) {

				//heat equation with Robin boundaries (based on Newton's law of cooling) and coupling to lattice
				heatEq_RHS[idx] = (pMesh->Temp.delsq_robin(idx, K) * K - G_el * (pMesh->Temp[idx] - pMesh->Temp_l[idx])) / cro_e;

				//add Joule heating if set
				if (pMesh->E.linear_size()) {

					double joule_eff = pMesh->joule_eff;
					pMesh->update_parameters_tcoarse(idx, pMesh->joule_eff, joule_eff);

					if (IsNZ(joule_eff)) {

						DBL3 position = pMesh->Temp.cellidx_to_position(idx);

						double elC_value = pMesh->elC.weighted_average(position, pMesh->Temp.h);
						DBL3 E_value = pMesh->E.weighted_average(position, pMesh->Temp.h);

						//add Joule heating source term
						heatEq_RHS[idx] += joule_eff * (elC_value * E_value * E_value) / cro_e;
					}
				}

				//add heat source contribution if set
				if (IsNZ(pMesh->Q.get0())) {

					double Q = pMesh->Q;
					pMesh->update_parameters_tcoarse(idx, pMesh->Q, Q);

					heatEq_RHS[idx] += Q / cro_e;
				}
			}

			//2. Lattice Temperature

			//lattice specific heat capacity + electron specific heat capacity gives the total specific heat capacity
			double cro_l = density * (shc - shc_e);

			pMesh->Temp_l[idx] += dT * G_el * (pMesh->Temp[idx] - pMesh->Temp_l[idx]) / cro_l;
		}
	}

	/////////////////////////////////////////
	// Q set using text equation
	/////////////////////////////////////////

	else {

		double time = pSMesh->GetStageTime();

		//1. First solve the RHS of the heat equation (centered space) : dT/dt = k del_sq T + j^2, where k = K/ c*ro , j^2 = Jc^2 / (c*ro*sigma)
		for (int j = 0; j < pMesh->n_t.y; j++) {
			for (int k = 0; k < pMesh->n_t.z; k++) {
				for (int i = 0; i < pMesh->n_t.x; i++) {

					int idx = i + j * pMesh->n_t.x + k * pMesh->n_t.x*pMesh->n_t.y;

					if (!pMesh->Temp.is_not_empty(idx)) continue;

					double density = pMesh->density;
					double shc = pMesh->shc;
					double shc_e = pMesh->shc_e;
					double G_el = pMesh->G_e;
					double thermCond = pMesh->thermCond;
					pMesh->update_parameters_tcoarse(idx, pMesh->density, density, pMesh->shc, shc, pMesh->shc_e, shc_e, pMesh->G_e, G_el, pMesh->thermCond, thermCond);

					double cro_e = density * shc_e;
					double K = thermCond;

					//1. Itinerant Electrons Temperature

					if (pMesh->Temp.is_not_cmbnd(idx)) {

						//heat equation with Robin boundaries (based on Newton's law of cooling) and coupling to lattice
						heatEq_RHS[idx] = (pMesh->Temp.delsq_robin(idx, K) * K - G_el * (pMesh->Temp[idx] - pMesh->Temp_l[idx])) / cro_e;

						//add Joule heating if set
						if (pMesh->E.linear_size()) {

							double joule_eff = pMesh->joule_eff;
							pMesh->update_parameters_tcoarse(idx, pMesh->joule_eff, joule_eff);

							if (IsNZ(joule_eff)) {

								DBL3 position = pMesh->Temp.cellidx_to_position(idx);

								double elC_value = pMesh->elC.weighted_average(position, pMesh->Temp.h);
								DBL3 E_value = pMesh->E.weighted_average(position, pMesh->Temp.h);

								//add Joule heating source term
								heatEq_RHS[idx] += joule_eff * (elC_value * E_value * E_value) / cro_e;
							}
						}

						//add heat source contribution
						DBL3 relpos = DBL3(i + 0.5, j + 0.5, k + 0.5) & pMesh->h_t;
						double Q = Q_equation.evaluate(relpos.x, relpos.y, relpos.z, time);

						heatEq_RHS[idx] += Q / cro_e;
					}

					//2. Lattice Temperature

					//lattice specific heat capacity + electron specific heat capacity gives the total specific heat capacity
					double cro_l = density * (shc - shc_e);

					pMesh->Temp_l[idx] += dT * G_el * (pMesh->Temp[idx] - pMesh->Temp_l[idx]) / cro_l;
				}
			}
		}
	}

	//2. Now use forward time to advance by dT:
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->Temp.linear_size(); idx++) {

		if (!pMesh->Temp.is_not_empty(idx) || !pMesh->Temp.is_not_cmbnd(idx)) continue;

		pMesh->Temp[idx] += dT * heatEq_RHS[idx];
	}
}

//-------------------CMBND computation methods

//CMBND values set based on continuity of temperature and heat flux

//heat flux, f(T) = -K * grad T = a + b * grad T -> a = 0, b = -K
double Heat::afunc_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	return 0.0;
}

double Heat::Heat::afunc_pri(int cell1_idx, int cell2_idx, DBL3 shift) const
{
	return 0.0;
}

double Heat::bfunc_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const
{
	double thermCond = pMesh->thermCond;
	pMesh->update_parameters_atposition(relpos_m1, pMesh->thermCond, thermCond);

	return -1.0 * thermCond;
}

double Heat::bfunc_pri(int cell1_idx, int cell2_idx) const
{
	double thermCond = pMesh->thermCond;
	pMesh->update_parameters_tcoarse(cell1_idx, pMesh->thermCond, thermCond);

	return -1.0 * thermCond;
}

//second order differential of T at cells either side of the boundary; delsq T = -Jc^2 / K * elC - Q / K - many-temperature model coupling terms / K
double Heat::diff2_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const
{
	double thermCond = pMesh->thermCond;

	if (pMesh->E.linear_size() || IsNZ(pMesh->Q.get0()) || Q_equation.is_set()) {

		pMesh->update_parameters_atposition(relpos_m1, pMesh->thermCond, thermCond);
	}
	else return 0.0;

	double value = 0.0;

	//Joule heating
	if (pMesh->E.linear_size()) {

		double joule_eff = pMesh->joule_eff;
		pMesh->update_parameters_atposition(relpos_m1, pMesh->joule_eff, joule_eff);

		if (IsNZ(joule_eff)) {

			double elCval = pMesh->elC.weighted_average(relpos_m1, stencil);
			DBL3 Eval = pMesh->E.weighted_average(relpos_m1, stencil);
			value = -joule_eff * (elCval * Eval * Eval) / thermCond;
		}
	}

	if (!Q_equation.is_set()) {

		//heat source contribution if set
		if (IsNZ(pMesh->Q.get0())) {

			double Q = pMesh->Q;
			pMesh->update_parameters_atposition(relpos_m1, pMesh->Q, Q);
			value -= Q / thermCond;
		}
	}
	else {

		double Q = Q_equation.evaluate(relpos_m1.x, relpos_m1.y, relpos_m1.z, pSMesh->GetStageTime());
		value -= Q / thermCond;
	}

	if (pMesh->Temp_l.linear_size()) {

		double G_el = pMesh->G_e;
		pMesh->update_parameters_atposition(relpos_m1, pMesh->G_e, G_el);

		value += G_el * (pMesh->Temp.weighted_average(relpos_m1, stencil) - pMesh->Temp_l.weighted_average(relpos_m1, stencil)) / thermCond;
	}

	return value;
}

double Heat::diff2_pri(int cell1_idx, DBL3 shift) const
{
	double thermCond = pMesh->thermCond;

	if (pMesh->E.linear_size() || IsNZ(pMesh->Q.get0()) || Q_equation.is_set()) {

		pMesh->update_parameters_tcoarse(cell1_idx, pMesh->thermCond, thermCond);
	}
	else return 0.0;

	double value = 0.0;

	//Joule heating
	if (pMesh->E.linear_size()) {

		double joule_eff = pMesh->joule_eff;
		pMesh->update_parameters_tcoarse(cell1_idx, pMesh->joule_eff, joule_eff);

		if (IsNZ(joule_eff)) {

			int idx1_E = pMesh->E.position_to_cellidx(pMesh->Temp.cellidx_to_position(cell1_idx));
			value = -joule_eff * (pMesh->elC[idx1_E] * pMesh->E[idx1_E] * pMesh->E[idx1_E]) / thermCond;
		}
	}

	if (!Q_equation.is_set()) {

		//heat source contribution if set
		if (IsNZ(pMesh->Q.get0())) {

			double Q = pMesh->Q;
			pMesh->update_parameters_tcoarse(cell1_idx, pMesh->Q, Q);
			value -= Q / thermCond;
		}
	}
	else {

		DBL3 relpos = pMesh->Temp.cellidx_to_position(cell1_idx);
		double Q = Q_equation.evaluate(relpos.x, relpos.y, relpos.z, pSMesh->GetStageTime());
		value -= Q / thermCond;
	}

	if (pMesh->Temp_l.linear_size()) {

		double G_el = pMesh->G_e;
		pMesh->update_parameters_tcoarse(cell1_idx, pMesh->G_e, G_el);

		value += G_el * (pMesh->Temp[cell1_idx] - pMesh->Temp_l[cell1_idx]) / thermCond;
	}

	return value;
}

#endif