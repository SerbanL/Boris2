#include "stdafx.h"
#include "Atom_DiffEqCubic.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC
#ifdef ODE_EVAL_COMPILATION_RK4

#include "Atom_Mesh_Cubic.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

//--------------------------------------------- RK4

void Atom_DifferentialEquationCubic::RunRK4_Step0_withReductions(void)
{
	bool stochastic = H_Thermal.linear_size() != 0;

	//RK4 can be used for stochastic equations - generate thermal VECs at the start
	if (H_Thermal.linear_size()) GenerateThermalField();

	if (stochastic) {

		mxh_av_reduction.new_average_reduction();

		//multiplicative conversion factor from atomic moment (units of muB) to A/m
		double conversion = MUB / paMesh->h.dim();

#pragma omp parallel for
		for (int idx = 0; idx < paMesh->n.dim(); idx++) {

			if (paMesh->M1.is_not_empty(idx)) {

				//Save current moment for later use
				sM1[idx] = paMesh->M1[idx];

				if (!paMesh->M1.is_skipcell(idx)) {

					//obtained maximum normalized torque term
					double Mnorm = paMesh->M1[idx].norm();
					mxh_av_reduction.reduce_average((paMesh->M1[idx] ^ paMesh->Heff1[idx]) / (conversion * Mnorm * Mnorm));

					//First evaluate RHS of set equation at the current time step
					sEval0[idx] = CALLFP(this, equation)(idx);

					//Now estimate moment using RK4 midle step
					paMesh->M1[idx] += sEval0[idx] * (dT / 2);
				}
			}
		}

		if (paMesh->grel.get0()) {

			mxh_reduction.max = GetMagnitude(mxh_av_reduction.average());
		}
		else {

			mxh_reduction.max = 0.0;
		}
	}
	else {

		mxh_reduction.new_minmax_reduction();

		//multiplicative conversion factor from atomic moment (units of muB) to A/m
		double conversion = MUB / paMesh->h.dim();

#pragma omp parallel for
		for (int idx = 0; idx < paMesh->n.dim(); idx++) {

			if (paMesh->M1.is_not_empty(idx)) {

				//Save current moment for later use
				sM1[idx] = paMesh->M1[idx];

				if (!paMesh->M1.is_skipcell(idx)) {

					//obtained maximum normalized torque term
					double Mnorm = paMesh->M1[idx].norm();
					double _mxh = GetMagnitude(paMesh->M1[idx] ^ paMesh->Heff1[idx]) / (conversion * Mnorm * Mnorm);
					mxh_reduction.reduce_max(_mxh);

					//First evaluate RHS of set equation at the current time step
					sEval0[idx] = CALLFP(this, equation)(idx);

					//Now estimate moment using RK4 midle step
					paMesh->M1[idx] += sEval0[idx] * (dT / 2);
				}
			}
		}

		if (paMesh->grel.get0()) {

			mxh_reduction.maximum();
		}
		else {

			mxh_reduction.max = 0.0;
		}
	}
}

void Atom_DifferentialEquationCubic::RunRK4_Step0(void)
{
	//RK4 can be used for stochastic equations - generate thermal VECs at the start
	if (H_Thermal.linear_size()) GenerateThermalField();

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			//Save current moment for later use
			sM1[idx] = paMesh->M1[idx];

			if (!paMesh->M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval0[idx] = CALLFP(this, equation)(idx);

				//Now estimate moment using RK4 midle step
				paMesh->M1[idx] += sEval0[idx] * (dT / 2);
			}
		}
	}
}

void Atom_DifferentialEquationCubic::RunRK4_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx) && !paMesh->M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval1[idx] = CALLFP(this, equation)(idx);

			//Now estimate moment using RK4 midle step
			paMesh->M1[idx] = sM1[idx] + sEval1[idx] * (dT / 2);
		}
	}
}

void Atom_DifferentialEquationCubic::RunRK4_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx) && !paMesh->M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval2[idx] = CALLFP(this, equation)(idx);

			//Now estimate moment using RK4 last step
			paMesh->M1[idx] = sM1[idx] + sEval2[idx] * dT;
		}
	}
}

void Atom_DifferentialEquationCubic::RunRK4_Step3_withReductions(void)
{
	bool stochastic = H_Thermal.linear_size() != 0;

	if (stochastic) {

		dmdt_av_reduction.new_average_reduction();

		//multiplicative conversion factor from atomic moment (units of muB) to A/m
		double conversion = MUB / paMesh->h.dim();

#pragma omp parallel for
		for (int idx = 0; idx < paMesh->n.dim(); idx++) {

			if (paMesh->M1.is_not_empty(idx)) {

				if (!paMesh->M1.is_skipcell(idx)) {

					//First evaluate RHS of set equation at the current time step
					DBL3 rhs = CALLFP(this, equation)(idx);

					//Now estimate moment using previous RK4 evaluations
					paMesh->M1[idx] = sM1[idx] + (sEval0[idx] + 2 * sEval1[idx] + 2 * sEval2[idx] + rhs) * (dT / 6);

					if (renormalize) {

						double mu_s = paMesh->mu_s;
						paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);
						paMesh->M1[idx].renormalize(mu_s);
					}

					//obtained maximum dmdt term
					double Mnorm = paMesh->M1[idx].norm();
					dmdt_av_reduction.reduce_average((paMesh->M1[idx] - sM1[idx]) / (dT * GAMMA * Mnorm * conversion * Mnorm));
				}
			}
		}

		if (paMesh->grel.get0()) {

			dmdt_reduction.max = GetMagnitude(dmdt_av_reduction.average());
		}
		else {

			dmdt_reduction.max = 0.0;
		}
	}
	else {

		dmdt_reduction.new_minmax_reduction();

		//multiplicative conversion factor from atomic moment (units of muB) to A/m
		double conversion = MUB / paMesh->h.dim();

#pragma omp parallel for
		for (int idx = 0; idx < paMesh->n.dim(); idx++) {

			if (paMesh->M1.is_not_empty(idx)) {

				if (!paMesh->M1.is_skipcell(idx)) {

					//First evaluate RHS of set equation at the current time step
					DBL3 rhs = CALLFP(this, equation)(idx);

					//Now estimate moment using previous RK4 evaluations
					paMesh->M1[idx] = sM1[idx] + (sEval0[idx] + 2 * sEval1[idx] + 2 * sEval2[idx] + rhs) * (dT / 6);

					if (renormalize) {

						double mu_s = paMesh->mu_s;
						paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);
						paMesh->M1[idx].renormalize(mu_s);
					}

					//obtained maximum dmdt term
					double Mnorm = paMesh->M1[idx].norm();
					double _dmdt = GetMagnitude(paMesh->M1[idx] - sM1[idx]) / (dT * GAMMA * Mnorm * conversion * Mnorm);
					dmdt_reduction.reduce_max(_dmdt);
				}
			}
		}

		if (paMesh->grel.get0()) {

			dmdt_reduction.maximum();
		}
		else {

			dmdt_reduction.max = 0.0;
		}
	}
}

void Atom_DifferentialEquationCubic::RunRK4_Step3(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			if (!paMesh->M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//Now estimate moment using previous RK4 evaluations
				paMesh->M1[idx] = sM1[idx] + (sEval0[idx] + 2 * sEval1[idx] + 2 * sEval2[idx] + rhs) * (dT / 6);

				if (renormalize) {

					double mu_s = paMesh->mu_s;
					paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);
					paMesh->M1[idx].renormalize(mu_s);
				}
			}
		}
	}
}

#endif
#endif