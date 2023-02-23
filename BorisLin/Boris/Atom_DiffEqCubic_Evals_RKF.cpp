#include "stdafx.h"
#include "Atom_DiffEqCubic.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC
#ifdef ODE_EVAL_COMPILATION_RKF45

#include "Atom_Mesh_Cubic.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

//--------------------------------------------- RUNGE KUTTA FEHLBERG (4th order solution, 5th order error)

void Atom_DifferentialEquationCubic::RunRKF45_Step0_withReductions(void)
{
	mxh_reduction.new_minmax_reduction();

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	double conversion = MUB / paMesh->h.dim();

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			//Save current moment
			sM1[idx] = paMesh->M1[idx];

			if (!paMesh->M1.is_skipcell(idx)) {

				//obtain maximum normalized torque term
				double Mnorm = paMesh->M1[idx].norm();
				double _mxh = GetMagnitude(paMesh->M1[idx] ^ paMesh->Heff1[idx]) / (conversion * Mnorm * Mnorm);
				mxh_reduction.reduce_max(_mxh);

				//First evaluate RHS of set equation at the current time step
				sEval0[idx] = CALLFP(this, equation)(idx);

				//Now estimate moment using RKF first step
				paMesh->M1[idx] += sEval0[idx] * (2 * dT / 9);
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

void Atom_DifferentialEquationCubic::RunRKF45_Step0(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			//Save current moment
			sM1[idx] = paMesh->M1[idx];

			if (!paMesh->M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval0[idx] = CALLFP(this, equation)(idx);

				//Now estimate moment using RKF first step
				paMesh->M1[idx] += sEval0[idx] * (2 * dT / 9);
			}
		}
	}
}

void Atom_DifferentialEquationCubic::RunRKF45_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx) && !paMesh->M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval1[idx] = CALLFP(this, equation)(idx);

			//Now estimate moment using RKF midle step 1
			paMesh->M1[idx] = sM1[idx] + (sEval0[idx] / 12 + sEval1[idx] / 4) * dT;
		}
	}
}

void Atom_DifferentialEquationCubic::RunRKF45_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx) && !paMesh->M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval2[idx] = CALLFP(this, equation)(idx);

			//Now estimate moment using RKF midle step 2
			paMesh->M1[idx] = sM1[idx] + (69 * sEval0[idx] / 128 - 243 * sEval1[idx] / 128 + 135 * sEval2[idx] / 64) * dT;
		}
	}
}

void Atom_DifferentialEquationCubic::RunRKF45_Step3(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx) && !paMesh->M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval3[idx] = CALLFP(this, equation)(idx);

			//Now estimate moment using RKF midle step 3
			paMesh->M1[idx] = sM1[idx] + (-17 * sEval0[idx] / 12 + 27 * sEval1[idx] / 4 - 27 * sEval2[idx] / 5 + 16 * sEval3[idx] / 15) * dT;
		}
	}
}

void Atom_DifferentialEquationCubic::RunRKF45_Step4(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx) && !paMesh->M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval4[idx] = CALLFP(this, equation)(idx);

			//Now estimate moment using RKF midle step 4
			paMesh->M1[idx] = sM1[idx] + (65 * sEval0[idx] / 432 - 5 * sEval1[idx] / 16 + 13 * sEval2[idx] / 16 + 4 * sEval3[idx] / 27 + 5 * sEval4[idx] / 144) * dT;
		}
	}
}

void Atom_DifferentialEquationCubic::RunRKF45_Step5_withReductions(void)
{
	dmdt_reduction.new_minmax_reduction();
	lte_reduction.new_minmax_reduction();

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	double conversion = MUB / paMesh->h.dim();

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			if (!paMesh->M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//4th order evaluation
				paMesh->M1[idx] = sM1[idx] + (sEval0[idx] / 9 + 9 * sEval2[idx] / 20 + 16 * sEval3[idx] / 45 + sEval4[idx] / 12) * dT;

				//5th order evaluation
				DBL3 prediction = sM1[idx] + (47 * sEval0[idx] / 450 + 12 * sEval2[idx] / 25 + 32 * sEval3[idx] / 225 + 1 * sEval4[idx] / 30 + 6 * rhs / 25) * dT;

				if (renormalize) {

					double mu_s = paMesh->mu_s;
					paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);
					paMesh->M1[idx].renormalize(mu_s);
				}

				//obtained maximum dmdt term
				double Mnorm = paMesh->M1[idx].norm();
				double _dmdt = GetMagnitude(paMesh->M1[idx] - sM1[idx]) / (dT * GAMMA * Mnorm * conversion * Mnorm);
				dmdt_reduction.reduce_max(_dmdt);

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(paMesh->M1[idx] - prediction) / paMesh->M1[idx].norm();
				lte_reduction.reduce_max(_lte);
			}
		}
	}

	lte_reduction.maximum();

	if (paMesh->grel.get0()) {

		dmdt_reduction.maximum();
	}
	else {

		dmdt_reduction.max = 0.0;
	}
}

void Atom_DifferentialEquationCubic::RunRKF45_Step5(void)
{
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			if (!paMesh->M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//4th order evaluation
				paMesh->M1[idx] = sM1[idx] + (sEval0[idx] / 9 + 9 * sEval2[idx] / 20 + 16 * sEval3[idx] / 45 + sEval4[idx] / 12) * dT;

				//5th order evaluation
				DBL3 prediction = sM1[idx] + (47 * sEval0[idx] / 450 + 12 * sEval2[idx] / 25 + 32 * sEval3[idx] / 225 + 1 * sEval4[idx] / 30 + 6 * rhs / 25) * dT;

				if (renormalize) {

					double mu_s = paMesh->mu_s;
					paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);
					paMesh->M1[idx].renormalize(mu_s);
				}

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(paMesh->M1[idx] - prediction) / paMesh->M1[idx].norm();
				lte_reduction.reduce_max(_lte);
			}
		}
	}

	lte_reduction.maximum();
}

#endif
#endif