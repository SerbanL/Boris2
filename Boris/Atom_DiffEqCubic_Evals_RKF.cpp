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
				paMesh->M1[idx] += sEval0[idx] * (dT / 4);
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
				paMesh->M1[idx] += sEval0[idx] * (dT / 4);
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
			paMesh->M1[idx] = sM1[idx] + (3 * sEval0[idx] + 9 * sEval1[idx]) * dT / 32;
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
			paMesh->M1[idx] = sM1[idx] + (1932 * sEval0[idx] - 7200 * sEval1[idx] + 7296 * sEval2[idx]) * dT / 2197;
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
			paMesh->M1[idx] = sM1[idx] + (439 * sEval0[idx] / 216 - 8 * sEval1[idx] + 3680 * sEval2[idx] / 513 - 845 * sEval3[idx] / 4104) * dT;
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
			paMesh->M1[idx] = sM1[idx] + (-8 * sEval0[idx] / 27 + 2 * sEval1[idx] - 3544 * sEval2[idx] / 2565 + 1859 * sEval3[idx] / 4104 - 11 * sEval4[idx] / 40) * dT;
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
				DBL3 prediction = sM1[idx] + (25 * sEval0[idx] / 216 + 1408 * sEval2[idx] / 2565 + 2197 * sEval3[idx] / 4101 - sEval4[idx] / 5) * dT;

				//5th order evaluation -> keep this as the new value, not the 4th order; relaxation doesn't work well the other way around.
				paMesh->M1[idx] = sM1[idx] + (16 * sEval0[idx] / 135 + 6656 * sEval2[idx] / 12825 + 28561 * sEval3[idx] / 56430 - 9 * sEval4[idx] / 50 + 2 * rhs / 55) * dT;

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
				DBL3 prediction = sM1[idx] + (25 * sEval0[idx] / 216 + 1408 * sEval2[idx] / 2565 + 2197 * sEval3[idx] / 4101 - sEval4[idx] / 5) * dT;

				//5th order evaluation -> keep this as the new value, not the 4th order; relaxation doesn't work well the other way around.
				paMesh->M1[idx] = sM1[idx] + (16 * sEval0[idx] / 135 + 6656 * sEval2[idx] / 12825 + 28561 * sEval3[idx] / 56430 - 9 * sEval4[idx] / 50 + 2 * rhs / 55) * dT;

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