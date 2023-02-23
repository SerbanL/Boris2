#include "stdafx.h"
#include "Atom_DiffEqCubic.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC
#ifdef ODE_EVAL_COMPILATION_RKF56

#include "Atom_Mesh_Cubic.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

//--------------------------------------------- RUNGE KUTTA FEHLBERG (5th order solution, 6th order error)

void Atom_DifferentialEquationCubic::RunRKF56_Step0_withReductions(void)
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
				paMesh->M1[idx] += sEval0[idx] * (dT / 6);
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

void Atom_DifferentialEquationCubic::RunRKF56_Step0(void)
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
				paMesh->M1[idx] += sEval0[idx] * (dT / 6);
			}
		}
	}
}

void Atom_DifferentialEquationCubic::RunRKF56_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx) && !paMesh->M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval1[idx] = CALLFP(this, equation)(idx);

			//Now estimate moment using RKF midle step 1
			paMesh->M1[idx] = sM1[idx] + (4 * sEval0[idx] + 16 * sEval1[idx]) * dT / 75;
		}
	}
}

void Atom_DifferentialEquationCubic::RunRKF56_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx) && !paMesh->M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval2[idx] = CALLFP(this, equation)(idx);

			//Now estimate moment using RKF midle step 2
			paMesh->M1[idx] = sM1[idx] + (5 * sEval0[idx] / 6 - 8 * sEval1[idx] / 3 + 5 * sEval2[idx] / 2) * dT;
		}
	}
}

void Atom_DifferentialEquationCubic::RunRKF56_Step3(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx) && !paMesh->M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval3[idx] = CALLFP(this, equation)(idx);

			//Now estimate moment using RKF midle step 3
			paMesh->M1[idx] = sM1[idx] + (-8 * sEval0[idx] / 5 + 144 * sEval1[idx] / 25 - 4 * sEval2[idx] + 16 * sEval3[idx] / 25) * dT;
		}
	}
}

void Atom_DifferentialEquationCubic::RunRKF56_Step4(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx) && !paMesh->M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval4[idx] = CALLFP(this, equation)(idx);

			//Now estimate moment using RKF midle step 4
			paMesh->M1[idx] = sM1[idx] + (361 * sEval0[idx] / 320 - 18 * sEval1[idx] / 5 + 407 * sEval2[idx] / 128 - 11 * sEval3[idx] / 80 + 55 * sEval4[idx] / 128) * dT;
		}
	}
}

void Atom_DifferentialEquationCubic::RunRKF56_Step5(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx) && !paMesh->M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval5[idx] = CALLFP(this, equation)(idx);

			paMesh->M1[idx] = sM1[idx] + (-11 * sEval0[idx] / 640 + 11 * sEval2[idx] / 256 - 11 * sEval3[idx] / 160 + 11 * sEval4[idx] / 256) * dT;
		}
	}
}

void Atom_DifferentialEquationCubic::RunRKF56_Step6(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx) && !paMesh->M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval6[idx] = CALLFP(this, equation)(idx);

			paMesh->M1[idx] = sM1[idx] + (93 * sEval0[idx] / 640 - 18 * sEval1[idx] / 5 + 803 * sEval2[idx] / 256 - 11 * sEval3[idx] / 160 + 99 * sEval4[idx] / 256 + sEval6[idx]) * dT;
		}
	}
}

void Atom_DifferentialEquationCubic::RunRKF56_Step7_withReductions(void)
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

				//5th order evaluation
				paMesh->M1[idx] = sM1[idx] + (31 * sEval0[idx] / 384 + 1125 * sEval2[idx] / 2816 + 9 * sEval3[idx] / 32 + 125 * sEval4[idx] / 768 + 5 * sEval5[idx] / 66) * dT;

				//local truncation error from 5th order evaluation and 6th order evaluation
				DBL3 lte_diff = 5 * (sEval0[idx] + sEval5[idx] - sEval6[idx] - rhs) * dT / 66;

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
				double _lte = GetMagnitude(lte_diff) / paMesh->M1[idx].norm();
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

void Atom_DifferentialEquationCubic::RunRKF56_Step7(void)
{
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			if (!paMesh->M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//5th order evaluation
				paMesh->M1[idx] = sM1[idx] + (31 * sEval0[idx] / 384 + 1125 * sEval2[idx] / 2816 + 9 * sEval3[idx] / 32 + 125 * sEval4[idx] / 768 + 5 * sEval5[idx] / 66) * dT;

				//local truncation error from 5th order evaluation and 6th order evaluation
				DBL3 lte_diff = 5 * (sEval0[idx] + sEval5[idx] - sEval6[idx] - rhs) * dT / 66;

				if (renormalize) {

					double mu_s = paMesh->mu_s;
					paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);
					paMesh->M1[idx].renormalize(mu_s);
				}

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(lte_diff) / paMesh->M1[idx].norm();
				lte_reduction.reduce_max(_lte);
			}
		}
	}

	lte_reduction.maximum();
}

#endif
#endif