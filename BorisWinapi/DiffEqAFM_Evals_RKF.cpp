#include "stdafx.h"
#include "DiffEqAFM.h"
#include "Mesh_AntiFerromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

#ifdef ODE_EVAL_RKF

//--------------------------------------------- RUNGE KUTTA FEHLBERG (4th order solution, 5th order error)

void DifferentialEquationAFM::RunRKF45_Step0_withReductions(void)
{
	mxh_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization for later use
			sM1[idx] = pMesh->M[idx];
			sM1_2[idx] = pMesh->M2[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//obtain maximum normalized torque term
				double Mnorm = pMesh->M[idx].norm();
				double _mxh = GetMagnitude(pMesh->M[idx] ^ pMesh->Heff[idx]) / (Mnorm * Mnorm);
				mxh_reduction.reduce_max(_mxh);

				//First evaluate RHS of set equation at the current time step
				sEval0[idx] = CALLFP(this, equation)(idx);
				sEval0_2[idx] = Equation_Eval_2[omp_get_thread_num()];

				//Now estimate magnetization using RKF first step
				pMesh->M[idx] += sEval0[idx] * (dT / 4);
				pMesh->M2[idx] += sEval0_2[idx] * (dT / 4);
			}
		}
	}

	if (pMesh->grel.get0()) {

		//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics are disabled in this mesh)
		mxh_reduction.maximum();
	}
	else {

		mxh_reduction.max = 0.0;
	}
}

void DifferentialEquationAFM::RunRKF45_Step0(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization for later use
			sM1[idx] = pMesh->M[idx];
			sM1_2[idx] = pMesh->M2[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval0[idx] = CALLFP(this, equation)(idx);
				sEval0_2[idx] = Equation_Eval_2[omp_get_thread_num()];

				//Now estimate magnetization using RKF first step
				pMesh->M[idx] += sEval0[idx] * (dT / 4);
				pMesh->M2[idx] += sEval0_2[idx] * (dT / 4);
			}
		}
	}
}

void DifferentialEquationAFM::RunRKF45_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval1[idx] = CALLFP(this, equation)(idx);
			sEval1_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RKF midle step 1
			pMesh->M[idx] = sM1[idx] + (3 * sEval0[idx] + 9 * sEval1[idx]) * dT / 32;
			pMesh->M2[idx] = sM1_2[idx] + (3 * sEval0_2[idx] + 9 * sEval1_2[idx]) * dT / 32;
		}
	}
}

void DifferentialEquationAFM::RunRKF45_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval2[idx] = CALLFP(this, equation)(idx);
			sEval2_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RKF midle step 2
			pMesh->M[idx] = sM1[idx] + (1932 * sEval0[idx] - 7200 * sEval1[idx] + 7296 * sEval2[idx]) * dT / 2197;
			pMesh->M2[idx] = sM1_2[idx] + (1932 * sEval0_2[idx] - 7200 * sEval1_2[idx] + 7296 * sEval2_2[idx]) * dT / 2197;
		}
	}
}

void DifferentialEquationAFM::RunRKF45_Step3(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval3[idx] = CALLFP(this, equation)(idx);
			sEval3_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RKF midle step 3
			pMesh->M[idx] = sM1[idx] + (439 * sEval0[idx] / 216 - 8 * sEval1[idx] + 3680 * sEval2[idx] / 513 - 845 * sEval3[idx] / 4104) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (439 * sEval0_2[idx] / 216 - 8 * sEval1_2[idx] + 3680 * sEval2_2[idx] / 513 - 845 * sEval3_2[idx] / 4104) * dT;
		}
	}
}

void DifferentialEquationAFM::RunRKF45_Step4(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval4[idx] = CALLFP(this, equation)(idx);
			sEval4_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RKF midle step 4
			pMesh->M[idx] = sM1[idx] + (-8 * sEval0[idx] / 27 + 2 * sEval1[idx] - 3544 * sEval2[idx] / 2565 + 1859 * sEval3[idx] / 4104 - 11 * sEval4[idx] / 40) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (-8 * sEval0_2[idx] / 27 + 2 * sEval1_2[idx] - 3544 * sEval2_2[idx] / 2565 + 1859 * sEval3_2[idx] / 4104 - 11 * sEval4_2[idx] / 40) * dT;
		}
	}
}

void DifferentialEquationAFM::RunRKF45_Step5_withReductions(void)
{
	dmdt_reduction.new_minmax_reduction();
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);
				DBL3 rhs_2 = Equation_Eval_2[omp_get_thread_num()];

				//4th order evaluation
				DBL3 prediction = sM1[idx] + (25 * sEval0[idx] / 216 + 1408 * sEval2[idx] / 2565 + 2197 * sEval3[idx] / 4101 - sEval4[idx] / 5) * dT;

				//5th order evaluation -> keep this as the new value, not the 4th order; relaxation doesn't work well the other way around.
				pMesh->M[idx] = sM1[idx] + (16 * sEval0[idx] / 135 + 6656 * sEval2[idx] / 12825 + 28561 * sEval3[idx] / 56430 - 9 * sEval4[idx] / 50 + 2 * rhs / 55) * dT;
				pMesh->M2[idx] = sM1_2[idx] + (16 * sEval0_2[idx] / 135 + 6656 * sEval2_2[idx] / 12825 + 28561 * sEval3_2[idx] / 56430 - 9 * sEval4_2[idx] / 50 + 2 * rhs_2 / 55) * dT;

				if (renormalize) {

					DBL2 Ms_AFM = pMesh->Ms_AFM;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
					pMesh->M[idx].renormalize(Ms_AFM.i);
					pMesh->M2[idx].renormalize(Ms_AFM.j);
				}

				//obtained maximum dmdt term
				double Mnorm = pMesh->M[idx].norm();
				double _dmdt = GetMagnitude(pMesh->M[idx] - sM1[idx]) / (dT * GAMMA * Mnorm * Mnorm);
				dmdt_reduction.reduce_max(_dmdt);

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(pMesh->M[idx] - prediction) / pMesh->M[idx].norm();
				lte_reduction.reduce_max(_lte);
			}
			else {

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
				pMesh->M[idx].renormalize(Ms_AFM.i);
				pMesh->M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	if (pMesh->grel.get0()) {

		//only reduce for dmdt if grel is not zero (if it's zero this means magnetisation dynamics are disabled in this mesh)
		dmdt_reduction.maximum();
	}
	else {

		dmdt_reduction.max = 0.0;
	}

	lte_reduction.maximum();
}

void DifferentialEquationAFM::RunRKF45_Step5(void)
{
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);
				DBL3 rhs_2 = Equation_Eval_2[omp_get_thread_num()];

				//4th order evaluation
				DBL3 prediction = sM1[idx] + (25 * sEval0[idx] / 216 + 1408 * sEval2[idx] / 2565 + 2197 * sEval3[idx] / 4101 - sEval4[idx] / 5) * dT;

				//5th order evaluation -> keep this as the new value, not the 4th order; relaxation doesn't work well the other way around.
				pMesh->M[idx] = sM1[idx] + (16 * sEval0[idx] / 135 + 6656 * sEval2[idx] / 12825 + 28561 * sEval3[idx] / 56430 - 9 * sEval4[idx] / 50 + 2 * rhs / 55) * dT;
				pMesh->M2[idx] = sM1_2[idx] + (16 * sEval0_2[idx] / 135 + 6656 * sEval2_2[idx] / 12825 + 28561 * sEval3_2[idx] / 56430 - 9 * sEval4_2[idx] / 50 + 2 * rhs_2 / 55) * dT;

				if (renormalize) {

					DBL2 Ms_AFM = pMesh->Ms_AFM;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
					pMesh->M[idx].renormalize(Ms_AFM.i);
					pMesh->M2[idx].renormalize(Ms_AFM.j);
				}

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(pMesh->M[idx] - prediction) / pMesh->M[idx].norm();
				lte_reduction.reduce_max(_lte);
			}
			else {

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
				pMesh->M[idx].renormalize(Ms_AFM.i);
				pMesh->M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}

	lte_reduction.maximum();
}

#endif
