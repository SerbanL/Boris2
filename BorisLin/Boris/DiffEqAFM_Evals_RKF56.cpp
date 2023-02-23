#include "stdafx.h"
#include "DiffEqAFM.h"

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC
#ifdef ODE_EVAL_COMPILATION_RKF56

#include "Mesh_AntiFerromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//--------------------------------------------- RUNGE KUTTA FEHLBERG (4th order solution, 5th order error)

void DifferentialEquationAFM::RunRKF56_Step0_withReductions(void)
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
				pMesh->M[idx] += sEval0[idx] * (dT / 6);
				pMesh->M2[idx] += sEval0_2[idx] * (dT / 6);
			}
		}
	}

	if (pMesh->grel.get0()) {

		//only reduce for mxh if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
		mxh_reduction.maximum();
	}
	else {

		mxh_reduction.max = 0.0;
	}
}

void DifferentialEquationAFM::RunRKF56_Step0(void)
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
				pMesh->M[idx] += sEval0[idx] * (dT / 6);
				pMesh->M2[idx] += sEval0_2[idx] * (dT / 6);
			}
		}
	}
}

void DifferentialEquationAFM::RunRKF56_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval1[idx] = CALLFP(this, equation)(idx);
			sEval1_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RKF midle step 1
			pMesh->M[idx] = sM1[idx] + (4 * sEval0[idx] + 16 * sEval1[idx]) * dT / 75;
			pMesh->M2[idx] = sM1_2[idx] + (4 * sEval0_2[idx] + 16 * sEval1_2[idx]) * dT / 75;
		}
	}
}

void DifferentialEquationAFM::RunRKF56_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval2[idx] = CALLFP(this, equation)(idx);
			sEval2_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RKF midle step 2
			pMesh->M[idx] = sM1[idx] + (5 * sEval0[idx] / 6 - 8 * sEval1[idx] / 3 + 5 * sEval2[idx] / 2) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (5 * sEval0_2[idx] / 6 - 8 * sEval1_2[idx] / 3 + 5 * sEval2_2[idx] / 2) * dT;
		}
	}
}

void DifferentialEquationAFM::RunRKF56_Step3(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval3[idx] = CALLFP(this, equation)(idx);
			sEval3_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RKF midle step 3
			pMesh->M[idx] = sM1[idx] + (-8 * sEval0[idx] / 5 + 144 * sEval1[idx] / 25 - 4 * sEval2[idx] + 16 * sEval3[idx] / 25) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (-8 * sEval0_2[idx] / 5 + 144 * sEval1_2[idx] / 25 - 4 * sEval2_2[idx] + 16 * sEval3_2[idx] / 25) * dT;
		}
	}
}

void DifferentialEquationAFM::RunRKF56_Step4(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval4[idx] = CALLFP(this, equation)(idx);
			sEval4_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RKF midle step 4
			pMesh->M[idx] = sM1[idx] + (361 * sEval0[idx] / 320 - 18 * sEval1[idx] / 5 + 407 * sEval2[idx] / 128 - 11 * sEval3[idx] / 80 + 55 * sEval4[idx] / 128) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (361 * sEval0_2[idx] / 320 - 18 * sEval1_2[idx] / 5 + 407 * sEval2_2[idx] / 128 - 11 * sEval3_2[idx] / 80 + 55 * sEval4_2[idx] / 128) * dT;
		}
	}
}

void DifferentialEquationAFM::RunRKF56_Step5(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval5[idx] = CALLFP(this, equation)(idx);
			sEval5_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RKF midle step 4
			pMesh->M[idx] = sM1[idx] + (-11 * sEval0[idx] / 640 + 11 * sEval2[idx] / 256 - 11 * sEval3[idx] / 160 + 11 * sEval4[idx] / 256) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (-11 * sEval0_2[idx] / 640 + 11 * sEval2_2[idx] / 256 - 11 * sEval3_2[idx] / 160 + 11 * sEval4_2[idx] / 256) * dT;
		}
	}
}

void DifferentialEquationAFM::RunRKF56_Step6(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval6[idx] = CALLFP(this, equation)(idx);
			sEval6_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RKF midle step 4
			pMesh->M[idx] = sM1[idx] + (93 * sEval0[idx] / 640 - 18 * sEval1[idx] / 5 + 803 * sEval2[idx] / 256 - 11 * sEval3[idx] / 160 + 99 * sEval4[idx] / 256 + sEval6[idx]) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (93 * sEval0_2[idx] / 640 - 18 * sEval1_2[idx] / 5 + 803 * sEval2_2[idx] / 256 - 11 * sEval3_2[idx] / 160 + 99 * sEval4_2[idx] / 256 + sEval6_2[idx]) * dT;
		}
	}
}

void DifferentialEquationAFM::RunRKF56_Step7_withReductions(void)
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

				//5th order evaluation
				pMesh->M[idx] = sM1[idx] + (31 * sEval0[idx] / 384 + 1125 * sEval2[idx] / 2816 + 9 * sEval3[idx] / 32 + 125 * sEval4[idx] / 768 + 5 * sEval5[idx] / 66) * dT;
				pMesh->M2[idx] = sM1_2[idx] + (31 * sEval0_2[idx] / 384 + 1125 * sEval2_2[idx] / 2816 + 9 * sEval3_2[idx] / 32 + 125 * sEval4_2[idx] / 768 + 5 * sEval5_2[idx] / 66) * dT;

				//local truncation error from 5th order evaluation and 6th order evaluation
				DBL3 lte_diff = 5 * (sEval0[idx] + sEval5[idx] - sEval6[idx] - rhs) * dT / 66;

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
				double _lte = GetMagnitude(lte_diff) / pMesh->M[idx].norm();
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

		//only reduce for dmdt if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
		dmdt_reduction.maximum();
	}
	else {

		dmdt_reduction.max = 0.0;
	}

	lte_reduction.maximum();
}

void DifferentialEquationAFM::RunRKF56_Step7(void)
{
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);
				DBL3 rhs_2 = Equation_Eval_2[omp_get_thread_num()];

				//5th order evaluation
				pMesh->M[idx] = sM1[idx] + (31 * sEval0[idx] / 384 + 1125 * sEval2[idx] / 2816 + 9 * sEval3[idx] / 32 + 125 * sEval4[idx] / 768 + 5 * sEval5[idx] / 66) * dT;
				pMesh->M2[idx] = sM1_2[idx] + (31 * sEval0_2[idx] / 384 + 1125 * sEval2_2[idx] / 2816 + 9 * sEval3_2[idx] / 32 + 125 * sEval4_2[idx] / 768 + 5 * sEval5_2[idx] / 66) * dT;

				//local truncation error from 5th order evaluation and 6th order evaluation
				DBL3 lte_diff = 5 * (sEval0[idx] + sEval5[idx] - sEval6[idx] - rhs) * dT / 66;

				if (renormalize) {

					DBL2 Ms_AFM = pMesh->Ms_AFM;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
					pMesh->M[idx].renormalize(Ms_AFM.i);
					pMesh->M2[idx].renormalize(Ms_AFM.j);
				}

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(lte_diff) / pMesh->M[idx].norm();
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
#endif