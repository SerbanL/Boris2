#include "stdafx.h"
#include "DiffEqAFM.h"

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC
#ifdef ODE_EVAL_COMPILATION_RKF45

#include "Mesh_AntiFerromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

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
				pMesh->M[idx] += sEval0[idx] * (2 * dT / 9);
				pMesh->M2[idx] += sEval0_2[idx] * (2 * dT / 9);
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
				pMesh->M[idx] += sEval0[idx] * (2 * dT / 9);
				pMesh->M2[idx] += sEval0_2[idx] * (2 * dT / 9);
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
			pMesh->M[idx] = sM1[idx] + (sEval0[idx] / 12 + sEval1[idx] / 4) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (sEval0_2[idx] / 12 + sEval1_2[idx] / 4) * dT;
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
			pMesh->M[idx] = sM1[idx] + (69 * sEval0[idx] / 128 - 243 * sEval1[idx] / 128 + 135 * sEval2[idx] / 64) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (69 * sEval0_2[idx] / 128 - 243 * sEval1_2[idx] / 128 + 135 * sEval2_2[idx] / 64) * dT;
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
			pMesh->M[idx] = sM1[idx] + (-17 * sEval0[idx] / 12 + 27 * sEval1[idx] / 4 - 27 * sEval2[idx] / 5 + 16 * sEval3[idx] / 15) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (-17 * sEval0_2[idx] / 12 + 27 * sEval1_2[idx] / 4 - 27 * sEval2_2[idx] / 5 + 16 * sEval3_2[idx] / 15) * dT;
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
			pMesh->M[idx] = sM1[idx] + (65 * sEval0[idx] / 432 - 5 * sEval1[idx] / 16 + 13 * sEval2[idx] / 16 + 4 * sEval3[idx] / 27 + 5 * sEval4[idx] / 144) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (65 * sEval0_2[idx] / 432 - 5 * sEval1_2[idx] / 16 + 13 * sEval2_2[idx] / 16 + 4 * sEval3_2[idx] / 27 + 5 * sEval4_2[idx] / 144) * dT;
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
				pMesh->M[idx] = sM1[idx] + (sEval0[idx] / 9 + 9 * sEval2[idx] / 20 + 16 * sEval3[idx] / 45 + sEval4[idx] / 12) * dT;
				pMesh->M2[idx] = sM1_2[idx] + (sEval0_2[idx] / 9 + 9 * sEval2_2[idx] / 20 + 16 * sEval3_2[idx] / 45 + sEval4_2[idx] / 12) * dT;

				//5th order evaluation
				DBL3 prediction = sM1[idx] + (47 * sEval0[idx] / 450 + 12 * sEval2[idx] / 25 + 32 * sEval3[idx] / 225 + 1 * sEval4[idx] / 30 + 6 * rhs / 25) * dT;

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

		//only reduce for dmdt if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
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
				pMesh->M[idx] = sM1[idx] + (sEval0[idx] / 9 + 9 * sEval2[idx] / 20 + 16 * sEval3[idx] / 45 + sEval4[idx] / 12) * dT;
				pMesh->M2[idx] = sM1_2[idx] + (sEval0_2[idx] / 9 + 9 * sEval2_2[idx] / 20 + 16 * sEval3_2[idx] / 45 + sEval4_2[idx] / 12) * dT;

				//5th order evaluation
				DBL3 prediction = sM1[idx] + (47 * sEval0[idx] / 450 + 12 * sEval2[idx] / 25 + 32 * sEval3[idx] / 225 + 1 * sEval4[idx] / 30 + 6 * rhs / 25) * dT;

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
#endif