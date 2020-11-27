#include "stdafx.h"
#include "DiffEqAFM.h"

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC
#ifdef ODE_EVAL_COMPILATION_RKDP

#include "Mesh_AntiFerromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//--------------------------------------------- RUNGE KUTTA DORMAND-PRINCE (4th order solution, 5th order error)

void DifferentialEquationAFM::RunRKDP54_Step0_withReductions(void)
{
	mxh_reduction.new_minmax_reduction();
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//obtain maximum normalized torque term
				double Mnorm = pMesh->M[idx].norm();
				double _mxh = GetMagnitude(pMesh->M[idx] ^ pMesh->Heff[idx]) / (Mnorm * Mnorm);
				mxh_reduction.reduce_max(_mxh);

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);
				DBL3 rhs_2 = Equation_Eval_2[omp_get_thread_num()];

				//Now calculate 5th order evaluation for adaptive time step -> FSAL property (a full pass required for this to be valid)
				DBL3 prediction = sM1[idx] + (5179 * sEval0[idx] / 57600 + 7571 * sEval2[idx] / 16695 + 393 * sEval3[idx] / 640 - 92097 * sEval4[idx] / 339200 + 187 * sEval5[idx] / 2100 + rhs / 40) * dT;

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(pMesh->M[idx] - prediction) / Mnorm;
				lte_reduction.reduce_max(_lte);

				//save evaluation for later use
				sEval0[idx] = rhs;
				sEval0_2[idx] = rhs_2;
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

	lte_reduction.maximum();
}

void DifferentialEquationAFM::RunRKDP54_Step0(void)
{
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);
				DBL3 rhs_2 = Equation_Eval_2[omp_get_thread_num()];

				//Now calculate 5th order evaluation for adaptive time step -> FSAL property (a full pass required for this to be valid)
				DBL3 prediction = sM1[idx] + (5179 * sEval0[idx] / 57600 + 7571 * sEval2[idx] / 16695 + 393 * sEval3[idx] / 640 - 92097 * sEval4[idx] / 339200 + 187 * sEval5[idx] / 2100 + rhs / 40) * dT;

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(pMesh->M[idx] - prediction) / pMesh->M[idx].norm();
				lte_reduction.reduce_max(_lte);

				//save evaluation for later use
				sEval0[idx] = rhs;
				sEval0_2[idx] = rhs_2;
			}
		}
	}

	lte_reduction.maximum();
}

void DifferentialEquationAFM::RunRKDP54_Step0_Advance(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//Save current magnetization for later use
				sM1[idx] = pMesh->M[idx];
				sM1_2[idx] = pMesh->M2[idx];

				//Now estimate magnetization using RKDP first step
				pMesh->M[idx] += sEval0[idx] * (dT / 5);
				pMesh->M2[idx] += sEval0_2[idx] * (dT / 5);
			}
		}
	}
}

void DifferentialEquationAFM::RunRKDP54_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval1[idx] = CALLFP(this, equation)(idx);
			sEval1_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RKDP midle step 1
			pMesh->M[idx] = sM1[idx] + (3 * sEval0[idx] / 40 + 9 * sEval1[idx] / 40) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (3 * sEval0_2[idx] / 40 + 9 * sEval1_2[idx] / 40) * dT;
		}
	}
}

void DifferentialEquationAFM::RunRKDP54_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval2[idx] = CALLFP(this, equation)(idx);
			sEval2_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RKDP midle step 2
			pMesh->M[idx] = sM1[idx] + (44 * sEval0[idx] / 45 - 56 * sEval1[idx] / 15 + 32 * sEval2[idx] / 9) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (44 * sEval0_2[idx] / 45 - 56 * sEval1_2[idx] / 15 + 32 * sEval2_2[idx] / 9) * dT;
		}
	}
}

void DifferentialEquationAFM::RunRKDP54_Step3(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval3[idx] = CALLFP(this, equation)(idx);
			sEval3_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RKDP midle step 3
			pMesh->M[idx] = sM1[idx] + (19372 * sEval0[idx] / 6561 - 25360 * sEval1[idx] / 2187 + 64448 * sEval2[idx] / 6561 - 212 * sEval3[idx] / 729) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (19372 * sEval0_2[idx] / 6561 - 25360 * sEval1_2[idx] / 2187 + 64448 * sEval2_2[idx] / 6561 - 212 * sEval3_2[idx] / 729) * dT;
		}
	}
}

void DifferentialEquationAFM::RunRKDP54_Step4(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval4[idx] = CALLFP(this, equation)(idx);
			sEval4_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RKDP midle step 4
			pMesh->M[idx] = sM1[idx] + (9017 * sEval0[idx] / 3168 - 355 * sEval1[idx] / 33 + 46732 * sEval2[idx] / 5247 + 49 * sEval3[idx] / 176 - 5103 * sEval4[idx] / 18656) * dT;
			pMesh->M2[idx] = sM1_2[idx] + (9017 * sEval0_2[idx] / 3168 - 355 * sEval1_2[idx] / 33 + 46732 * sEval2_2[idx] / 5247 + 49 * sEval3_2[idx] / 176 - 5103 * sEval4_2[idx] / 18656) * dT;
		}
	}
}

void DifferentialEquationAFM::RunRKDP54_Step5_withReductions(void)
{
	dmdt_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval5[idx] = CALLFP(this, equation)(idx);
				sEval5_2[idx] = Equation_Eval_2[omp_get_thread_num()];

				//RKDP54 : 5th order evaluation
				pMesh->M[idx] = sM1[idx] + (35 * sEval0[idx] / 384 + 500 * sEval2[idx] / 1113 + 125 * sEval3[idx] / 192 - 2187 * sEval4[idx] / 6784 + 11 * sEval5[idx] / 84) * dT;
				pMesh->M2[idx] = sM1_2[idx] + (35 * sEval0_2[idx] / 384 + 500 * sEval2_2[idx] / 1113 + 125 * sEval3_2[idx] / 192 - 2187 * sEval4_2[idx] / 6784 + 11 * sEval5_2[idx] / 84) * dT;

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
}

void DifferentialEquationAFM::RunRKDP54_Step5(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval5[idx] = CALLFP(this, equation)(idx);
				sEval5_2[idx] = Equation_Eval_2[omp_get_thread_num()];

				//RKDP54 : 5th order evaluation
				pMesh->M[idx] = sM1[idx] + (35 * sEval0[idx] / 384 + 500 * sEval2[idx] / 1113 + 125 * sEval3[idx] / 192 - 2187 * sEval4[idx] / 6784 + 11 * sEval5[idx] / 84) * dT;
				pMesh->M2[idx] = sM1_2[idx] + (35 * sEval0_2[idx] / 384 + 500 * sEval2_2[idx] / 1113 + 125 * sEval3_2[idx] / 192 - 2187 * sEval4_2[idx] / 6784 + 11 * sEval5_2[idx] / 84) * dT;

				if (renormalize) {

					DBL2 Ms_AFM = pMesh->Ms_AFM;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
					pMesh->M[idx].renormalize(Ms_AFM.i);
					pMesh->M2[idx].renormalize(Ms_AFM.j);
				}
			}
			else {

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
				pMesh->M[idx].renormalize(Ms_AFM.i);
				pMesh->M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}
}

#endif
#endif