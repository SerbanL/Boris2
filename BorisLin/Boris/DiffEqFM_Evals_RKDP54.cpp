#include "stdafx.h"
#include "DiffEqFM.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC
#ifdef ODE_EVAL_COMPILATION_RKDP

#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//--------------------------------------------- RUNGE KUTTA DORMAND-PRINCE (4th order solution, 5th order error)

void DifferentialEquationFM::RunRKDP54_Step0_withReductions(void)
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

				//Now calculate 5th order evaluation for adaptive time step -> FSAL property (a full pass required for this to be valid)
				DBL3 prediction = sM1[idx] + (5179 * sEval0[idx] / 57600 + 7571 * sEval2[idx] / 16695 + 393 * sEval3[idx] / 640 - 92097 * sEval4[idx] / 339200 + 187 * sEval5[idx] / 2100 + rhs / 40) * dT;

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(pMesh->M[idx] - prediction) / Mnorm;
				lte_reduction.reduce_max(_lte);

				//save evaluation for later use
				sEval0[idx] = rhs;
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

void DifferentialEquationFM::RunRKDP54_Step0(void)
{
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//Now calculate 5th order evaluation for adaptive time step -> FSAL property (a full pass required for this to be valid)
				DBL3 prediction = sM1[idx] + (5179 * sEval0[idx] / 57600 + 7571 * sEval2[idx] / 16695 + 393 * sEval3[idx] / 640 - 92097 * sEval4[idx] / 339200 + 187 * sEval5[idx] / 2100 + rhs / 40) * dT;

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(pMesh->M[idx] - prediction) / pMesh->M[idx].norm();
				lte_reduction.reduce_max(_lte);

				//save evaluation for later use
				sEval0[idx] = rhs;
			}
		}
	}

	lte_reduction.maximum();
}

void DifferentialEquationFM::RunRKDP54_Step0_Advance(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//Save current magnetization for later use
				sM1[idx] = pMesh->M[idx];

				//Now estimate magnetization using RKDP first step
				pMesh->M[idx] += sEval0[idx] * (dT / 5);
			}
		}
	}
}

void DifferentialEquationFM::RunRKDP54_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval1[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKDP midle step 1
			pMesh->M[idx] = sM1[idx] + (3 * sEval0[idx] / 40 + 9 * sEval1[idx] / 40) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKDP54_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval2[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKDP midle step 2
			pMesh->M[idx] = sM1[idx] + (44 * sEval0[idx] / 45 - 56 * sEval1[idx] / 15 + 32 * sEval2[idx] / 9) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKDP54_Step3(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval3[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKDP midle step 3
			pMesh->M[idx] = sM1[idx] + (19372 * sEval0[idx] / 6561 - 25360 * sEval1[idx] / 2187 + 64448 * sEval2[idx] / 6561 - 212 * sEval3[idx] / 729) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKDP54_Step4(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval4[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKDP midle step 4
			pMesh->M[idx] = sM1[idx] + (9017 * sEval0[idx] / 3168 - 355 * sEval1[idx] / 33 + 46732 * sEval2[idx] / 5247 + 49 * sEval3[idx] / 176 - 5103 * sEval4[idx] / 18656) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKDP54_Step5_withReductions(void)
{
	dmdt_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval5[idx] = CALLFP(this, equation)(idx);

				//RKDP54 : 5th order evaluation
				pMesh->M[idx] = sM1[idx] + (35 * sEval0[idx] / 384 + 500 * sEval2[idx] / 1113 + 125 * sEval3[idx] / 192 - 2187 * sEval4[idx] / 6784 + 11 * sEval5[idx] / 84) * dT;

				if (renormalize) {

					double Ms = pMesh->Ms;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
					pMesh->M[idx].renormalize(Ms);
				}

				//obtained maximum dmdt term
				double Mnorm = pMesh->M[idx].norm();
				double _dmdt = GetMagnitude(pMesh->M[idx] - sM1[idx]) / (dT * GAMMA * Mnorm * Mnorm);
				dmdt_reduction.reduce_max(_dmdt);
			}
			else {

				double Ms = pMesh->Ms;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
				pMesh->M[idx].renormalize(Ms);
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

void DifferentialEquationFM::RunRKDP54_Step5(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval5[idx] = CALLFP(this, equation)(idx);

				//RKDP54 : 5th order evaluation
				pMesh->M[idx] = sM1[idx] + (35 * sEval0[idx] / 384 + 500 * sEval2[idx] / 1113 + 125 * sEval3[idx] / 192 - 2187 * sEval4[idx] / 6784 + 11 * sEval5[idx] / 84) * dT;

				if (renormalize) {

					double Ms = pMesh->Ms;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
					pMesh->M[idx].renormalize(Ms);
				}
			}
			else {

				double Ms = pMesh->Ms;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
				pMesh->M[idx].renormalize(Ms);
			}
		}
	}
}

#endif
#endif