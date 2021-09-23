#include "stdafx.h"
#include "DiffEqFM.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC
#ifdef ODE_EVAL_COMPILATION_RKF56

#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//--------------------------------------------- RUNGE KUTTA FEHLBERG (4th order solution, 5th order error)

void DifferentialEquationFM::RunRKF56_Step0_withReductions(void)
{
	mxh_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization for later use
			sM1[idx] = pMesh->M[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//obtain maximum normalized torque term
				double Mnorm = pMesh->M[idx].norm();
				double _mxh = GetMagnitude(pMesh->M[idx] ^ pMesh->Heff[idx]) / (Mnorm * Mnorm);
				mxh_reduction.reduce_max(_mxh);

				//First evaluate RHS of set equation at the current time step
				sEval0[idx] = CALLFP(this, equation)(idx);

				//Now estimate magnetization using RKF first step
				pMesh->M[idx] += sEval0[idx] * (dT / 6);
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

void DifferentialEquationFM::RunRKF56_Step0(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization for later use
			sM1[idx] = pMesh->M[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval0[idx] = CALLFP(this, equation)(idx);

				//Now estimate magnetization using RKF first step
				pMesh->M[idx] += sEval0[idx] * (dT / 6);
			}
		}
	}
}

void DifferentialEquationFM::RunRKF56_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval1[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKF midle step 1
			pMesh->M[idx] = sM1[idx] + (4 * sEval0[idx] + 16 * sEval1[idx]) * dT / 75;
		}
	}
}

void DifferentialEquationFM::RunRKF56_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval2[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKF midle step 2
			pMesh->M[idx] = sM1[idx] + (5 * sEval0[idx] / 6 - 8 * sEval1[idx] / 3 + 5 * sEval2[idx] / 2) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKF56_Step3(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval3[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKF midle step 3
			pMesh->M[idx] = sM1[idx] + (-8 * sEval0[idx] / 5 + 144 * sEval1[idx] / 25 - 4 * sEval2[idx] + 16 * sEval3[idx] / 25) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKF56_Step4(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval4[idx] = CALLFP(this, equation)(idx);

			pMesh->M[idx] = sM1[idx] + (361 * sEval0[idx] / 320 - 18 * sEval1[idx] / 5 + 407 * sEval2[idx] / 128 - 11 * sEval3[idx] / 80 + 55 * sEval4[idx] / 128) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKF56_Step5(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval5[idx] = CALLFP(this, equation)(idx);

			pMesh->M[idx] = sM1[idx] + (-11 * sEval0[idx] / 640 + 11 * sEval2[idx] / 256 - 11 * sEval3[idx] / 160 + 11 * sEval4[idx] / 256) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKF56_Step6(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval6[idx] = CALLFP(this, equation)(idx);

			pMesh->M[idx] = sM1[idx] + (93 * sEval0[idx] / 640 - 18 * sEval1[idx] / 5 + 803 * sEval2[idx] / 256 - 11 * sEval3[idx] / 160 + 99 * sEval4[idx] / 256 + sEval6[idx]) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKF56_Step7_withReductions(void)
{
	dmdt_reduction.new_minmax_reduction();
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//5th order evaluation
				pMesh->M[idx] = sM1[idx] + (31 * sEval0[idx] / 384 + 1125 * sEval2[idx] / 2816 + 9 * sEval3[idx] / 32 + 125 * sEval4[idx] / 768 + 5 * sEval5[idx] / 66) * dT;

				//local truncation error from 5th order evaluation and 6th order evaluation
				DBL3 lte_diff = 5 * (sEval0[idx] + sEval5[idx] - sEval6[idx] - rhs) * dT / 66;

				if (renormalize) {

					double Ms = pMesh->Ms;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
					pMesh->M[idx].renormalize(Ms);
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

	lte_reduction.maximum();
}

void DifferentialEquationFM::RunRKF56_Step7(void)
{
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//5th order evaluation
				pMesh->M[idx] = sM1[idx] + (31 * sEval0[idx] / 384 + 1125 * sEval2[idx] / 2816 + 9 * sEval3[idx] / 32 + 125 * sEval4[idx] / 768 + 5 * sEval5[idx] / 66) * dT;

				//local truncation error from 5th order evaluation and 6th order evaluation
				DBL3 lte_diff = 5 * (sEval0[idx] + sEval5[idx] - sEval6[idx] - rhs) * dT / 66;

				if (renormalize) {

					double Ms = pMesh->Ms;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
					pMesh->M[idx].renormalize(Ms);
				}

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(lte_diff) / pMesh->M[idx].norm();
				lte_reduction.reduce_max(_lte);
			}
			else {

				double Ms = pMesh->Ms;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
				pMesh->M[idx].renormalize(Ms);
			}
		}
	}

	lte_reduction.maximum();
}

#endif
#endif