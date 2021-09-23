#include "stdafx.h"
#include "DiffEqFM.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC
#ifdef ODE_EVAL_COMPILATION_RKF45

#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//--------------------------------------------- RUNGE KUTTA FEHLBERG (4th order solution, 5th order error)

void DifferentialEquationFM::RunRKF45_Step0_withReductions(void)
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
				pMesh->M[idx] += sEval0[idx] * (2 * dT / 9);
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

void DifferentialEquationFM::RunRKF45_Step0(void)
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
				pMesh->M[idx] += sEval0[idx] * (2 * dT / 9);
			}
		}
	}
}

void DifferentialEquationFM::RunRKF45_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval1[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKF midle step 1
			pMesh->M[idx] = sM1[idx] + (sEval0[idx] / 12 + sEval1[idx] / 4) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKF45_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval2[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKF midle step 2
			pMesh->M[idx] = sM1[idx] + (69 * sEval0[idx] / 128 - 243 * sEval1[idx] / 128 + 135 * sEval2[idx] / 64) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKF45_Step3(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval3[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKF midle step 3
			pMesh->M[idx] = sM1[idx] + (-17 * sEval0[idx] / 12 + 27 * sEval1[idx] / 4 - 27 * sEval2[idx] / 5 + 16 * sEval3[idx] / 15) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKF45_Step4(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval4[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKF midle step 4
			pMesh->M[idx] = sM1[idx] + (65 * sEval0[idx] / 432 - 5 * sEval1[idx] / 16 + 13 * sEval2[idx] / 16 + 4 * sEval3[idx] / 27 + 5 * sEval4[idx] / 144) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKF45_Step5_withReductions(void)
{
	dmdt_reduction.new_minmax_reduction();
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//4th order evaluation
				pMesh->M[idx] = sM1[idx] + (sEval0[idx] / 9 + 9 * sEval2[idx] / 20 + 16 * sEval3[idx] / 45 + sEval4[idx] / 12) * dT;

				//5th order evaluation
				DBL3 prediction = sM1[idx] + (47 * sEval0[idx] / 450 + 12 * sEval2[idx] / 25 + 32 * sEval3[idx] / 225 + 1 * sEval4[idx] / 30 + 6 * rhs / 25) * dT;

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
				double _lte = GetMagnitude(pMesh->M[idx] - prediction) / pMesh->M[idx].norm();
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

void DifferentialEquationFM::RunRKF45_Step5(void)
{
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//4th order evaluation
				pMesh->M[idx] = sM1[idx] + (sEval0[idx] / 9 + 9 * sEval2[idx] / 20 + 16 * sEval3[idx] / 45 + sEval4[idx] / 12) * dT;

				//5th order evaluation
				DBL3 prediction = sM1[idx] + (47 * sEval0[idx] / 450 + 12 * sEval2[idx] / 25 + 32 * sEval3[idx] / 225 + 1 * sEval4[idx] / 30 + 6 * rhs / 25) * dT;

				if (renormalize) {

					double Ms = pMesh->Ms;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
					pMesh->M[idx].renormalize(Ms);
				}

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(pMesh->M[idx] - prediction) / pMesh->M[idx].norm();
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