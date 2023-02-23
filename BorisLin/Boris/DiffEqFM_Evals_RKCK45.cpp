#include "stdafx.h"
#include "DiffEqFM.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC
#ifdef ODE_EVAL_COMPILATION_RKCK

#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//--------------------------------------------- RUNGE KUTTA CASH-KARP (4th order solution, 5th order error)

void DifferentialEquationFM::RunRKCK45_Step0_withReductions(void)
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

				//Now estimate magnetization using RKCK first step
				pMesh->M[idx] += sEval0[idx] * (dT / 5);
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

void DifferentialEquationFM::RunRKCK45_Step0(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization for later use
			sM1[idx] = pMesh->M[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval0[idx] = CALLFP(this, equation)(idx);

				//Now estimate magnetization using RKCK first step
				pMesh->M[idx] += sEval0[idx] * (dT / 5);
			}
		}
	}
}

void DifferentialEquationFM::RunRKCK45_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval1[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKCK midle step 1
			pMesh->M[idx] = sM1[idx] + (3 * sEval0[idx] + 9 * sEval1[idx]) * dT / 40;
		}
	}
}

void DifferentialEquationFM::RunRKCK45_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval2[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKCK midle step 2
			pMesh->M[idx] = sM1[idx] + (3 * sEval0[idx] / 10 - 9 * sEval1[idx] / 10 + 6 * sEval2[idx] / 5) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKCK45_Step3(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval3[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKCK midle step 3
			pMesh->M[idx] = sM1[idx] + (-11 * sEval0[idx] / 54 + 5 * sEval1[idx] / 2 - 70 * sEval2[idx] / 27 + 35 * sEval3[idx] / 27) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKCK45_Step4(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval4[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKCK midle step 4
			pMesh->M[idx] = sM1[idx] + (1631 * sEval0[idx] / 55296 + 175 * sEval1[idx] / 512 + 575 * sEval2[idx] / 13824 + 44275 * sEval3[idx] / 110592 + 253 * sEval4[idx] / 4096) * dT;
		}
	}
}

void DifferentialEquationFM::RunRKCK45_Step5_withReductions(void)
{
	dmdt_reduction.new_minmax_reduction();
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//RKCK45 : 4th order evaluation
				pMesh->M[idx] = sM1[idx] + (2825 * sEval0[idx] / 27648 + 18575 * sEval2[idx] / 48384 + 13525 * sEval3[idx] / 55296 + 277 * sEval4[idx] / 14336 + rhs / 4) * dT;

				//Now calculate 5th order evaluation for adaptive time step
				DBL3 prediction = sM1[idx] + (37 * sEval0[idx] / 378 + 250 * sEval2[idx] / 621 + 125 * sEval3[idx] / 594 + 512 * rhs / 1771) * dT;

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

void DifferentialEquationFM::RunRKCK45_Step5(void)
{
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//RKCK45 : 4th order evaluation
				pMesh->M[idx] = sM1[idx] + (2825 * sEval0[idx] / 27648 + 18575 * sEval2[idx] / 48384 + 13525 * sEval3[idx] / 55296 + 277 * sEval4[idx] / 14336 + rhs / 4) * dT;

				//Now calculate 5th order evaluation for adaptive time step
				DBL3 prediction = sM1[idx] + (37 * sEval0[idx] / 378 + 250 * sEval2[idx] / 621 + 125 * sEval3[idx] / 594 + 512 * rhs / 1771) * dT;

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