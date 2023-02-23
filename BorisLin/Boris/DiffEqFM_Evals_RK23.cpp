#include "stdafx.h"
#include "DiffEqFM.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC
#ifdef ODE_EVAL_COMPILATION_RK23

#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//--------------------------------------------- RUNGE KUTTA 23 (Bogacki-Shampine) (2nd order adaptive step with FSAL, 3rd order evaluation)

void DifferentialEquationFM::RunRK23_Step0_withReductions(void)
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

				//2nd order evaluation for adaptive step
				DBL3 prediction = sM1[idx] + (7 * sEval0[idx] / 24 + 1 * sEval1[idx] / 4 + 1 * sEval2[idx] / 3 + 1 * rhs / 8) * dT;

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(pMesh->M[idx] - prediction) / Mnorm;
				lte_reduction.reduce_max(_lte);

				//save evaluation for later use
				sEval0[idx] = rhs;
			}
		}
	}

	lte_reduction.maximum();

	if (pMesh->grel.get0()) {

		//only reduce for mxh if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
		mxh_reduction.maximum();
	}
	else {

		mxh_reduction.max = 0.0;
	}
}

void DifferentialEquationFM::RunRK23_Step0(void)
{
	//lte reductions needed for adaptive time step
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//2nd order evaluation for adaptive step
				DBL3 prediction = sM1[idx] + (7 * sEval0[idx] / 24 + 1 * sEval1[idx] / 4 + 1 * sEval2[idx] / 3 + 1 * rhs / 8) * dT;

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

void DifferentialEquationFM::RunRK23_Step0_Advance(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//Save current magnetization for later use
				sM1[idx] = pMesh->M[idx];

				//Now estimate magnetization using RK23 first step
				pMesh->M[idx] += sEval0[idx] * (dT / 2);
			}
		}
	}
}

void DifferentialEquationFM::RunRK23_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval1[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RK23 midle step 1
			pMesh->M[idx] = sM1[idx] + 3 * sEval1[idx] * dT / 4;
		}
	}
}

void DifferentialEquationFM::RunRK23_Step2_withReductions(void)
{
	dmdt_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval2[idx] = CALLFP(this, equation)(idx);

				//Now calculate 3rd order evaluation
				pMesh->M[idx] = sM1[idx] + (2 * sEval0[idx] / 9 + 1 * sEval1[idx] / 3 + 4 * sEval2[idx] / 9) * dT;

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

void DifferentialEquationFM::RunRK23_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval2[idx] = CALLFP(this, equation)(idx);

				//Now calculate 3rd order evaluation
				pMesh->M[idx] = sM1[idx] + (2 * sEval0[idx] / 9 + 1 * sEval1[idx] / 3 + 4 * sEval2[idx] / 9) * dT;

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