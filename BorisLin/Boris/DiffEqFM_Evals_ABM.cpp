#include "stdafx.h"
#include "DiffEqFM.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC
#ifdef ODE_EVAL_COMPILATION_ABM

#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//--------------------------------------------- ADAMS-BASHFORTH-MOULTON

void DifferentialEquationFM::RunABM_Predictor_withReductions(void)
{
	mxh_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			sM1[idx] = pMesh->M[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//obtained maximum normalized torque term
				double Mnorm = pMesh->M[idx].norm();
				double _mxh = GetMagnitude(pMesh->M[idx] ^ pMesh->Heff[idx]) / (Mnorm * Mnorm);
				mxh_reduction.reduce_max(_mxh);

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//ABM predictor : pk+1 = mk + (dt/2) * (3*fk - fk-1)
				if (alternator) {

					pMesh->M[idx] += dT * (3 * rhs - sEval0[idx]) / 2;
					sEval1[idx] = rhs;
				}
				else {

					pMesh->M[idx] += dT * (3 * rhs - sEval1[idx]) / 2;
					sEval0[idx] = rhs;
				}
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

void DifferentialEquationFM::RunABM_Predictor(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			sM1[idx] = pMesh->M[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//ABM predictor : pk+1 = mk + (dt/2) * (3*fk - fk-1)
				if (alternator) {

					pMesh->M[idx] += dT * (3 * rhs - sEval0[idx]) / 2;
					sEval1[idx] = rhs;
				}
				else {

					pMesh->M[idx] += dT * (3 * rhs - sEval1[idx]) / 2;
					sEval0[idx] = rhs;
				}
			}
		}
	}
}

void DifferentialEquationFM::RunABM_Corrector_withReductions(void)
{
	dmdt_reduction.new_minmax_reduction();
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//First save predicted magnetization for lte calculation
				DBL3 saveM = pMesh->M[idx];

				//ABM corrector : mk+1 = mk + (dt/2) * (fk+1 + fk)
				if (alternator) {

					pMesh->M[idx] = sM1[idx] + dT * (rhs + sEval1[idx]) / 2;
				}
				else {

					pMesh->M[idx] = sM1[idx] + dT * (rhs + sEval0[idx]) / 2;
				}

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
				double _lte = GetMagnitude(pMesh->M[idx] - saveM) / pMesh->M[idx].norm();
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

	if (pMesh->grel.get0()) {

		//only reduce for dmdt if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
		dmdt_reduction.maximum();
	}
	else {

		dmdt_reduction.max = 0.0;
	}
}

void DifferentialEquationFM::RunABM_Corrector(void)
{
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//First save predicted magnetization for lte calculation
				DBL3 saveM = pMesh->M[idx];

				//ABM corrector : mk+1 = mk + (dt/2) * (fk+1 + fk)
				if (alternator) {

					pMesh->M[idx] = sM1[idx] + dT * (rhs + sEval1[idx]) / 2;
				}
				else {

					pMesh->M[idx] = sM1[idx] + dT * (rhs + sEval0[idx]) / 2;
				}

				if (renormalize) {

					double Ms = pMesh->Ms;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
					pMesh->M[idx].renormalize(Ms);
				}

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(pMesh->M[idx] - saveM) / pMesh->M[idx].norm();
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

void DifferentialEquationFM::RunABM_TEuler0(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			sM1[idx] = pMesh->M[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval0[idx] = CALLFP(this, equation)(idx);

				//Now estimate magnetization for the next time step
				pMesh->M[idx] += sEval0[idx] * dT;
			}
		}
	}
}

void DifferentialEquationFM::RunABM_TEuler1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			DBL3 rhs = CALLFP(this, equation)(idx);

			//Now estimate magnetization using the second trapezoidal Euler step equation
			pMesh->M[idx] = (sM1[idx] + pMesh->M[idx] + rhs * dT) / 2;
		}
	}
}

#endif
#endif