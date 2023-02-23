#include "stdafx.h"
#include "Atom_DiffEqCubic.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC
#ifdef ODE_EVAL_COMPILATION_ABM

#include "Atom_Mesh_Cubic.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

//--------------------------------------------- ADAMS-BASHFORTH-MOULTON

void Atom_DifferentialEquationCubic::RunABM_Predictor_withReductions(void)
{
	mxh_reduction.new_minmax_reduction();

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	double conversion = MUB / paMesh->h.dim();

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			//Save current moment for the next step
			sM1[idx] = paMesh->M1[idx];

			if (!paMesh->M1.is_skipcell(idx)) {

				//obtained maximum normalized torque term
				double Mnorm = paMesh->M1[idx].norm();
				double _mxh = GetMagnitude(paMesh->M1[idx] ^ paMesh->Heff1[idx]) / (conversion * Mnorm * Mnorm);
				mxh_reduction.reduce_max(_mxh);

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//ABM predictor : pk+1 = mk + (dt/2) * (3*fk - fk-1)
				if (alternator) {

					paMesh->M1[idx] += dT * (3 * rhs - sEval0[idx]) / 2;
					sEval1[idx] = rhs;
				}
				else {

					paMesh->M1[idx] += dT * (3 * rhs - sEval1[idx]) / 2;
					sEval0[idx] = rhs;
				}
			}
		}
	}

	if (paMesh->grel.get0()) {

		//only reduce if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
		mxh_reduction.maximum();
	}
	else {

		mxh_reduction.max = 0.0;
	}
}

void Atom_DifferentialEquationCubic::RunABM_Predictor(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			//Save current moment for the next step
			sM1[idx] = paMesh->M1[idx];

			if (!paMesh->M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//ABM predictor : pk+1 = mk + (dt/2) * (3*fk - fk-1)
				if (alternator) {

					paMesh->M1[idx] += dT * (3 * rhs - sEval0[idx]) / 2;
					sEval1[idx] = rhs;
				}
				else {

					paMesh->M1[idx] += dT * (3 * rhs - sEval1[idx]) / 2;
					sEval0[idx] = rhs;
				}
			}
		}
	}
}

void Atom_DifferentialEquationCubic::RunABM_Corrector_withReductions(void)
{
	dmdt_reduction.new_minmax_reduction();
	lte_reduction.new_minmax_reduction();

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	double conversion = MUB / paMesh->h.dim();

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			if (!paMesh->M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//First save predicted moment for lte calculation
				DBL3 saveM = paMesh->M1[idx];

				//ABM corrector : mk+1 = mk + (dt/2) * (fk+1 + fk)
				if (alternator) {

					paMesh->M1[idx] = sM1[idx] + dT * (rhs + sEval1[idx]) / 2;
				}
				else {

					paMesh->M1[idx] = sM1[idx] + dT * (rhs + sEval0[idx]) / 2;
				}

				if (renormalize) {

					double mu_s = paMesh->mu_s;
					paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);
					paMesh->M1[idx].renormalize(mu_s);
				}

				//obtained maximum dmdt term
				double Mnorm = paMesh->M1[idx].norm();
				double _dmdt = GetMagnitude(paMesh->M1[idx] - sM1[idx]) / (dT * GAMMA * Mnorm * conversion * Mnorm);
				dmdt_reduction.reduce_max(_dmdt);

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(paMesh->M1[idx] - saveM) / paMesh->M1[idx].norm();
				lte_reduction.reduce_max(_lte);
			}
		}
	}

	lte_reduction.maximum();

	if (paMesh->grel.get0()) {

		//only reduce for dmdt if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
		dmdt_reduction.maximum();
	}
	else {

		dmdt_reduction.max = 0.0;
	}
}

void Atom_DifferentialEquationCubic::RunABM_Corrector(void)
{
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			if (!paMesh->M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//First save predicted moment for lte calculation
				DBL3 saveM = paMesh->M1[idx];

				//ABM corrector : mk+1 = mk + (dt/2) * (fk+1 + fk)
				if (alternator) {

					paMesh->M1[idx] = sM1[idx] + dT * (rhs + sEval1[idx]) / 2;
				}
				else {

					paMesh->M1[idx] = sM1[idx] + dT * (rhs + sEval0[idx]) / 2;
				}

				if (renormalize) {

					double mu_s = paMesh->mu_s;
					paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);
					paMesh->M1[idx].renormalize(mu_s);
				}

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(paMesh->M1[idx] - saveM) / paMesh->M1[idx].norm();
				lte_reduction.reduce_max(_lte);
			}
		}
	}

	lte_reduction.maximum();
}

void Atom_DifferentialEquationCubic::RunABM_TEuler0(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			//Save current moment for the next step
			sM1[idx] = paMesh->M1[idx];

			if (!paMesh->M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval0[idx] = CALLFP(this, equation)(idx);

				//Now estimate moment for the next time step
				paMesh->M1[idx] += sEval0[idx] * dT;
			}
		}
	}
}

void Atom_DifferentialEquationCubic::RunABM_TEuler1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx) && !paMesh->M1.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			DBL3 rhs = CALLFP(this, equation)(idx);

			//Now estimate moment using the second trapezoidal Euler step equation
			paMesh->M1[idx] = (sM1[idx] + paMesh->M1[idx] + rhs * dT) / 2;
		}
	}
}

#endif
#endif