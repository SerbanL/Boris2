#include "stdafx.h"
#include "Atom_DiffEqCubic.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC
#ifdef ODE_EVAL_COMPILATION_AHEUN

#include "Atom_Mesh_Cubic.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

//--------------------------------------------- TRAPEZOIDAL EULER

void Atom_DifferentialEquationCubic::RunAHeun_Step0_withReductions(void)
{
	mxh_av_reduction.new_average_reduction();

	//Trapezoidal Euler can be used for stochastic equations - generate thermal VECs at the start
	if (H_Thermal.linear_size()) GenerateThermalField();

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	double conversion = MUB / paMesh->h.dim();

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			//Save current moment for the next step
			sM1[idx] = paMesh->M1[idx];

			if (!paMesh->M1.is_skipcell(idx)) {

				//obtained average normalized torque term
				double Mnorm = paMesh->M1[idx].norm();
				mxh_av_reduction.reduce_average((paMesh->M1[idx] ^ paMesh->Heff1[idx]) / (conversion * Mnorm * Mnorm));

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//Now estimate moment for the next time step
				paMesh->M1[idx] += rhs * dT;
			}
		}
	}

	//magnitude of average mxh torque, set in mxh_reduction.max as this will be used to set the mxh value in ODECommon
	if (paMesh->grel.get0()) {

		//only reduce for mxh if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
		mxh_reduction.max = GetMagnitude(mxh_av_reduction.average());
	}
	else mxh_reduction.max = 0.0;
}

void Atom_DifferentialEquationCubic::RunAHeun_Step0(void)
{
	//Trapezoidal Euler can be used for stochastic equations - generate thermal VECs at the start
	if (H_Thermal.linear_size()) GenerateThermalField();

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			//Save current moment for the next step
			sM1[idx] = paMesh->M1[idx];

			if (!paMesh->M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//Now estimate moment for the next time step
				paMesh->M1[idx] += rhs * dT;
			}
		}
	}
}

void Atom_DifferentialEquationCubic::RunAHeun_Step1_withReductions(void)
{
	dmdt_av_reduction.new_average_reduction();
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

				//Now estimate moment using the second trapezoidal Euler step equation
				paMesh->M1[idx] = (sM1[idx] + paMesh->M1[idx] + rhs * dT) / 2;

				if (renormalize) {

					double mu_s = paMesh->mu_s;
					paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);
					paMesh->M1[idx].renormalize(mu_s);
				}

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(paMesh->M1[idx] - saveM) / paMesh->M1[idx].norm();
				lte_reduction.reduce_max(_lte);

				//obtained average dmdt term
				double Mnorm = paMesh->M1[idx].norm();
				dmdt_av_reduction.reduce_average((paMesh->M1[idx] - sM1[idx]) / (dT * GAMMA * Mnorm * conversion * Mnorm));
			}
		}
	}

	lte_reduction.maximum();

	if (paMesh->grel.get0()) {

		//only reduce for dmdt if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
		dmdt_reduction.max = GetMagnitude(dmdt_av_reduction.average());
	}
	else {

		dmdt_reduction.max = 0.0;
	}
}

void Atom_DifferentialEquationCubic::RunAHeun_Step1(void)
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

				//Now estimate moment using the second trapezoidal Euler step equation
				paMesh->M1[idx] = (sM1[idx] + paMesh->M1[idx] + rhs * dT) / 2;

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

#endif
#endif