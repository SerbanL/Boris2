#include "stdafx.h"
#include "Atom_DiffEqCubic.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC
#ifdef ODE_EVAL_COMPILATION_EULER

#include "Atom_Mesh_Cubic.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

//--------------------------------------------- EULER

void Atom_DifferentialEquationCubic::RunEuler_withReductions(void)
{
	mxh_av_reduction.new_average_reduction();
	dmdt_av_reduction.new_average_reduction();

	//Euler can be used for stochastic equations
	if (H_Thermal.linear_size()) GenerateThermalField();

	//multiplicative conversion factor from atomic moment (units of muB) to A/m
	double conversion = MUB / paMesh->h.dim();

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			//Save current moment
			sM1[idx] = paMesh->M1[idx];

			if (!paMesh->M1.is_skipcell(idx)) {

				//obtained average normalized torque term
				double Mnorm = paMesh->M1[idx].norm();
				mxh_av_reduction.reduce_average((paMesh->M1[idx] ^ paMesh->Heff1[idx]) / (conversion * Mnorm * Mnorm));

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//Now estimate moment for the next time step
				paMesh->M1[idx] += rhs * dT;

				if (renormalize) {

					double mu_s = paMesh->mu_s;
					paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);
					paMesh->M1[idx].renormalize(mu_s);
				}

				dmdt_av_reduction.reduce_average((paMesh->M1[idx] - sM1[idx]) / (dT * GAMMA * Mnorm * conversion * Mnorm));
			}
		}
	}

	if (paMesh->grel.get0()) {

		mxh_reduction.max = GetMagnitude(mxh_av_reduction.average());
		dmdt_reduction.max = GetMagnitude(dmdt_av_reduction.average());
	}
	else {

		mxh_reduction.max = 0.0;
		dmdt_reduction.max = 0.0;
	}
}

void Atom_DifferentialEquationCubic::RunEuler(void)
{
	//Euler can be used for stochastic equations
	if (H_Thermal.linear_size()) GenerateThermalField();

#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++) {

		if (paMesh->M1.is_not_empty(idx)) {

			//Save current moment
			sM1[idx] = paMesh->M1[idx];

			if (!paMesh->M1.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//Now estimate moment for the next time step
				paMesh->M1[idx] += rhs * dT;

				if (renormalize) {

					double mu_s = paMesh->mu_s;
					paMesh->update_parameters_mcoarse(idx, paMesh->mu_s, mu_s);
					paMesh->M1[idx].renormalize(mu_s);
				}
			}
		}
	}
}

#endif
#endif