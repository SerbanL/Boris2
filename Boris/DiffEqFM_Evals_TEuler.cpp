#include "stdafx.h"
#include "DiffEqFM.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC
#ifdef ODE_EVAL_COMPILATION_TEULER

#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//--------------------------------------------- TRAPEZOIDAL EULER

void DifferentialEquationFM::RunTEuler_Step0_withReductions(void)
{
	mxh_av_reduction.new_average_reduction();

	//Trapezoidal Euler can be used for stochastic equations - generate thermal VECs at the start
	if (H_Thermal.linear_size() && Torque_Thermal.linear_size()) GenerateThermalField_and_Torque();
	else if (H_Thermal.linear_size()) GenerateThermalField();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			sM1[idx] = pMesh->M[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//obtained average normalized torque term
				double Mnorm = pMesh->M[idx].norm();
				mxh_av_reduction.reduce_average((pMesh->M[idx] ^ pMesh->Heff[idx]) / (Mnorm * Mnorm));

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//Now estimate magnetization for the next time step
				pMesh->M[idx] += rhs * dT;
			}
		}
	}

	//magnitude of average mxh torque, set in mxh_reduction.max as this will be used to set the mxh value in ODECommon
	if (pMesh->grel.get0()) {

		//only reduce for mxh if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
		mxh_reduction.max = GetMagnitude(mxh_av_reduction.average());
	}
	else mxh_reduction.max = 0.0;
}

void DifferentialEquationFM::RunTEuler_Step0(void)
{
	//Trapezoidal Euler can be used for stochastic equations - generate thermal VECs at the start
	if (H_Thermal.linear_size() && Torque_Thermal.linear_size()) GenerateThermalField_and_Torque();
	else if (H_Thermal.linear_size()) GenerateThermalField();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization for the next step
			sM1[idx] = pMesh->M[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//Now estimate magnetization for the next time step
				pMesh->M[idx] += rhs * dT;
			}
		}
	}
}

void DifferentialEquationFM::RunTEuler_Step1_withReductions(void)
{
	dmdt_av_reduction.new_average_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//Now estimate magnetization using the second trapezoidal Euler step equation
				pMesh->M[idx] = (sM1[idx] + pMesh->M[idx] + rhs * dT) / 2;

				if (renormalize) {

					double Ms = pMesh->Ms;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
					pMesh->M[idx].renormalize(Ms);
				}

				//obtained average dmdt term
				double Mnorm = pMesh->M[idx].norm();
				dmdt_av_reduction.reduce_average((pMesh->M[idx] - sM1[idx]) / (dT * GAMMA * Mnorm * Mnorm));
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
		dmdt_reduction.max = GetMagnitude(dmdt_av_reduction.average());
	}
	else {

		dmdt_reduction.max = 0.0;
	}
}

void DifferentialEquationFM::RunTEuler_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//Now estimate magnetization using the second trapezoidal Euler step equation
				pMesh->M[idx] = (sM1[idx] + pMesh->M[idx] + rhs * dT) / 2;

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