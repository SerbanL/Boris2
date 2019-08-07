#include "stdafx.h"
#include "DiffEq.h"
#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

#ifdef ODE_EVAL_TEULER

//--------------------------------------------- TRAPEZOIDAL EULER

void DifferentialEquation::RunAHeun_Step0_withReductions(void)
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
				if (H_Thermal.linear_size()) {

					mxh_av_reduction.reduce_average((pMesh->M[idx] ^ (pMesh->Heff[idx] + H_Thermal[idx])) / (Mnorm * Mnorm));
				}
				else mxh_av_reduction.reduce_average((pMesh->M[idx] ^ pMesh->Heff[idx]) / (Mnorm * Mnorm));

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//Now estimate magnetization for the next time step
				pMesh->M[idx] += rhs * dT;
			}
		}
	}

	//magnitude of average mxh torque, set in mxh_reduction.max as this will be used to set the mxh value in ODECommon
	if (pMesh->grel.get0()) {

		//only reduce for mxh if grel is not zero (if it's zero this means magnetisation dynamics are disabled in this mesh)
		mxh_reduction.max = GetMagnitude(mxh_av_reduction.average());
	}
	else mxh_reduction.max = 0.0;
}

void DifferentialEquation::RunAHeun_Step0(void)
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

void DifferentialEquation::RunAHeun_Step1_withReductions(void)
{
	dmdt_av_reduction.new_average_reduction();
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//First save predicted magnetization for lte calculation
				DBL3 saveM = pMesh->M[idx];

				//Now estimate magnetization using the second trapezoidal Euler step formula
				pMesh->M[idx] = (sM1[idx] + pMesh->M[idx] + rhs * dT) / 2;

				if (renormalize) {

					double Ms = pMesh->Ms;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
					pMesh->M[idx].renormalize(Ms);
				}

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(pMesh->M[idx] - saveM) / pMesh->M[idx].norm();
				lte_reduction.reduce_max(_lte);

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

	lte_reduction.maximum();

	if (pMesh->grel.get0()) {

		//only reduce for dmdt if grel is not zero (if it's zero this means magnetisation dynamics are disabled in this mesh)
		dmdt_reduction.max = GetMagnitude(dmdt_av_reduction.average());
	}
	else {

		dmdt_reduction.max = 0.0;
	}
}

void DifferentialEquation::RunAHeun_Step1(void)
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

				//Now estimate magnetization using the second trapezoidal Euler step formula
				pMesh->M[idx] = (sM1[idx] + pMesh->M[idx] + rhs * dT) / 2;

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

#endif