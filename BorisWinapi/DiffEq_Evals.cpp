#include "stdafx.h"
#include "DiffEq.h"
#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//---------------------------------------- OTHERS

//Restore magnetisation after a failed step for adaptive time-step methods
void DifferentialEquation::RestoreMagnetisation(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++)
		pMesh->M[idx] = sM1[idx];
}

//--------------------------------------------- EULER

void DifferentialEquation::RunEuler(void)
{
	//mxh_av_reduction.new_average_reduction();

	//Euler can be used for stochastic equations
	if (H_Thermal.linear_size() && Torque_Thermal.linear_size()) GenerateThermalField_and_Torque();
	else if (H_Thermal.linear_size()) GenerateThermalField();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization
			sM1[idx] = pMesh->M[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//Now estimate magnetization for the next time step
				pMesh->M[idx] += rhs * dT;

				if (renormalize) {

					double Ms = pMesh->Ms;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
					pMesh->M[idx].renormalize(Ms);
				}

				//obtained average normalized torque term
				mxh_av_reduction.reduce_average((pMesh->M[idx] ^ pMesh->Heff[idx]) / (pMesh->Ms * pMesh->Ms));
			}
			else {

				double Ms = pMesh->Ms;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
				pMesh->M[idx].renormalize(Ms);		//re-normalize the skipped cells no matter what - temperature can change
			}
		}
	}

	//magnitude of average mxh torque, set in mxh_reduction.max as this will be used to set the mxh value in ODECommon
	if (pMesh->grel.get0()) {

		//only reduce for mxh if grel is not zero (if it's zero this means magnetisaiton dynamics is disabled in this mesh)
		mxh_reduction.max = GetMagnitude(mxh_av_reduction.average());
	}
	else mxh_reduction.max = 0.0;
}

//--------------------------------------------- TRAPEZOIDAL EULER

void DifferentialEquation::RunTEuler_Step0(void)
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

void DifferentialEquation::RunTEuler_Step1(void)
{
	mxh_av_reduction.new_average_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//Now estimate magnetization using the second trapezoidal Euler step formula
				pMesh->M[idx] = (sM1[idx] + pMesh->M[idx] + rhs * dT) / 2;

				if (renormalize) {

					double Ms = pMesh->Ms;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
					pMesh->M[idx].renormalize(Ms);
				}

				//obtained average normalized torque term
				mxh_av_reduction.reduce_average((pMesh->M[idx] ^ pMesh->Heff[idx]) / (pMesh->Ms * pMesh->Ms));
			}
			else {

				double Ms = pMesh->Ms;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
				pMesh->M[idx].renormalize(Ms);
			}
		}
	}

	//magnitude of average mxh torque, set in mxh_reduction.max as this will be used to set the mxh value in ODECommon
	if (pMesh->grel.get0()) {

		mxh_reduction.max = GetMagnitude(mxh_av_reduction.average());
	}
	else mxh_reduction.max = 0.0;
}

//--------------------------------------------- RK4

void DifferentialEquation::RunRK4_Step0(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization for later use
			sM1[idx] = pMesh->M[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval0[idx] = CALLFP(this, equation)(idx);

				//Now estimate magnetization using RK4 midle step
				pMesh->M[idx] += sEval0[idx] * (dT / 2);
			}
		}
	}
}

void DifferentialEquation::RunRK4_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval1[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RK4 midle step
			pMesh->M[idx] = sM1[idx] + sEval1[idx] * (dT / 2);
		}
	}
}

void DifferentialEquation::RunRK4_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval2[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RK4 last step
			pMesh->M[idx] = sM1[idx] + sEval2[idx] * dT;
		}
	}
}

void DifferentialEquation::RunRK4_Step3(void)
{
	mxh_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//Now estimate magnetization using previous RK4 evaluations
				pMesh->M[idx] = sM1[idx] + (sEval0[idx] + 2 * sEval1[idx] + 2 * sEval2[idx] + rhs) * (dT / 6);

				if (renormalize) {

					double Ms = pMesh->Ms;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
					pMesh->M[idx].renormalize(Ms);
				}

				//obtained maximum normalized torque term
				double _mxh = GetMagnitude(pMesh->M[idx] ^ pMesh->Heff[idx]) / (pMesh->Ms * pMesh->Ms);

				mxh_reduction.reduce_max(_mxh);
			}
			else {

				double Ms = pMesh->Ms;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
				pMesh->M[idx].renormalize(Ms);
			}
		}
	}

	if (pMesh->grel.get0()) {

		mxh_reduction.maximum();
	}
	else mxh_reduction.max = 0.0;
}

//--------------------------------------------- ADAMS-BASHFORTH-MOULTON

void DifferentialEquation::RunABM_Predictor(void)
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

void DifferentialEquation::RunABM_Corrector(void)
{
	mxh_reduction.new_minmax_reduction();
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
				double _lte = GetMagnitude(pMesh->M[idx] - saveM) / pMesh->Ms;

				//obtained maximum normalized torque term
				double _mxh = GetMagnitude(pMesh->M[idx] ^ pMesh->Heff[idx]) / (pMesh->Ms * pMesh->Ms);

				mxh_reduction.reduce_max(_mxh);
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

		mxh_reduction.maximum();
		lte_reduction.maximum();
	}
	else {

		mxh_reduction.max = 0.0;
		lte_reduction.max = 0.0;
	}
}

void DifferentialEquation::RunABM_TEuler0(void)
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

void DifferentialEquation::RunABM_TEuler1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			DBL3 rhs = CALLFP(this, equation)(idx);

			//Now estimate magnetization using the second trapezoidal Euler step formula
			pMesh->M[idx] = (sM1[idx] + pMesh->M[idx] + rhs * dT) / 2;
		}
	}
}

//--------------------------------------------- RUNGE KUTTA FEHLBERG (4th order adaptive step, 5th order evaluation)

void DifferentialEquation::RunRKF45_Step0(void)
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
				pMesh->M[idx] += sEval0[idx] * (dT / 4);
			}
		}
	}
}

void DifferentialEquation::RunRKF45_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval1[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKF midle step 1
			pMesh->M[idx] = sM1[idx] + (3 * sEval0[idx] + 9 * sEval1[idx]) * dT / 32;
		}
	}
}

void DifferentialEquation::RunRKF45_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval2[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKF midle step 2
			pMesh->M[idx] = sM1[idx] + (1932 * sEval0[idx] - 7200 * sEval1[idx] + 7296 * sEval2[idx]) * dT / 2197;
		}
	}
}

void DifferentialEquation::RunRKF45_Step3(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval3[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKF midle step 3
			pMesh->M[idx] = sM1[idx] + (439 * sEval0[idx] / 216 - 8 * sEval1[idx] + 3680 * sEval2[idx] / 513 - 845 * sEval3[idx] / 4104) * dT;
		}
	}
}

void DifferentialEquation::RunRKF45_Step4(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval4[idx] = CALLFP(this, equation)(idx);

			//Now estimate magnetization using RKF midle step 4
			pMesh->M[idx] = sM1[idx] + (-8 * sEval0[idx] / 27 + 2 * sEval1[idx] - 3544 * sEval2[idx] / 2565 + 1859 * sEval3[idx] / 4104 - 11 * sEval4[idx] / 40) * dT;
		}
	}
}

void DifferentialEquation::RunRKF45_Step5(void)
{
	mxh_reduction.new_minmax_reduction();
	lte_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);

				//RKF45 : 4th order predictor for adaptive time step
				DBL3 prediction = sM1[idx] + (25 * sEval0[idx] / 216 + 1408 * sEval2[idx] / 2565 + 2197 * sEval3[idx] / 4101 - sEval4[idx] / 5) * dT;

				//Now calculate 5th order evaluation
				pMesh->M[idx] = sM1[idx] + (16 * sEval0[idx] / 135 + 6656 * sEval2[idx] / 12825 + 28561 * sEval3[idx] / 56430 - 9 * sEval4[idx] / 50 + 2 * rhs / 55) * dT;

				if (renormalize) {

					double Ms = pMesh->Ms;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
					pMesh->M[idx].renormalize(Ms);
				}

				//obtained maximum normalized torque term
				double _mxh = GetMagnitude(pMesh->M[idx] ^ pMesh->Heff[idx]) / (pMesh->Ms * pMesh->Ms);

				//local truncation error (between predicted and corrected)
				double _lte = GetMagnitude(pMesh->M[idx] - prediction) / pMesh->Ms;

				mxh_reduction.reduce_max(_mxh);
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

		mxh_reduction.maximum();
		lte_reduction.maximum();
	}
	else {

		mxh_reduction.max = 0.0;
		lte_reduction.max = 0.0;
	}
}
