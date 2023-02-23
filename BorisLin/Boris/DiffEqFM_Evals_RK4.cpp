#include "stdafx.h"
#include "DiffEqFM.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC
#ifdef ODE_EVAL_COMPILATION_RK4

#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//--------------------------------------------- RK4

void DifferentialEquationFM::RunRK4_Step0_withReductions(void)
{
	bool stochastic = H_Thermal.linear_size() != 0;

	//RK4 can be used for stochastic equations - generate thermal VECs at the start
	if (H_Thermal.linear_size() && Torque_Thermal.linear_size()) GenerateThermalField_and_Torque();
	else if (H_Thermal.linear_size()) GenerateThermalField();

	if (stochastic) {

		//Stochastic : reduce for average mxh

		mxh_av_reduction.new_average_reduction();

#pragma omp parallel for
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				//Save current magnetization for later use
				sM1[idx] = pMesh->M[idx];

				if (!pMesh->M.is_skipcell(idx)) {

					//obtained maximum normalized torque term
					double Mnorm = pMesh->M[idx].norm();
					mxh_av_reduction.reduce_average((pMesh->M[idx] ^ pMesh->Heff[idx]) / (Mnorm * Mnorm));

					//First evaluate RHS of set equation at the current time step
					sEval0[idx] = CALLFP(this, equation)(idx);

					//Now estimate magnetization using RK4 midle step
					pMesh->M[idx] += sEval0[idx] * (dT / 2);
				}
			}
		}

		//magnitude of average mxh torque, set in mxh_reduction.max as this will be used to set the mxh value in ODECommon
		if (pMesh->grel.get0()) {

			//only reduce for mxh if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
			mxh_reduction.max = GetMagnitude(mxh_av_reduction.average());
		}
		else {

			mxh_reduction.max = 0.0;
		}
	}
	else {

		//Not stochastic : reduce for maximum mxh

		mxh_reduction.new_minmax_reduction();

#pragma omp parallel for
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				//Save current magnetization for later use
				sM1[idx] = pMesh->M[idx];

				if (!pMesh->M.is_skipcell(idx)) {

					//obtained maximum normalized torque term
					double Mnorm = pMesh->M[idx].norm();
					double _mxh = GetMagnitude(pMesh->M[idx] ^ pMesh->Heff[idx]) / (Mnorm * Mnorm);
					mxh_reduction.reduce_max(_mxh);

					//First evaluate RHS of set equation at the current time step
					sEval0[idx] = CALLFP(this, equation)(idx);

					//Now estimate magnetization using RK4 midle step
					pMesh->M[idx] += sEval0[idx] * (dT / 2);
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
}

void DifferentialEquationFM::RunRK4_Step0(void)
{
	//RK4 can be used for stochastic equations - generate thermal VECs at the start
	if (H_Thermal.linear_size() && Torque_Thermal.linear_size()) GenerateThermalField_and_Torque();
	else if (H_Thermal.linear_size()) GenerateThermalField();

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

void DifferentialEquationFM::RunRK4_Step1(void)
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

void DifferentialEquationFM::RunRK4_Step2(void)
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

void DifferentialEquationFM::RunRK4_Step3_withReductions(void)
{
	bool stochastic = H_Thermal.linear_size() != 0;

	if (stochastic) {

		//Stochastic : reduce for average dmdt

		dmdt_av_reduction.new_average_reduction();

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

					//obtained maximum dmdt term
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
	else {

		//Not stochastic : reduce for maximum dmdt

		dmdt_reduction.new_minmax_reduction();

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
}

void DifferentialEquationFM::RunRK4_Step3(void)
{
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