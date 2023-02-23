#include "stdafx.h"
#include "DiffEqAFM.h"

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC
#ifdef ODE_EVAL_COMPILATION_RK4

#include "Mesh_AntiFerromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//--------------------------------------------- RK4

void DifferentialEquationAFM::RunRK4_Step0_withReductions(void)
{
	bool stochastic = H_Thermal.linear_size() != 0;

	//RK4 can be used for stochastic equations - generate thermal VECs at the start
	if (H_Thermal.linear_size() && Torque_Thermal.linear_size()) GenerateThermalField_and_Torque();
	else if (H_Thermal.linear_size()) GenerateThermalField();

	if (stochastic) {

		mxh_av_reduction.new_average_reduction();

#pragma omp parallel for
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				//Save current magnetization for later use
				sM1[idx] = pMesh->M[idx];
				sM1_2[idx] = pMesh->M2[idx];

				if (!pMesh->M.is_skipcell(idx)) {

					//obtained maximum normalized torque term
					double Mnorm = pMesh->M[idx].norm();
					mxh_av_reduction.reduce_average((pMesh->M[idx] ^ pMesh->Heff[idx]) / (Mnorm * Mnorm));

					//First evaluate RHS of set equation at the current time step
					sEval0[idx] = CALLFP(this, equation)(idx);
					sEval0_2[idx] = Equation_Eval_2[omp_get_thread_num()];

					//Now estimate magnetization using RK4 midle step
					pMesh->M[idx] += sEval0[idx] * (dT / 2);
					pMesh->M2[idx] += sEval0_2[idx] * (dT / 2);
				}
			}
		}

		if (pMesh->grel_AFM.get0().i) {

			//only reduce for mxh if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
			mxh_reduction.max = GetMagnitude(mxh_av_reduction.average());
		}
		else {

			mxh_reduction.max = 0.0;
		}
	}
	else {

		mxh_reduction.new_minmax_reduction();

#pragma omp parallel for
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				//Save current magnetization for later use
				sM1[idx] = pMesh->M[idx];
				sM1_2[idx] = pMesh->M2[idx];

				if (!pMesh->M.is_skipcell(idx)) {

					//obtained maximum normalized torque term
					double Mnorm = pMesh->M[idx].norm();
					double _mxh = GetMagnitude(pMesh->M[idx] ^ pMesh->Heff[idx]) / (Mnorm * Mnorm);
					mxh_reduction.reduce_max(_mxh);

					//First evaluate RHS of set equation at the current time step
					sEval0[idx] = CALLFP(this, equation)(idx);
					sEval0_2[idx] = Equation_Eval_2[omp_get_thread_num()];

					//Now estimate magnetization using RK4 midle step
					pMesh->M[idx] += sEval0[idx] * (dT / 2);
					pMesh->M2[idx] += sEval0_2[idx] * (dT / 2);
				}
			}
		}

		if (pMesh->grel_AFM.get0().i) {

			//only reduce for mxh if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
			mxh_reduction.maximum();
		}
		else {

			mxh_reduction.max = 0.0;
		}
	}
}

void DifferentialEquationAFM::RunRK4_Step0(void)
{
	//RK4 can be used for stochastic equations - generate thermal VECs at the start
	if (H_Thermal.linear_size() && Torque_Thermal.linear_size()) GenerateThermalField_and_Torque();
	else if (H_Thermal.linear_size()) GenerateThermalField();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization for later use
			sM1[idx] = pMesh->M[idx];
			sM1_2[idx] = pMesh->M2[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				sEval0[idx] = CALLFP(this, equation)(idx);
				sEval0_2[idx] = Equation_Eval_2[omp_get_thread_num()];

				//Now estimate magnetization using RK4 midle step
				pMesh->M[idx] += sEval0[idx] * (dT / 2);
				pMesh->M2[idx] += sEval0_2[idx] * (dT / 2);
			}
		}
	}
}

void DifferentialEquationAFM::RunRK4_Step1(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval1[idx] = CALLFP(this, equation)(idx);
			sEval1_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RK4 midle step
			pMesh->M[idx] = sM1[idx] + sEval1[idx] * (dT / 2);
			pMesh->M2[idx] = sM1_2[idx] + sEval1_2[idx] * (dT / 2);
		}
	}
}

void DifferentialEquationAFM::RunRK4_Step2(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx) && !pMesh->M.is_skipcell(idx)) {

			//First evaluate RHS of set equation at the current time step
			sEval2[idx] = CALLFP(this, equation)(idx);
			sEval2_2[idx] = Equation_Eval_2[omp_get_thread_num()];

			//Now estimate magnetization using RK4 last step
			pMesh->M[idx] = sM1[idx] + sEval2[idx] * dT;
			pMesh->M2[idx] = sM1_2[idx] + sEval2_2[idx] * dT;
		}
	}
}

void DifferentialEquationAFM::RunRK4_Step3_withReductions(void)
{
	bool stochastic = H_Thermal.linear_size() != 0;

	if (stochastic) {

		dmdt_av_reduction.new_average_reduction();

#pragma omp parallel for
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				if (!pMesh->M.is_skipcell(idx)) {

					//First evaluate RHS of set equation at the current time step
					DBL3 rhs = CALLFP(this, equation)(idx);
					DBL3 rhs_2 = Equation_Eval_2[omp_get_thread_num()];

					//Now estimate magnetization using previous RK4 evaluations
					pMesh->M[idx] = sM1[idx] + (sEval0[idx] + 2 * sEval1[idx] + 2 * sEval2[idx] + rhs) * (dT / 6);
					pMesh->M2[idx] = sM1_2[idx] + (sEval0_2[idx] + 2 * sEval1_2[idx] + 2 * sEval2_2[idx] + rhs_2) * (dT / 6);

					if (renormalize) {

						DBL2 Ms_AFM = pMesh->Ms_AFM;
						pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
						pMesh->M[idx].renormalize(Ms_AFM.i);
						pMesh->M2[idx].renormalize(Ms_AFM.j);
					}

					//obtained maximum dmdt term
					double Mnorm = pMesh->M[idx].norm();
					dmdt_av_reduction.reduce_average((pMesh->M[idx] - sM1[idx]) / (dT * GAMMA * Mnorm * Mnorm));
				}
				else {

					DBL2 Ms_AFM = pMesh->Ms_AFM;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
					pMesh->M[idx].renormalize(Ms_AFM.i);
					pMesh->M2[idx].renormalize(Ms_AFM.j);
				}
			}
		}

		if (pMesh->grel_AFM.get0().i) {

			//only reduce for dmdt if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
			dmdt_reduction.max = GetMagnitude(dmdt_av_reduction.average());
		}
		else {

			dmdt_reduction.max = 0.0;
		}
	}
	else {

		dmdt_reduction.new_minmax_reduction();

#pragma omp parallel for
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				if (!pMesh->M.is_skipcell(idx)) {

					//First evaluate RHS of set equation at the current time step
					DBL3 rhs = CALLFP(this, equation)(idx);
					DBL3 rhs_2 = Equation_Eval_2[omp_get_thread_num()];

					//Now estimate magnetization using previous RK4 evaluations
					pMesh->M[idx] = sM1[idx] + (sEval0[idx] + 2 * sEval1[idx] + 2 * sEval2[idx] + rhs) * (dT / 6);
					pMesh->M2[idx] = sM1_2[idx] + (sEval0_2[idx] + 2 * sEval1_2[idx] + 2 * sEval2_2[idx] + rhs_2) * (dT / 6);

					if (renormalize) {

						DBL2 Ms_AFM = pMesh->Ms_AFM;
						pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
						pMesh->M[idx].renormalize(Ms_AFM.i);
						pMesh->M2[idx].renormalize(Ms_AFM.j);
					}

					//obtained maximum dmdt term
					double Mnorm = pMesh->M[idx].norm();
					double _dmdt = GetMagnitude(pMesh->M[idx] - sM1[idx]) / (dT * GAMMA * Mnorm * Mnorm);
					dmdt_reduction.reduce_max(_dmdt);
				}
				else {

					DBL2 Ms_AFM = pMesh->Ms_AFM;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
					pMesh->M[idx].renormalize(Ms_AFM.i);
					pMesh->M2[idx].renormalize(Ms_AFM.j);
				}
			}
		}

		if (pMesh->grel_AFM.get0().i) {

			//only reduce for dmdt if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
			dmdt_reduction.maximum();
		}
		else {

			dmdt_reduction.max = 0.0;
		}
	}
}

void DifferentialEquationAFM::RunRK4_Step3(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);
				DBL3 rhs_2 = Equation_Eval_2[omp_get_thread_num()];

				//Now estimate magnetization using previous RK4 evaluations
				pMesh->M[idx] = sM1[idx] + (sEval0[idx] + 2 * sEval1[idx] + 2 * sEval2[idx] + rhs) * (dT / 6);
				pMesh->M2[idx] = sM1_2[idx] + (sEval0_2[idx] + 2 * sEval1_2[idx] + 2 * sEval2_2[idx] + rhs_2) * (dT / 6);

				if (renormalize) {

					DBL2 Ms_AFM = pMesh->Ms_AFM;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
					pMesh->M[idx].renormalize(Ms_AFM.i);
					pMesh->M2[idx].renormalize(Ms_AFM.j);
				}
			}
			else {

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
				pMesh->M[idx].renormalize(Ms_AFM.i);
				pMesh->M2[idx].renormalize(Ms_AFM.j);
			}
		}
	}
}

#endif
#endif