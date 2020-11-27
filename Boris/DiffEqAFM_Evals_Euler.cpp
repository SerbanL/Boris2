#include "stdafx.h"
#include "DiffEqAFM.h"

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC
#ifdef ODE_EVAL_COMPILATION_EULER

#include "Mesh_AntiFerromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//--------------------------------------------- EULER

void DifferentialEquationAFM::RunEuler_withReductions(void)
{
	mxh_av_reduction.new_average_reduction();
	dmdt_av_reduction.new_average_reduction();

	//Euler can be used for stochastic equations
	if (H_Thermal.linear_size() && Torque_Thermal.linear_size()) GenerateThermalField_and_Torque();
	else if (H_Thermal.linear_size()) GenerateThermalField();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization
			sM1[idx] = pMesh->M[idx];
			sM1_2[idx] = pMesh->M2[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//obtained average normalized torque term
				double Mnorm = pMesh->M[idx].norm();
				mxh_av_reduction.reduce_average((pMesh->M[idx] ^ pMesh->Heff[idx]) / (Mnorm * Mnorm));

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);
				DBL3 rhs_2 = Equation_Eval_2[omp_get_thread_num()];

				//Now estimate magnetization for the next time step
				pMesh->M[idx] += rhs * dT;
				pMesh->M2[idx] += rhs_2 * dT;

				if (renormalize) {

					DBL2 Ms_AFM = pMesh->Ms_AFM;
					pMesh->update_parameters_mcoarse(idx, pMesh->Ms_AFM, Ms_AFM);
					pMesh->M[idx].renormalize(Ms_AFM.i);
					pMesh->M2[idx].renormalize(Ms_AFM.j);
				}

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

	//magnitude of average mxh torque, set in mxh_reduction.max as this will be used to set the mxh value in ODECommon
	if (pMesh->grel.get0()) {

		//only reduce for mxh if grel is not zero (if it's zero this means magnetization dynamics are disabled in this mesh)
		mxh_reduction.max = GetMagnitude(mxh_av_reduction.average());
		dmdt_reduction.max = GetMagnitude(dmdt_av_reduction.average());
	}
	else {

		mxh_reduction.max = 0.0;
		dmdt_reduction.max = 0.0;
	}
}

void DifferentialEquationAFM::RunEuler(void)
{
	//Euler can be used for stochastic equations
	if (H_Thermal.linear_size() && Torque_Thermal.linear_size()) GenerateThermalField_and_Torque();
	else if (H_Thermal.linear_size()) GenerateThermalField();

#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			//Save current magnetization
			sM1[idx] = pMesh->M[idx];
			sM1_2[idx] = pMesh->M2[idx];

			if (!pMesh->M.is_skipcell(idx)) {

				//First evaluate RHS of set equation at the current time step
				DBL3 rhs = CALLFP(this, equation)(idx);
				DBL3 rhs_2 = Equation_Eval_2[omp_get_thread_num()];

				//Now estimate magnetization for the next time step
				pMesh->M[idx] += rhs * dT;
				pMesh->M2[idx] += rhs_2 * dT;

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