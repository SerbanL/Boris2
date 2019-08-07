#include "stdafx.h"
#include "DiffEq.h"
#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"

//---------------------------------------- ITERATE METHODS

//calculate time step for adaptive methods based on current error values and evaluation method settings - in DiffEq_Iterate.cpp
bool ODECommon::SetAdaptiveTimeStep(void)
{
	double lte = Get_lte();

	//adaptive time step based on lte - is lte over acceptable relative error?
	if (lte > err_high) {

		//failed - try again. Restore state to start of evaluation and try again with smaller time step - output not available yet.
		//do not redo if time step is at or below the minimum value allowed - this is used as a "timeout" to stop the solver from getting stuck on the same iteration; in this case the method is not good enough for the problem.
		if (lte > err_high_fail && dT > dT_min) {

			time -= dT;
			stagetime -= dT;

			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RestoreMagnetisation();

			//reduce time step based on error ratio; use a 0.99 back-off factor so the solver doesn't get stuck
			dT *= 0.99 * pow(err_high_fail / lte, 0.25);

			//failed - must repeat
			return false;
		}

		//not failed but still need to reduce time step
		dT *= 0.99 * pow(err_high / lte, 0.25);

		//must not go below minimum time step
		if (dT < dT_min) dT = dT_min;
	}

	//is lte below minimum relative error ? If yes we can go quicker
	if (lte < err_low) {

		//increase by a small constant factor
		dT *= dT_increase;

		//must not go above maximum time step
		if (dT > dT_max) dT = dT_max;
	}

	//good, next step
	return true;
}

bool ODECommon::SetAdaptiveTimeStep_SingleThreshold(void)
{
	double lte = Get_lte();

	//failed - try again. Restore state to start of evaluation and try again with smaller time step - output not available yet.
	//do not redo if time step is at or below the minimum value allowed - this is used as a "timeout" to stop the solver from getting stuck on the same iteration; in this case the method is not good enough for the problem.
	if (lte > err_high_fail && dT > dT_min) {

		time -= dT;
		stagetime -= dT;

		for (int idx = 0; idx < (int)pODE.size(); idx++)
			pODE[idx]->RestoreMagnetisation();

		//reduce time step based on error ratio; use a 0.99 back-off factor so the solver doesn't get stuck
		dT *= 0.99 * pow(err_high_fail / lte, 1.0/3.0);

		//failed - must repeat
		return false;
	}

	//not failed, adjust time step
	dT *= sqrt(err_high_fail / lte);

	//must not go below minimum time step
	if (dT < dT_min) dT = dT_min;

	//must not go above maximum time step
	if (dT > dT_max) dT = dT_max;

	//good, next step
	return true;
}

//advance time using the set ODE solver
void ODECommon::Iterate(void)
{
	//save current dT value in case it changes (adaptive time step methods)
	dT_last = dT;

	switch (evalMethod) {

	case EVAL_EULER:
	{
#ifdef ODE_EVAL_EULER
		if (calculate_mxh || calculate_dmdt) {

			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunEuler_withReductions();

			if (calculate_mxh) Set_mxh();
			if (calculate_dmdt) Set_dmdt();

			calculate_mxh = false;
			calculate_dmdt = false;
		}
		else {

			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunEuler();
		}

		time += dT;
		stagetime += dT;
		iteration++;
		stageiteration++;
#endif
	}
	break;

	case EVAL_TEULER:
	{
#ifdef ODE_EVAL_TEULER
		if (evalStep == 0) {

			if (calculate_mxh) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunTEuler_Step0_withReductions();

				calculate_mxh = false;

				Set_mxh();
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunTEuler_Step0();
			}

			evalStep = 1;
			available = false;
		}
		else {

			if (calculate_dmdt) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunTEuler_Step1_withReductions();

				calculate_dmdt = false;

				Set_dmdt();
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunTEuler_Step1();
			}

			evalStep = 0;
			available = true;

			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;
		}
#endif
	}
	break;

	case EVAL_AHEUN:
	{
#ifdef ODE_EVAL_AHEUN
		if (evalStep == 0) {

			if (calculate_mxh) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunAHeun_Step0_withReductions();

				calculate_mxh = false;

				Set_mxh();
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunAHeun_Step0();
			}

			evalStep = 1;
			available = false;
		}
		else {

			if (calculate_dmdt) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunAHeun_Step1_withReductions();

				calculate_dmdt = false;

				Set_dmdt();
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunAHeun_Step1();
			}

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;

			//if you use the SingleThreshold Method with a stochastic equation you can end up backtracking very often -> slow
			//if (!SetAdaptiveTimeStep_SingleThreshold()) break;
			//this works better
			if (!SetAdaptiveTimeStep()) break;

			available = true;
		}
#endif
	}
	break;

	case EVAL_ABM:
	{
#ifdef ODE_EVAL_ABM
		if (primed) {

			if (evalStep == 0) {

				if (calculate_mxh) {

					for (int idx = 0; idx < (int)pODE.size(); idx++)
						pODE[idx]->RunABM_Predictor_withReductions();

					calculate_mxh = false;

					Set_mxh();
				}
				else {

					for (int idx = 0; idx < (int)pODE.size(); idx++)
						pODE[idx]->RunABM_Predictor();
				}

				evalStep = 1;
				available = false;
			}
			else {

				if (calculate_dmdt) {

					for (int idx = 0; idx < (int)pODE.size(); idx++)
						pODE[idx]->RunABM_Corrector_withReductions();

					calculate_dmdt = false;

					Set_dmdt();
				}
				else {

					for (int idx = 0; idx < (int)pODE.size(); idx++)
						pODE[idx]->RunABM_Corrector();
				}

				evalStep = 0;
				time += dT;
				stagetime += dT;
				iteration++;
				stageiteration++;

				if (!SetAdaptiveTimeStep()) break;

				//done for this step. Make it available.
				available = true;

				//switch between saved equation evaluations (use latest)
				if (alternator) alternator = false;
				else alternator = true;
			}
		}
		else {

			if (evalStep == 0) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunABM_TEuler0();

				evalStep = 1;
				available = false;
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunABM_TEuler1();

				evalStep = 0;
				available = true;

				time += dT;
				stagetime += dT;
				iteration++;
				stageiteration++;

				primed = true;
				alternator = false;		//evaluated equation is now in sEval0
			}
		}
#endif
	}
	break;

	case EVAL_RK23:
	{
#ifdef ODE_EVAL_RK23
		switch (evalStep) {

		case 0:
		{
			if (calculate_mxh) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRK23_Step0_withReductions();

				calculate_mxh = false;

				Set_mxh();
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRK23_Step0();
			}

			evalStep++;
			available = false;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRK23_Step1();

			evalStep++;
		}
		break;

		case 2:
		{
			if (calculate_dmdt) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRK23_Step2_withReductions();

				calculate_dmdt = false;

				Set_dmdt();
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRK23_Step2();
			}

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;

			if (primed) {

				if (!SetAdaptiveTimeStep()) break;
			}
			//need one pass of rk23 before we can use it for adaptive time step (FSAL)
			else primed = true;

			//done for this step. Make it available.
			available = true;
		}
		break;
		}
#endif
	}
	break;

	case EVAL_RK4:
	{
#ifdef ODE_EVAL_RK4
		switch (evalStep) {

		case 0:
		{
			if (calculate_mxh) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRK4_Step0_withReductions();

				calculate_mxh = false;

				Set_mxh();
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRK4_Step0();
			}

			evalStep++;
			available = false;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRK4_Step1();

			evalStep++;
		}
		break;

		case 2:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRK4_Step2();

			evalStep++;
		}
		break;

		case 3:
		{
			if (calculate_dmdt) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRK4_Step3_withReductions();

				calculate_dmdt = false;

				Set_dmdt();
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRK4_Step3();
			}

			evalStep = 0;
			available = true;

			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;
		}
		break;
		}
#endif
	}
	break;

	case EVAL_RKF:
	{
#ifdef ODE_EVAL_RKF
		switch (evalStep) {

		case 0:
		{
			if (calculate_mxh) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRKF45_Step0_withReductions();

				calculate_mxh = false;

				Set_mxh();
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRKF45_Step0();
			}

			evalStep++;
			available = false;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKF45_Step1();

			evalStep++;
		}
		break;

		case 2:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKF45_Step2();

			evalStep++;
		}
		break;

		case 3:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKF45_Step3();

			evalStep++;
		}
		break;

		case 4:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKF45_Step4();

			evalStep++;
		}
		break;

		case 5:
		{
			if (calculate_dmdt) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRKF45_Step5_withReductions();

				calculate_dmdt = false;

				Set_dmdt();
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRKF45_Step5();
			}

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;

			if (!SetAdaptiveTimeStep()) break;

			//done for this step. Make it available.
			available = true;
		}
		break;
		}
#endif
	}
	break;

	case EVAL_SD:
	{
#ifdef ODE_EVAL_SD
		if (primed) {

			//1. calculate parameters for Barzilai-Borwein stepsizes -> solver must be primed with step 0 (after it is primed next loop starts from step 1)
			//must reset the static delta_... quantities before running these across all meshes
			delta_M_sq = 0.0;
			delta_G_sq = 0.0;
			delta_M_dot_delta_G = 0.0;

			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunSD_BB();

			//2. Set stepsize - alternate between BB values
			if (iteration % 2) {

				if (delta_M_dot_delta_G) {

					dT = delta_M_sq / delta_M_dot_delta_G;
				}
				else dT = SD_DEFAULT_DT;
			}
			else {

				if (delta_G_sq) {

					dT = delta_M_dot_delta_G / delta_G_sq;
				}
				else dT = SD_DEFAULT_DT;
			}
			//don't enforce dT limits, even if dT becomes negative.

			//3. set new magnetization vectors
			if (calculate_mxh || calculate_dmdt) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunSD_Advance_withReductions();

				if (calculate_mxh) Set_mxh();
				if (calculate_dmdt) Set_dmdt();

				calculate_mxh = false;
				calculate_dmdt = false;
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunSD_Advance();
			}

			iteration++;
			stageiteration++;
			time += dT;
			stagetime += dT;
		}
		else {

			dT = SD_DEFAULT_DT;

			//0. prime the SD solver
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunSD_Start();

			evalStep = 0;
			iteration++;
			stageiteration++;
			time += dT;
			stagetime += dT;

			primed = true;
		}
#endif
	}
	break;
	}
}

void ODECommon::MovingMeshAlgorithm(SuperMesh* pSMesh)
{
	//should check moving_mesh flag before calling this, so another check is not done here
	for (int idx = 0; idx < pODE.size(); idx++) {

		double mesh_shift = pODE[idx]->pMesh->CheckMoveMesh();

		if (IsNZ(mesh_shift)) {

			moving_mesh_dwshift -= mesh_shift;

			//do the actual mesh shifts
			for (int idx_mesh = 0; idx_mesh < pSMesh->size(); idx_mesh++) {

				(*pSMesh)[idx_mesh]->MoveMesh(mesh_shift);
			}

			//only one trigger used
			break;
		}
	}
}