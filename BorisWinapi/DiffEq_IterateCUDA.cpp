#include "stdafx.h"
#include "DiffEq.h"
#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"

#if COMPILECUDA == 1

//---------------------------------------- ITERATE METHODS : CUDA


bool ODECommon::SetAdaptiveTimeStepCUDA(void)
{
	double lte = pODECUDA->Get_lte();

	//adaptive time step based on lte - is lte over acceptable relative error?
	if (lte > err_high) {

		//failed - try again. Restore state to start of evaluation and try again with smaller time step - output not available yet.
		//do not redo if time step is at or below the minimum value allowed - this is used as a "timeout" to stop the solver from getting stuck on the same iteration; in this case the method is not good enough for the problem.
		if (lte > err_high_fail && dT > dT_min) {

			time -= dT;
			stagetime -= dT;

			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->pmeshODECUDA->RestoreMagnetisation();

			//reduce time step based on error ratio
			dT *= sqrt(err_high_fail / (2 * lte));

			pODECUDA->Sync_dT();

			//failed - must repeat
			return false;
		}

		//not failed but still need to reduce time step
		dT *= pow(err_high / (2 * lte), 0.25);

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

	pODECUDA->Sync_dT();

	//good, next step
	return true;
}

bool ODECommon::SetAdaptiveTimeStepCUDA_SingleThreshold(void)
{
	double lte = pODECUDA->Get_lte();

	//failed - try again. Restore state to start of evaluation and try again with smaller time step - output not available yet.
	//do not redo if time step is at or below the minimum value allowed - this is used as a "timeout" to stop the solver from getting stuck on the same iteration; in this case the method is not good enough for the problem.
	if (lte > err_high_fail && dT > dT_min) {

		time -= dT;
		stagetime -= dT;

		for (int idx = 0; idx < (int)pODE.size(); idx++)
			pODE[idx]->pmeshODECUDA->RestoreMagnetisation();

		//reduce time step based on error ratio; use a 0.99 back-off factor so the solver doesn't get stuck
		dT *= 0.99 * pow(err_high_fail / lte, 1.0/3.0);

		pODECUDA->Sync_dT();

		//failed - must repeat
		return false;
	}

	//not failed but still need to reduce time step
	dT *= sqrt(err_high_fail / lte);

	//must not go below minimum time step
	if (dT < dT_min) dT = dT_min;

	//must not go above maximum time step
	if (dT > dT_max) dT = dT_max;

	pODECUDA->Sync_dT();

	//good, next step
	return true;
}

bool ODECommon::SetAdaptiveTimeStepCUDA_SingleThreshold2(void)
{
	double lte = pODECUDA->Get_lte();

	//failed - try again. Restore state to start of evaluation and try again with smaller time step - output not available yet.
	//do not redo if time step is at or below the minimum value allowed - this is used as a "timeout" to stop the solver from getting stuck on the same iteration; in this case the method is not good enough for the problem.
	if (lte > err_high_fail && dT > dT_min) {

		time -= dT;
		stagetime -= dT;

		for (int idx = 0; idx < (int)pODE.size(); idx++)
			pODE[idx]->pmeshODECUDA->RestoreMagnetisation();

		//reduce time step based on error ratio
		dT *= pow(err_high_fail / (2 * lte), 0.25);

		pODECUDA->Sync_dT();

		//failed - must repeat
		return false;
	}

	//not failed but still need to reduce time step
	dT *= pow(err_high_fail / (2 * lte), 0.25);

	//must not go below minimum time step
	if (dT < dT_min) dT = dT_min;

	//must not go above maximum time step
	if (dT > dT_max) dT = dT_max;

	pODECUDA->Sync_dT();

	//good, next step
	return true;
}

void ODECommon::IterateCUDA(void)
{
	//save current dT value in case it changes (adaptive time step methods)
	dT_last = dT;
	pODECUDA->Sync_dT_last();

	switch (evalMethod) {

	case EVAL_EULER:
	{
#ifdef ODE_EVAL_EULER
		if (calculate_mxh && calculate_dmdt) pODECUDA->Zero_reduction_values();
		else if (calculate_mxh) pODECUDA->Zero_mxh_lte_values();
		else if (calculate_dmdt) pODECUDA->Zero_dmdt_lte_values();

		for (int idx = 0; idx < (int)pODE.size(); idx++) {

			//Euler can be used for stochastic equations
			if (pODE[idx]->H_Thermal.linear_size() && pODE[idx]->Torque_Thermal.linear_size()) pODE[idx]->pmeshODECUDA->GenerateThermalField_and_Torque();
			else if (pODE[idx]->H_Thermal.linear_size()) pODE[idx]->pmeshODECUDA->GenerateThermalField();

			if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunEuler_LLG(calculate_mxh, calculate_dmdt);
			else pODE[idx]->pmeshODECUDA->RunEuler(calculate_mxh, calculate_dmdt);
		}

		if (calculate_mxh) pODECUDA->mxhav_to_mxh();
		if (calculate_dmdt) pODECUDA->dmdtav_to_dmdt();

		calculate_mxh = false;
		calculate_dmdt = false;

		time += dT;
		stagetime += dT;
		iteration++;
		stageiteration++;
		available = true;
#endif
	}
	break;

	case EVAL_TEULER:
	{
#ifdef ODE_EVAL_TEULER
		if (evalStep == 0) {

			if (calculate_mxh) pODECUDA->Zero_mxh_lte_values();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				//Trapezoidal Euler can be used for stochastic equations
				if (pODE[idx]->H_Thermal.linear_size() && pODE[idx]->Torque_Thermal.linear_size()) pODE[idx]->pmeshODECUDA->GenerateThermalField_and_Torque();
				else if (pODE[idx]->H_Thermal.linear_size()) pODE[idx]->pmeshODECUDA->GenerateThermalField();

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunTEuler_LLG(0, calculate_mxh, calculate_dmdt);
				else pODE[idx]->pmeshODECUDA->RunTEuler(0, calculate_mxh, calculate_dmdt);
			}

			if (calculate_mxh) pODECUDA->mxhav_to_mxh();

			calculate_mxh = false;

			evalStep = 1;
			available = false;
		}
		else {

			if (calculate_dmdt) pODECUDA->Zero_dmdt_lte_values();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunTEuler_LLG(1, calculate_mxh, calculate_dmdt);
				else pODE[idx]->pmeshODECUDA->RunTEuler(1, calculate_mxh, calculate_dmdt);
			}

			if (calculate_dmdt) pODECUDA->dmdtav_to_dmdt();

			calculate_dmdt = false;

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

			if (calculate_mxh) pODECUDA->Zero_mxh_lte_values();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				//Adaptive Heun can be used for stochastic equations
				if (pODE[idx]->H_Thermal.linear_size() && pODE[idx]->Torque_Thermal.linear_size()) pODE[idx]->pmeshODECUDA->GenerateThermalField_and_Torque();
				else if (pODE[idx]->H_Thermal.linear_size()) pODE[idx]->pmeshODECUDA->GenerateThermalField();

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunAHeun_LLG(0, calculate_mxh, calculate_dmdt);
				else pODE[idx]->pmeshODECUDA->RunAHeun(0, calculate_mxh, calculate_dmdt);
			}

			if (calculate_mxh) pODECUDA->mxhav_to_mxh();

			calculate_mxh = false;

			evalStep = 1;
			available = false;
		}
		else {

			if (calculate_dmdt) pODECUDA->Zero_dmdt_lte_values();
			//need this for lte reduction
			else pODECUDA->Zero_lte_value();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunAHeun_LLG(1, calculate_mxh, calculate_dmdt);
				else pODE[idx]->pmeshODECUDA->RunAHeun(1, calculate_mxh, calculate_dmdt);
			}

			if (calculate_dmdt) pODECUDA->dmdtav_to_dmdt();

			calculate_dmdt = false;

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;

			//if you use the SingleThreshold Method with a stochastic equation you can end up backtracking very often -> slow
			//if (!SetAdaptiveTimeStepCUDA_SingleThreshold()) break;
			if (!SetAdaptiveTimeStepCUDA()) break;

			available = true;
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
			if (calculate_mxh) pODECUDA->Zero_mxh_lte_values();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunRK4_LLG(0, calculate_mxh, calculate_dmdt);
				else pODE[idx]->pmeshODECUDA->RunRK4(0, calculate_mxh, calculate_dmdt);
			}

			calculate_mxh = false;

			evalStep++;
			available = false;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunRK4_LLG(1);
				else pODE[idx]->pmeshODECUDA->RunRK4(1);
			}

			evalStep++;
		}
		break;

		case 2:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunRK4_LLG(2);
				else pODE[idx]->pmeshODECUDA->RunRK4(2);
			}

			evalStep++;
		}
		break;

		case 3:
		{
			if (calculate_dmdt) pODECUDA->Zero_dmdt_lte_values();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunRK4_LLG(3, calculate_mxh, calculate_dmdt);
				else pODE[idx]->pmeshODECUDA->RunRK4(3, calculate_mxh, calculate_dmdt);
			}

			calculate_dmdt = false;

			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;

			evalStep = 0;
			available = true;
		}
		break;
		}
#endif
	}
	break;

	case EVAL_ABM:
	{
#ifdef ODE_EVAL_ABM
		if (primed) {

			if (evalStep == 0) {

				if (calculate_mxh) pODECUDA->Zero_mxh_lte_values();

				for (int idx = 0; idx < (int)pODE.size(); idx++) {

					if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunABM_LLG(0, calculate_mxh, calculate_dmdt);
					else pODE[idx]->pmeshODECUDA->RunABM(0, calculate_mxh, calculate_dmdt);
				}

				calculate_mxh = false;

				evalStep = 1;
				available = false;
			}
			else {

				if (calculate_dmdt) pODECUDA->Zero_dmdt_lte_values();
				//need this for lte reduction
				else pODECUDA->Zero_lte_value();

				for (int idx = 0; idx < (int)pODE.size(); idx++) {

					if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunABM_LLG(1, calculate_mxh, calculate_dmdt);
					else pODE[idx]->pmeshODECUDA->RunABM(1, calculate_mxh, calculate_dmdt);
				}

				calculate_dmdt = false;

				evalStep = 0;
				time += dT;
				stagetime += dT;
				iteration++;
				stageiteration++;

				if (!SetAdaptiveTimeStepCUDA()) break;
				
				//done for this step. Make it available.
				available = true;

				//switch between saved equation evaluations (use latest)
				if (alternator) alternator = false;
				else alternator = true;
				
				pODECUDA->Sync_alternator();
			}
		}
		else {

			if (evalStep == 0) {

				for (int idx = 0; idx < (int)pODE.size(); idx++) {

					if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunABMTEuler_LLG(0);
					pODE[idx]->pmeshODECUDA->RunABMTEuler(0);
				}

				evalStep = 1;
				available = false;
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++) {

					if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunABMTEuler_LLG(1);
					pODE[idx]->pmeshODECUDA->RunABMTEuler(1);
				}

				evalStep = 0;
				available = true;

				time += dT;
				stagetime += dT;
				iteration++;
				stageiteration++;

				primed = true;
				alternator = false;		//evaluated equation is now in sEval0

				pODECUDA->Sync_alternator();
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
			if (calculate_mxh) pODECUDA->Zero_mxh_lte_values();
			//need this for lte reduction
			else pODECUDA->Zero_lte_value();

			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->pmeshODECUDA->RunRK23_Step0_NoAdvance(calculate_mxh);

			calculate_mxh = false;

			available = false;

			//RK23 has the FSAL property, thus the error is calculated on this first stage -> must be primed by a full pass first
			//It also means the above methods do not advance the magnetization yet as we need to calculate the new stepsize first
			//Magnetization is advanced below if the step has not failed
			//If the step has failed then restore magnetization, and also set primed to false -> we must now take a full pass again before we can calculate a new stepsize.
			if (primed) {

				if (!SetAdaptiveTimeStepCUDA()) {

					primed = false;
					break;
				}
			}
			else primed = true;

			//Advance magnetization with new stepsize
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->pmeshODECUDA->RunRK23(0);

			evalStep++;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->pmeshODECUDA->RunRK23(1);

			evalStep++;
		}
		break;

		case 2:
		{
			if (calculate_dmdt) pODECUDA->Zero_dmdt_lte_values();

			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->pmeshODECUDA->RunRK23(2, calculate_mxh, calculate_dmdt);

			calculate_dmdt = false;

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;

			//done for this step. Make it available.
			available = true;
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
			if (calculate_mxh) pODECUDA->Zero_mxh_lte_values();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunRKF45(0, calculate_mxh, calculate_dmdt);
				else pODE[idx]->pmeshODECUDA->RunRKF45(0, calculate_mxh, calculate_dmdt);
			}

			calculate_mxh = false;

			evalStep++;
			available = false;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunRKF45_LLG(1);
				else pODE[idx]->pmeshODECUDA->RunRKF45(1);
			}

			evalStep++;
		}
		break;

		case 2:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunRKF45_LLG(2);
				else pODE[idx]->pmeshODECUDA->RunRKF45(2);
			}

			evalStep++;
		}
		break;

		case 3:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunRKF45_LLG(3);
				else pODE[idx]->pmeshODECUDA->RunRKF45(3);
			}

			evalStep++;
		}
		break;

		case 4:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunRKF45_LLG(4);
				else pODE[idx]->pmeshODECUDA->RunRKF45(4);
			}

			evalStep++;
		}
		break;

		case 5:
		{
			if (calculate_dmdt) pODECUDA->Zero_dmdt_lte_values();
			//need this for lte reduction
			else pODECUDA->Zero_lte_value();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunRKF45_LLG(5, calculate_mxh, calculate_dmdt);
				else pODE[idx]->pmeshODECUDA->RunRKF45(5, calculate_mxh, calculate_dmdt);
			}

			calculate_dmdt = false;

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;

			if (!SetAdaptiveTimeStepCUDA()) break;

			//done for this step. Make it available.
			available = true;
		}
		break;
		}
#endif
	}
	break;

	case EVAL_RKCK:
	{
#ifdef ODE_EVAL_RKCK
		switch (evalStep) {

		case 0:
		{
			if (calculate_mxh) pODECUDA->Zero_mxh_lte_values();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				pODE[idx]->pmeshODECUDA->RunRKCK45(0, calculate_mxh, calculate_dmdt);
			}

			calculate_mxh = false;

			evalStep++;
			available = false;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				pODE[idx]->pmeshODECUDA->RunRKCK45(1);
			}

			evalStep++;
		}
		break;

		case 2:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				pODE[idx]->pmeshODECUDA->RunRKCK45(2);
			}

			evalStep++;
		}
		break;

		case 3:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				pODE[idx]->pmeshODECUDA->RunRKCK45(3);
			}

			evalStep++;
		}
		break;

		case 4:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				pODE[idx]->pmeshODECUDA->RunRKCK45(4);
			}

			evalStep++;
		}
		break;

		case 5:
		{
			if (calculate_dmdt) pODECUDA->Zero_dmdt_lte_values();
			//need this for lte reduction
			else pODECUDA->Zero_lte_value();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				pODE[idx]->pmeshODECUDA->RunRKCK45(5, calculate_mxh, calculate_dmdt);
			}

			calculate_dmdt = false;

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;

			if (!SetAdaptiveTimeStepCUDA()) break;

			//done for this step. Make it available.
			available = true;
		}
		break;
		}
#endif
	}
	break;

	case EVAL_RKDP:
	{
#ifdef ODE_EVAL_RKDP
		switch (evalStep) {

		case 0:
		{
			if (calculate_mxh) pODECUDA->Zero_mxh_lte_values();
			//need this for lte reduction
			else pODECUDA->Zero_lte_value();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				pODE[idx]->pmeshODECUDA->RunRKDP54_Step0_NoAdvance(calculate_mxh);
			}

			calculate_mxh = false;

			available = false;

			//RKDP has the FSAL property, thus the error is calculated on this first stage -> must be primed by a full pass first
			//It also means the above methods do not advance the magnetization yet as we need to calculate the new stepsize first
			//Magnetization is advanced below if the step has not failed
			//If the step has failed then restore magnetization, and also set primed to false -> we must now take a full pass again before we can calculate a new stepsize.
			if (primed) {

				if (!SetAdaptiveTimeStepCUDA()) {

					primed = false;
					break;
				}
			}
			else primed = true;

			//Advance magnetization with new stepsize
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->pmeshODECUDA->RunRKDP54(0);

			evalStep++;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				pODE[idx]->pmeshODECUDA->RunRKDP54(1);
			}

			evalStep++;
		}
		break;

		case 2:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				pODE[idx]->pmeshODECUDA->RunRKDP54(2);
			}

			evalStep++;
		}
		break;

		case 3:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				pODE[idx]->pmeshODECUDA->RunRKDP54(3);
			}

			evalStep++;
		}
		break;

		case 4:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				pODE[idx]->pmeshODECUDA->RunRKDP54(4);
			}

			evalStep++;
		}
		break;

		case 5:
		{
			if (calculate_dmdt) pODECUDA->Zero_dmdt_lte_values();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				pODE[idx]->pmeshODECUDA->RunRKDP54(5, calculate_mxh, calculate_dmdt);
			}

			calculate_dmdt = false;

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;

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
			
			//zero all quantities used for Barzilai-Borwein stepsize calculations
			pODECUDA->Zero_SD_Solver_BB_Values();

			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->pmeshODECUDA->RunSD_BB();

			//2. Set stepsize - alternate between BB values
			//first transfer BB values to cpu memory
			pODECUDA->Get_SD_Solver_BB_Values(&delta_M_sq, &delta_G_sq, &delta_M_dot_delta_G, &delta_M2_sq, &delta_G2_sq, &delta_M2_dot_delta_G2);
			
			if (iteration % 2) {

				if (delta_M_dot_delta_G) {

					dT = delta_M_sq / delta_M_dot_delta_G;

					//for antiferromagnetic meshes also consider sub-lattice B, take smallest dT value
					if (delta_M2_dot_delta_G2) {

						double dT_2 = delta_M2_sq / delta_M2_dot_delta_G2;
						dT = (dT_2 < dT ? dT_2 : dT);
					}
				}
				else dT = SD_DEFAULT_DT;
			}
			else {

				if (delta_G_sq) {

					dT = delta_M_dot_delta_G / delta_G_sq;

					//for antiferromagnetic meshes also consider sub-lattice B, take smallest dT value
					if (delta_G2_sq) {

						double dT_2 = delta_M2_dot_delta_G2 / delta_G2_sq;
						dT = (dT_2 < dT ? dT_2 : dT);
					}
				}
				else dT = SD_DEFAULT_DT;
			}
			//don't enforce dT limits, even if dT becomes negative.

			//make sure to transfer dT value to GPU
			pODECUDA->Sync_dT();

			//3. set new magnetization vectors
			if (calculate_mxh && calculate_dmdt) pODECUDA->Zero_reduction_values();
			else if (calculate_mxh) pODECUDA->Zero_mxh_lte_values();
			else if (calculate_dmdt) pODECUDA->Zero_dmdt_lte_values();

			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->pmeshODECUDA->RunSD_Advance(calculate_mxh, calculate_dmdt);

			calculate_mxh = false;
			calculate_dmdt = false;

			iteration++;
			stageiteration++;
			time += dT;
			stagetime += dT;
		}
		else {

			dT = SD_DEFAULT_DT;

			//make sure to transfer dT value to GPU
			pODECUDA->Sync_dT();

			//0. prime the SD solver
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->pmeshODECUDA->RunSD_Start();

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

#endif