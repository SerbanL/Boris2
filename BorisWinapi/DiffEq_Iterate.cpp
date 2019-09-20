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

			//reduce time step based on error ratio
			dT *= sqrt(err_high_fail / (2 * lte));

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

bool ODECommon::SetAdaptiveTimeStep_SingleThreshold2(void)
{
	double lte = Get_lte();

	//failed - try again. Restore state to start of evaluation and try again with smaller time step - output not available yet.
	//do not redo if time step is at or below the minimum value allowed - this is used as a "timeout" to stop the solver from getting stuck on the same iteration; in this case the method is not good enough for the problem.
	if (lte > err_high_fail && dT > dT_min) {

		time -= dT;
		stagetime -= dT;

		for (int idx = 0; idx < (int)pODE.size(); idx++)
			pODE[idx]->RestoreMagnetisation();

		//reduce time step based on error ratio;
		dT *= pow(err_high_fail / (2*lte), 0.25);

		//failed - must repeat
		return false;
	}

	//not failed, adjust time step
	dT *= pow(err_high_fail / (2*lte), 0.25);

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

			available = false;

			//RK23 has the FSAL property, thus the error is calculated on this first stage -> must be primed by a full pass first
			//It also means the above methods do not advance the magnetization yet as we need to calculate the new stepsize first
			//Magnetization is advanced below if the step has not failed
			//If the step has failed then restore magnetization, and also set primed to false -> we must now take a full pass again before we can calculate a new stepsize.
			if (primed) {

				if (!SetAdaptiveTimeStep()) {

					primed = false;
					break;
				}
			}
			else primed = true;

			//Advance magnetization with new stepsize
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRK23_Step0_Advance();

			evalStep++;
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

	case EVAL_RKCK:
	{
#ifdef ODE_EVAL_RKCK
		switch (evalStep) {

		case 0:
		{
			if (calculate_mxh) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRKCK45_Step0_withReductions();

				calculate_mxh = false;

				Set_mxh();
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRKCK45_Step0();
			}

			evalStep++;
			available = false;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKCK45_Step1();

			evalStep++;
		}
		break;

		case 2:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKCK45_Step2();

			evalStep++;
		}
		break;

		case 3:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKCK45_Step3();

			evalStep++;
		}
		break;

		case 4:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKCK45_Step4();

			evalStep++;
		}
		break;

		case 5:
		{
			if (calculate_dmdt) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRKCK45_Step5_withReductions();

				calculate_dmdt = false;

				Set_dmdt();
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRKCK45_Step5();
			}

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;

			//if (!SetAdaptiveTimeStep_SingleThreshold2()) break;
			if (!SetAdaptiveTimeStep()) break;

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
			if (calculate_mxh) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRKDP54_Step0_withReductions();

				calculate_mxh = false;

				Set_mxh();
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRKDP54_Step0();
			}

			available = false;

			//RKDP has the FSAL property, thus the error is calculated on this first stage -> must be primed by a full pass first
			//It also means the above methods do not advance the magnetization yet as we need to calculate the new stepsize first
			//Magnetization is advanced below if the step has not failed
			//If the step has failed then restore magnetization, and also set primed to false -> we must now take a full pass again before we can calculate a new stepsize.
			if (primed) {

				if (!SetAdaptiveTimeStep()) {

					primed = false;
					break;
				}
			}
			else primed = true;

			//Advance magnetization with new stepsize
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKDP54_Step0_Advance();

			evalStep++;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKDP54_Step1();

			evalStep++;
		}
		break;

		case 2:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKDP54_Step2();

			evalStep++;
		}
		break;

		case 3:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKDP54_Step3();

			evalStep++;
		}
		break;

		case 4:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKDP54_Step4();

			evalStep++;
		}
		break;

		case 5:
		{
			if (calculate_dmdt) {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRKDP54_Step5_withReductions();

				calculate_dmdt = false;

				Set_dmdt();
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunRKDP54_Step5();
			}

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

//depending on the set evaluation method, some effective fields (demag field in particular) might not need to be updated at certain steps in the evaluation - keep the previous calculated field contribution
//Modules can call this method to check if they should be updating the field (true) or not (false).
//EVALSPEEDUPSTEP_SKIP : do not update field, use previous calculation if available
//EVALSPEEDUPSTEP_COMPUTE_NO_SAVE : update field and do not save calculation for next step (since at the next step we'll have to calculate field again so no point saving it)
//EVALSPEEDUPSTEP_COMPUTE_AND_SAVE : update field and save calculation for next time (since at the next step we'll need to re-use calculation)
//To enable this mode you need to set use_evaluation_speedup != EVALSPEEDUP_NONE
int ODECommon::Check_Step_Update(void)
{
	//must enable by setting use_evaluation_speedup != EVALSPEEDUP_NONE
	if (use_evaluation_speedup == EVALSPEEDUP_NONE) return EVALSPEEDUPSTEP_COMPUTE_NO_SAVE;

	//it is possible to indicate skipping on step 0, but the caller will have to ensure a previous field evaluation exists which can be used instead - we do not keep track of that here.
	
	//TEuler
	static int evalspeedup_teuler[3][2] = {
		//accurate
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE },
		//aggressive
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE },
		//extreme
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE }
	};

	//AHeun
	static int evalspeedup_aheun[3][2] = {
		//accurate
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE },
		//aggressive
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE },
		//extreme
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE }
	};

	//ABM
	static int evalspeedup_abm[3][2] = {
		//accurate
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE },
		//aggressive
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE },
		//extreme
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE }
	};

	//RK23
	static int evalspeedup_rk23[3][3] = {
		//accurate
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_NO_SAVE , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE },
		//aggressive
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE , EVALSPEEDUPSTEP_SKIP },
		//extreme
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE , EVALSPEEDUPSTEP_SKIP }
	};

	//RK4
	static int evalspeedup_rk4[3][4] = {
		//accurate
		{ EVALSPEEDUPSTEP_COMPUTE_NO_SAVE , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_NO_SAVE },
		//aggressive
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE },
		//extreme
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE }
	};

	//RKF
	static int evalspeedup_rkf[3][6] = {
		//accurate
		{ EVALSPEEDUPSTEP_COMPUTE_NO_SAVE , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE, EVALSPEEDUPSTEP_SKIP, EVALSPEEDUPSTEP_SKIP },
		//aggressive
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE, EVALSPEEDUPSTEP_SKIP, EVALSPEEDUPSTEP_SKIP },
		//extreme
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE, EVALSPEEDUPSTEP_SKIP, EVALSPEEDUPSTEP_SKIP }
	};

	//RKCK
	static int evalspeedup_rkck[3][6] = {
		//accurate
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_NO_SAVE, EVALSPEEDUPSTEP_COMPUTE_AND_SAVE, EVALSPEEDUPSTEP_SKIP },
		//aggressive
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_NO_SAVE, EVALSPEEDUPSTEP_COMPUTE_AND_SAVE, EVALSPEEDUPSTEP_SKIP },
		//extreme
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_SKIP, EVALSPEEDUPSTEP_COMPUTE_AND_SAVE, EVALSPEEDUPSTEP_SKIP }
	};

	//RKDP
	static int evalspeedup_rkdp[3][6] = {
		//accurate
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_NO_SAVE , EVALSPEEDUPSTEP_COMPUTE_NO_SAVE , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE, EVALSPEEDUPSTEP_SKIP, EVALSPEEDUPSTEP_SKIP },
		//aggressive
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE, EVALSPEEDUPSTEP_SKIP, EVALSPEEDUPSTEP_SKIP },
		//extreme
		{ EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_SKIP , EVALSPEEDUPSTEP_COMPUTE_AND_SAVE, EVALSPEEDUPSTEP_SKIP, EVALSPEEDUPSTEP_SKIP }
	};

	switch (evalMethod) {

	case EVAL_EULER:
	{
		return EVALSPEEDUPSTEP_COMPUTE_NO_SAVE;
	}
	break;

	case EVAL_TEULER:
	{
		//0: 0, 1: 1/2

		//Accurate, Aggresive, Extreme : skip on 0 (accurate)
		return evalspeedup_teuler[use_evaluation_speedup - 1][evalStep];
	}
	break;

	case EVAL_AHEUN:
	{
		//0: 0, 1: 1/2

		//Accurate, Aggresive, Extreme : skip on 0 (accurate)
		return evalspeedup_aheun[use_evaluation_speedup - 1][evalStep];
	}
	break;

	case EVAL_RK4:
	{
		//0: 0, 1: 1/2, 2: 1/2, 3: 1

		//Accurate:  skip on 2 (accurate)
		//Aggressive : skip on 0 and 2 (still accurate)
		//Extreme : evaluate last only (inaccurate)

		return evalspeedup_rk4[use_evaluation_speedup - 1][evalStep];
	}
	break;

	case EVAL_ABM:
	{
		//Accurate, Aggresive, Extreme : skip on 0 (accurate)

		if (primed) return evalspeedup_abm[use_evaluation_speedup - 1][evalStep];
	}
	break;

	case EVAL_RK23:
	{
		//0: 0, 1: 1/2, 2 : 3/4

		//Accurate: skip on 0 (accurate)
		//Aggressive, Extreme (some loss of accuracy): skip on 0 and 2

		return evalspeedup_rk23[use_evaluation_speedup - 1][evalStep];
	}
	break;

	case EVAL_SD:
	{
		return EVALSPEEDUPSTEP_COMPUTE_NO_SAVE;
	}
	break;

	case EVAL_RKF:
	{
		//0: 0, 1: 1/4, 2: 3/8, 3: 12/13, 4: 1, 5: 1/2

		//Accurate: skip on 2, 4, 5 (accurate)
		//Aggressive: skip on 0, 2, 4, 5 (still accurate)
		//Extreme : evaluate at 3 only (inaccurate but surprisingly close)

		return evalspeedup_rkf[use_evaluation_speedup - 1][evalStep];
	}
	break;

	case EVAL_RKCK:
	{
		//0: 0, 1: 1/5, 2: 3/10, 3: 3/5, 4: 1, 5: 7/8

		//Accurate: skip on 0, 2, 5 (accurate)
		//Aggressive: skip on 0, 1, 2, 5 (some loss of accuracy)
		//Extreme: evaluate at 4 only (inaccurate)

		return evalspeedup_rkck[use_evaluation_speedup - 1][evalStep];
	}
	break;

	case EVAL_RKDP:
	{
		//0: 0, 1: 1/5, 2: 3/10, 3: 4/5, 4: 8/9, 5: 1

		//Accurate: skip on 0, 4, 5 (accurate)
		//Aggressive: skip on 0, 2, 4, 5 (a little loss of accuracy)
		//Extreme: evaluate at 3 only (inaccurate)

		return evalspeedup_rkdp[use_evaluation_speedup - 1][evalStep];
	}
	break;
	}

	return EVALSPEEDUPSTEP_COMPUTE_AND_SAVE;
}