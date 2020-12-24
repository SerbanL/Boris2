#include "stdafx.h"
#include "DiffEq_CommonBase.h"

#include "DiffEq_Common.h"
#include "Atom_DiffEq_Common.h"

#include "DiffEq.h"
#include "Atom_DiffEq.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//calculate time step for adaptive methods based on current error values and evaluation method settings
//return true if current time step is good, else false if we must repeat it
//this uses a 2 level error threshold -> above the high threshold fail, adjust step based on max_error / error ratio. Below the low error threshold increase step by a small constant factor.
bool ODECommon_Base::SetAdaptiveTimeStep(void)
{
	//adaptive time step based on lte - is lte over acceptable relative error?
	if (lte > err_high) {

		//failed - try again. Restore state to start of evaluation and try again with smaller time step - output not available yet.
		//do not redo if time step is at or below the minimum value allowed - this is used as a "timeout" to stop the solver from getting stuck on the same iteration; in this case the method is not good enough for the problem.
		if (lte > err_high_fail && dT > dT_min) {

			time -= dT;
			stagetime -= dT;

			//reduce time step based on error ratio
			dT *= sqrt(err_high_fail / (2 * lte));

			//failed - must repeat
			available = false;
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ODECommon_Base::Iterate(void)
{
	//save current dT value in case it changes (adaptive time step methods)
	dT_last = dT;

	switch (evalMethod) {

	case EVAL_EULER:
	{
#ifdef ODE_EVAL_COMPILATION_EULER
		if (calculate_mxh || calculate_dmdt) {

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunEuler_withReductions();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunEuler_withReductions();
			}

			if (calculate_mxh) {

				calculate_mxh = false;
				mxh = 0.0;
				podeSolver->Set_mxh();
				patom_odeSolver->Set_mxh();
			}

			if (calculate_dmdt) {

				calculate_dmdt = false;
				dmdt = 0.0;
				podeSolver->Set_dmdt();
				patom_odeSolver->Set_dmdt();
			}
		}
		else {

			for (int idx = 0; idx < (int)podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunEuler();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunEuler();
			}
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
#ifdef ODE_EVAL_COMPILATION_TEULER
		if (evalStep == 0) {

			if (calculate_mxh) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunTEuler_Step0_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunTEuler_Step0_withReductions();
				}

				calculate_mxh = false;
				mxh = 0.0;
				podeSolver->Set_mxh();
				patom_odeSolver->Set_mxh();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunTEuler_Step0();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunTEuler_Step0();
				}
			}

			evalStep = 1;
			available = false;
		}
		else {

			if (calculate_dmdt) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunTEuler_Step1_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunTEuler_Step1_withReductions();
				}

				calculate_dmdt = false;
				dmdt = 0.0;
				podeSolver->Set_dmdt();
				patom_odeSolver->Set_dmdt();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunTEuler_Step1();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunTEuler_Step1();
				}
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
#ifdef ODE_EVAL_COMPILATION_AHEUN
		if (evalStep == 0) {

			if (calculate_mxh) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunAHeun_Step0_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunAHeun_Step0_withReductions();
				}

				calculate_mxh = false;
				mxh = 0.0;
				podeSolver->Set_mxh();
				patom_odeSolver->Set_mxh();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunAHeun_Step0();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunAHeun_Step0();
				}
			}

			evalStep = 1;
			available = false;
		}
		else {

			if (calculate_dmdt) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunAHeun_Step1_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunAHeun_Step1_withReductions();
				}

				calculate_dmdt = false;
				dmdt = 0.0;
				podeSolver->Set_dmdt();
				patom_odeSolver->Set_dmdt();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunAHeun_Step1();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunAHeun_Step1();
				}
			}

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;
			available = true;

			dT_last = dT;
			//accumulate lte value so we can use it to adjust the time step
			lte = 0.0;
			podeSolver->Set_lte();
			patom_odeSolver->Set_lte();

			if (!SetAdaptiveTimeStep()) {

				podeSolver->Restore();
				patom_odeSolver->Restore();
			}
		}
#endif
	}
	break;

	case EVAL_ABM:
	{
#ifdef ODE_EVAL_COMPILATION_ABM
		if (primed) {

			if (evalStep == 0) {

				if (calculate_mxh) {

					for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

						podeSolver->pODE[idx]->RunABM_Predictor_withReductions();
					}

					for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

						patom_odeSolver->pODE[idx]->RunABM_Predictor_withReductions();
					}

					calculate_mxh = false;

					mxh = 0.0;
					podeSolver->Set_mxh();
					patom_odeSolver->Set_mxh();
				}
				else {

					for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

						podeSolver->pODE[idx]->RunABM_Predictor();
					}

					for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

						patom_odeSolver->pODE[idx]->RunABM_Predictor();
					}
				}

				evalStep = 1;
				available = false;
			}
			else {

				if (calculate_dmdt) {

					for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

						podeSolver->pODE[idx]->RunABM_Corrector_withReductions();
					}

					for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

						patom_odeSolver->pODE[idx]->RunABM_Corrector_withReductions();
					}

					calculate_dmdt = false;

					dmdt = 0.0;
					podeSolver->Set_dmdt();
					patom_odeSolver->Set_dmdt();
				}
				else {

					for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

						podeSolver->pODE[idx]->RunABM_Corrector();
					}

					for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

						patom_odeSolver->pODE[idx]->RunABM_Corrector();
					}
				}

				evalStep = 0;
				time += dT;
				stagetime += dT;
				iteration++;
				stageiteration++;
				available = true;

				dT_last = dT;
				lte = 0.0;
				podeSolver->Set_lte();
				patom_odeSolver->Set_lte();

				if (!SetAdaptiveTimeStep()) {

					podeSolver->Restore();
					patom_odeSolver->Restore();
				}
				else {

					//switch between saved equation evaluations (use latest)
					alternator = !alternator;
				}
			}
		}
		else {

			if (evalStep == 0) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunABM_TEuler0();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunABM_TEuler0();
				}

				evalStep = 1;
				available = false;
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunABM_TEuler1();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunABM_TEuler1();
				}

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
#ifdef ODE_EVAL_COMPILATION_RK23
		switch (evalStep) {

		case 0:
		{
			if (calculate_mxh) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRK23_Step0_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRK23_Step0_withReductions();
				}

				calculate_mxh = false;
				mxh = 0.0;
				podeSolver->Set_mxh();
				patom_odeSolver->Set_mxh();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRK23_Step0();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRK23_Step0();
				}
			}

			available = false;

			//RK23 has the FSAL property, thus the error is calculated on this first stage -> must be primed by a full pass first
			//It also means the above methods do not advance the magnetization yet as we need to calculate the new stepsize first
			//Magnetization is advanced below if the step has not failed
			//If the step has failed then restore magnetization, and also set primed to false -> we must now take a full pass again before we can calculate a new stepsize.
			if (primed) {

				dT_last = dT;
				lte = 0.0;
				podeSolver->Set_lte();
				patom_odeSolver->Set_lte();

				if (!SetAdaptiveTimeStep()) {

					primed = false;
					podeSolver->Restore();
					patom_odeSolver->Restore();
					break;
				}
			}
			else primed = true;

			//Advance with new stepsize
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRK23_Step0_Advance();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRK23_Step0_Advance();
			}

			evalStep++;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRK23_Step1();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRK23_Step1();
			}

			evalStep++;
		}
		break;

		case 2:
		{
			if (calculate_dmdt) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRK23_Step2_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRK23_Step2_withReductions();
				}

				calculate_dmdt = false;
				dmdt = 0.0;
				podeSolver->Set_dmdt();
				patom_odeSolver->Set_dmdt();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRK23_Step2();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRK23_Step2();
				}
			}

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;
			available = true;
		}
		break;
		}
#endif
	}
	break;

	case EVAL_RK4:
	{
#ifdef ODE_EVAL_COMPILATION_RK4
		switch (evalStep) {

		case 0:
		{
			if (calculate_mxh) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRK4_Step0_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRK4_Step0_withReductions();
				}

				calculate_mxh = false;
				mxh = 0.0;
				podeSolver->Set_mxh();
				patom_odeSolver->Set_mxh();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRK4_Step0();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRK4_Step0();
				}
			}

			evalStep++;
			available = false;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRK4_Step1();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRK4_Step1();
			}

			evalStep++;
		}
		break;

		case 2:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRK4_Step2();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRK4_Step2();
			}

			evalStep++;
		}
		break;

		case 3:
		{
			if (calculate_dmdt) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRK4_Step3_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRK4_Step3_withReductions();
				}

				calculate_dmdt = false;
				dmdt = 0.0;
				podeSolver->Set_dmdt();
				patom_odeSolver->Set_dmdt();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRK4_Step3();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRK4_Step3();
				}
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
#ifdef ODE_EVAL_COMPILATION_RKF
		switch (evalStep) {

		case 0:
		{
			if (calculate_mxh) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKF45_Step0_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKF45_Step0_withReductions();
				}

				calculate_mxh = false;
				mxh = 0.0;
				podeSolver->Set_mxh();
				patom_odeSolver->Set_mxh();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKF45_Step0();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKF45_Step0();
				}
			}

			evalStep++;
			available = false;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKF45_Step1();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKF45_Step1();
			}

			evalStep++;
		}
		break;

		case 2:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKF45_Step2();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKF45_Step2();
			}

			evalStep++;
		}
		break;

		case 3:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKF45_Step3();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKF45_Step3();
			}

			evalStep++;
		}
		break;

		case 4:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKF45_Step4();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKF45_Step4();
			}

			evalStep++;
		}
		break;

		case 5:
		{
			if (calculate_dmdt) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKF45_Step5_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKF45_Step5_withReductions();
				}

				calculate_dmdt = false;
				dmdt = 0.0;
				podeSolver->Set_dmdt();
				patom_odeSolver->Set_dmdt();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKF45_Step5();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKF45_Step5();
				}
			}

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;
			available = true;

			dT_last = dT;
			lte = 0.0;
			podeSolver->Set_lte();
			patom_odeSolver->Set_lte();

			if (!SetAdaptiveTimeStep()) {

				podeSolver->Restore();
				patom_odeSolver->Restore();
			}
		}
		break;
		}
#endif
	}
	break;

	case EVAL_RKCK:
	{
#ifdef ODE_EVAL_COMPILATION_RKCK
		switch (evalStep) {

		case 0:
		{
			if (calculate_mxh) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKCK45_Step0_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKCK45_Step0_withReductions();
				}

				calculate_mxh = false;
				mxh = 0.0;
				podeSolver->Set_mxh();
				patom_odeSolver->Set_mxh();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKCK45_Step0();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKCK45_Step0();
				}
			}

			evalStep++;
			available = false;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKCK45_Step1();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKCK45_Step1();
			}

			evalStep++;
		}
		break;

		case 2:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKCK45_Step2();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKCK45_Step2();
			}

			evalStep++;
		}
		break;

		case 3:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKCK45_Step3();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKCK45_Step3();
			}

			evalStep++;
		}
		break;

		case 4:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKCK45_Step4();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKCK45_Step4();
			}

			evalStep++;
		}
		break;

		case 5:
		{
			if (calculate_dmdt) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKCK45_Step5_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKCK45_Step5_withReductions();
				}

				calculate_dmdt = false;
				dmdt = 0.0;
				podeSolver->Set_dmdt();
				patom_odeSolver->Set_dmdt();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKCK45_Step5();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKCK45_Step5();
				}
			}

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;
			available = true;

			dT_last = dT;
			lte = 0.0;
			podeSolver->Set_lte();
			patom_odeSolver->Set_lte();

			if (!SetAdaptiveTimeStep()) {

				podeSolver->Restore();
				patom_odeSolver->Restore();
			}
		}
		break;
		}
#endif
	}
	break;

	case EVAL_RKDP:
	{
#ifdef ODE_EVAL_COMPILATION_RKDP
		switch (evalStep) {

		case 0:
		{
			if (calculate_mxh) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKDP54_Step0_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKDP54_Step0_withReductions();
				}

				calculate_mxh = false;
				mxh = 0.0;
				podeSolver->Set_mxh();
				patom_odeSolver->Set_mxh();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKDP54_Step0();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKDP54_Step0();
				}
			}

			available = false;

			//RKDP has the FSAL property, thus the error is calculated on this first stage -> must be primed by a full pass first
			//It also means the above methods do not advance the magnetization yet as we need to calculate the new stepsize first
			//Magnetization is advanced below if the step has not failed
			//If the step has failed then restore magnetization, and also set primed to false -> we must now take a full pass again before we can calculate a new stepsize.
			if (primed) {

				dT_last = dT;
				lte = 0.0;
				podeSolver->Set_lte();
				patom_odeSolver->Set_lte();

				if (!SetAdaptiveTimeStep()) {

					podeSolver->Restore();
					patom_odeSolver->Restore();
					primed = false;
					break;
				}
			}
			else primed = true;

			//Advance magnetization with new stepsize
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKDP54_Step0_Advance();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKDP54_Step0_Advance();
			}

			evalStep++;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKDP54_Step1();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKDP54_Step1();
			}

			evalStep++;
		}
		break;

		case 2:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKDP54_Step2();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKDP54_Step2();
			}

			evalStep++;
		}
		break;

		case 3:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKDP54_Step3();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKDP54_Step3();
			}

			evalStep++;
		}
		break;

		case 4:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKDP54_Step4();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKDP54_Step4();
			}

			evalStep++;
		}
		break;

		case 5:
		{
			if (calculate_dmdt) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKDP54_Step5_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKDP54_Step5_withReductions();
				}

				calculate_dmdt = false;
				dmdt = 0.0;
				podeSolver->Set_dmdt();
				patom_odeSolver->Set_dmdt();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKDP54_Step5();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKDP54_Step5();
				}
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
#ifdef ODE_EVAL_COMPILATION_SD
		if (primed) {

			bool do_sd_reset = false;

			//1. calculate parameters for Barzilai-Borwein stepsizes -> solver must be primed with step 0 (after it is primed next loop starts from step 1)
			//must reset the static delta_... quantities before running these across all meshes
			podeSolver->delta_M_sq = 0.0;
			podeSolver->delta_G_sq = 0.0;
			podeSolver->delta_M_dot_delta_G = 0.0;

			podeSolver->delta_M2_sq = 0.0;
			podeSolver->delta_G2_sq = 0.0;
			podeSolver->delta_M2_dot_delta_G2 = 0.0;

			patom_odeSolver->delta_M_sq = 0.0;
			patom_odeSolver->delta_G_sq = 0.0;
			patom_odeSolver->delta_M_dot_delta_G = 0.0;

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunSD_BB();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunSD_BB();
			}

			//2. Set stepsize - alternate between BB values
			if (iteration % 2) {

				if (podeSolver->pODE.size()) {

					if (podeSolver->delta_M_dot_delta_G) {

						dT = podeSolver->delta_M_sq / podeSolver->delta_M_dot_delta_G;

						//for antiferromagnetic meshes also consider sub-lattice B, take smallest dT value
						if (podeSolver->delta_M2_dot_delta_G2) {

							double dT_2 = podeSolver->delta_M2_sq / podeSolver->delta_M2_dot_delta_G2;
							dT = (dT_2 < dT ? dT_2 : dT);
						}
					}
					else do_sd_reset = true;
				}

				if (patom_odeSolver->pODE.size()) {

					double atom_dT;

					if (patom_odeSolver->delta_M_dot_delta_G && !do_sd_reset) {

						atom_dT = patom_odeSolver->delta_M_sq / patom_odeSolver->delta_M_dot_delta_G;

						if (podeSolver->pODE.size()) dT = (atom_dT < dT ? atom_dT : dT);
						else dT = atom_dT;
					}
					else do_sd_reset = true;
				}
			}
			else {

				if (podeSolver->pODE.size()) {

					if (podeSolver->delta_G_sq) {

						dT = podeSolver->delta_M_dot_delta_G / podeSolver->delta_G_sq;

						//for antiferromagnetic meshes also consider sub-lattice B, take smallest dT value
						if (podeSolver->delta_G2_sq) {

							double dT_2 = podeSolver->delta_M2_dot_delta_G2 / podeSolver->delta_G2_sq;
							dT = (dT_2 < dT ? dT_2 : dT);
						}
					}
					else do_sd_reset = true;
				}

				if (patom_odeSolver->pODE.size()) {

					double atom_dT;

					if (patom_odeSolver->delta_G_sq && !do_sd_reset) {

						atom_dT = patom_odeSolver->delta_M_dot_delta_G / patom_odeSolver->delta_G_sq;

						if (podeSolver->pODE.size()) dT = (atom_dT < dT ? atom_dT : dT);
						else dT = atom_dT;
					}
					else do_sd_reset = true;
				}
			}
			
			if (dT < 0) do_sd_reset = true;

			if (do_sd_reset) {

				dT = (++sd_reset_consecutive_iters) * dT_min;
				if (dT > dT_max) dT = dT_max;
			}
			else sd_reset_consecutive_iters = 0;

			//3. set new magnetization vectors
			if (calculate_mxh || calculate_dmdt) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunSD_Advance_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunSD_Advance_withReductions();
				}

				if (calculate_mxh) {

					calculate_mxh = false;
					mxh = 0.0;
					podeSolver->Set_mxh();
					patom_odeSolver->Set_mxh();
				}

				if (calculate_dmdt) {

					calculate_dmdt = false;
					dmdt = 0.0;
					podeSolver->Set_dmdt();
					patom_odeSolver->Set_dmdt();
				}
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunSD_Advance();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunSD_Advance();
				}
			}

			iteration++;
			stageiteration++;
			time += dT;
			stagetime += dT;
		}
		else {

			dT = dT_min;

			//0. prime the SD solver
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunSD_Start();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunSD_Start();
			}

			evalStep = 0;
			iteration++;
			stageiteration++;
			time += dT;
			stagetime += dT;

			primed = true;
			sd_reset_consecutive_iters = 0;
		}
#endif
	}
	break;
	}
}

//depending on the set evaluation method, some effective fields (demag field in particular) might not need to be updated at certain steps in the evaluation - keep the previous calculated field contribution
//Modules can call this method to check if they should be updating the field (true) or not (false).
//EVALSPEEDUPSTEP_SKIP : do not update field, use previous calculation if available
//EVALSPEEDUPSTEP_COMPUTE_NO_SAVE : update field and do not save calculation for next step (since at the next step we'll have to calculate field again so no point saving it)
//EVALSPEEDUPSTEP_COMPUTE_AND_SAVE : update field and save calculation for next time (since at the next step we'll need to re-use calculation)
//To enable this mode you need to set use_evaluation_speedup != EVALSPEEDUP_NONE
int ODECommon_Base::Check_Step_Update(void)
{
	//must enable by setting use_evaluation_speedup != EVALSPEEDUP_NONE
	if (use_evaluation_speedup == EVALSPEEDUP_NONE) return EVALSPEEDUPSTEP_COMPUTE_NO_SAVE;

	if (use_evaluation_speedup == EVALSPEEDUP_EXTREME && !link_dTspeedup) {

		//extreme mode with time step not linked to dT (the intention is for the updating time step to be larger than dT, but in any case in this mode the updating will be done at most once per dT time step)
		//Note, in this mode the update will be done at the first evaluation step, whereas in extreme mode with link_dTspeedup true, it will be done as indicated below for each method (not necessarily at the first evaluation)
		
		//recommend skipping (and re-using previous value) if time step from previous evaluation is too small
		if (time < time_speedup + dTspeedup) return EVALSPEEDUPSTEP_SKIP;
		else {

			//update evaluation time and recommend updating
			time_speedup = time;
			return EVALSPEEDUPSTEP_COMPUTE_AND_SAVE;
		}
	}

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