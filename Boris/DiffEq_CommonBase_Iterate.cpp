#include "stdafx.h"
#include "DiffEq_CommonBase.h"

#include "DiffEq_Common.h"
#include "Atom_DiffEq_Common.h"

#include "DiffEq.h"
#include "Atom_DiffEq.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//calculate time step for adaptive methods based on current error values and evaluation method settings
//return true if current time step is good, else false if we must repeat it

bool ODECommon_Base::SetAdaptiveTimeStep(void)
{
	//fixed time step if limits same
	if (dT_min == dT_max) return true;

	//Integral controller adaptive time step.

	if (lte > err_high_fail && dT > dT_min) {
		
		//reject. The dT > dT_min check needed to stop solver getting stuck.
		time -= dT;
		stagetime -= dT;
		
		//Use I controller only when rejecting
		double c = pow(err_high_fail * 0.8 / lte, 1.0 / (eval_method_order + 1));

		if (c < 0.01) c = 0.01;
		dT *= c;

		//must not go below minimum time step
		if (dT < dT_min) dT = dT_min;

		//failed - must repeat
		available = false;
		return false;
	}
	else if (lte > 0) {

		double c = pow(err_high_fail * 0.8 / lte, 1.0 / (eval_method_order + 1));
		if (c > dT_increase) c = dT_increase;
		if (c < 0.01) c = 0.01;
		dT *= c;

		//must not go below minimum time step
		if (dT < dT_min) dT = dT_min;

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

	case EVAL_RKF45:
	{
#ifdef ODE_EVAL_COMPILATION_RKF45
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

	case EVAL_RKF56:
	{
#ifdef ODE_EVAL_COMPILATION_RKF56
		switch (evalStep) {

		case 0:
		{
			if (calculate_mxh) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKF56_Step0_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKF56_Step0_withReductions();
				}

				calculate_mxh = false;
				mxh = 0.0;
				podeSolver->Set_mxh();
				patom_odeSolver->Set_mxh();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKF56_Step0();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKF56_Step0();
				}
			}

			evalStep++;
			available = false;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKF56_Step1();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKF56_Step1();
			}

			evalStep++;
		}
		break;

		case 2:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKF56_Step2();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKF56_Step2();
			}

			evalStep++;
		}
		break;

		case 3:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKF56_Step3();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKF56_Step3();
			}

			evalStep++;
		}
		break;

		case 4:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKF56_Step4();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKF56_Step4();
			}

			evalStep++;
		}
		break;

		case 5:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKF56_Step5();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKF56_Step5();
			}

			evalStep++;
		}
		break;

		case 6:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->RunRKF56_Step6();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->RunRKF56_Step6();
			}

			evalStep++;
		}
		break;

		case 7:
		{
			if (calculate_dmdt) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKF56_Step7_withReductions();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKF56_Step7_withReductions();
				}

				calculate_dmdt = false;
				dmdt = 0.0;
				podeSolver->Set_dmdt();
				patom_odeSolver->Set_dmdt();
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					podeSolver->pODE[idx]->RunRKF56_Step7();
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					patom_odeSolver->pODE[idx]->RunRKF56_Step7();
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

	case EVAL_RKCK45:
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

	case EVAL_RKDP54:
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

//get total time with evaluation step resolution level
double ODECommon_Base::Get_EvalStep_Time(void)
{
	//Define c values in Butcher tableaux

	//Trapezoidal Euler
	static double evaltime_teuler[4] = { 0.0, 1.0 };

	//RK4
	static double evaltime_rk4[4] = { 0.0, 0.5, 0.5, 1.0 };

	//ABM
	static double evaltime_abm[2] = { 0.0, 1.0 };

	//RK23
	static double evaltime_rk23[4] = { 0.0, 0.5, 0.75 };

	//RKF45
	//static double evaltime_rkf45[6] = { 0.0, 0.25, 0.375, 12.0/13, 1.0, 0.5 };
	static double evaltime_rkf45[6] = { 0.0, 2.0/9, 1.0/3, 3.0/4, 1.0, 5.0/6 };

	//RKF56
	static double evaltime_rkf56[8] = { 0.0, 1.0/6, 4.0/15, 2.0/3, 4.0/5, 1.0, 0.0, 1.0 };

	//RKCK45
	static double evaltime_rkck45[6] = { 0.0, 0.2, 0.3, 0.6, 1.0, 0.875 };

	//RKDP54
	static double evaltime_rkdp54[6] = { 0.0, 0.2, 0.3, 0.8, 8.0/9, 1.0 };

	switch (evalMethod) {

	case EVAL_EULER:
	{
		return time;
	}
	break;

	case EVAL_TEULER:
	{
		return time + dT * evaltime_teuler[evalStep];
	}
	break;

	case EVAL_AHEUN:
	{
		return time + dT * evaltime_teuler[evalStep];
	}
	break;

	case EVAL_RK4:
	{
		return time + dT * evaltime_rk4[evalStep];
	}
	break;

	case EVAL_ABM:
	{
		return time + dT * evaltime_abm[evalStep];
	}
	break;

	case EVAL_RK23:
	{
		return time + dT * evaltime_rk23[evalStep];
	}
	break;

	case EVAL_SD:
	{
		return time;
	}
	break;

	case EVAL_RKF45:
	{
		return time + dT * evaltime_rkf45[evalStep];
	}
	break;

	case EVAL_RKF56:
	{
		return time + dT * evaltime_rkf56[evalStep];
	}
	break;

	case EVAL_RKCK45:
	{
		return time + dT * evaltime_rkck45[evalStep];
	}
	break;

	case EVAL_RKDP54:
	{
		return time + dT * evaltime_rkdp54[evalStep];
	}
	break;
	}

	return time;
}

//check if we should update the demag field if using evaluation speedup mode; 
//when in this mode any demag field updates will happen only when ODE step has finished (available == true), but if link_dTspeedup == false this is not necessarily every time step
bool ODECommon_Base::Check_Step_Update(void)
{
	if (use_evaluation_speedup == EVALSPEEDUP_NONE) return true;

	else if (!link_dTspeedup) {
		
		double time_eval = Get_EvalStep_Time();

		//recommend skipping (and re-using previous value) if time step from previous evaluation is too small
		if (time_eval < time_speedup + dTspeedup) return false;
		else {

			//update evaluation time and recommend updating
			time_speedup = time_eval;
			return true;
		}
	}
	else return available == true;
}