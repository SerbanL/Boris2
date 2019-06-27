#include "stdafx.h"
#include "DiffEq.h"
#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"

#if COMPILECUDA == 1

//---------------------------------------- ITERATE METHODS : CUDA

void ODECommon::IterateCUDA(void)
{
	//save current dT value in case it changes (adaptive time step methods)
	dT_last = dT;
	pODECUDA->Sync_dT_last();

	switch (evalMethod) {

	case EVAL_EULER:
	{
		pODECUDA->Zero_reduction_values();
		
		for (int idx = 0; idx < (int)pODE.size(); idx++) {

			//Euler can be used for stochastic equations
			if (pODE[idx]->H_Thermal.linear_size() && pODE[idx]->Torque_Thermal.linear_size()) pODE[idx]->pmeshODECUDA->GenerateThermalField_and_Torque();
			else if (pODE[idx]->H_Thermal.linear_size()) pODE[idx]->pmeshODECUDA->GenerateThermalField();

			if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunEuler_LLG();
			else pODE[idx]->pmeshODECUDA->RunEuler();
		}

		pODECUDA->mxhav_to_mxh();

		time += dT;
		stagetime += dT;
		iteration++;
		stageiteration++;
	}
	break;

	case EVAL_TEULER:
	{
		if (evalStep == 0) {

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				//Trapezoidal Euler can be used for stochastic equations
				if (pODE[idx]->H_Thermal.linear_size() && pODE[idx]->Torque_Thermal.linear_size()) pODE[idx]->pmeshODECUDA->GenerateThermalField_and_Torque();
				else if (pODE[idx]->H_Thermal.linear_size()) pODE[idx]->pmeshODECUDA->GenerateThermalField();

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunTEuler_LLG(0);
				else pODE[idx]->pmeshODECUDA->RunTEuler(0);
			}

			evalStep = 1;
			available = false;
		}
		else {

			pODECUDA->Zero_reduction_values();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunTEuler_LLG(1);
				else pODE[idx]->pmeshODECUDA->RunTEuler(1);
			}

			pODECUDA->mxhav_to_mxh();

			evalStep = 0;
			available = true;

			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;
		}
	}
	break;

	case EVAL_RK4:
	{
		switch (evalStep) {

		case 0:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunRK4_LLG(0);
				else pODE[idx]->pmeshODECUDA->RunRK4(0);
			}

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
			pODECUDA->Zero_reduction_values();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunRK4_LLG(3);
				else pODE[idx]->pmeshODECUDA->RunRK4(3);
			}

			evalStep = 0;
			available = true;

			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;

			Set_mxh();
		}
		break;
		}
	}
	break;

	case EVAL_ABM:
	{
		if (primed) {

			if (evalStep == 0) {

				for (int idx = 0; idx < (int)pODE.size(); idx++) {

					if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunABM_LLG(0);
					else pODE[idx]->pmeshODECUDA->RunABM(0);
				}

				evalStep = 1;
				available = false;
			}
			else {

				pODECUDA->Zero_reduction_values();

				for (int idx = 0; idx < (int)pODE.size(); idx++) {

					if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunABM_LLG(1);
					else pODE[idx]->pmeshODECUDA->RunABM(1);
				}

				evalStep = 0;
				time += dT;
				stagetime += dT;
				
				double lte = pODECUDA->Get_lte();

				//adaptive time step based on lte - is lte over acceptable relative error?
				if (lte > ABM_RELERRMAX) {

					//failed - try again. Restore state to start of evaluation and try again with smaller time step - output not available yet.
					//do not redo if time step is at or below the minimum value allowed - this is used as a "timeout" to stop the solver from getting stuck on the same iteration; in this case the method is not good enough for the problem.
					if (lte > ABM_RELERRFAIL && dT > ABM_MINDT) {

						time -= dT;
						stagetime -= dT;

						for (int idx = 0; idx < (int)pODE.size(); idx++)
							pODE[idx]->pmeshODECUDA->RestoreMagnetisation();

						dT *= ABM_DTREDUCE;

						pODECUDA->Sync_dT();
						break;
					}

					//not failed but still need to reduce time step
					dT *= ABM_DTREDUCE;

					//must not go below minimum time step
					if (dT < ABM_MINDT) dT = ABM_MINDT;
				}

				//is lte below minimum relative error ? If yes we can go quicker
				if (lte < ABM_RELERRMIN) {

					dT *= ABM_DTINCREASE;

					//must not go above maximum time step
					if (dT > ABM_MAXDT) dT = ABM_MAXDT;
				}
				
				//done for this step. Make it available.
				available = true;

				iteration++;
				stageiteration++;

				//switch between saved equation evaluations (use latest)
				if (alternator) alternator = false;
				else alternator = true;
				
				pODECUDA->Sync_dT();
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
	}
	break;

	case EVAL_RKF:
	{
		switch (evalStep) {

		case 0:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunRKF45_LLG(0);
				else pODE[idx]->pmeshODECUDA->RunRKF45(0);
			}

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
			pODECUDA->Zero_reduction_values();

			for (int idx = 0; idx < (int)pODE.size(); idx++) {

				if (setODE == ODE_LLG) pODE[idx]->pmeshODECUDA->RunRKF45_LLG(5);
				else pODE[idx]->pmeshODECUDA->RunRKF45(5);
			}

			evalStep = 0;
			time += dT;
			stagetime += dT;

			double lte = pODECUDA->Get_lte();

			//adaptive time step based on lte - is lte over acceptable relative error?
			if (lte > RKF_RELERRMAX) {

				//failed - try again. Restore state to start of evaluation and try again with smaller time step - output not available yet.
				//do not redo if time step is at or below the minimum value allowed - this is used as a "timeout" to stop the solver from getting stuck on the same iteration; in this case the method is not good enough for the problem.
				if (lte > RKF_RELERRFAIL && dT > RKF_MINDT) {

					time -= dT;
					stagetime -= dT;

					for (int idx = 0; idx < (int)pODE.size(); idx++)
						pODE[idx]->pmeshODECUDA->RestoreMagnetisation();

					dT *= RKF_DTREDUCE;
					
					pODECUDA->Sync_dT();
					break;
				}

				//not failed but still need to reduce time step
				dT *= RKF_DTREDUCE;

				//must not go below minimum time step
				if (dT < RKF_MINDT) dT = RKF_MINDT;

				pODECUDA->Sync_dT();
			}

			//is lte below minimum relative error ? If yes we can go quicker
			if (lte < RKF_RELERRMIN) {

				dT *= RKF_DTINCREASE;

				//must not go above maximum time step
				if (dT > RKF_MAXDT) dT = RKF_MAXDT;

				pODECUDA->Sync_dT();
			}

			//done for this step. Make it available.
			available = true;

			iteration++;
			stageiteration++;
		}
		break;
		}
	}
	break;
	}
}

#endif