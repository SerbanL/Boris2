#include "stdafx.h"
#include "DiffEq.h"
#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"

//---------------------------------------- ITERATE METHODS

void ODECommon::Iterate(void)
{
	//save current dT value in case it changes (adaptive time step methods)
	dT_last = dT;

	switch (evalMethod) {

	case EVAL_EULER:
	{
		for (int idx = 0; idx < (int)pODE.size(); idx++)
			pODE[idx]->RunEuler();

		time += dT;
		stagetime += dT;
		iteration++;
		stageiteration++;

		Set_mxh();
	}
	break;

	case EVAL_TEULER:
	{
		if (evalStep == 0) {

			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunTEuler_Step0();

			evalStep = 1;
			available = false;
		}
		else {

			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunTEuler_Step1();

			evalStep = 0;
			available = true;

			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;

			Set_mxh();
		}
	}
	break;

	case EVAL_RK4:
	{
		switch (evalStep) {

		case 0:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRK4_Step0();

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
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRK4_Step3();

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

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunABM_Predictor();

				evalStep = 1;
				available = false;
			}
			else {

				for (int idx = 0; idx < (int)pODE.size(); idx++)
					pODE[idx]->RunABM_Corrector();

				evalStep = 0;
				time += dT;
				stagetime += dT;

				double lte = Get_lte();

				//adaptive time step based on lte - is lte over acceptable relative error?
				if (lte > ABM_RELERRMAX) {

					//failed - try again. Restore state to start of evaluation and try again with smaller time step - output not available yet.
					//do not redo if time step is at or below the minimum value allowed - this is used as a "timeout" to stop the solver from getting stuck on the same iteration; in this case the method is not good enough for the problem.
					if (lte > ABM_RELERRFAIL && dT > ABM_MINDT) {

						time -= dT;
						stagetime -= dT;

						for (int idx = 0; idx < (int)pODE.size(); idx++)
							pODE[idx]->RestoreMagnetisation();

						dT *= ABM_DTREDUCE;
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

				Set_mxh();

				//done for this step. Make it available.
				available = true;

				iteration++;
				stageiteration++;

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
	}
	break;

	case EVAL_RKF:
	{
		switch (evalStep) {

		case 0:
		{
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKF45_Step0();

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
			for (int idx = 0; idx < (int)pODE.size(); idx++)
				pODE[idx]->RunRKF45_Step5();

			evalStep = 0;
			time += dT;
			stagetime += dT;

			double lte = Get_lte();

			//adaptive time step based on lte - is lte over acceptable relative error?
			if (lte > RKF_RELERRMAX) {

				//failed - try again. Restore state to start of evaluation and try again with smaller time step - output not available yet.
				//do not redo if time step is at or below the minimum value allowed - this is used as a "timeout" to stop the solver from getting stuck on the same iteration; in this case the method is not good enough for the problem.
				if (lte > RKF_RELERRFAIL && dT > RKF_MINDT) {

					time -= dT;
					stagetime -= dT;

					for (int idx = 0; idx < (int)pODE.size(); idx++)
						pODE[idx]->RestoreMagnetisation();

					dT *= RKF_DTREDUCE;
					break;
				}

				//not failed but still need to reduce time step
				dT *= RKF_DTREDUCE;

				//must not go below minimum time step
				if (dT < RKF_MINDT) dT = RKF_MINDT;
			}

			//is lte below minimum relative error ? If yes we can go quicker
			if (lte < RKF_RELERRMIN) {

				dT *= RKF_DTINCREASE;

				//must not go above maximum time step
				if (dT > RKF_MAXDT) dT = RKF_MAXDT;
			}

			Set_mxh();

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