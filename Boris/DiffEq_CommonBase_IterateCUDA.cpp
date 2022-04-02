#include "stdafx.h"
#include "DiffEq_CommonBase.h"
#include "DiffEq_Common.h"
#include "Atom_DiffEq_Common.h"

#include "DiffEq.h"
#include "Atom_DiffEq.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if COMPILECUDA == 1

void ODECommon_Base::IterateCUDA(void)
{
	switch (evalMethod) {

	case EVAL_EULER:
	{
#ifdef ODE_EVAL_COMPILATION_EULER
		if (calculate_mxh && calculate_dmdt) podeSolver->pODECUDA->Zero_reduction_values();
		else if (calculate_mxh) podeSolver->pODECUDA->Zero_mxh_lte_values();
		else if (calculate_dmdt) podeSolver->pODECUDA->Zero_dmdt_lte_values();

		for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

			//Euler can be used for stochastic equations
			if (podeSolver->pODE[idx]->H_Thermal.linear_size() && podeSolver->pODE[idx]->Torque_Thermal.linear_size()) podeSolver->pODE[idx]->pmeshODECUDA->GenerateThermalField_and_Torque();
			else if (podeSolver->pODE[idx]->H_Thermal.linear_size()) podeSolver->pODE[idx]->pmeshODECUDA->GenerateThermalField();

			if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunEuler_LLG(calculate_mxh, calculate_dmdt);
			else podeSolver->pODE[idx]->pmeshODECUDA->RunEuler(calculate_mxh, calculate_dmdt);
		}

		for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

			//Euler can be used for stochastic equations
			if (patom_odeSolver->pODE[idx]->H_Thermal.linear_size()) patom_odeSolver->pODE[idx]->pameshODECUDA->GenerateThermalField();

			if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunEuler_LLG(calculate_mxh, calculate_dmdt);
			else patom_odeSolver->pODE[idx]->pameshODECUDA->RunEuler(calculate_mxh, calculate_dmdt);
		}

		if (calculate_mxh) {

			calculate_mxh = false;
			podeSolver->pODECUDA->mxhav_to_mxh();
		}

		if (calculate_dmdt) {

			calculate_dmdt = false;
			podeSolver->pODECUDA->dmdtav_to_dmdt();
		}

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
#ifdef ODE_EVAL_COMPILATION_TEULER
		if (evalStep == 0) {

			if (calculate_mxh) podeSolver->pODECUDA->Zero_mxh_lte_values();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				//Trapezoidal Euler can be used for stochastic equations
				if (podeSolver->pODE[idx]->H_Thermal.linear_size() && podeSolver->pODE[idx]->Torque_Thermal.linear_size()) podeSolver->pODE[idx]->pmeshODECUDA->GenerateThermalField_and_Torque();
				else if (podeSolver->pODE[idx]->H_Thermal.linear_size()) podeSolver->pODE[idx]->pmeshODECUDA->GenerateThermalField();

				if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunTEuler_LLG(0, calculate_mxh, calculate_dmdt);
				else podeSolver->pODE[idx]->pmeshODECUDA->RunTEuler(0, calculate_mxh, calculate_dmdt);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				//Trapezoidal Euler can be used for stochastic equations
				if (patom_odeSolver->pODE[idx]->H_Thermal.linear_size()) patom_odeSolver->pODE[idx]->pameshODECUDA->GenerateThermalField();

				if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunTEuler_LLG(0, calculate_mxh, calculate_dmdt);
				else patom_odeSolver->pODE[idx]->pameshODECUDA->RunTEuler(0, calculate_mxh, calculate_dmdt);
			}

			if (calculate_mxh) {

				calculate_mxh = false;
				podeSolver->pODECUDA->mxhav_to_mxh();
			}

			evalStep = 1;
			available = false;
		}
		else {

			if (calculate_dmdt) podeSolver->pODECUDA->Zero_dmdt_lte_values();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunTEuler_LLG(1, calculate_mxh, calculate_dmdt);
				else podeSolver->pODE[idx]->pmeshODECUDA->RunTEuler(1, calculate_mxh, calculate_dmdt);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunTEuler_LLG(1, calculate_mxh, calculate_dmdt);
				else patom_odeSolver->pODE[idx]->pameshODECUDA->RunTEuler(1, calculate_mxh, calculate_dmdt);
			}

			if (calculate_dmdt) {

				calculate_dmdt = false;
				podeSolver->pODECUDA->dmdtav_to_dmdt();
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

			if (calculate_mxh) podeSolver->pODECUDA->Zero_mxh_lte_values();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				//Adaptive Heun can be used for stochastic equations
				if (podeSolver->pODE[idx]->H_Thermal.linear_size() && podeSolver->pODE[idx]->Torque_Thermal.linear_size()) podeSolver->pODE[idx]->pmeshODECUDA->GenerateThermalField_and_Torque();
				else if (podeSolver->pODE[idx]->H_Thermal.linear_size()) podeSolver->pODE[idx]->pmeshODECUDA->GenerateThermalField();

				if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunAHeun_LLG(0, calculate_mxh, calculate_dmdt);
				else podeSolver->pODE[idx]->pmeshODECUDA->RunAHeun(0, calculate_mxh, calculate_dmdt);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				//Adaptive Heun can be used for stochastic equations
				if (patom_odeSolver->pODE[idx]->H_Thermal.linear_size()) patom_odeSolver->pODE[idx]->pameshODECUDA->GenerateThermalField();

				if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunAHeun_LLG(0, calculate_mxh, calculate_dmdt);
				else patom_odeSolver->pODE[idx]->pameshODECUDA->RunAHeun(0, calculate_mxh, calculate_dmdt);
			}

			if (calculate_mxh) {

				calculate_mxh = false;
				podeSolver->pODECUDA->mxhav_to_mxh();
			}

			evalStep = 1;
			available = false;
		}
		else {

			if (calculate_dmdt) podeSolver->pODECUDA->Zero_dmdt_lte_values();
			//need this for lte reduction
			else podeSolver->pODECUDA->Zero_lte_value();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunAHeun_LLG(1, calculate_mxh, calculate_dmdt);
				else podeSolver->pODE[idx]->pmeshODECUDA->RunAHeun(1, calculate_mxh, calculate_dmdt);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunAHeun_LLG(1, calculate_mxh, calculate_dmdt);
				else patom_odeSolver->pODE[idx]->pameshODECUDA->RunAHeun(1, calculate_mxh, calculate_dmdt);
			}

			if (calculate_dmdt) {

				calculate_dmdt = false;
				podeSolver->pODECUDA->dmdtav_to_dmdt();
			}

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;
			available = true;

			dT_last = dT;
			podeSolver->pODECUDA->Sync_dT_last();
			lte = podeSolver->pODECUDA->Get_lte();

			if (!SetAdaptiveTimeStep()) {

				podeSolver->RestoreCUDA();
				patom_odeSolver->RestoreCUDA();
			}

			podeSolver->pODECUDA->Sync_dT();
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
			if (calculate_mxh) podeSolver->pODECUDA->Zero_mxh_lte_values();

			bool stochastic = false;

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				//RK4 can be used for stochastic equations
				if (podeSolver->pODE[idx]->H_Thermal.linear_size() && podeSolver->pODE[idx]->Torque_Thermal.linear_size()) { podeSolver->pODE[idx]->pmeshODECUDA->GenerateThermalField_and_Torque(); stochastic = true; }
				else if (podeSolver->pODE[idx]->H_Thermal.linear_size()) { podeSolver->pODE[idx]->pmeshODECUDA->GenerateThermalField(); stochastic = true; }
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				//RK4 can be used for stochastic equations
				if (patom_odeSolver->pODE[idx]->H_Thermal.linear_size()) { patom_odeSolver->pODE[idx]->pameshODECUDA->GenerateThermalField(); stochastic = true; }
			}

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunRK4_LLG(0, calculate_mxh, calculate_dmdt);
				else podeSolver->pODE[idx]->pmeshODECUDA->RunRK4(0, calculate_mxh, calculate_dmdt, stochastic);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunRK4_LLG(0, calculate_mxh, calculate_dmdt);
				else patom_odeSolver->pODE[idx]->pameshODECUDA->RunRK4(0, calculate_mxh, calculate_dmdt, stochastic);
			}

			if (calculate_mxh) {

				calculate_mxh = false;
				if (stochastic) podeSolver->pODECUDA->mxhav_to_mxh();
			}

			evalStep++;
			available = false;
		}
		break;

		case 1:
		case 2:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunRK4_LLG(evalStep);
				else podeSolver->pODE[idx]->pmeshODECUDA->RunRK4(evalStep);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunRK4_LLG(evalStep);
				else patom_odeSolver->pODE[idx]->pameshODECUDA->RunRK4(evalStep);
			}

			evalStep++;
		}
		break;

		case 3:
		{
			if (calculate_dmdt) podeSolver->pODECUDA->Zero_dmdt_lte_values();

			bool stochastic = false;

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) stochastic |= podeSolver->pODE[idx]->H_Thermal.linear_size() != 0;
			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) stochastic |= patom_odeSolver->pODE[idx]->H_Thermal.linear_size() != 0;

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunRK4_LLG(3, calculate_mxh, calculate_dmdt);
				else podeSolver->pODE[idx]->pmeshODECUDA->RunRK4(3, calculate_mxh, calculate_dmdt, stochastic);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunRK4_LLG(3, calculate_mxh, calculate_dmdt);
				else patom_odeSolver->pODE[idx]->pameshODECUDA->RunRK4(3, calculate_mxh, calculate_dmdt, stochastic);
			}

			if (calculate_dmdt) {

				calculate_dmdt = false;
				if (stochastic) podeSolver->pODECUDA->dmdtav_to_dmdt();
			}

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
#ifdef ODE_EVAL_COMPILATION_ABM
		if (primed) {

			if (evalStep == 0) {

				if (calculate_mxh) podeSolver->pODECUDA->Zero_mxh_lte_values();

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunABM_LLG(0, calculate_mxh, calculate_dmdt);
					else podeSolver->pODE[idx]->pmeshODECUDA->RunABM(0, calculate_mxh, calculate_dmdt);
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunABM_LLG(0, calculate_mxh, calculate_dmdt);
					else patom_odeSolver->pODE[idx]->pameshODECUDA->RunABM(0, calculate_mxh, calculate_dmdt);
				}

				calculate_mxh = false;

				evalStep = 1;
				available = false;
			}
			else {

				if (calculate_dmdt) podeSolver->pODECUDA->Zero_dmdt_lte_values();
				//need this for lte reduction
				else podeSolver->pODECUDA->Zero_lte_value();

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunABM_LLG(1, calculate_mxh, calculate_dmdt);
					else podeSolver->pODE[idx]->pmeshODECUDA->RunABM(1, calculate_mxh, calculate_dmdt);
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunABM_LLG(1, calculate_mxh, calculate_dmdt);
					else patom_odeSolver->pODE[idx]->pameshODECUDA->RunABM(1, calculate_mxh, calculate_dmdt);
				}

				calculate_dmdt = false;

				evalStep = 0;
				time += dT;
				stagetime += dT;
				iteration++;
				stageiteration++;
				available = true;

				dT_last = dT;
				podeSolver->pODECUDA->Sync_dT_last();
				lte = podeSolver->pODECUDA->Get_lte();

				if (!SetAdaptiveTimeStep()) {

					podeSolver->RestoreCUDA();
					patom_odeSolver->RestoreCUDA();
				}
				else {

					//switch between saved equation evaluations (use latest)
					alternator = !alternator;
					podeSolver->pODECUDA->Sync_alternator();
				}

				podeSolver->pODECUDA->Sync_dT();
			}
		}
		else {

			if (evalStep == 0) {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunABMTEuler_LLG(0);
					else podeSolver->pODE[idx]->pmeshODECUDA->RunABMTEuler(0);
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunABMTEuler_LLG(0);
					else patom_odeSolver->pODE[idx]->pameshODECUDA->RunABMTEuler(0);
				}

				evalStep = 1;
				available = false;
			}
			else {

				for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

					if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunABMTEuler_LLG(1);
					else podeSolver->pODE[idx]->pmeshODECUDA->RunABMTEuler(1);
				}

				for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

					if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunABMTEuler_LLG(1);
					else patom_odeSolver->pODE[idx]->pameshODECUDA->RunABMTEuler(1);
				}

				evalStep = 0;
				available = true;

				time += dT;
				stagetime += dT;
				iteration++;
				stageiteration++;

				primed = true;
				alternator = false;		//evaluated equation is now in sEval0
				podeSolver->pODECUDA->Sync_alternator();
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
			if (calculate_mxh) podeSolver->pODECUDA->Zero_mxh_lte_values();
			//need this for lte reduction
			else podeSolver->pODECUDA->Zero_lte_value();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunRK23_Step0_NoAdvance(calculate_mxh);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunRK23_Step0_NoAdvance(calculate_mxh);
			}

			calculate_mxh = false;

			available = false;

			//RK23 has the FSAL property, thus the error is calculated on this first stage -> must be primed by a full pass first
			//It also means the above methods do not advance the magnetization yet as we need to calculate the new stepsize first
			//Magnetization is advanced below if the step has not failed
			//If the step has failed then restore magnetization, and also set primed to false -> we must now take a full pass again before we can calculate a new stepsize.
			if (primed) {

				dT_last = dT;
				podeSolver->pODECUDA->Sync_dT_last();
				lte = podeSolver->pODECUDA->Get_lte();

				if (!SetAdaptiveTimeStep()) {

					podeSolver->RestoreCUDA();
					patom_odeSolver->RestoreCUDA();
					
					primed = false;
					podeSolver->pODECUDA->Sync_dT();
					break;
				}

				podeSolver->pODECUDA->Sync_dT();
			}
			else primed = true;

			//Advance with new stepsize
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunRK23(0);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunRK23(0);
			}

			evalStep++;
		}
		break;

		case 1:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunRK23(1);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunRK23(1);
			}

			evalStep++;
		}
		break;

		case 2:
		{
			if (calculate_dmdt) podeSolver->pODECUDA->Zero_dmdt_lte_values();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunRK23(2, calculate_mxh, calculate_dmdt);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunRK23(2, calculate_mxh, calculate_dmdt);
			}

			calculate_dmdt = false;

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

	case EVAL_RKF45:
	{
#ifdef ODE_EVAL_COMPILATION_RKF45
		switch (evalStep) {

		case 0:
		{
			if (calculate_mxh) podeSolver->pODECUDA->Zero_mxh_lte_values();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunRKF45_LLG(0, calculate_mxh, calculate_dmdt);
				else podeSolver->pODE[idx]->pmeshODECUDA->RunRKF45(0, calculate_mxh, calculate_dmdt);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKF45_LLG(0, calculate_mxh, calculate_dmdt);
				else patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKF45(0, calculate_mxh, calculate_dmdt);
			}

			calculate_mxh = false;

			evalStep++;
			available = false;
		}
		break;

		case 1:
		case 2:
		case 3:
		case 4:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunRKF45_LLG(evalStep);
				else podeSolver->pODE[idx]->pmeshODECUDA->RunRKF45(evalStep);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKF45_LLG(evalStep);
				else patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKF45(evalStep);
			}

			evalStep++;
		}
		break;

		case 5:
		{
			if (calculate_dmdt) podeSolver->pODECUDA->Zero_dmdt_lte_values();
			//need this for lte reduction
			else podeSolver->pODECUDA->Zero_lte_value();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				if (podeSolver->setODE == ODE_LLG) podeSolver->pODE[idx]->pmeshODECUDA->RunRKF45_LLG(5, calculate_mxh, calculate_dmdt);
				else podeSolver->pODE[idx]->pmeshODECUDA->RunRKF45(5, calculate_mxh, calculate_dmdt);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				if (patom_odeSolver->setODE == ODE_LLG) patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKF45_LLG(5, calculate_mxh, calculate_dmdt);
				else patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKF45(5, calculate_mxh, calculate_dmdt);
			}

			calculate_dmdt = false;

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;
			available = true;

			dT_last = dT;
			podeSolver->pODECUDA->Sync_dT_last();
			lte = podeSolver->pODECUDA->Get_lte();

			if (!SetAdaptiveTimeStep()) {

				podeSolver->RestoreCUDA();
				patom_odeSolver->RestoreCUDA();
			}

			podeSolver->pODECUDA->Sync_dT();
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
			if (calculate_mxh) podeSolver->pODECUDA->Zero_mxh_lte_values();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunRKF56(0, calculate_mxh, calculate_dmdt);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKF56(0, calculate_mxh, calculate_dmdt);
			}

			calculate_mxh = false;

			evalStep++;
			available = false;
		}
		break;

		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
		case 6:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunRKF56(evalStep);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKF56(evalStep);
			}

			evalStep++;
		}
		break;

		case 7:
		{
			if (calculate_dmdt) podeSolver->pODECUDA->Zero_dmdt_lte_values();
			//need this for lte reduction
			else podeSolver->pODECUDA->Zero_lte_value();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunRKF56(7, calculate_mxh, calculate_dmdt);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKF56(7, calculate_mxh, calculate_dmdt);
			}

			calculate_dmdt = false;

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;
			available = true;

			dT_last = dT;
			podeSolver->pODECUDA->Sync_dT_last();
			lte = podeSolver->pODECUDA->Get_lte();

			if (!SetAdaptiveTimeStep()) {

				podeSolver->RestoreCUDA();
				patom_odeSolver->RestoreCUDA();
			}

			podeSolver->pODECUDA->Sync_dT();
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
			if (calculate_mxh) podeSolver->pODECUDA->Zero_mxh_lte_values();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunRKCK45(0, calculate_mxh, calculate_dmdt);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKCK45(0, calculate_mxh, calculate_dmdt);
			}

			calculate_mxh = false;

			evalStep++;
			available = false;
		}
		break;

		case 1:
		case 2:
		case 3:
		case 4:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunRKCK45(evalStep);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKCK45(evalStep);
			}

			evalStep++;
		}
		break;

		case 5:
		{
			if (calculate_dmdt) podeSolver->pODECUDA->Zero_dmdt_lte_values();
			//need this for lte reduction
			else podeSolver->pODECUDA->Zero_lte_value();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunRKCK45(5, calculate_mxh, calculate_dmdt);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKCK45(5, calculate_mxh, calculate_dmdt);
			}

			calculate_dmdt = false;

			evalStep = 0;
			time += dT;
			stagetime += dT;
			iteration++;
			stageiteration++;
			available = true;

			dT_last = dT;
			podeSolver->pODECUDA->Sync_dT_last();
			lte = podeSolver->pODECUDA->Get_lte();

			if (!SetAdaptiveTimeStep()) {

				podeSolver->RestoreCUDA();
				patom_odeSolver->RestoreCUDA();
			}

			podeSolver->pODECUDA->Sync_dT();
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
			if (calculate_mxh) podeSolver->pODECUDA->Zero_mxh_lte_values();
			//need this for lte reduction
			else podeSolver->pODECUDA->Zero_lte_value();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunRKDP54_Step0_NoAdvance(calculate_mxh);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKDP54_Step0_NoAdvance(calculate_mxh);
			}

			calculate_mxh = false;

			available = false;

			//RKDP has the FSAL property, thus the error is calculated on this first stage -> must be primed by a full pass first
			//It also means the above methods do not advance the magnetization yet as we need to calculate the new stepsize first
			//Magnetization is advanced below if the step has not failed
			//If the step has failed then restore magnetization, and also set primed to false -> we must now take a full pass again before we can calculate a new stepsize.
			if (primed) {

				dT_last = dT;
				podeSolver->pODECUDA->Sync_dT_last();
				lte = podeSolver->pODECUDA->Get_lte();

				if (!SetAdaptiveTimeStep()) {

					podeSolver->RestoreCUDA();
					patom_odeSolver->RestoreCUDA();

					primed = false;
					podeSolver->pODECUDA->Sync_dT();
					break;
				}

				podeSolver->pODECUDA->Sync_dT();
			}
			else primed = true;

			//Advance with new stepsize
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunRKDP54(0);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKDP54(0);
			}

			evalStep++;
		}
		break;

		case 1:
		case 2:
		case 3:
		case 4:
		{
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunRKDP54(evalStep);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKDP54(evalStep);
			}

			evalStep++;
		}
		break;

		case 5:
		{
			if (calculate_dmdt) podeSolver->pODECUDA->Zero_dmdt_lte_values();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunRKDP54(5, calculate_mxh, calculate_dmdt);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunRKDP54(5, calculate_mxh, calculate_dmdt);
			}

			calculate_dmdt = false;

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

	case EVAL_SD:
	{
#ifdef ODE_EVAL_COMPILATION_SD
		if (primed) {

			bool do_sd_reset = false;

			//1. calculate parameters for Barzilai-Borwein stepsizes -> solver must be primed with step 0 (after it is primed next loop starts from step 1)
			//must reset the static delta_... quantities before running these across all meshes

			//zero all quantities used for Barzilai-Borwein stepsize calculations
			if (podeSolver->pODE.size()) podeSolver->pODECUDA->Zero_SD_Solver_BB_Values();
			if (patom_odeSolver->pODE.size()) patom_odeSolver->pODECUDA->Zero_SD_Solver_BB_Values();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunSD_BB();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunSD_BB();
			}

			//2. Set stepsize - alternate between BB values
			//first transfer BB values to cpu memory
			if (podeSolver->pODE.size()) podeSolver->pODECUDA->Get_SD_Solver_BB_Values(&podeSolver->delta_m_sq, &podeSolver->delta_G_sq, &podeSolver->delta_m_dot_delta_G, &podeSolver->delta_m2_sq, &podeSolver->delta_G2_sq, &podeSolver->delta_m2_dot_delta_G2);
			if (patom_odeSolver->pODE.size()) patom_odeSolver->pODECUDA->Get_SD_Solver_BB_Values(&patom_odeSolver->delta_m_sq, &patom_odeSolver->delta_G_sq, &patom_odeSolver->delta_m_dot_delta_G);

			double delta_m_sq = podeSolver->delta_m_sq + podeSolver->delta_m2_sq + patom_odeSolver->delta_m_sq;
			double delta_m_G = podeSolver->delta_m_dot_delta_G + podeSolver->delta_m2_dot_delta_G2 + patom_odeSolver->delta_m_dot_delta_G;
			double delta_G_sq = podeSolver->delta_G_sq + podeSolver->delta_G2_sq + patom_odeSolver->delta_G_sq;

			//2. Set stepsize - alternate between BB values
			if (iteration % 2) {

				if (delta_m_G && delta_m_sq * delta_m_G > 0.0) dT = delta_m_sq / delta_m_G;
				else if (delta_G_sq && delta_m_G * delta_G_sq > 0.0) dT = delta_m_G / delta_G_sq;
				else do_sd_reset = true;
			}
			else {

				if (delta_G_sq && delta_m_G * delta_G_sq > 0.0) dT = delta_m_G / delta_G_sq;
				else if (delta_m_G && delta_m_sq * delta_m_G > 0.0) dT = delta_m_sq / delta_m_G;
				else do_sd_reset = true;
			}

			if (dT < dT_min) dT = dT_min;
			if (dT > dT_max) dT = dT_max;

			if (do_sd_reset) {

				dT = (++sd_reset_consecutive_iters) * dT_min;
				if (dT > dT_max) dT = dT_max;
			}
			else sd_reset_consecutive_iters = 0;

			//make sure to transfer dT value to GPU
			podeSolver->pODECUDA->Sync_dT();

			//3. set new magnetization vectors
			if (calculate_mxh && calculate_dmdt) podeSolver->pODECUDA->Zero_reduction_values();
			else if (calculate_mxh) podeSolver->pODECUDA->Zero_mxh_lte_values();
			else if (calculate_dmdt) podeSolver->pODECUDA->Zero_dmdt_lte_values();

			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunSD_Advance(calculate_mxh, calculate_dmdt);
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunSD_Advance(calculate_mxh, calculate_dmdt);
			}

			calculate_mxh = false;
			calculate_dmdt = false;

			iteration++;
			stageiteration++;
			time += dT;
			stagetime += dT;
		}
		else {

			dT = dT_min;

			//make sure to transfer dT value to GPU
			podeSolver->pODECUDA->Sync_dT();

			//0. prime the SD solver
			for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

				podeSolver->pODE[idx]->pmeshODECUDA->RunSD_Start();
			}

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				patom_odeSolver->pODE[idx]->pameshODECUDA->RunSD_Start();
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

	podeSolver->pODECUDA->Sync_time();
}

#endif

