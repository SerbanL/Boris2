#include "stdafx.h"
#include "DiffEq_CommonBase.h"

#include "DiffEq_Common.h"
#include "Atom_DiffEq_Common.h"

//----------------------------------- Evaluation Method and Control

BError ODECommon_Base::SetEvaluationMethod(EVAL_ evalMethod_)
{
	BError error(__FUNCTION__);

	evalMethod = evalMethod_;

	//set default parameters for given evaluation method
	switch (evalMethod) {

	case EVAL_EULER:
	{
		dT = EULER_DEFAULT_DT;
		eval_method_order = 1;
	}
	break;

	case EVAL_TEULER:
	{
		dT = TEULER_DEFAULT_DT;
		eval_method_order = 2;
	}
	break;

	case EVAL_AHEUN:
	{
		dT = AHEUN_DEFAULT_DT;

		err_high_fail = AHEUN_RELERRFAIL;
		dT_increase = AHEUN_DTINCREASE;
		dT_max = AHEUN_MAXDT;
		dT_min = AHEUN_MINDT;
		eval_method_order = 2;
	}
	break;

	case EVAL_RK4:
	{
		dT = RK4_DEFAULT_DT;
		eval_method_order = 4;
	}
	break;

	case EVAL_ABM:
	{
		dT = ABM_DEFAULT_DT;

		err_high_fail = ABM_RELERRFAIL;
		dT_increase = ABM_DTINCREASE;
		dT_max = ABM_MAXDT;
		dT_min = ABM_MINDT;
		eval_method_order = 2;
	}
	break;

	case EVAL_RK23:
	{
		dT = RK23_DEFAULT_DT;

		err_high_fail = RK23_RELERRFAIL;
		dT_increase = RK23_DTINCREASE;
		dT_max = RK23_MAXDT;
		dT_min = RK23_MINDT;
		eval_method_order = 3;
	}
	break;

	case EVAL_SD:
	{
		//set starting dT - set a very conservative initial stepsize otherwise solver priming can be bad.
		dT = SD_DEFAULT_DT;
		//reset to dT_min when needed
		dT_min = SD_MINDT;
		dT_max = SD_MAXDT;
		eval_method_order = 1;
	}
	break;

	default:
	case EVAL_RKF:
	{
		dT = RKF_DEFAULT_DT;

		err_high_fail = RKF_RELERRFAIL;
		dT_increase = RKF_DTINCREASE;
		dT_max = RKF_MAXDT;
		dT_min = RKF_MINDT;
		eval_method_order = 5;
	}
	break;

	case EVAL_RKCK:
	{
		dT = RKCK_DEFAULT_DT;

		err_high_fail = RKCK_RELERRFAIL;
		dT_increase = RKCK_DTINCREASE;
		dT_max = RKCK_MAXDT;
		dT_min = RKCK_MINDT;
		eval_method_order = 4;
	}
	break;

	case EVAL_RKDP:
	{
		dT = RKDP_DEFAULT_DT;

		err_high_fail = RKDP_RELERRFAIL;
		dT_increase = RKDP_DTINCREASE;
		dT_max = RKDP_MAXDT;
		dT_min = RKDP_MINDT;
		eval_method_order = 5;
	}
	break;
	}

	//initial settings
	dT_last = dT;

	mxh = 1;
	dmdt = 1;

	available = true;
	evalStep = 0;

	alternator = false;
	primed = false;

	calculate_mxh = true;
	calculate_dmdt = true;

	if (link_dTspeedup) dTspeedup = dT;

#if COMPILECUDA == 1
	if (podeSolver->pODECUDA) podeSolver->pODECUDA->SyncODEValues();
	if (patom_odeSolver->pODECUDA) patom_odeSolver->pODECUDA->SyncODEValues();
#endif

	return error;
}

void ODECommon_Base::Reset(void)
{
	iteration = 0;
	stageiteration = 0;
	time = 0.0;
	stagetime = 0.0;

	time_stoch = 0.0;

	time_speedup = 0.0;

	mxh = 1.0;
	dmdt = 1.0;

	available = true;
	evalStep = 0;

	alternator = false;
	primed = false;

	calculate_mxh = true;
	calculate_dmdt = true;

	moving_mesh_dwshift = 0.0;

#if COMPILECUDA == 1
	if (podeSolver->pODECUDA) podeSolver->pODECUDA->SyncODEValues();
	if (patom_odeSolver->pODECUDA) patom_odeSolver->pODECUDA->SyncODEValues();
#endif
}

void ODECommon_Base::NewStage(void)
{
	stagetime = 0;
	stageiteration = 0;

	available = true;
	evalStep = 0;

	mxh = 1;
	dmdt = 1;

	alternator = false;
	primed = false;

	calculate_mxh = true;
	calculate_dmdt = true;

#if COMPILECUDA == 1
	if (podeSolver->pODECUDA) podeSolver->pODECUDA->SyncODEValues();
	if (patom_odeSolver->pODECUDA) patom_odeSolver->pODECUDA->SyncODEValues();
#endif
}

void ODECommon_Base::SetdT(double dT)
{
	this->dT = dT;

	if (link_dTspeedup) dTspeedup = dT;

#if COMPILECUDA == 1
	if (podeSolver->pODECUDA) podeSolver->pODECUDA->SyncODEValues();
	if (patom_odeSolver->pODECUDA) patom_odeSolver->pODECUDA->SyncODEValues();
#endif
}

void ODECommon_Base::SetStochTimeStep(double dTstoch_)
{
	link_dTstoch = false;

	time_stoch = time;
	dTstoch = dTstoch_;
}

void ODECommon_Base::SetLink_dTstoch(bool link_dTstoch_)
{
	link_dTstoch = link_dTstoch_;

	if (!link_dTstoch) time_stoch = time;
	else dTstoch = dT;
}

void ODECommon_Base::SetSpeedupTimeStep(double dTspeedup_)
{
	link_dTspeedup = false;

	time_speedup = time;
	dTspeedup = dTspeedup_;
}

void ODECommon_Base::SetLink_dTspeedup(bool link_dTspeedup_)
{
	link_dTspeedup = link_dTspeedup_;

	if (!link_dTspeedup) time_speedup = time;
	else dTspeedup = dT;
}

void ODECommon_Base::SetAdaptiveTimeStepCtrl(double err_high_fail, double dT_increase, double dT_min, double dT_max)
{
	this->err_high_fail = err_high_fail;
	this->dT_increase = dT_increase;
	this->dT_min = dT_min;
	this->dT_max = dT_max;
}
