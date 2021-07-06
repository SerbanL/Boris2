#include "stdafx.h"
#include "DiffEqFM.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC

#include "Mesh_Ferromagnetic.h"
#include "MeshParamsControl.h"

DifferentialEquationFM::DifferentialEquationFM(FMesh *pMesh):
	DifferentialEquation(pMesh)
{
	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
}

DifferentialEquationFM::~DifferentialEquationFM()
{
}

//---------------------------------------- OTHERS

//Restore magnetization after a failed step for adaptive time-step methods
void DifferentialEquationFM::RestoreMagnetization(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++)
		pMesh->M[idx] = sM1[idx];
}

//renormalize vectors to set magnetization length value (which could have a spatial variation)
void DifferentialEquationFM::RenormalizeMagnetization(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->M.linear_size(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			double Ms = pMesh->Ms;
			pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms);
 
			pMesh->M[idx].renormalize(Ms);
		}
	}
}

//---------------------------------------- SET-UP METHODS

BError DifferentialEquationFM::AllocateMemory(void)
{
	BError error(CLASS_STR(DifferentialEquationFM));

	//first make sure everything not needed is cleared
	CleanupMemory();

	if (!sM1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);

	switch (evalMethod) {

	case EVAL_EULER:
		break;

	case EVAL_TEULER:
	case EVAL_AHEUN:
		break;

	case EVAL_ABM:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RK4:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RK23:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RKF45:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval3.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval4.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RKCK45:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval3.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval4.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RKDP54:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval3.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval4.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval5.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_SD:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;
	}

	//For stochastic equations must also allocate memory for thermal VECs

	switch (setODE) {

	case ODE_SLLG:
	case ODE_SLLGSTT:
	case ODE_SLLGSA:
		if (!H_Thermal.resize(pMesh->h_s, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case ODE_SLLB:
	case ODE_SLLBSTT:
	case ODE_SLLBSA:
		if (!H_Thermal.resize(pMesh->h_s, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!Torque_Thermal.resize(pMesh->h_s, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;
	}

	//----------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pmeshODECUDA) {

		if (!error) error = pmeshODECUDA->AllocateMemory();
	}
#endif

	return error;
}

void DifferentialEquationFM::CleanupMemory(void)
{
	//Only clear vectors not used for current evaluation method
	sM1.clear();

	if (evalMethod != EVAL_RK4 &&
		evalMethod != EVAL_ABM &&
		evalMethod != EVAL_RK23 &&
		evalMethod != EVAL_RKF45 &&
		evalMethod != EVAL_RKCK45 &&
		evalMethod != EVAL_RKDP54 &&
		evalMethod != EVAL_SD) {

		sEval0.clear();
	}

	if (evalMethod != EVAL_RK4 &&
		evalMethod != EVAL_ABM &&
		evalMethod != EVAL_RK23 &&
		evalMethod != EVAL_RKF45 &&
		evalMethod != EVAL_RKCK45 &&
		evalMethod != EVAL_RKDP54) {

		sEval1.clear();
	}

	if (evalMethod != EVAL_RK4 &&
		evalMethod != EVAL_RK23 &&
		evalMethod != EVAL_RKF45 &&
		evalMethod != EVAL_RKCK45 &&
		evalMethod != EVAL_RKDP54) {

		sEval2.clear();
	}

	if (evalMethod != EVAL_RK4 &&
		evalMethod != EVAL_RKF45 &&
		evalMethod != EVAL_RKCK45 &&
		evalMethod != EVAL_RKDP54) {

		sEval3.clear();
		sEval4.clear();
	}

	if (evalMethod != EVAL_RKDP54) {

		sEval5.clear();
	}

	//For thermal vecs only clear if not used for current set ODE
	if (setODE != ODE_SLLG &&
		setODE != ODE_SLLGSTT &&
		setODE != ODE_SLLB &&
		setODE != ODE_SLLBSTT &&
		setODE != ODE_SLLGSA &&
		setODE != ODE_SLLBSA) {

		H_Thermal.clear();
	}

	if (setODE != ODE_SLLB &&
		setODE != ODE_SLLBSTT &&
		setODE != ODE_SLLBSA) {

		Torque_Thermal.clear();
	}

	//----------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pmeshODECUDA) { pmeshODECUDA->CleanupMemory(); }
#endif
}

BError DifferentialEquationFM::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(DifferentialEquationFM));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_ODE_SOLVER)) {

		if (pMesh->link_stochastic) {

			pMesh->h_s = pMesh->h;
			pMesh->n_s = pMesh->n;
		}
		else {

			pMesh->n_s = round(pMesh->meshRect / pMesh->h_s);
			if (pMesh->n_s.x == 0) pMesh->n_s.x = 1;
			if (pMesh->n_s.y == 0) pMesh->n_s.y = 1;
			if (pMesh->n_s.z == 0) pMesh->n_s.z = 1;
			pMesh->h_s = pMesh->meshRect / pMesh->n_s;
		}

		error = AllocateMemory();
	}

	if (cfgMessage == UPDATECONFIG_PARAMVALUECHANGED_MLENGTH) RenormalizeMagnetization();

	//----------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pmeshODECUDA) {

		if (!error) error = pmeshODECUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

//switch CUDA state on/off
BError DifferentialEquationFM::SwitchCUDAState(bool cudaState)
{
	BError error(CLASS_STR(DifferentialEquationFM));

#if COMPILECUDA == 1

	//are we switching to cuda?
	if (cudaState) {

		if (!pmeshODECUDA) {

			pmeshODECUDA = new DifferentialEquationFMCUDA(this);
			error = pmeshODECUDA->Error_On_Create();
		}
	}
	else {

		//cuda switched off so delete cuda module object
		if (pmeshODECUDA) delete pmeshODECUDA;
		pmeshODECUDA = nullptr;
	}

#endif

	return error;
}

//---------------------------------------- GETTERS

//return dM by dT - should only be used when evaluation sequence has ended (TimeStepSolved() == true)
DBL3 DifferentialEquationFM::dMdt(int idx)
{
	return (pMesh->M[idx] - sM1[idx]) / dT_last;
}

#endif
