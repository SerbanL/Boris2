#include "stdafx.h"
#include "Atom_DiffEqCubic.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "Atom_Mesh_Cubic.h"

Atom_DifferentialEquationCubic::Atom_DifferentialEquationCubic(Atom_Mesh_Cubic *paMesh):
	Atom_DifferentialEquation(paMesh)
{
	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
}

Atom_DifferentialEquationCubic::~Atom_DifferentialEquationCubic()
{
}

//---------------------------------------- OTHERS

//Restore magnetisation after a failed step for adaptive time-step methods
void Atom_DifferentialEquationCubic::RestoreMoments(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < paMesh->n.dim(); idx++)
		paMesh->M1[idx] = sM1[idx];
}

//---------------------------------------- SET-UP METHODS

BError Atom_DifferentialEquationCubic::AllocateMemory(void)
{
	BError error(CLASS_STR(Atom_DifferentialEquationCubic));

	//first make sure everything not needed is cleared
	CleanupMemory();

	if (!sM1.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);

	switch (evalMethod) {

	case EVAL_EULER:
		break;

	case EVAL_TEULER:
	case EVAL_AHEUN:
		break;

	case EVAL_ABM:
		if (!sEval0.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RK4:
		if (!sEval0.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RK23:
		if (!sEval0.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RKF:
		if (!sEval0.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval3.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval4.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RKCK:
		if (!sEval0.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval3.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval4.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RKDP:
		if (!sEval0.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval3.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval4.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval5.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_SD:
		if (!sEval0.resize(paMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;
	}

	//For stochastic equations must also allocate memory for thermal VECs

	switch (setODE) {

	case ODE_SLLG:
	case ODE_SLLGSTT:
	case ODE_SLLGSA:
		if (!H_Thermal.resize(paMesh->h, paMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;
	}

	//----------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pameshODECUDA) {

		if (!error) error = pameshODECUDA->AllocateMemory();
	}
#endif

	return error;
}

void Atom_DifferentialEquationCubic::CleanupMemory(void)
{
	//Only clear vectors not used for current evaluation method
	sM1.clear();

	if (evalMethod != EVAL_RK4 &&
		evalMethod != EVAL_ABM &&
		evalMethod != EVAL_RK23 &&
		evalMethod != EVAL_RKF &&
		evalMethod != EVAL_RKCK &&
		evalMethod != EVAL_RKDP &&
		evalMethod != EVAL_SD) {

		sEval0.clear();
	}

	if (evalMethod != EVAL_RK4 &&
		evalMethod != EVAL_ABM &&
		evalMethod != EVAL_RK23 &&
		evalMethod != EVAL_RKF &&
		evalMethod != EVAL_RKCK &&
		evalMethod != EVAL_RKDP) {

		sEval1.clear();
	}

	if (evalMethod != EVAL_RK4 &&
		evalMethod != EVAL_RK23 &&
		evalMethod != EVAL_RKF &&
		evalMethod != EVAL_RKCK &&
		evalMethod != EVAL_RKDP) {

		sEval2.clear();
	}

	if (evalMethod != EVAL_RK4 &&
		evalMethod != EVAL_RKF &&
		evalMethod != EVAL_RKCK &&
		evalMethod != EVAL_RKDP) {

		sEval3.clear();
		sEval4.clear();
	}

	if (evalMethod != EVAL_RKDP) {

		sEval5.clear();
	}

	//For thermal vecs only clear if not used for current set ODE
	if (setODE != ODE_SLLG && 
		setODE != ODE_SLLGSTT &&
		setODE != ODE_SLLGSA) {

		H_Thermal.clear();
	}

	//----------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pameshODECUDA) { pameshODECUDA->CleanupMemory(); }
#endif
}


BError Atom_DifferentialEquationCubic::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_DifferentialEquationCubic));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_ODE_SOLVER)) {

		error = AllocateMemory();
	}

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_ODE_MOVEMESH)) {

		if (!error) {

			//set skip cells flags for moving mesh if enabled
			if (moving_mesh) {

				Rect mesh_rect = paMesh->GetMeshRect();

				DBL3 end_size = mesh_rect.size() & DBL3(MOVEMESH_ENDRATIO, 1, 1);

				Rect end_rect_left = Rect(mesh_rect.s, mesh_rect.s + end_size);
				Rect end_rect_right = Rect(mesh_rect.e - end_size, mesh_rect.e);

				paMesh->M1.set_skipcells(end_rect_left);
				paMesh->M1.set_skipcells(end_rect_right);
			}
			else {

				paMesh->M1.clear_skipcells();
			}
		}
	}

	//----------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pameshODECUDA) {

		if (!error) error = pameshODECUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

//switch CUDA state on/off
BError Atom_DifferentialEquationCubic::SwitchCUDAState(bool cudaState)
{
	BError error(CLASS_STR(Atom_DifferentialEquationCubic));

#if COMPILECUDA == 1

	//are we switching to cuda?
	if (cudaState) {

		if (!pameshODECUDA) {

			pameshODECUDA = new Atom_DifferentialEquationCubicCUDA(this);
			error = pameshODECUDA->Error_On_Create();
		}
	}
	else {

		//cuda switched off so delete cuda module object
		if (pameshODECUDA) delete pameshODECUDA;
		pameshODECUDA = nullptr;
	}

#endif

	return error;
}

//---------------------------------------- GETTERS

//return dM by dT - should only be used when evaluation sequence has ended (TimeStepSolved() == true)
DBL3 Atom_DifferentialEquationCubic::dMdt(int idx)
{
	return (paMesh->M1[idx] - sM1[idx]) / dT_last;
}

#endif
