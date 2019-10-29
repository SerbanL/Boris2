#include "stdafx.h"
#include "DiffEqAFM.h"
#include "Mesh_AntiFerromagnetic.h"

DifferentialEquationAFM::DifferentialEquationAFM(AFMesh *pMesh) :
	DifferentialEquation(pMesh)
{
	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	OmpThreads = omp_get_num_procs();

	Equation_Eval_2.resize(OmpThreads);
}

DifferentialEquationAFM::~DifferentialEquationAFM()
{

}

//---------------------------------------- OTHERS

//Restore magnetisation after a failed step for adaptive time-step methods
void DifferentialEquationAFM::RestoreMagnetisation(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		pMesh->M[idx] = sM1[idx];
		pMesh->M2[idx] = sM1_2[idx];
	}
}

//---------------------------------------- SET-UP METHODS

BError DifferentialEquationAFM::AllocateMemory(void)
{
	BError error(CLASS_STR(DifferentialEquationAFM));

	//first make sure everything not needed is cleared
	CleanupMemory();

	if (!sM1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
	if (!sM1_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);

	switch (evalMethod) {

	case EVAL_EULER:
		break;

	case EVAL_TEULER:
	case EVAL_AHEUN:
		break;

	case EVAL_ABM:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval0_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RK4:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval0_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RK23:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval0_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RKF:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval3.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval4.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval0_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval3_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval4_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RKCK:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval3.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval4.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval0_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval3_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval4_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RKDP:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval3.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval4.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval5.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval0_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval3_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval4_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval5_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_SD:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval0_2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;
	}

	//For stochastic equations must also allocate memory for thermal VECs

	switch (setODE) {

	case ODE_SLLG:
	case ODE_SLLGSTT:
	case ODE_SLLGSA:
		if (!H_Thermal.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case ODE_SLLB:
	case ODE_SLLBSTT:
	case ODE_SLLBSA:
		if (!H_Thermal.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!Torque_Thermal.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
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

void DifferentialEquationAFM::CleanupMemory(void)
{
	//Only clear vectors not used for current evaluation method
	sM1.clear();
	sM1_2.clear();

	if (evalMethod != EVAL_RK4 &&
		evalMethod != EVAL_ABM &&
		evalMethod != EVAL_RK23 &&
		evalMethod != EVAL_RKF &&
		evalMethod != EVAL_RKCK &&
		evalMethod != EVAL_RKDP &&
		evalMethod != EVAL_SD) {

		sEval0.clear();
		sEval0_2.clear();
	}

	if (evalMethod != EVAL_RK4 &&
		evalMethod != EVAL_ABM &&
		evalMethod != EVAL_RK23 &&
		evalMethod != EVAL_RKF &&
		evalMethod != EVAL_RKCK &&
		evalMethod != EVAL_RKDP) {

		sEval1.clear();
		sEval1_2.clear();
	}

	if (evalMethod != EVAL_RK4 &&
		evalMethod != EVAL_RK23 &&
		evalMethod != EVAL_RKF &&
		evalMethod != EVAL_RKCK &&
		evalMethod != EVAL_RKDP) {

		sEval2.clear();
		sEval2_2.clear();
	}

	if (evalMethod != EVAL_RK4 &&
		evalMethod != EVAL_RKF &&
		evalMethod != EVAL_RKCK &&
		evalMethod != EVAL_RKDP) {

		sEval3.clear();
		sEval4.clear();
		sEval3_2.clear();
		sEval4_2.clear();
	}

	if (evalMethod != EVAL_RKDP) {

		sEval5.clear();
		sEval5_2.clear();
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


BError DifferentialEquationAFM::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(DifferentialEquationAFM));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_ODE_SOLVER)) {

		error = AllocateMemory();
	}

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_ODE_MOVEMESH)) {

		if (!error) {

			//set skip cells flags for moving mesh if enabled
			if (moving_mesh) {

				Rect mesh_rect = pMesh->GetMeshRect();

				DBL3 end_size = mesh_rect.size() & DBL3(MOVEMESH_ENDRATIO, 1, 1);

				Rect end_rect_left = Rect(mesh_rect.s, mesh_rect.s + end_size);
				Rect end_rect_right = Rect(mesh_rect.e - end_size, mesh_rect.e);

				pMesh->M.set_skipcells(end_rect_left);
				pMesh->M.set_skipcells(end_rect_right);
				pMesh->M2.set_skipcells(end_rect_left);
				pMesh->M2.set_skipcells(end_rect_right);
			}
			else {

				pMesh->M.clear_skipcells();
				pMesh->M2.clear_skipcells();
			}
		}
	}

	//----------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pmeshODECUDA) {

		if (!error) error = pmeshODECUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

//switch CUDA state on/off
BError DifferentialEquationAFM::SwitchCUDAState(bool cudaState)
{
	BError error(CLASS_STR(DifferentialEquationAFM));

#if COMPILECUDA == 1

	//are we switching to cuda?
	if (cudaState) {

		///delete cuda module object and null (just in case)
		if (pmeshODECUDA) delete pmeshODECUDA;
		pmeshODECUDA = nullptr;

		pmeshODECUDA = new DifferentialEquationAFMCUDA(this);
		error = pmeshODECUDA->Error_On_Create();
		//setup the managed cuDiffEq object with pointers to all required data for differential equations calculations
		if (!error) error = (reinterpret_cast<DifferentialEquationAFMCUDA*>(pmeshODECUDA)->Get_ManagedDiffEqCUDA())()->set_pointers(reinterpret_cast<DifferentialEquationAFMCUDA*>(pmeshODECUDA));

		//based on setODE value, setup device method pointers in cuDiffEq held in each DifferentialEquationCUDA object
		if (!error) pmeshODECUDA->SetODEMethodPointers();
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
DBL3 DifferentialEquationAFM::dMdt(int idx)
{
	return (pMesh->M[idx] - sM1[idx]) / dT_last;
}