#include "stdafx.h"
#include "DiffEq.h"
#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"

///////////////////////////////////////////////////////////////////////////////

DifferentialEquation::DifferentialEquation(FMesh *pMesh) :
	ODECommon(true),
	prng(GetTickCount())
{
	//store this pointer in the pODE vector and get a unique id so it can be erased later in the destructor
	odeId.minor = pODE.push_back(this, odeId.major);

	this->pMesh = pMesh;

	error_on_create = UpdateConfiguration();
}

DifferentialEquation::~DifferentialEquation() 
{
	//erase the pODE entry
	if(pODE.is_id_set(odeId))
		pODE.erase(odeId);

	//if this mesh (ferromagnetic mesh) had the moving mesh trigger set, then it is now gone so must set the moving_mesh trigger to false in the ODE common too.
	if (pMesh->GetMoveMeshTrigger()) moving_mesh = false;

#if COMPILECUDA == 1
	if (pmeshODECUDA) delete pmeshODECUDA;
	pmeshODECUDA = nullptr;
#endif
}

BError DifferentialEquation::AllocateMemory(void)
{
	BError error(CLASS_STR(DifferentialEquation));
	
	//first make sure everything not needed is cleared
	CleanupMemory();
	
	if (!sM1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);

	switch(evalMethod) {

	case EVAL_EULER:
		break;

	case EVAL_TEULER:
		break;

	case EVAL_RK4:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_ABM:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;

	case EVAL_RKF:
		if (!sEval0.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval2.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval3.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		if (!sEval4.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
		break;
	}

	//For stochastic equations must also allocate memory for thermal VECs

	switch (setODE) {

	case ODE_SLLG:
	case ODE_SLLGSTT:
	case ODE_SLLGSA:
		if(!H_Thermal.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);
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

void DifferentialEquation::CleanupMemory(void) 
{
	//Only clear vectors not used for current evaluation method
	sM1.clear();
	
	if (evalMethod != EVAL_RK4 && evalMethod != EVAL_ABM && evalMethod != EVAL_RKF) {

		sEval0.clear();
	}

	if (evalMethod != EVAL_RK4 && evalMethod != EVAL_ABM && evalMethod != EVAL_RKF) {

		sEval1.clear();
	}

	if (evalMethod != EVAL_RK4 && evalMethod != EVAL_RKF) {

		sEval2.clear();
		sEval3.clear();
		sEval4.clear();
	}

	//For thermal vecs only clear if not used for current set ODE
	if (setODE != ODE_SLLG && setODE != ODE_SLLGSTT && setODE != ODE_SLLB && setODE != ODE_SLLBSTT && setODE != ODE_SLLGSA && setODE != ODE_SLLBSA) {

		H_Thermal.clear();
	}

	if (setODE != ODE_SLLB && setODE != ODE_SLLBSTT && setODE != ODE_SLLBSA) {

		Torque_Thermal.clear();
	}

	//----------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pmeshODECUDA) { pmeshODECUDA->CleanupMemory(); }
#endif
}


BError DifferentialEquation::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(DifferentialEquation));
	
	error = AllocateMemory();
	
	if (!error) {

		//set skip cells flags for moving mesh if enabled
		if (moving_mesh) {

			Rect mesh_rect = pMesh->GetMeshRect();

			DBL3 end_size = mesh_rect.size() & DBL3(MOVEMESH_ENDRATIO, 1, 1);

			Rect end_rect_left = Rect(mesh_rect.s, mesh_rect.s + end_size);
			Rect end_rect_right = Rect(mesh_rect.e - end_size, mesh_rect.e);

			pMesh->M.set_skipcells(end_rect_left);
			pMesh->M.set_skipcells(end_rect_right);
		}
		else {

			pMesh->M.clear_skipcells();
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
BError DifferentialEquation::SwitchCUDAState(bool cudaState)
{
	BError error(CLASS_STR(DifferentialEquation));

#if COMPILECUDA == 1

	//are we switching to cuda?
	if (cudaState) {

		///delete cuda module object and null (just in case)
		if (pmeshODECUDA) delete pmeshODECUDA;
		pmeshODECUDA = nullptr;

		pmeshODECUDA = new DifferentialEquationCUDA(this);
		error = pmeshODECUDA->Error_On_Create();
		//setup the managed cuDiffEq object with pointers to all required data for differential equations calculations
		if (!error) error = pmeshODECUDA->cuDiffEq()->set_pointers(pmeshODECUDA);
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
DBL3 DifferentialEquation::dMdt(int idx) 
{ 
	return (pMesh->M[idx] - sM1[idx]) / dT_last; 
}