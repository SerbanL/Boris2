#include "stdafx.h"
#include "DiffEqDM.h"

#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "Mesh_Diamagnetic.h"

DifferentialEquationDM::DifferentialEquationDM(DiaMesh *pMesh):
	DifferentialEquation(pMesh)
{
	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
}

DifferentialEquationDM::~DifferentialEquationDM()
{

}

//---------------------------------------- OTHERS

//Restore magnetisation after a failed step for adaptive time-step methods
void DifferentialEquationDM::RestoreMagnetisation(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n.dim(); idx++)
		pMesh->M[idx] = sM1[idx];
}

//---------------------------------------- SET-UP METHODS

BError DifferentialEquationDM::AllocateMemory(void)
{
	BError error(CLASS_STR(DifferentialEquationDM));

	//first make sure everything not needed is cleared
	CleanupMemory();

	if (!sM1.resize(pMesh->n)) return error(BERROR_OUTOFMEMORY_CRIT);

	//----------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pmeshODECUDA) {

		if (!error) error = pmeshODECUDA->AllocateMemory();
	}
#endif

	return error;
}

void DifferentialEquationDM::CleanupMemory(void)
{
	//Only clear vectors not used for current evaluation method
	sM1.clear();

	//----------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pmeshODECUDA) { pmeshODECUDA->CleanupMemory(); }
#endif
}


BError DifferentialEquationDM::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(DifferentialEquationDM));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_ODE_SOLVER)) {

		error = AllocateMemory();
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
BError DifferentialEquationDM::SwitchCUDAState(bool cudaState)
{
	BError error(CLASS_STR(DifferentialEquationDM));

#if COMPILECUDA == 1

	//are we switching to cuda?
	if (cudaState) {

		if (!pmeshODECUDA) {

			pmeshODECUDA = new DifferentialEquationDMCUDA(this);
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
DBL3 DifferentialEquationDM::dMdt(int idx)
{
	return DBL3();
}

#endif
