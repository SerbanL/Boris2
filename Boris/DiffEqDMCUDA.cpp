#include "stdafx.h"
#include "DiffEqDMCUDA.h"

#if COMPILECUDA == 1
#ifdef MESH_COMPILATION_DIAMAGNETIC

#include "DiffEqDM.h"

#include "Mesh_Diamagnetic.h"
#include "Mesh_DiamagneticCUDA.h"

DifferentialEquationDMCUDA::DifferentialEquationDMCUDA(DifferentialEquation *pmeshODE) :
	DifferentialEquationCUDA(pmeshODE)
{
	error_on_create = AllocateMemory(true);

	cuDiffEq()->set_pointers(this);
	SetODEMethodPointers();
}

DifferentialEquationDMCUDA::~DifferentialEquationDMCUDA()
{
	//If called_from_destructor is true then do not attempt to transfer data to cpu where this is held in the derived class of DifferentialEquation
	//This derived class will have already destructed so attempting to copy over data to it is an invalid action and can crash the program.
	//This can happen when a mesh is deleted with CUDA switched on
	//It's also possible this destructor was called simply due to switching CUDA off, in which case called_from_destructor should be false
	CleanupMemory(!pmeshODE->called_from_destructor);
}

//---------------------------------------- SET-UP METHODS  : DiffEqCUDA.cpp and DiffEqCUDA.cu

//allocate memory depending on set evaluation method - also cleans up previously allocated memory by calling CleanupMemory();
BError DifferentialEquationDMCUDA::AllocateMemory(bool copy_from_cpu)
{
	BError error(CLASS_STR(DifferentialEquationDMCUDA));

	//first make sure everything not needed is cleared
	CleanupMemory();

	if (!sM1()->resize((cuSZ3)pMesh->n)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	//copy values from cpu : it's possible the user switches to CUDA during a simulation and expects it to continue as normal
	//some evaluation methods need saved data (sM1, sEval... etc.) to carry over.
	else if (copy_from_cpu) sM1()->copy_from_cpuvec(pmeshODE->sM1);

	return error;
}

//deallocate memory before re-allocating it (depending on evaluation method previously allocated memory might not be used again, so need clean-up before)
void DifferentialEquationDMCUDA::CleanupMemory(bool copy_to_cpu)
{
	//copy values to cpu before erasing : it's possible the user switches CUDA off during a simulation and expects it to continue as normal
	//some evaluation methods need saved data (sM1, sEval... etc.) to carry over.
	//CleanupMemory will be called by destructor in this case, so before destroying gpu data copy it over to cpu if possible
	//CleanupMemory may also be called in other circumstances, in particular from the cpu version of CleanupMemory, after having cleaned cpu vecs, thus in this case the copy methods will not run.

	//Only clear vectors not used for current evaluation method
	if (copy_to_cpu && sM1()->size_cpu() == pmeshODE->sM1.size()) sM1()->copy_to_cpuvec(pmeshODE->sM1);
	sM1()->clear();
}


BError DifferentialEquationDMCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(DifferentialEquationDMCUDA));

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHDELETED)) {

		//if a mesh is deleted then a DifferentialEquationCUDA object can be deleted
		//this results in deletion of static data in ODECommonCUDA
		//whilst the static data is remade by UpdateConfiguration in ODECommonCUDA following this, our ManagedDiffEqDMCUDA object now has pointers which are not linked correctly, so need to update them
		cuDiffEq()->set_pointers(this);
	}

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_ODE_SOLVER)) {

		error = AllocateMemory();
	}

	return error;
}

#endif
#endif