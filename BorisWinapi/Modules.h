#pragma once

#include <omp.h>

#include "ErrorHandler.h"

#include "BorisLib.h"
#include "Boris_Enums_Defs.h"

#if COMPILECUDA == 1
#include "ModulesCUDA.h"
#endif

using namespace std;

//Modules (MODS_ entries are super-mesh versions and not available as normal module handles
//Add new entries at the end to keep older simulation files compatible
//If you need to take something out you'll need to figure a way to keep older simulation files compatible
enum MOD_ {
	MOD_ALL = -1,
	MOD_ERROR = 0,
	MOD_DEMAG_N, MOD_DEMAG, MODS_SDEMAG,
	MOD_EXCHANGE6NGBR, MOD_DMEXCHANGE, MOD_IDMEXCHANGE, MOD_SURFEXCHANGE,
	MOD_ZEEMAN,
	MOD_ANIUNI, MOD_ANICUBI,
	MOD_TRANSPORT, MODS_STRANSPORT, MODS_OERSTED,
	MODS_STRAYFIELD,
	MOD_HEAT, MODS_SHEAT,
	MOD_SOTFIELD,
	MOD_ROUGHNESS,
	MOD_SDEMAG_DEMAG
};

class Mesh;

class Modules {

private:

protected:

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

	bool initialized = false;

	//energy value for this effective field term
	double energy = 0.0;

	//The CUDA version of this module (ModulesCUDA is the interface and when CUDA is switched on an implementation is created depending on module type)
#if COMPILECUDA == 1
	ModulesCUDA* pModuleCUDA = nullptr;
#endif

private: //Private methods

public:

	//-------------------------- Constructor and Destructor

	Modules(void) {}
	virtual ~Modules()
	{
#if COMPILECUDA == 1
		if (pModuleCUDA) {

			//at this point the Modules implementation is destroyed, so when the pModuleCUDA destructor is called it cannot do anything it might normally do to Modules implementation data. Use flag to check.
			pModuleCUDA->Holder_Module_Destroyed();

			delete pModuleCUDA;
			pModuleCUDA = nullptr;
		}
#endif
	}
	
	//-------------------------- Error report

	BError Error_On_Create(void) { return error_on_create; }

	//-------------------------- Initialization / Uninitialization

	//This is called in order to pre-compute required parameters before the simulation can run.
	//Should be executed only if initialized flag is false - thus if it definitely needs to (re-)initialize set initialized = false before calling it
	//1. Call it from constructor.
	//2. Call it from MeshChange method.
	//3. Call it after run command entered in console (might be initialized already, initialized flag will be checked)
	//4. For some console commands a module will need to be re-initialized (e.g. an electrode changed for the Transport module). Deal with these special cases by setting a special flag in the Mesh module.
	//This flag will be checked in UpdateField in those special modules, calling Initialize.
	virtual BError Initialize(void) = 0;

#if COMPILECUDA == 1
	//Only call this if cuda is switched on : if cuda on call InitializeCUDA chain, if cuda off call Initialize chain.
	BError InitializeCUDA(void) 
	{ 
		if (pModuleCUDA) return pModuleCUDA->Initialize(); 
		else return BError();
	}
#endif

	//To implement special behaviour implement this.
	virtual void Uninitialize(void) = 0;

#if COMPILECUDA == 1
	//Only call this if cuda is switched on
	void UninitializeCUDA(void) { if (pModuleCUDA) pModuleCUDA->Uninitialize(); }
#endif

	//-------------------------- UpdateConfiguration

	//Call this when the mesh size has changed or the body inside the mesh has changed, or CUDA flag has changed.
	virtual BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC) = 0;

	//-------------------------- CUDA Switch

	//switch CUDA state on/off -> this is the logical outline of the switch
	BError SwitchCUDAState(bool cudaState);

	//The switch method calls this to make the implementation of ModulesCUDA in *pModuleCUDA
	virtual BError MakeCUDAModule(void) = 0;

	//-------------------------- UpdateField

	//Simulation run-time method used to do calculations.
	//return total volume energy density -> each module will have a contribution, so sum it
	virtual double UpdateField(void) = 0;

#if COMPILECUDA == 1
	//Only call this if cuda is switched on : if cuda on call UpdateFieldCUDA chain, if cuda off call UpdateField chain.
	void UpdateFieldCUDA(void) { if (pModuleCUDA) pModuleCUDA->UpdateField(); }
#endif

	//-------------------------- Getters

#if COMPILECUDA == 1
	ModulesCUDA* GetCUDAModule(void) { return pModuleCUDA; }
#endif

	bool IsInitialized(void) 
	{
#if COMPILECUDA == 1
		if (pModuleCUDA) return pModuleCUDA->IsInitialized();
#endif
		return initialized;
	}

	double GetEnergy(void) 
	{ 
#if COMPILECUDA == 1
		if (pModuleCUDA) return pModuleCUDA->GetEnergy();
#endif
		return energy;
	}
};