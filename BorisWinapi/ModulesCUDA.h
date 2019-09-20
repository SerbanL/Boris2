#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "ErrorHandler.h"

#include "BorisCUDALib.h"

class ModulesCUDA {

private:

protected:

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

	bool initialized = false;

	//energy value for this effective field term
	cu_obj<cuBReal> energy;

private:

	bool holder_module_destroyed = false;

protected:

	//-------------------------- Kernel Launchers

	void ZeroEnergy(void);

public:

	//-------------------------- Constructor and Destructor

	ModulesCUDA(void) { ZeroEnergy(); }

	virtual ~ModulesCUDA() {}

	//-------------------------- Error report / Management

	BError Error_On_Create(void) { return error_on_create; }

	//this is called from the Module base in the destruction process just before this object is deleted. Remember destruction is always in reverse order of construction, top down.
	//When a module is deleted first the derived class is destroyed, then the Module base is destroyyed. The Module base calls for ModulesCUDA pointer to be deleted (if allocated).
	//If ModulesCUDA implementations do something in the destructor using the Module implementation data, we must stop it otherwise illegal memory access will result.
	//Use holder_module_destroyed flag to check before doing anything.
	void Holder_Module_Destroyed(void) { holder_module_destroyed = true; }
	bool Holder_Module_Available(void) { return !holder_module_destroyed; }

	//-------------------------- Initialization / Uninitialization

	//This is called in order to pre-compute required parameters before the simulation can run.
	//Should be executed only if initialized flag is false - thus if it definitely needs to (re-)initialize set initialized = false before calling it
	//1. Call it from constructor.
	//2. Call it from MeshChange method.
	//3. Call it after run command entered in console (might be initialized already, initialized flag will be checked)
	//4. For some console commands a module will need to be re-initialized (e.g. an electrode changed for the Transport module). Deal with these special cases by setting a special flag in the Mesh module.
	//This flag will be checked in UpdateField in those special modules, calling Initialize.
	virtual BError Initialize(void) = 0;

	//To implement special behaviour implement this.
	virtual void Uninitialize(void) = 0;

	//-------------------------- UpdateConfiguration

	//Call this when the mesh size has changed or the body inside the mesh has changed, or CUDA flag has changed.
	virtual BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC) = 0;

	//-------------------------- UpdateField

	//Simulation run-time method used to do calculations. This will launch a kernel.
	virtual void UpdateField(void) = 0;

	//-------------------------- Getters

	bool IsInitialized(void) { return initialized; }

	cuBReal GetEnergy(void) { return energy.to_cpu(); }
};

#endif
