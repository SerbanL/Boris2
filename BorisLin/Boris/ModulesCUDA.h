#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "ErrorHandler.h"

#include "BorisCUDALib.h"

#include "ManagedModulesCUDA.h"

#include "ModulesDefs.h"

class ModulesCUDA {

	friend ManagedModulesCUDA;

private:

	bool holder_module_destroyed = false;

protected:

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

	bool initialized = false;

	//energy value for this effective field term
	cu_obj<cuBReal> energy;

	//used to obtain average torque value
	cu_obj<cuReal3> torque;

	//auxiliary for obtaining average energy in a custom rectangle : count non-zero points for average.
	cu_obj<size_t> points_count;

	//effective field in this module : if sized then module should be updating this if appropriate, else skip saving effective field
	//2nd VEC used for 2-sublattice modules
	cu_obj<cuVEC<cuReal3>> Module_Heff, Module_Heff2;
	size_t Module_Heff_size = 0, Module_Heff2_size = 0;

	//energy (density) spatial variation : if sized then module should be updating this is appropriate, else skip saving energy density
	//2nd VEC used for 2-sublattice modules
	cu_obj<cuVEC<cuBReal>> Module_energy, Module_energy2;
	size_t Module_energy_size = 0, Module_energy2_size = 0;

	//Managed Module : just pass this into cuda kernels so you can access the above data through cuModule
	cu_obj<ManagedModulesCUDA> cuModule;

protected:

	//-------------------------- Kernel Launchers

	//zero the auxiliary values (energy and points_count)
	void ZeroEnergy(void);

	//zero the cuVECs if not empty
	void ZeroModuleVECs(void);

	//return cross product of M with Module_Heff, averaged in given rect (relative)
	cuReal3 CalculateTorque(cu_obj<cuVEC_VC<cuReal3>>& M, cuRect& avRect);

public:

	//-------------------------- Constructor and Destructor

	ModulesCUDA(void) 
	{ 
		ZeroEnergy(); 
		//setup ManagedModulesCUDA object
		cuModule()->set_pointers(this);
	}

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

	//Call this when the object configuration has changed
	virtual BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) = 0;

	//This is a "softer" version of UpdateConfiguration, which can be used any time and doesn't require the object to be Uninitialized; 
	//this will typically involve changing a value across multiple objects, thus better to call this method rather than try to remember which objects need the value changed.
	virtual void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) = 0;

	//-------------------------- UpdateField

	//Simulation run-time method used to do calculations. This will launch a kernel.
	virtual void UpdateField(void) = 0;

	//-------------------------- Effective field and energy VECs

	//Make sure memory is allocated correctly for display data if used, else free memory
	BError Update_Module_Display_VECs(cuReal3 h, cuRect meshRect, bool Module_Heff_used, bool Module_Energy_used, bool twosublattice = false);

	//Get VECs for display
	cu_obj<cuVEC<cuReal3>>& Get_Module_Heff(void) { return Module_Heff; }
	cu_obj<cuVEC<cuReal3>>& Get_Module_Heff2(void) { return Module_Heff2; }

	//Get VECs for display
	cu_obj<cuVEC<cuBReal>>& Get_Module_Energy(void) { return Module_energy; }
	cu_obj<cuVEC<cuBReal>>& Get_Module_Energy2(void) { return Module_energy2; }

	//-------------------------- Getters

	bool IsInitialized(void) { return initialized; }

	//-------------------------- Energies and Torques

	cuBReal GetEnergyDensity(void) { return energy.to_cpu(); }
	
	//Calculate the energy density in the given rect only
	cuBReal GetEnergyDensity(cuRect avRect);

	//implement as needed : this is normally the torque obtained as M cross Module_Heff
	virtual cuReal3 GetTorque(cuRect avRect) { return cuReal3(); }
};

#endif
