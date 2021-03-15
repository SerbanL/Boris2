#pragma once

#include <omp.h>

#include "ErrorHandler.h"

#include "BorisLib.h"
#include "Boris_Enums_Defs.h"

#include "ModulesDefs.h"

#if COMPILECUDA == 1
#include "ModulesCUDA.h"
#endif

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

	//effective field in this module : if sized then module should be updating this if appropriate, else skip saving effective field
	//2nd VEC used for 2-sublattice modules
	VEC<DBL3> Module_Heff, Module_Heff2;

	//energy (density) spatial variation : if sized then module should be updating this is appropriate, else skip saving energy density
	//2nd VEC used for 2-sublattice modules
	VEC<double> Module_energy, Module_energy2;

protected:

	//zero the VECs if not empty
	void ZeroModuleVECs(void);

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

	//Call this when the object configuration has changed.
	virtual BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) = 0;

	//This is a "softer" version of UpdateConfiguration, which can be used any time and doesn't require the object to be Uninitialized; 
	//this will typically involve changing a value across multiple objects, thus better to call this method rather than try to remember which objects need the value changed.
	virtual void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) = 0;

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

	//-------------------------- Effective field and energy VECs

	//Make sure memory is allocated correctly for display data if used, else free memory
	BError Update_Module_Display_VECs(DBL3 h, Rect meshRect, bool Module_Heff_used, bool Module_Energy_used, bool twosublattice = false);

	//Get VECs for display
	VEC<DBL3>& Get_Module_Heff(void) { return Module_Heff; }
	VEC<DBL3>& Get_Module_Heff2(void) { return Module_Heff2; }

	//Get VECs for display
	VEC<double>& Get_Module_Energy(void) { return Module_energy; }
	VEC<double>& Get_Module_Energy2(void) { return Module_energy2; }

#if COMPILECUDA == 1
	//IMPORTANT : only use these if pModuleCUDA is not nullptr (so only before checking if CUDA switched on)
	//Get VECs for display
	cu_obj<cuVEC<cuReal3>>& Get_Module_HeffCUDA(void) { return pModuleCUDA->Get_Module_Heff(); }
	cu_obj<cuVEC<cuReal3>>& Get_Module_Heff2CUDA(void) { return pModuleCUDA->Get_Module_Heff2(); }

	//Get VECs for display
	cu_obj<cuVEC<cuBReal>>& Get_Module_EnergyCUDA(void) { return pModuleCUDA->Get_Module_Energy(); }
	cu_obj<cuVEC<cuBReal>>& Get_Module_Energy2CUDA(void) { return pModuleCUDA->Get_Module_Energy2(); }
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

	//-------------------------- Energies

	//Get energy density averaged over the entire mesh during the UpdateField call
	double GetEnergyDensity(void);

	//Calculate the energy density in the given rect only
	double GetEnergyDensity(Rect& avRect);

	//return atomistic energy change from current spin to new spin (for a simple cubic mesh current spin coincides with the index in M1, but for other crystal structures it does not - individual mesh modules will know what to do).
	//Implement as appropriate.
	//Return new - old energy to be used with monte Carlo methods
	virtual double Get_Atomistic_EnergyChange(int spin_index, DBL3 Mnew) { return 0.0; }
};