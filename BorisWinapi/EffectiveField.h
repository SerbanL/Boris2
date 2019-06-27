#pragma once

#include <omp.h>

#include "Types.h"
#include "Boris_Enums_Defs.h"

using namespace std;

class EffectiveField {

protected: //Protected data

	bool initialized;

	//energy value for this effective field term
	double energy;

public:  //Public data

private: //Private methods

public:

	EffectiveField(void) { initialized = false; energy = 0; }
	virtual ~EffectiveField() { }
	
	//This is called in order to pre-compute required parameters before the simulation can run.
	//Should be executed only if initialized flag is false - thus if it definitely needs to (re-)initialize set initialized = false before calling it
	//1. Call it from constructor.
	//2. Call it from MeshChange method.
	//3. Call it after run command entered in console (might be initialized already, initialized flag will be checked)
	//4. For some console commands a module will need to be re-initialized (e.g. an electrode changed for the Transport module). Deal with these special cases by setting a special flag in the Mesh module.
	//This flag will be checked in UpdateField in those special modules, calling Initialize.
	virtual bool Initialize(void) = 0;

	//To implement special behaviour implement this.
	virtual void Uninitialize(void) = 0;

	//Call this when the mesh size has changed or the body inside the mesh has changed, or CUDA flag has changed.
	virtual void MeshSizeChange(void) = 0;

	//Simulation run-time method used to do calculations.
	virtual void UpdateField(void) = 0;

	bool IsInitialized(void) {return initialized;}
};