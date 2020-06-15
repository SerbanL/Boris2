#pragma once

#include "BorisLib.h"
#include "Modules.h"

#if COMPILECUDA == 1
#include "StrayFieldCUDA.h"
#endif

using namespace std;

class SuperMesh;
class DipoleMesh;

#ifdef MODULE_STRAYFIELD

class StrayField :
	public Modules,
	public ProgramState<StrayField, tuple<>, tuple<>>
{

#if COMPILECUDA == 1
	friend StrayFieldCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module (supermesh)
	SuperMesh *pSMesh;

	//strayField VEC taking values on the supermesh - transfer to and from ferromagnetic meshes
	VEC<DBL3> strayField;

	//all currently set dipole meshes
	vector<DipoleMesh*> dipoles;

private:

	//calculate stray field from all dipoles
	void CalculateStrayField(void);

	//check if stray field needs to be recalculated
	bool CheckRecalculateStrayField(void);

public:

	StrayField(SuperMesh *pSMesh_);
	~StrayField() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Getters

	VEC<DBL3>& GetStrayField(void) { return strayField; }

#if COMPILECUDA == 1
	cu_obj<cuVEC<cuReal3>>& GetStrayFieldCUDA(void) { return dynamic_cast<StrayFieldCUDA*>(pModuleCUDA)->GetStrayField(); }
#endif
};

#else

class StrayField :
	public Modules
{

private:

	//strayField VEC taking values on the supermesh - transfer to and from ferromagnetic meshes
	VEC<DBL3> strayField;

#if COMPILECUDA == 1
	cu_obj<cuVEC<cuReal3>> strayField_CUDA;
#endif

private:

public:

	StrayField(SuperMesh *pSMesh_) {}
	~StrayField() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Getters

	VEC<DBL3>& GetStrayField(void) { return strayField; }

#if COMPILECUDA == 1
	cu_obj<cuVEC<cuReal3>>& GetStrayFieldCUDA(void) { return strayField_CUDA; }
#endif
};

#endif
