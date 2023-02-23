#pragma once

#include "BorisLib.h"
#include "Modules.h"

#if COMPILECUDA == 1
#include "StrayFieldCUDA.h"
#endif

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "StrayField_Base.h"

class StrayField :
	public StrayField_Base,
	public Modules,
	public ProgramState<StrayField, std::tuple<>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend StrayFieldCUDA;
#endif

private:

public:

	StrayField(SuperMesh *pSMesh_);
	~StrayField();

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

class SuperMesh;

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
