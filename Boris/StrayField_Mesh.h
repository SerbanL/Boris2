#pragma once

#include "BorisLib.h"
#include "Modules.h"

#if COMPILECUDA == 1
#include "StrayField_MeshCUDA.h"
#endif

class Mesh;

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "StrayField_Base.h"

class StrayField_Mesh :
	public StrayField_Base,
	public Modules,
	public ProgramState<StrayField_Mesh, std::tuple<>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend StrayField_MeshCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Mesh *pMesh;

public:

	StrayField_Mesh(Mesh *pMesh_);
	~StrayField_Mesh();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Energy methods

	double Get_EnergyChange(int spin_index, DBL3 Mnew);

	//AFM mesh
	DBL2 Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B);

	//-------------------Getters

	VEC<DBL3>& GetStrayField(void) { return strayField; }

#if COMPILECUDA == 1
	cu_obj<cuVEC<cuReal3>>& GetStrayFieldCUDA(void) { return dynamic_cast<StrayField_MeshCUDA*>(pModuleCUDA)->GetStrayField(); }
#endif
};

#else

class StrayField_Mesh :
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

	StrayField(Mesh *pMesh_) {}
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
