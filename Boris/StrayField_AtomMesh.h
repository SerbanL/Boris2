#pragma once

#include "BorisLib.h"
#include "Modules.h"

#if COMPILECUDA == 1
#include "StrayField_AtomMeshCUDA.h"
#endif

class Atom_Mesh;

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "StrayField_Base.h"

class StrayField_AtomMesh :
	public StrayField_Base,
	public Modules,
	public ProgramState<StrayField_AtomMesh, std::tuple<>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend StrayField_AtomMeshCUDA;
#endif

private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh *paMesh;

	//divide energy by this to obtain energy density : this is the energy density in the entire mesh, which may not be rectangular.
	double non_empty_volume = 0.0;

public:

	StrayField_AtomMesh(Atom_Mesh *paMesh_);
	~StrayField_AtomMesh();

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

	//For simple cubic mesh spin_index coincides with index in M1
	double Get_EnergyChange(int spin_index, DBL3 Mnew);

	//-------------------Getters

	VEC<DBL3>& GetStrayField(void) { return strayField; }

#if COMPILECUDA == 1
	cu_obj<cuVEC<cuReal3>>& GetStrayFieldCUDA(void) { return dynamic_cast<StrayField_AtomMeshCUDA*>(pModuleCUDA)->GetStrayField(); }
#endif
};

#else

class StrayField_AtomMesh :
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

	StrayField_AtomMesh(Atom_Mesh *paMesh_) {}
	~StrayField_AtomMesh() {}

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