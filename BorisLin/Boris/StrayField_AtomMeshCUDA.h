#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "ModulesCUDA.h"

#include "StrayField_BaseCUDA.h"

class Atom_Mesh;
class Atom_MeshCUDA;
class StrayField_AtomMesh;

class StrayField_AtomMeshCUDA :
	public StrayField_BaseCUDA,
	public ModulesCUDA
{

private:

	//pointer to mesh object holding this effective field module
	Atom_Mesh *paMesh;

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	Atom_MeshCUDA* paMeshCUDA;

	//the cpu version module holding this CUDA version
	StrayField_AtomMesh* pStrayField;

private:

	void UpdateFieldCUDA(void);

	void set_StrayField_AtomMeshCUDA_pointers(void);

public:

	StrayField_AtomMeshCUDA(Atom_Mesh* paMesh_, StrayField_AtomMesh* pStrayField_);
	~StrayField_AtomMeshCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);

	//-------------------Getters

	cu_obj<cuVEC<cuReal3>>& GetStrayField(void) { return strayField; }
};

#else

class StrayField_AtomMeshCUDA
{
};

#endif

#endif




