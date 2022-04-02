#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "ModulesCUDA.h"

#include "StrayField_BaseCUDA.h"

class Mesh;
class MeshCUDA;
class StrayField_Mesh;

class StrayField_MeshCUDA :
	public StrayField_BaseCUDA,
	public ModulesCUDA
{

private:

	//pointer to mesh object holding this effective field module
	Mesh *pMesh;

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

	//the cpu version module holding this CUDA version
	StrayField_Mesh* pStrayField;

private:

	void UpdateFieldCUDA(void);

	void set_StrayField_MeshCUDA_pointers(void);

public:

	StrayField_MeshCUDA(Mesh* pMesh_, StrayField_Mesh* pStrayField_);
	~StrayField_MeshCUDA();

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

class StrayField_MeshCUDA
{
};

#endif

#endif




