#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_STRAYFIELD

#include "ModulesCUDA.h"

class SuperMesh;
class DipoleMeshCUDA;

class StrayField;

class StrayFieldCUDA :
	public ModulesCUDA
{

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	SuperMesh* pSMesh;

	//the SDemag module (cpu version) holding this CUDA version
	StrayField* pStrayField;

	//strayField cuVEC taking values on the supermesh - transfer to and from ferromagnetic meshes
	cu_obj<cuVEC<cuReal3>> strayField;
	size_t strayField_size;

	//all currently set dipole meshes - get pointers to cuda M so we can pass them to a cuda kernel 
	cu_arr<cuVEC_VC<cuReal3>> Mdipoles;

private:

	//calculate stray field from all dipoles
	void CalculateStrayField(void);

public:

	StrayFieldCUDA(SuperMesh* pSMesh_, StrayField* pStrayField_);
	~StrayFieldCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	void UpdateField(void);

	//-------------------Getters

	cu_obj<cuVEC<cuReal3>>& GetStrayField(void) { return strayField; }
};

#else

class StrayFieldCUDA
{
};

#endif

#endif




