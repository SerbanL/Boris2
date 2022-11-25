#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STFIELD

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class MeshCUDA;
class STField;
class ManagedMeshCUDA;
class ManagedAtom_MeshCUDA;

class STFieldCUDA :
	public ModulesCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

	//pointer to cpu version of STField
	STField* pSTField;

	//use fixed polarization version or not? STp not zero => fixed polarization. Set at initialization time
	bool fixed_polarization = true;

	//cu arrays with pointers to other meshes in coupling with the mesh holding this module, top and bottom, (ferromagnetic)
	cu_arr<ManagedMeshCUDA> pMeshFM_Bot;
	cu_arr<ManagedMeshCUDA> pMeshFM_Top;

	//cu arrays with pointers to other meshes in coupling with the mesh holding this module, top and bottom, (atomistic)
	cu_arr<ManagedAtom_MeshCUDA> pMeshAtom_Bot;
	cu_arr<ManagedAtom_MeshCUDA> pMeshAtom_Top;

public:

	STFieldCUDA(MeshCUDA* pMeshCUDA_, STField* pSTField_);
	~STFieldCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);
};

#else

class STFieldCUDA
{
};

#endif

#endif
