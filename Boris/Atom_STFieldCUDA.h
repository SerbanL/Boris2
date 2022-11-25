#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_STFIELD) && ATOMISTIC == 1

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class Atom_STField;

class Atom_MeshCUDA;
class ManagedAtom_MeshCUDA;
class ManagedMeshCUDA;

class Atom_STFieldCUDA :
	public ModulesCUDA
{

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	Atom_MeshCUDA* paMeshCUDA;

	//pointer to cpu version of STField
	Atom_STField* pSTField;

	//use fixed polarization version or not? STp not zero => fixed polarization. Set at initialization time
	bool fixed_polarization = true;

	//cu arrays with pointers to other meshes in coupling with the mesh holding this module, top and bottom, (atomistic)
	cu_arr<ManagedAtom_MeshCUDA> paMesh_Bot;
	cu_arr<ManagedAtom_MeshCUDA> paMesh_Top;

	//cu arrays with pointers to other meshes in coupling with the mesh holding this module, top and bottom, (ferromagnetic)
	cu_arr<ManagedMeshCUDA> pMeshFM_Bot;
	cu_arr<ManagedMeshCUDA> pMeshFM_Top;

public:

	Atom_STFieldCUDA(Atom_MeshCUDA* paMeshCUDA_, Atom_STField* pSTField_);
	~Atom_STFieldCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);
};

#else

class Atom_STFieldCUDA
{
};

#endif

#endif
