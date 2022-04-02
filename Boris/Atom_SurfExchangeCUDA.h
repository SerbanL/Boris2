#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_SURFEXCHANGE) && ATOMISTIC == 1

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

class Atom_SurfExchange;

class Atom_MeshCUDA;
class ManagedAtom_MeshCUDA;

class ManagedMeshCUDA;

class Atom_SurfExchangeCUDA :
	public ModulesCUDA
{

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	Atom_MeshCUDA* paMeshCUDA;

	//pointer to cpu version of SurfExchange
	Atom_SurfExchange* paSurfExch;

	//cu arrays with pointers to other atomistic meshes in surface exchange coupling with the mesh holding this module, top and bottom
	cu_arr<ManagedAtom_MeshCUDA> paMesh_Bot;
	cu_arr<ManagedAtom_MeshCUDA> paMesh_Top;

	//cu arrays with pointers to micromagnetic ferromagnetic meshes in surface exchange coupling with the mesh holding this module, top and bottom
	cu_arr<ManagedMeshCUDA> pMeshFM_Bot;
	cu_arr<ManagedMeshCUDA> pMeshFM_Top;

	//cu arrays with pointers to micromagnetic antiferromagnetic meshes in surface exchange coupling with the mesh holding this module, top and bottom
	cu_arr<ManagedMeshCUDA> pMeshAFM_Bot;
	cu_arr<ManagedMeshCUDA> pMeshAFM_Top;

private:

	//Set pointers in ManagedAtom_MeshCUDA so we can access them in device code. This is used by MonteCarlo algorithm.
	void set_Atom_SurfExchangeCUDA_pointers(void);

public:

	Atom_SurfExchangeCUDA(Atom_MeshCUDA* paMeshCUDA_, Atom_SurfExchange* paSurfExch_);
	~Atom_SurfExchangeCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);

	//-------------------Torque methods

	cuReal3 GetTorque(cuRect avRect);

};

#else

class Atom_SurfExchangeCUDA
{
};

#endif

#endif