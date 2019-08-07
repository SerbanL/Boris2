#pragma once

#include "BorisLib.h"
#include "Modules.h"

using namespace std;

class Mesh;
class FMesh;
class SDemag;

#ifdef MODULE_SDEMAG

#include "Convolution.h"
#include "DemagKernelCollection.h"

#if COMPILECUDA == 1
#include "SDemagCUDA_Demag.h"
#endif

#if COMPILECUDA == 1
class SDemagCUDA;
#endif

////////////////////////////////////////////////////////////////////////////////////////////////
//
// Demag module used by super-mesh demag module, used to calculate contributions to super-mesh demag calculations from individual meshes

class SDemag_Demag :
	public Modules,
	public Convolution<DemagKernelCollection>,
	public ProgramState<SDemag_Demag, tuple<SZ3, DBL3, Rect>, tuple<>>
{

	friend SDemag;

#if COMPILECUDA == 1
	friend SDemagCUDA;
	friend SDemagCUDA_Demag;
#endif

private:

	//pointer to mesh object holding this effective field module
	FMesh* pMesh;

	//pointer to managing SDemag module
	SDemag *pSDemag = nullptr;

	//save ferromagnetic mesh values so we know if we need to uninitialize
	SZ3 n;
	DBL3 h;
	Rect meshRect;

	//transfer values from M of this mesh to a VEC with fixed number of cells -> use same meshing for all layers.
	VEC<DBL3> transfer;

	//do transfer as M -> transfer -> convolution -> transfer -> Heff if true
	//if false then don't use the transfer VEC but can do M -> convolution -> Heff directly
	//this flag will be set false only if the convolution rect matches that of M and n_common matches the discretisation of M (i.e. in this case a mesh transfer would be pointless).
	bool do_transfer = true;

	//number of non-empty cells in transfer
	int non_empty_cells;

	//different meshes have different weights when contributing to the total energy density -> ratio of their non-empty volume to total non-empty volume
	double energy_density_weight = 0.0;

private:

	//initialize transfer object
	BError Initialize_Mesh_Transfer(void);

public:

	SDemag_Demag(Mesh *pMesh_);
	~SDemag_Demag() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) { Uninitialize(); }

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	BError MakeCUDAModule(void);

	double UpdateField(void) { return 0.0; }

	//-------------------Setters

	void Set_SDemag_Pointer(SDemag *pSDemag_);
};

#else

class SDemag_Demag :
	public Modules
{

private:

public:

	SDemag_Demag(Mesh *pMesh_) {}
	~SDemag_Demag() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC) { return BError(); }

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Setters

	void Set_SDemag_Pointer(SDemag *pSDemag_) {}
};

#endif