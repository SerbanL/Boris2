#pragma once

#include "BorisLib.h"
#include "Modules.h"

class MeshBase;
class Mesh;
class Atom_Mesh;
class SDemag;

#ifdef MODULE_COMPILATION_SDEMAG

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
	public Convolution<SDemag_Demag, DemagKernelCollection>,
	public ProgramState<SDemag_Demag, std::tuple<SZ3, DBL3, Rect, int>, std::tuple<>>
{

	friend SDemag;

#if COMPILECUDA == 1
	friend SDemagCUDA;
	friend SDemagCUDA_Demag;
#endif

private:

	MeshBase *pMeshBase = nullptr;

	//pointer to mesh object holding this effective field module (either micromagnetic or atomistic - only one will be not nullptr, so check)
	Mesh* pMesh = nullptr;
	Atom_Mesh* paMesh = nullptr;

	//pointer to managing SDemag module
	SDemag *pSDemag = nullptr;

	//save ferromagnetic mesh values so we know if we need to uninitialize
	SZ3 n;
	DBL3 h;
	Rect meshRect;

	//when used with 2d layering mode this will identify the layer number along z direction. Must have layer_number_2d < n.z.
	//if layer_number_2d = -1 then 2d layering is disabled.
	//can use this layer number to set convolution data, as well as detect if this module needs to be deleted (layer_number_2d >= n.z)
	//If n.z changes then we'll either need to delete some layers (SDemag_Demag modules) - done here - or add some more SDemag_Demag modules - done in SDemag
	int layer_number_2d = -1;

	//transfer values from M of this mesh to a VEC with fixed number of cells -> use same meshing for all layers.
	VEC<DBL3> transfer;

	//if displaying module Heff or energy, and mesh transfer is required, then use these (capture output field and energy, then do transfer)
	VEC<DBL3> transfer_Module_Heff;
	VEC<double> transfer_Module_energy;

	//do transfer as M -> transfer -> convolution -> transfer -> Heff if true
	//if false then don't use the transfer VEC but can do M -> convolution -> Heff directly
	//this flag will be set false only if the convolution rect matches that of M and n_common matches the discretisation of M (i.e. in this case a mesh transfer would be pointless).
	bool do_transfer = true;

	//number of non-empty cells in transfer
	int non_empty_cells;

	//different meshes have different weights when contributing to the total energy density -> ratio of their non-empty volume to total non-empty volume
	double energy_density_weight = 0.0;

	//Evaluation speedup mode data

	//vec for demagnetizing field polynomial extrapolation
	VEC<DBL3> Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, Hdemag6;

	//-Nxx, -Nyy, -Nzz values at r = r0
	DBL3 selfDemagCoeff = DBL3();

private:

	//initialize transfer object
	BError Initialize_Mesh_Transfer(void);

	//allocate memory and initialize mesh transfer for module Heff and energy display data
	BError Initialize_Module_Display(void);

public:

	SDemag_Demag(MeshBase *pMeshBase_);
	~SDemag_Demag() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) { Uninitialize(); }

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Energy methods

	//FM mesh
	double Get_EnergyChange(int spin_index, DBL3 Mnew);

	//AFM mesh
	DBL2 Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B);

	//-------------------Setters

	void Set_SDemag_Pointer(SDemag *pSDemag_);

	//Force the SDemag_Demag modules to calculate the layer_number_2d value, used for 2D layered convolution
	//this value is obtained from the minor id of this module held in *pMesh : the minor ids are guaranteed to be sequential and start from zero
	//thus if we've added enough of these SDemag_Demag modules, all layers will be covered
	void Set_2D_Layering(void);
};

#else

class SDemag_Demag :
	public Modules
{

private:

public:

	SDemag_Demag(MeshBase *pMeshBase_) {}
	~SDemag_Demag() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Setters

	void Set_SDemag_Pointer(SDemag *pSDemag_) {}
};

#endif