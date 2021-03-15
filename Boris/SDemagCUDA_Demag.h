#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SDEMAG

#include "ModulesCUDA.h"
#include "BorisCUDALib.h"

#include "ConvolutionCUDA.h"
#include "DemagKernelCollectionCUDA.h"

class MeshCUDA;
class ManagedMeshCUDA;
class SDemag_Demag;
class SDemagCUDA;

class SDemagCUDA_Demag :
	public ModulesCUDA,
	public ConvolutionCUDA<SDemagCUDA_Demag, DemagKernelCollectionCUDA>
{
	friend SDemagCUDA;

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA * pMeshCUDA;

	//pointer to cpu version of this module
	SDemag_Demag *pSDemag_Demag;

	//transfer values from M of this mesh to a cuVEC with fixed number of cells -> use same meshing for all layers.
	cu_obj<cuVEC<cuReal3>> transfer;

	//if displaying module Heff or energy, and mesh transfer is required, then use these (capture output field and energy, then do transfer)
	cu_obj<cuVEC<cuReal3>> transfer_Module_Heff;
	cu_obj<cuVEC<cuBReal>> transfer_Module_energy;

	//do transfer as M -> transfer -> convolution -> transfer -> Heff if true
	//if false then don't use the transfer VEC but can do M -> convolution -> Heff directly
	//this flag will be set false only if the convolution rect matches that of M and n_common matches the discretisation of M (i.e. in this case a mesh transfer would be pointless).
	bool do_transfer = true;

	//different meshes have different weights when contributing to the total energy density -> ratio of their non-empty volume to total non-empty volume
	cu_obj<cuBReal> energy_density_weight;

	//Evaluation speedup mode data

	//1 Hdemag: no extrapolation, just save evaluation and reuse
	//2 Hdemag: linear extrapolation, need 2
	//3 Hdemag: quadratic extrapolation, need 3
	cu_obj<cuVEC<cuReal3>> Hdemag, Hdemag2, Hdemag3;

	//-Nxx, -Nyy, -Nzz values at r = r0
	cu_obj<cuReal3> selfDemagCoeff;

public:

	SDemagCUDA_Demag(MeshCUDA* pMeshCUDA_, SDemag_Demag *pSDemag_Demag_);
	~SDemagCUDA_Demag();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);
	
	//-------------------Getters

	//add energy in this module to a running total
	void Add_Energy(cu_obj<cuBReal>& total_energy);

	//-------------------Setters
};

#else

class SDemagCUDA_Demag
{
};

#endif

#endif


