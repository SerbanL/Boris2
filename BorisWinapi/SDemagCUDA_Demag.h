#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_SDEMAG

#include "ModulesCUDA.h"
#include "BorisCUDALib.h"

#include "ConvolutionCUDA.h"
#include "DemagKernelCollectionCUDA.h"

class FMeshCUDA;
class ManagedMeshCUDA;
class SDemag_Demag;
class SDemagCUDA;

class SDemagCUDA_Demag :
	public ModulesCUDA,
	public ConvolutionCUDA<DemagKernelCollectionCUDA>
{
	friend SDemagCUDA;

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	FMeshCUDA * pMeshCUDA;

	//pointer to cpu version of this module
	SDemag_Demag *pSDemag_Demag;

	//transfer values from M of this mesh to a cuVEC with fixed number of cells -> use same meshing for all layers.
	cu_obj<cuVEC<cuReal3>> transfer;

	//do transfer as M -> transfer -> convolution -> transfer -> Heff if true
	//if false then don't use the transfer VEC but can do M -> convolution -> Heff directly
	//this flag will be set false only if the convolution rect matches that of M and n_common matches the discretisation of M (i.e. in this case a mesh transfer would be pointless).
	bool do_transfer = true;

	//different meshes have different weights when contributing to the total energy density -> ratio of their non-empty volume to total non-empty volume
	cu_obj<cuReal> energy_density_weight;

public:

	SDemagCUDA_Demag(FMeshCUDA* pMeshCUDA_, SDemag_Demag *pSDemag_Demag_);
	~SDemagCUDA_Demag();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	void UpdateField(void) {}
	
	//-------------------Getters

	//add energy in this module to a running total
	void Add_Energy(cu_obj<cuReal>& total_energy);

	//-------------------Setters
};

#else

class SDemagCUDA_Demag
{
};

#endif

#endif


