#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_SDEMAG

#include "ModulesCUDA.h"

#include "ConvolutionCUDA.h"
#include "DemagKernelCUDA.h"
#include "DemagKernelCollectionCUDA.h"

#include "SDemagCUDA_Demag.h"

class SuperMesh;
class SDemag;
class SDemagCUDA_Demag;
class ManagedMeshCUDA;

class SDemagCUDA :
	public ModulesCUDA,
	public ConvolutionCUDA<DemagKernelCUDA>
{

	friend SDemagCUDA_Demag;

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	SuperMesh* pSMesh;

	//the SDemag module (cpu version) holding this CUDA version
	SDemag* pSDemag;

	//--------

	//super-mesh magnetization values used for computing demag field on the super-mesh
	cu_obj<cuVEC<cuReal3>> sm_Vals;

	//total non-empty volume from all meshes participating in convolution
	double total_nonempty_volume = 0.0;

	//--------

	//collection of all SDemagCUDA_Demag modules in individual ferromagnetic meshes
	vector<SDemagCUDA_Demag*> pSDemagCUDA_Demag;

	//collect FFT input spaces : after Forward FFT the ffts of M from the individual meshes will be found here
	//These are used as inputs to kernel multiplications. Same order as pSDemag_Demag.
	std::vector<cu_arr<cuBComplex>*> FFT_Spaces_x_Input, FFT_Spaces_y_Input, FFT_Spaces_z_Input;

	//collection of rectangles of meshes, same ordering as for pSDemag_Demag and FFT_Spaces, used in multi-layered convolution
	//these are not necessarily the rectangles of the input M meshes, but are the rectangles of the transfer meshes (M -> transfer -> convolution)
	//for instance in 3D mode, all rectangles in multi-layered convolution must have same size
	//in 2D mode the rectangles can differ in thickness but must have the same xy size
	//thus in 3D mode find largest one and extend all the other rectangles to match (if possible try to have them overlapping in xy-plane projections so we can use kernel symmetries)
	//in 2D mode find largest xy dimension and extend all xy dimensions -> again try to overlap their xy-plane projections
	vector<cuRect> Rect_collection;

	//demag kernels used for multilayered convolution, one collection per mesh/SDemag_Demag module. Don't recalculate redundant kernels in the collection.
	vector<DemagKernelCollectionCUDA*> kernel_collection;

	//The demag field computed separately (supermesh convolution): at certain steps in the ODE evaluation method we don't need to recalculate the demag field but can use a previous evaluation with an acceptable impact on the numerical error.
	//This mode needs to be enabled by the user, and can be much faster than the default mode. The default mode is to re-evaluate the demag field at every step.
	cu_obj<cuVEC<cuReal3>> Hdemag;

	//when using the evaluation speedup method we must ensure we have a previous Hdemag evaluation available; this flag applies to both supermesh and multilayered convolution
	bool Hdemag_calculated = false;

public:

	SDemagCUDA(SuperMesh* pSMesh_, SDemag* pSDemag_);
	~SDemagCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	void UninitializeAll(void);

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);

	//-------------------

	cu_obj<cuVEC<cuReal3>>& GetDemagField(void) { return sm_Vals; }
};

#else

class SDemagCUDA
{
};

#endif

#endif


