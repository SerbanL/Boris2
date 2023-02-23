#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SDEMAG

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
	public ConvolutionCUDA<SDemagCUDA, DemagKernelCUDA>
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
	std::vector<SDemagCUDA_Demag*> pSDemagCUDA_Demag;

	//collect FFT input spaces : after Forward FFT the ffts of M from the individual meshes will be found here
	//These are used as inputs to kernel multiplications. Same order as pSDemag_Demag.
	std::vector<cu_arr<cuBComplex>*> FFT_Spaces_x_Input, FFT_Spaces_y_Input, FFT_Spaces_z_Input;

	//collection of rectangles of meshes, same ordering as for pSDemag_Demag and FFT_Spaces, used in multi-layered convolution
	//these are not necessarily the rectangles of the input M meshes, but are the rectangles of the transfer meshes (M -> transfer -> convolution)
	//for instance in 3D mode, all rectangles in multi-layered convolution must have same size
	//in 2D mode the rectangles can differ in thickness but must have the same xy size
	//thus in 3D mode find largest one and extend all the other rectangles to match (if possible try to have them overlapping in xy-plane projections so we can use kernel symmetries)
	//in 2D mode find largest xy dimension and extend all xy dimensions -> again try to overlap their xy-plane projections
	std::vector<cuRect> Rect_collection;

	//demag kernels used for multilayered convolution, one collection per mesh/SDemag_Demag module. Don't recalculate redundant kernels in the collection.
	std::vector<DemagKernelCollectionCUDA*> kernel_collection;

	//Evaluation speedup mode data

	//times at which evaluations were done, used for extrapolation
	double time_demag1 = 0.0, time_demag2 = 0.0, time_demag3 = 0.0, time_demag4 = 0.0, time_demag5 = 0.0, time_demag6 = 0.0;

	int num_Hdemag_saved = 0;

private:

	//Subtract self demag contribution from H; M is cuVEC, not cuVEC_VC, since it's a transfer mesh : transfer mode
	void SDemag_EvalSpeedup_SubSelf(size_t size, cu_obj<cuVEC<cuReal3>>& H, cu_obj<cuVEC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff);

	//Add H to Heff and Heff2, then subtract self demag contribution : AFM, non-transfer mode
	void SDemag_EvalSpeedup_AddField_SubSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
		cu_obj<cuVEC<cuReal3>>& H,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2, cu_obj<cuReal3>& selfDemagCoeff);

	//Add H to Heff, then subtract self demag contribution : FM, non-transfer mode
	void SDemag_EvalSpeedup_AddField_SubSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& Heff,
		cu_obj<cuVEC<cuReal3>>& H,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, QUINTIC
	void SDemag_EvalSpeedup_SetExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& H_M,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4, cu_obj<cuVEC<cuReal3>>& Hdemag5, cu_obj<cuVEC<cuReal3>>& Hdemag6,
		cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6,
		cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, QUINTIC
	void SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4, cu_obj<cuVEC<cuReal3>>& Hdemag5, cu_obj<cuVEC<cuReal3>>& Hdemag6,
		cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2, cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, QUINTIC
	void SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& Heff,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4, cu_obj<cuVEC<cuReal3>>& Hdemag5, cu_obj<cuVEC<cuReal3>>& Hdemag6,
		cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, QUARTIC
	void SDemag_EvalSpeedup_SetExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& H_M,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4, cu_obj<cuVEC<cuReal3>>& Hdemag5,
		cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5,
		cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, QUARTIC
	void SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4, cu_obj<cuVEC<cuReal3>>& Hdemag5,
		cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2, cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, QUARTIC
	void SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& Heff,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4, cu_obj<cuVEC<cuReal3>>& Hdemag5,
		cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, CUBIC
	void SDemag_EvalSpeedup_SetExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& H_M,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4,
		cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4,
		cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, CUBIC
	void SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4,
		cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2, cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, CUBIC
	void SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& Heff,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4,
		cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, QUADRATIC
	void SDemag_EvalSpeedup_SetExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& H_M,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3,
		cuBReal a1, cuBReal a2, cuBReal a3,
		cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, QUADRATIC
	void SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3,
		cuBReal a1, cuBReal a2, cuBReal a3,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2, cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, QUADRATIC
	void SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& Heff,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3,
		cuBReal a1, cuBReal a2, cuBReal a3,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, LINEAR
	void SDemag_EvalSpeedup_SetExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& H_M,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2,
		cuBReal a1, cuBReal a2,
		cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, LINEAR
	void SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2,
		cuBReal a1, cuBReal a2,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2, cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, LINEAR
	void SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& Heff,
		cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2,
		cuBReal a1, cuBReal a2,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, STEP
	void SDemag_EvalSpeedup_SetExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& H_M,
		cu_obj<cuVEC<cuReal3>>& Hdemag,
		cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, STEP
	void SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
		cu_obj<cuVEC<cuReal3>>& Hdemag,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2, cu_obj<cuReal3>& selfDemagCoeff);

	//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, STEP
	void SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
		cu_obj<cuVEC<cuReal3>>& Heff,
		cu_obj<cuVEC<cuReal3>>& Hdemag,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff);

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


