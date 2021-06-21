#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_DEMAG

#include "ModulesCUDA.h"

#include "ConvolutionCUDA.h"
#include "DemagKernelCUDA.h"

class MeshCUDA;
class Demag;

class DemagCUDA :
	public ModulesCUDA,	
	public ConvolutionCUDA<DemagCUDA, DemagKernelCUDA>
{

private:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

	//point to cpu version of this module
	Demag *pDemag;

	//Evaluation speedup mode data

	//1 Hdemag: no extrapolation, just save evaluation and reuse
	//2 Hdemag: linear extrapolation, need 2
	//3 Hdemag: quadratic extrapolation, need 3
	cu_obj<cuVEC<cuReal3>> Hdemag, Hdemag2, Hdemag3;

	//times at which evaluations were done, used for extrapolation
	double time_demag1 = 0.0, time_demag2 = 0.0, time_demag3 = 0.0;

	int num_Hdemag_saved = 0;

	//-Nxx, -Nyy, -Nzz values at r = r0
	cu_obj<cuReal3> selfDemagCoeff;

private:

	void set_DemagCUDA_pointers(void);

	//Add newly computed field to Heff and Heff2, then subtract self demag contribution from it : AFM
	void Demag_EvalSpeedup_AddField_SubSelf(
		cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2, 
		cu_obj<cuVEC<cuReal3>>& HField, 
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2);

	//Add newly computed field to Heff and Heff2, then subtract self demag contribution from it : FM
	void Demag_EvalSpeedup_AddField_SubSelf(
		cu_obj<cuVEC<cuReal3>>& Heff,
		cu_obj<cuVEC<cuReal3>>& HField,
		cu_obj<cuVEC_VC<cuReal3>>& M);

	//Add extrapolated field together with self demag contribution : AFM, QUADRATIC
	void Demag_EvalSpeedup_AddExtrapField_AddSelf(
		cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2, 
		cuBReal a1, cuBReal a2, cuBReal a3, 
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2);

	//Add extrapolated field together with self demag contribution : FM, QUADRATIC
	void Demag_EvalSpeedup_AddExtrapField_AddSelf(
		cu_obj<cuVEC<cuReal3>>& Heff,
		cuBReal a1, cuBReal a2, cuBReal a3,
		cu_obj<cuVEC_VC<cuReal3>>& M);

	//Add extrapolated field together with self demag contribution : AFM, LINEAR
	void Demag_EvalSpeedup_AddExtrapField_AddSelf(
		cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
		cuBReal a1, cuBReal a2,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2);

	//Add extrapolated field together with self demag contribution : FM, LINEAR
	void Demag_EvalSpeedup_AddExtrapField_AddSelf(
		cu_obj<cuVEC<cuReal3>>& Heff,
		cuBReal a1, cuBReal a2,
		cu_obj<cuVEC_VC<cuReal3>>& M);

	//Add extrapolated field together with self demag contribution : AFM, STEP
	void Demag_EvalSpeedup_AddExtrapField_AddSelf(
		cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
		cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2);

	//Add extrapolated field together with self demag contribution : FM, STEP
	void Demag_EvalSpeedup_AddExtrapField_AddSelf(
		cu_obj<cuVEC<cuReal3>>& Heff,
		cu_obj<cuVEC_VC<cuReal3>>& M);

public:

	DemagCUDA(MeshCUDA* pMeshCUDA_, Demag *pDemag_);
	~DemagCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);
};

#else

class DemagCUDA
{
};

#endif

#endif


