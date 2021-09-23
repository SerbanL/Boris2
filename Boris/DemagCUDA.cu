#include "DemagCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_DEMAG

#include "BorisCUDALib.cuh"
#include "MeshCUDA.h"

//----------------------- Initialization

__global__ void set_DemagCUDA_pointers_kernel(
	ManagedMeshCUDA& cuMesh, cuVEC<cuReal3>& Module_Heff)
{
	if (threadIdx.x == 0) cuMesh.pDemag_Heff = &Module_Heff;
}

void DemagCUDA::set_DemagCUDA_pointers(void)
{
	set_DemagCUDA_pointers_kernel <<< 1, CUDATHREADS >>>
		(pMeshCUDA->cuMesh, Module_Heff);
}

//----------------------- LAUNCHERS

//SUBTRACT SELF DEMAG

//Add newly computed field to Heff and Heff2, then subtract self demag contribution from it : AFM
__global__ void Demag_EvalSpeedup_AddField_SubSelf_Kernel(
	cuVEC<cuReal3>& Heff, cuVEC<cuReal3>& Heff2,
	cuVEC<cuReal3>& HField,
	cuVEC_VC<cuReal3>& M, cuVEC_VC<cuReal3>& M2, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < HField.linear_size()) {

		Heff[idx] += HField[idx];
		Heff2[idx] += HField[idx];

		HField[idx] -= (selfDemagCoeff & (M[idx] + M2[idx]) / 2);
	}
}

//Add newly computed field to Heff and Heff2, then subtract self demag contribution from it : FM
__global__ void Demag_EvalSpeedup_AddField_SubSelf_Kernel(
	cuVEC<cuReal3>& Heff,
	cuVEC<cuReal3>& HField,
	cuVEC_VC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < HField.linear_size()) {

		Heff[idx] += HField[idx];

		HField[idx] -= (selfDemagCoeff & M[idx]);
	}
}

//QUINTIC

//Add extrapolated field together with self demag contribution : AFM, QUINTIC
__global__ void Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff, cuVEC<cuReal3>& Heff2,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4, cuVEC<cuReal3>& Hdemag5, cuVEC<cuReal3>& Hdemag6,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6,
	cuVEC_VC<cuReal3>& M, cuVEC_VC<cuReal3>& M2, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hdemag_value = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + Hdemag6[idx] * a6 + (selfDemagCoeff & (M[idx] + M2[idx]) / 2);
		Heff[idx] += Hdemag_value;
		Heff2[idx] += Hdemag_value;
	}
}

//Add extrapolated field together with self demag contribution : FM, QUINTIC
__global__ void Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4, cuVEC<cuReal3>& Hdemag5, cuVEC<cuReal3>& Hdemag6,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6,
	cuVEC_VC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		Heff[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + Hdemag6[idx] * a6 + (selfDemagCoeff & M[idx]);
	}
}

//QUARTIC

//Add extrapolated field together with self demag contribution : AFM, QUARTIC
__global__ void Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff, cuVEC<cuReal3>& Heff2,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4, cuVEC<cuReal3>& Hdemag5,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5,
	cuVEC_VC<cuReal3>& M, cuVEC_VC<cuReal3>& M2, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hdemag_value = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + (selfDemagCoeff & (M[idx] + M2[idx]) / 2);
		Heff[idx] += Hdemag_value;
		Heff2[idx] += Hdemag_value;
	}
}

//Add extrapolated field together with self demag contribution : FM, QUARTIC
__global__ void Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4, cuVEC<cuReal3>& Hdemag5,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5,
	cuVEC_VC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		Heff[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + (selfDemagCoeff & M[idx]);
	}
}

//CUBIC

//Add extrapolated field together with self demag contribution : AFM, CUBIC
__global__ void Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff, cuVEC<cuReal3>& Heff2,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4,
	cuVEC_VC<cuReal3>& M, cuVEC_VC<cuReal3>& M2, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hdemag_value = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + (selfDemagCoeff & (M[idx] + M2[idx]) / 2);
		Heff[idx] += Hdemag_value;
		Heff2[idx] += Hdemag_value;
	}
}

//Add extrapolated field together with self demag contribution : FM, CUBIC
__global__ void Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4,
	cuVEC_VC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		Heff[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + (selfDemagCoeff & M[idx]);
	}
}

//QUADRATIC

//Add extrapolated field together with self demag contribution : AFM, QUADRATIC
__global__ void Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff, cuVEC<cuReal3>& Heff2,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3,
	cuBReal a1, cuBReal a2, cuBReal a3,
	cuVEC_VC<cuReal3>& M, cuVEC_VC<cuReal3>& M2, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hdemag_value = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + (selfDemagCoeff & (M[idx] + M2[idx]) / 2);
		Heff[idx] += Hdemag_value;
		Heff2[idx] += Hdemag_value;
	}
}

//Add extrapolated field together with self demag contribution : FM, QUADRATIC
__global__ void Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3,
	cuBReal a1, cuBReal a2, cuBReal a3,
	cuVEC_VC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		Heff[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + (selfDemagCoeff & M[idx]);
	}
}

//LINEAR

//Add extrapolated field together with self demag contribution : AFM, LINEAR
__global__ void Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff, cuVEC<cuReal3>& Heff2,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2,
	cuBReal a1, cuBReal a2,
	cuVEC_VC<cuReal3>& M, cuVEC_VC<cuReal3>& M2, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hdemag_value = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + (selfDemagCoeff & (M[idx] + M2[idx]) / 2);
		Heff[idx] += Hdemag_value;
		Heff2[idx] += Hdemag_value;
	}
}

//Add extrapolated field together with self demag contribution : FM, LINEAR
__global__ void Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2,
	cuBReal a1, cuBReal a2,
	cuVEC_VC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		Heff[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + (selfDemagCoeff & M[idx]);
	}
}

//STEP

//Add extrapolated field together with self demag contribution : AFM, STEP
__global__ void Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff, cuVEC<cuReal3>& Heff2,
	cuVEC<cuReal3>& Hdemag,
	cuVEC_VC<cuReal3>& M, cuVEC_VC<cuReal3>& M2, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hdemag_value = Hdemag[idx] + (selfDemagCoeff & (M[idx] + M2[idx]) / 2);
		Heff[idx] += Hdemag_value;
		Heff2[idx] += Hdemag_value;
	}
}

//Add extrapolated field together with self demag contribution : FM, STEP
__global__ void Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff,
	cuVEC<cuReal3>& Hdemag,
	cuVEC_VC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		Heff[idx] += Hdemag[idx] + (selfDemagCoeff & M[idx]);
	}
}

//----------------------- LAUNCHERS

//Add newly computed field to Heff and Heff2, then subtract self demag contribution from it : AFM
void DemagCUDA::Demag_EvalSpeedup_AddField_SubSelf(
	cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
	cu_obj<cuVEC<cuReal3>>& HField,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2)
{
	Demag_EvalSpeedup_AddField_SubSelf_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Heff2, HField, M, M2, selfDemagCoeff);
}

//Add newly computed field to Heff and Heff2, then subtract self demag contribution from it : FM
void DemagCUDA::Demag_EvalSpeedup_AddField_SubSelf(
	cu_obj<cuVEC<cuReal3>>& Heff,
	cu_obj<cuVEC<cuReal3>>& HField,
	cu_obj<cuVEC_VC<cuReal3>>& M)
{
	Demag_EvalSpeedup_AddField_SubSelf_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, HField, M, selfDemagCoeff);
}

//QUINTIC

//Add extrapolated field together with self demag contribution : AFM, QUINTIC
void DemagCUDA::Demag_EvalSpeedup_AddExtrapField_AddSelf(
	cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2)
{
	Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Heff2, Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, Hdemag6, a1, a2, a3, a4, a5, a6, M, M2, selfDemagCoeff);
}

//Add extrapolated field together with self demag contribution : FM, QUINTIC
void DemagCUDA::Demag_EvalSpeedup_AddExtrapField_AddSelf(
	cu_obj<cuVEC<cuReal3>>& Heff,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6,
	cu_obj<cuVEC_VC<cuReal3>>& M)
{
	Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel << < (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >
		(Heff, Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, Hdemag6, a1, a2, a3, a4, a5, a6, M, selfDemagCoeff);
}

//QUARTIC

//Add extrapolated field together with self demag contribution : AFM, QUARTIC
void DemagCUDA::Demag_EvalSpeedup_AddExtrapField_AddSelf(
	cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2)
{
	Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Heff2, Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, a1, a2, a3, a4, a5, M, M2, selfDemagCoeff);
}

//Add extrapolated field together with self demag contribution : FM, QUARTIC
void DemagCUDA::Demag_EvalSpeedup_AddExtrapField_AddSelf(
	cu_obj<cuVEC<cuReal3>>& Heff,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5,
	cu_obj<cuVEC_VC<cuReal3>>& M)
{
	Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, a1, a2, a3, a4, a5, M, selfDemagCoeff);
}

//CUBIC

//Add extrapolated field together with self demag contribution : AFM, CUBIC
void DemagCUDA::Demag_EvalSpeedup_AddExtrapField_AddSelf(
	cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2)
{
	Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Heff2, Hdemag, Hdemag2, Hdemag3, Hdemag4, a1, a2, a3, a4, M, M2, selfDemagCoeff);
}

//Add extrapolated field together with self demag contribution : FM, CUBIC
void DemagCUDA::Demag_EvalSpeedup_AddExtrapField_AddSelf(
	cu_obj<cuVEC<cuReal3>>& Heff,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4,
	cu_obj<cuVEC_VC<cuReal3>>& M)
{
	Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Hdemag, Hdemag2, Hdemag3, Hdemag4, a1, a2, a3, a4, M, selfDemagCoeff);
}

//QUADRATIC

//Add extrapolated field together with self demag contribution : AFM, QUADRATIC
void DemagCUDA::Demag_EvalSpeedup_AddExtrapField_AddSelf(
	cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
	cuBReal a1, cuBReal a2, cuBReal a3,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2)
{
	Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Heff2, Hdemag, Hdemag2, Hdemag3, a1, a2, a3, M, M2, selfDemagCoeff);
}

//Add extrapolated field together with self demag contribution : FM, QUADRATIC
void DemagCUDA::Demag_EvalSpeedup_AddExtrapField_AddSelf(
	cu_obj<cuVEC<cuReal3>>& Heff,
	cuBReal a1, cuBReal a2, cuBReal a3,
	cu_obj<cuVEC_VC<cuReal3>>& M)
{
	Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Hdemag, Hdemag2, Hdemag3, a1, a2, a3, M, selfDemagCoeff);
}

//LINEAR

//Add extrapolated field together with self demag contribution : AFM, LINEAR
void DemagCUDA::Demag_EvalSpeedup_AddExtrapField_AddSelf(
	cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
	cuBReal a1, cuBReal a2,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2)
{
	Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Heff2, Hdemag, Hdemag2, a1, a2, M, M2, selfDemagCoeff);
}

//Add extrapolated field together with self demag contribution : FM, LINEAR
void DemagCUDA::Demag_EvalSpeedup_AddExtrapField_AddSelf(
	cu_obj<cuVEC<cuReal3>>& Heff,
	cuBReal a1, cuBReal a2,
	cu_obj<cuVEC_VC<cuReal3>>& M)
{
	Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Hdemag, Hdemag2, a1, a2, M, selfDemagCoeff);
}

//STEP

//Add extrapolated field together with self demag contribution : AFM, STEP
void DemagCUDA::Demag_EvalSpeedup_AddExtrapField_AddSelf(
	cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2)
{
	Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Heff2, Hdemag, M, M2, selfDemagCoeff);
}

//Add extrapolated field together with self demag contribution : FM, STEP
void DemagCUDA::Demag_EvalSpeedup_AddExtrapField_AddSelf(
	cu_obj<cuVEC<cuReal3>>& Heff,
	cu_obj<cuVEC_VC<cuReal3>>& M)
{
	Demag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Hdemag, M, selfDemagCoeff);
}

#endif

#endif