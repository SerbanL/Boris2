#include "SDemagCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SDEMAG

#include "BorisCUDALib.cuh"

//Subtract self demag contribution from H; M is cuVEC, not cuVEC_VC, since it's a transfer mesh : transfer mode
__global__ void SDemag_EvalSpeedup_SubSelf_Kernel(cuVEC<cuReal3>& H, cuVEC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		H[idx] -= (M[idx] & selfDemagCoeff);
	}
}

//Add H to Heff and Heff2, then subtract self demag contribution : AFM, non-transfer mode
__global__ void SDemag_EvalSpeedup_AddField_SubSelf_Kernel(
	cuVEC<cuReal3>& Heff, cuVEC<cuReal3>& Heff2,
	cuVEC<cuReal3>& H,
	cuVEC_VC<cuReal3>& M, cuVEC_VC<cuReal3>& M2, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		Heff[idx] += H[idx];
		Heff2[idx] += H[idx];

		H[idx] -= (((M[idx] + M2[idx]) / 2) & selfDemagCoeff);
	}
}

//Add H to Heff, then subtract self demag contribution : FM, non-transfer mode
__global__ void SDemag_EvalSpeedup_AddField_SubSelf_Kernel(
	cuVEC<cuReal3>& Heff,
	cuVEC<cuReal3>& H,
	cuVEC_VC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		Heff[idx] += H[idx];

		H[idx] -= (M[idx] & selfDemagCoeff);
	}
}

//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, QUINTIC
__global__ void SDemag_EvalSpeedup_SetExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& H_M,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4, cuVEC<cuReal3>& Hdemag5, cuVEC<cuReal3>& Hdemag6,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6,
	cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		H_M[idx] = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + Hdemag6[idx] * a6 + (H_M[idx] & selfDemagCoeff);
	}
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, QUINTIC
__global__ void SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff, cuVEC<cuReal3>& Heff2,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4, cuVEC<cuReal3>& Hdemag5, cuVEC<cuReal3>& Hdemag6,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6,
	cuVEC_VC<cuReal3>& M, cuVEC_VC<cuReal3>& M2, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hvalue = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + Hdemag6[idx] * a6 + (((M[idx] + M2[idx]) / 2) & selfDemagCoeff);

		Heff[idx] += Hvalue;
		Heff2[idx] += Hvalue;
	}
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, QUINTIC
__global__ void SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4, cuVEC<cuReal3>& Hdemag5, cuVEC<cuReal3>& Hdemag6,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6,
	cuVEC_VC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hvalue = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + Hdemag6[idx] * a6 + (M[idx] & selfDemagCoeff);

		Heff[idx] += Hvalue;
	}
}


//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, QUARTIC
__global__ void SDemag_EvalSpeedup_SetExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& H_M,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4, cuVEC<cuReal3>& Hdemag5,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5,
	cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		H_M[idx] = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + (H_M[idx] & selfDemagCoeff);
	}
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, QUARTIC
__global__ void SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff, cuVEC<cuReal3>& Heff2,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4, cuVEC<cuReal3>& Hdemag5,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5,
	cuVEC_VC<cuReal3>& M, cuVEC_VC<cuReal3>& M2, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hvalue = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + (((M[idx] + M2[idx]) / 2) & selfDemagCoeff);

		Heff[idx] += Hvalue;
		Heff2[idx] += Hvalue;
	}
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, QUARTIC
__global__ void SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4, cuVEC<cuReal3>& Hdemag5,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5,
	cuVEC_VC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hvalue = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + (M[idx] & selfDemagCoeff);

		Heff[idx] += Hvalue;
	}
}


//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, CUBIC
__global__ void SDemag_EvalSpeedup_SetExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& H_M,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4,
	cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		H_M[idx] = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + (H_M[idx] & selfDemagCoeff);
	}
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, CUBIC
__global__ void SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff, cuVEC<cuReal3>& Heff2,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4,
	cuVEC_VC<cuReal3>& M, cuVEC_VC<cuReal3>& M2, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hvalue = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + (((M[idx] + M2[idx]) / 2) & selfDemagCoeff);

		Heff[idx] += Hvalue;
		Heff2[idx] += Hvalue;
	}
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, CUBIC
__global__ void SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4,
	cuVEC_VC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hvalue = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + (M[idx] & selfDemagCoeff);

		Heff[idx] += Hvalue;
	}
}


//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, QUADRATIC
__global__ void SDemag_EvalSpeedup_SetExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& H_M,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3,
	cuBReal a1, cuBReal a2, cuBReal a3,
	cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		H_M[idx] = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + (H_M[idx] & selfDemagCoeff);
	}
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, QUADRATIC
__global__ void SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff, cuVEC<cuReal3>& Heff2,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3,
	cuBReal a1, cuBReal a2, cuBReal a3,
	cuVEC_VC<cuReal3>& M, cuVEC_VC<cuReal3>& M2, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hvalue = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + (((M[idx] + M2[idx]) / 2) & selfDemagCoeff);

		Heff[idx] += Hvalue;
		Heff2[idx] += Hvalue;
	}
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, QUADRATIC
__global__ void SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3,
	cuBReal a1, cuBReal a2, cuBReal a3,
	cuVEC_VC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hvalue = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + (M[idx] & selfDemagCoeff);

		Heff[idx] += Hvalue;
	}
}

//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, LINEAR
__global__ void SDemag_EvalSpeedup_SetExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& H_M,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2,
	cuBReal a1, cuBReal a2,
	cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		H_M[idx] = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + (H_M[idx] & selfDemagCoeff);
	}
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, LINEAR
__global__ void SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff, cuVEC<cuReal3>& Heff2,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2,
	cuBReal a1, cuBReal a2,
	cuVEC_VC<cuReal3>& M, cuVEC_VC<cuReal3>& M2, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hvalue = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + (((M[idx] + M2[idx]) / 2) & selfDemagCoeff);

		Heff[idx] += Hvalue;
		Heff2[idx] += Hvalue;
	}
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, LINEAR
__global__ void SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2,
	cuBReal a1, cuBReal a2,
	cuVEC_VC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hvalue = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + (M[idx] & selfDemagCoeff);

		Heff[idx] += Hvalue;
	}
}

//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, STEP
__global__ void SDemag_EvalSpeedup_SetExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& H_M,
	cuVEC<cuReal3>& Hdemag,
	cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		H_M[idx] = Hdemag[idx] + (H_M[idx] & selfDemagCoeff);
	}
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, STEP
__global__ void SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff, cuVEC<cuReal3>& Heff2,
	cuVEC<cuReal3>& Hdemag,
	cuVEC_VC<cuReal3>& M, cuVEC_VC<cuReal3>& M2, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hvalue = Hdemag[idx] + (((M[idx] + M2[idx]) / 2) & selfDemagCoeff);

		Heff[idx] += Hvalue;
		Heff2[idx] += Hvalue;
	}
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, STEP
__global__ void SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& Heff,
	cuVEC<cuReal3>& Hdemag,
	cuVEC_VC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Hdemag.linear_size()) {

		cuReal3 Hvalue = Hdemag[idx] + (M[idx] & selfDemagCoeff);

		Heff[idx] += Hvalue;
	}
}

//----------------------- LAUNCHERS

//Subtract self demag contribution from H; M is cuVEC, not cuVEC_VC, since it's a transfer mesh : transfer mode
void SDemagCUDA::SDemag_EvalSpeedup_SubSelf(size_t size, cu_obj<cuVEC<cuReal3>>& H, cu_obj<cuVEC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_SubSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H, M, selfDemagCoeff);
}

//Add H to Heff and Heff2, then subtract self demag contribution : AFM, non-transfer mode
void SDemagCUDA::SDemag_EvalSpeedup_AddField_SubSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
	cu_obj<cuVEC<cuReal3>>& H,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_AddField_SubSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Heff2, H, M, M2, selfDemagCoeff);
}

//Add H to Heff, then subtract self demag contribution : FM, non-transfer mode
void SDemagCUDA::SDemag_EvalSpeedup_AddField_SubSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& Heff,
	cu_obj<cuVEC<cuReal3>>& H,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_AddField_SubSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, H, M, selfDemagCoeff);
}

//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, QUINTIC
void SDemagCUDA::SDemag_EvalSpeedup_SetExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& H_M,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4, cu_obj<cuVEC<cuReal3>>& Hdemag5, cu_obj<cuVEC<cuReal3>>& Hdemag6,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6,
	cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_SetExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H_M, Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, Hdemag6, a1, a2, a3, a4, a5, a6, selfDemagCoeff);
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, QUINTIC
void SDemagCUDA::SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4, cu_obj<cuVEC<cuReal3>>& Hdemag5, cu_obj<cuVEC<cuReal3>>& Hdemag6,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Heff2, Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, Hdemag6, a1, a2, a3, a4, a5, a6, M, M2, selfDemagCoeff);
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff and : FM, non-transfer mode, QUINTIC
void SDemagCUDA::SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& Heff,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4, cu_obj<cuVEC<cuReal3>>& Hdemag5, cu_obj<cuVEC<cuReal3>>& Hdemag6,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, Hdemag6, a1, a2, a3, a4, a5, a6, M, selfDemagCoeff);
}

//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, QUARTIC
void SDemagCUDA::SDemag_EvalSpeedup_SetExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& H_M,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4, cu_obj<cuVEC<cuReal3>>& Hdemag5,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5,
	cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_SetExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H_M, Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, a1, a2, a3, a4, a5, selfDemagCoeff);
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, QUARTIC
void SDemagCUDA::SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4, cu_obj<cuVEC<cuReal3>>& Hdemag5,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Heff2, Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, a1, a2, a3, a4, a5, M, M2, selfDemagCoeff);
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, QUARTIC
void SDemagCUDA::SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& Heff,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4, cu_obj<cuVEC<cuReal3>>& Hdemag5,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, a1, a2, a3, a4, a5, M, selfDemagCoeff);
}


//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, CUBIC
void SDemagCUDA::SDemag_EvalSpeedup_SetExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& H_M,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4,
	cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_SetExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H_M, Hdemag, Hdemag2, Hdemag3, Hdemag4, a1, a2, a3, a4, selfDemagCoeff);
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, CUBIC
void SDemagCUDA::SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Heff2, Hdemag, Hdemag2, Hdemag3, Hdemag4, a1, a2, a3, a4, M, M2, selfDemagCoeff);
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, CUBIC
void SDemagCUDA::SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& Heff,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3, cu_obj<cuVEC<cuReal3>>& Hdemag4,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Hdemag, Hdemag2, Hdemag3, Hdemag4, a1, a2, a3, a4, M, selfDemagCoeff);
}


//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, QUADRATIC
void SDemagCUDA::SDemag_EvalSpeedup_SetExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& H_M,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3,
	cuBReal a1, cuBReal a2, cuBReal a3,
	cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_SetExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H_M, Hdemag, Hdemag2, Hdemag3, a1, a2, a3, selfDemagCoeff);
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, QUADRATIC
void SDemagCUDA::SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3,
	cuBReal a1, cuBReal a2, cuBReal a3,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Heff2, Hdemag, Hdemag2, Hdemag3, a1, a2, a3, M, M2, selfDemagCoeff);
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, QUADRATIC
void SDemagCUDA::SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& Heff,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2, cu_obj<cuVEC<cuReal3>>& Hdemag3,
	cuBReal a1, cuBReal a2, cuBReal a3,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Hdemag, Hdemag2, Hdemag3, a1, a2, a3, M, selfDemagCoeff);
}

//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, LINEAR
void SDemagCUDA::SDemag_EvalSpeedup_SetExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& H_M,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2,
	cuBReal a1, cuBReal a2,
	cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_SetExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H_M, Hdemag, Hdemag2, a1, a2, selfDemagCoeff);
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, LINEAR
void SDemagCUDA::SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2,
	cuBReal a1, cuBReal a2,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Heff2, Hdemag, Hdemag2, a1, a2, M, M2, selfDemagCoeff);
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, LINEAR
void SDemagCUDA::SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& Heff,
	cu_obj<cuVEC<cuReal3>>& Hdemag, cu_obj<cuVEC<cuReal3>>& Hdemag2,
	cuBReal a1, cuBReal a2,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Hdemag, Hdemag2, a1, a2, M, selfDemagCoeff);
}

//Calculate extrapolated field in H_M, also adding in the self contribution given that H_M initially contains M : AFM or FM transfer mode, STEP
void SDemagCUDA::SDemag_EvalSpeedup_SetExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& H_M,
	cu_obj<cuVEC<cuReal3>>& Hdemag,
	cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_SetExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H_M, Hdemag, selfDemagCoeff);
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff and Heff2 : AFM, non-transfer mode, STEP
void SDemagCUDA::SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& Heff, cu_obj<cuVEC<cuReal3>>& Heff2,
	cu_obj<cuVEC<cuReal3>>& Hdemag,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuVEC_VC<cuReal3>>& M2, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Heff2, Hdemag, M, M2, selfDemagCoeff);
}

//Calculate extrapolated field, together with self demag contribution, and set it in Heff : FM, non-transfer mode, STEP
void SDemagCUDA::SDemag_EvalSpeedup_AddExtrapField_AddSelf(size_t size,
	cu_obj<cuVEC<cuReal3>>& Heff,
	cu_obj<cuVEC<cuReal3>>& Hdemag,
	cu_obj<cuVEC_VC<cuReal3>>& M, cu_obj<cuReal3>& selfDemagCoeff)
{
	SDemag_EvalSpeedup_AddExtrapField_AddSelf_Kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(Heff, Hdemag, M, selfDemagCoeff);
}

#endif

#endif