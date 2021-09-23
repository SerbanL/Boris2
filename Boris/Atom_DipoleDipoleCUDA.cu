#include "Atom_DipoleDipoleCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ATOM_DIPOLEDIPOLE) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"
#include "Atom_MeshCUDA.h"

//----------------------- Initialization

__global__ void set_Atom_DipoleDipoleCUDA_pointers_kernel(
	ManagedAtom_MeshCUDA& cuaMesh, cuVEC<cuReal3>& Module_Heff)
{
	if (threadIdx.x == 0) cuaMesh.pAtom_Demag_Heff = &Module_Heff;
}

void Atom_DipoleDipoleCUDA::set_Atom_DipoleDipoleCUDA_pointers(void)
{
	set_Atom_DipoleDipoleCUDA_pointers_kernel <<< 1, CUDATHREADS >>>
		(paMeshCUDA->cuaMesh, Module_Heff);
}

//----------------------- Auxiliary

__global__ void Energy_to_EnergyDensity_Kernel(cuBReal& energy, cuVEC<cuReal3>& V)
{
	if (threadIdx.x == 0) energy *= (cuBReal)MUB / V.h.dim();
}

//convert value in energy to energy density by dividing by cellsize volume of V
void Atom_DipoleDipoleCUDA::Energy_to_EnergyDensity(cu_obj<cuVEC<cuReal3>>& V)
{
	Energy_to_EnergyDensity_Kernel <<< 1, CUDATHREADS >>> (energy, V);
}

//----------------------- KERNELS

__global__ void Atom_DipoleDipole_EvalSpeedup_SubSelf_Kernel(
	cuVEC<cuReal3>& H,
	cuVEC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		H[idx] -= (selfDemagCoeff & M[idx]);
	}
}

//Macrocell mode, QUINTIC
__global__ void Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& H,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4, cuVEC<cuReal3>& Hdemag5, cuVEC<cuReal3>& Hdemag6,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6,
	cuVEC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		H[idx] = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + Hdemag6[idx] * a6 + (selfDemagCoeff & M[idx]);
	}
}

//Macrocell mode, QUARTIC
__global__ void Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& H,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4, cuVEC<cuReal3>& Hdemag5,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5,
	cuVEC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		H[idx] = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + (selfDemagCoeff & M[idx]);
	}
}

//Macrocell mode, CUBIC
__global__ void Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& H,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4,
	cuVEC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		H[idx] = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + (selfDemagCoeff & M[idx]);
	}
}

//Macrocell mode, QUADRATIC
__global__ void Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& H,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3,
	cuBReal a1, cuBReal a2, cuBReal a3,
	cuVEC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		H[idx] = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + (selfDemagCoeff & M[idx]);
	}
}

//Macrocell mode, LINEAR
__global__ void Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& H,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2,
	cuBReal a1, cuBReal a2,
	cuVEC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		H[idx] = Hdemag[idx] * a1 + Hdemag2[idx] * a2 + (selfDemagCoeff & M[idx]);
	}
}

//Macrocell mode, STEP
__global__ void Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf_Kernel(
	cuVEC<cuReal3>& H,
	cuVEC<cuReal3>& Hdemag,
	cuVEC<cuReal3>& M, cuReal3& selfDemagCoeff)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		H[idx] = Hdemag[idx] + (selfDemagCoeff & M[idx]);
	}
}

//QUINTIC
__global__ void Atom_DipoleDipole_EvalSpeedup_AddExtrapField_Kernel(
	cuVEC<cuReal3>& H,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4, cuVEC<cuReal3>& Hdemag5, cuVEC<cuReal3>& Hdemag6,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		H[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5 + Hdemag6[idx] * a6;
	}
}

//QUARTIC
__global__ void Atom_DipoleDipole_EvalSpeedup_AddExtrapField_Kernel(
	cuVEC<cuReal3>& H,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4, cuVEC<cuReal3>& Hdemag5,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		H[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4 + Hdemag5[idx] * a5;
	}
}

//CUBIC
__global__ void Atom_DipoleDipole_EvalSpeedup_AddExtrapField_Kernel(
	cuVEC<cuReal3>& H,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3, cuVEC<cuReal3>& Hdemag4,
	cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		H[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3 + Hdemag4[idx] * a4;
	}
}

//QUADRATIC
__global__ void Atom_DipoleDipole_EvalSpeedup_AddExtrapField_Kernel(
	cuVEC<cuReal3>& H,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2, cuVEC<cuReal3>& Hdemag3,
	cuBReal a1, cuBReal a2, cuBReal a3)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		H[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2 + Hdemag3[idx] * a3;
	}
}

//LINEAR
__global__ void Atom_DipoleDipole_EvalSpeedup_AddExtrapField_Kernel(
	cuVEC<cuReal3>& H,
	cuVEC<cuReal3>& Hdemag, cuVEC<cuReal3>& Hdemag2,
	cuBReal a1, cuBReal a2)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		H[idx] += Hdemag[idx] * a1 + Hdemag2[idx] * a2;
	}
}

//STEP
__global__ void Atom_DipoleDipole_EvalSpeedup_AddExtrapField_Kernel(
	cuVEC<cuReal3>& H,
	cuVEC<cuReal3>& Hdemag)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < H.linear_size()) {

		H[idx] += Hdemag[idx];
	}
}

//----------------------- LAUNCHERS

//Macrocell mode: subtract self contribution from calculated field
void Atom_DipoleDipoleCUDA::Atom_DipoleDipole_EvalSpeedup_SubSelf(cu_obj<cuVEC<cuReal3>>& H)
{
	Atom_DipoleDipole_EvalSpeedup_SubSelf_Kernel <<< (paMeshCUDA->n_dm.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H, M, selfDemagCoeff);
}

//Macrocell mode, QUINTIC: extrapolate field and add self contribution
void Atom_DipoleDipoleCUDA::Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6)
{
	Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf_Kernel <<< (paMeshCUDA->n_dm.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H, Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, Hdemag6, a1, a2, a3, a4, a5, a6, M, selfDemagCoeff);
}


//Macrocell mode, QUARTIC: extrapolate field and add self contribution
void Atom_DipoleDipoleCUDA::Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5)
{
	Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf_Kernel <<< (paMeshCUDA->n_dm.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H, Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, a1, a2, a3, a4, a5, M, selfDemagCoeff);
}

//Macrocell mode, CUBIC: extrapolate field and add self contribution
void Atom_DipoleDipoleCUDA::Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4)
{
	Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf_Kernel <<< (paMeshCUDA->n_dm.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H, Hdemag, Hdemag2, Hdemag3, Hdemag4, a1, a2, a3, a4, M, selfDemagCoeff);
}

//Macrocell mode, QUADRATIC: extrapolate field and add self contribution
void Atom_DipoleDipoleCUDA::Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2, cuBReal a3)
{
	Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf_Kernel <<< (paMeshCUDA->n_dm.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H, Hdemag, Hdemag2, Hdemag3, a1, a2, a3, M, selfDemagCoeff);
}

//Macrocell mode, LINEAR: extrapolate field and add self contribution
void Atom_DipoleDipoleCUDA::Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2)
{
	Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf_Kernel <<< (paMeshCUDA->n_dm.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H, Hdemag, Hdemag2, a1, a2, M, selfDemagCoeff);
}

//Macrocell mode, STEP: extrapolate field and add self contribution
void Atom_DipoleDipoleCUDA::Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf(cu_obj<cuVEC<cuReal3>>& H)
{
	Atom_DipoleDipole_EvalSpeedup_SetExtrapField_AddSelf_Kernel <<< (paMeshCUDA->n_dm.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H, Hdemag, M, selfDemagCoeff);
}

//QUINTIC: extrapolate field
void Atom_DipoleDipoleCUDA::Atom_DipoleDipole_EvalSpeedup_AddExtrapField(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5, cuBReal a6)
{
	Atom_DipoleDipole_EvalSpeedup_AddExtrapField_Kernel <<< (paMeshCUDA->n_dm.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H, Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, Hdemag6, a1, a2, a3, a4, a5, a6);
}

//QUARTIC: extrapolate field
void Atom_DipoleDipoleCUDA::Atom_DipoleDipole_EvalSpeedup_AddExtrapField(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4, cuBReal a5)
{
	Atom_DipoleDipole_EvalSpeedup_AddExtrapField_Kernel <<< (paMeshCUDA->n_dm.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H, Hdemag, Hdemag2, Hdemag3, Hdemag4, Hdemag5, a1, a2, a3, a4, a5);
}

//CUBIC: extrapolate field
void Atom_DipoleDipoleCUDA::Atom_DipoleDipole_EvalSpeedup_AddExtrapField(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2, cuBReal a3, cuBReal a4)
{
	Atom_DipoleDipole_EvalSpeedup_AddExtrapField_Kernel <<< (paMeshCUDA->n_dm.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H, Hdemag, Hdemag2, Hdemag3, Hdemag4, a1, a2, a3, a4);
}

//QUADRATIC: extrapolate field
void Atom_DipoleDipoleCUDA::Atom_DipoleDipole_EvalSpeedup_AddExtrapField(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2, cuBReal a3)
{
	Atom_DipoleDipole_EvalSpeedup_AddExtrapField_Kernel <<< (paMeshCUDA->n_dm.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H, Hdemag, Hdemag2, Hdemag3, a1, a2, a3);
}

//LINEAR: extrapolate field
void Atom_DipoleDipoleCUDA::Atom_DipoleDipole_EvalSpeedup_AddExtrapField(cu_obj<cuVEC<cuReal3>>& H, cuBReal a1, cuBReal a2)
{
	Atom_DipoleDipole_EvalSpeedup_AddExtrapField_Kernel <<< (paMeshCUDA->n_dm.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H, Hdemag, Hdemag2, a1, a2);
}

//STEP: extrapolate field
void Atom_DipoleDipoleCUDA::Atom_DipoleDipole_EvalSpeedup_AddExtrapField(cu_obj<cuVEC<cuReal3>>& H)
{
	Atom_DipoleDipole_EvalSpeedup_AddExtrapField_Kernel <<< (paMeshCUDA->n_dm.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(H, Hdemag);
}

#endif

#endif