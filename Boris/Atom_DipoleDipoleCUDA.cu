#include "Atom_DipoleDipoleCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_ATOM_DIPOLEDIPOLE) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"
#include "Atom_MeshCUDA.h"

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