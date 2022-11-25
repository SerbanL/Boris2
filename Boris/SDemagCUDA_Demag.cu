#include "SDemagCUDA_Demag.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SDEMAG

#include "BorisCUDALib.cuh"
#include "MeshCUDA.h"
#include "Atom_MeshCUDA.h"

//----------------------- Initialization

__global__ void set_SDemag_DemagCUDA_pointers_kernel(
	ManagedMeshCUDA& cuMesh, cuVEC<cuReal3>& Module_Heff)
{
	if (threadIdx.x == 0) cuMesh.pDemag_Heff = &Module_Heff;
}

__global__ void set_SDemag_DemagCUDA_pointers_atomistic_kernel(
	ManagedAtom_MeshCUDA& cuaMesh, cuVEC<cuReal3>& Module_Heff)
{
	if (threadIdx.x == 0) cuaMesh.pAtom_Demag_Heff = &Module_Heff;
}

void SDemagCUDA_Demag::set_SDemag_DemagCUDA_pointers(void)
{
	if (pMeshCUDA) {

		set_SDemag_DemagCUDA_pointers_kernel <<< 1, CUDATHREADS >>>
			(pMeshCUDA->cuMesh, Module_Heff);
	}
	else if (paMeshCUDA) {

		set_SDemag_DemagCUDA_pointers_atomistic_kernel <<< 1, CUDATHREADS >>>
			(paMeshCUDA->cuaMesh, Module_Heff);
	}
}

//-------------------Getters

__global__ void Add_Energy_Kernel(cuBReal& energy, cuBReal& total_energy)
{
	if (threadIdx.x == 0) total_energy += energy;
}

//add energy in this module to a running total
void SDemagCUDA_Demag::Add_Energy(cu_obj<cuBReal>& total_energy)
{
	Add_Energy_Kernel <<< (1 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (energy, total_energy);
}

#endif

#endif