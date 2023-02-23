#include "StrayField_AtomMeshCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "BorisCUDALib.cuh"

#include "Atom_MeshCUDA.h"
#include "MeshDefs.h"

//----------------------- Initialization

__global__ void set_StrayField_AtomMeshCUDA_pointers_kernel(
	ManagedAtom_MeshCUDA& cuaMesh, cuVEC<cuReal3>& strayField)
{
	if (threadIdx.x == 0) cuaMesh.pstrayField = &strayField;
}

void StrayField_AtomMeshCUDA::set_StrayField_AtomMeshCUDA_pointers(void)
{
	set_StrayField_AtomMeshCUDA_pointers_kernel <<< 1, CUDATHREADS >>>
		(paMeshCUDA->cuaMesh, strayField);
}

//----------------------- Computation

__global__ void UpdateStrayField_ASC_kernel(ManagedAtom_MeshCUDA& cuaMesh, ManagedModulesCUDA& cuModule, cuVEC<cuReal3>& strayField, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff1.linear_size()) {

		cuReal3 Hstray = strayField[idx];

		Heff1[idx] += Hstray;

		if (do_reduction) {

			int non_empty_cells = M1.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MUB_MU0 * M1[idx] * Hstray / (non_empty_cells * M1.h.dim());
		}

		if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = Hstray;
		if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = -(cuBReal)MUB_MU0 * M1[idx] * Hstray / M1.h.dim();
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

void StrayField_AtomMeshCUDA::UpdateFieldCUDA(void)
{
	if (paMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();

		UpdateStrayField_ASC_kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, cuModule, strayField, true);
	}
	else {

		UpdateStrayField_ASC_kernel <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, cuModule, strayField, false);
	}
}

#endif

#endif