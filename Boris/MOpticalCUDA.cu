#include "MOpticalCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_MOPTICAL

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"
#include "MeshCUDA.h"
#include "MeshParamsControlCUDA.h"

__global__ void MOpticalCUDA_UpdateField_FM(ManagedMeshCUDA& cuMesh, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuBReal cHmo = *cuMesh.pcHmo;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHmo, cHmo);

		//magneto-optical field along z direction only : spatial and time dependence set through the usual material parameter mechanism
		Heff[idx] += cuReal3(0, 0, cHmo);

		if (do_reduction) {

			int non_empty_cells = M.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MU0 * M[idx] * cuReal3(0, 0, cHmo) / non_empty_cells;
		}

		if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = cuReal3(0, 0, cHmo);
		if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = -(cuBReal)MU0 * M[idx] * cuReal3(0, 0, cHmo);
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

__global__ void MOpticalCUDA_UpdateField_AFM(ManagedMeshCUDA& cuMesh, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff.linear_size()) {

		cuBReal cHmo = *cuMesh.pcHmo;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHmo, cHmo);

		//magneto-optical field along z direction only : spatial and time dependence set through the usual material parameter mechanism
		Heff[idx] += cuReal3(0, 0, cHmo);
		Heff2[idx] += cuReal3(0, 0, cHmo);

		if (do_reduction) {

			int non_empty_cells = M.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MU0 * (M[idx] + M2[idx]) * cuReal3(0, 0, cHmo) / (2 * non_empty_cells);
		}

		if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = cuReal3(0, 0, cHmo);
		if (do_reduction && cuModule.pModule_Heff2->linear_size()) (*cuModule.pModule_Heff2)[idx] = cuReal3(0, 0, cHmo);
		if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = -MU0 * M[idx] * cuReal3(0, 0, cHmo);
		if (do_reduction && cuModule.pModule_energy2->linear_size()) (*cuModule.pModule_energy2)[idx] = -MU0 * M2[idx] * cuReal3(0, 0, cHmo);
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//----------------------- UpdateField LAUNCHER

void MOpticalCUDA::UpdateField(void)
{
	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			MOpticalCUDA_UpdateField_AFM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, true);
		}
		else MOpticalCUDA_UpdateField_AFM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, false);
	}

	else {

		if (pMeshCUDA->CurrentTimeStepSolved()) {

			ZeroEnergy();

			MOpticalCUDA_UpdateField_FM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, true);
		}
		else MOpticalCUDA_UpdateField_FM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule, false);
	}
}

#endif

#endif