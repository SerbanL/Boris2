#include "MOpticalCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_MOPTICAL

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

#include "Mesh_FerromagneticCUDA.h"
#include "Mesh_AntiFerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

__global__ void MOpticalCUDA_UpdateField_FMDM(ManagedMeshCUDA& cuMesh)
{
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Heff.linear_size()) {

		cuBReal cHmo = *cuMesh.pcHmo;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHmo, cHmo);

		//magneto-optical field along z direction only : spatial and time dependence set through the usual material parameter mechanism
		Heff[idx] += cuReal3(0, 0, cHmo);
	}
}

__global__ void MOpticalCUDA_UpdateField_AFM(ManagedMeshCUDA& cuMesh)
{
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Heff.linear_size()) {

		cuBReal cHmo = *cuMesh.pcHmo;
		cuMesh.update_parameters_mcoarse(idx, *cuMesh.pcHmo, cHmo);

		//magneto-optical field along z direction only : spatial and time dependence set through the usual material parameter mechanism
		Heff[idx] += cuReal3(0, 0, cHmo);
		Heff2[idx] += cuReal3(0, 0, cHmo);
	}
}

//----------------------- UpdateField LAUNCHER

void MOpticalCUDA::UpdateField(void)
{
	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		MOpticalCUDA_UpdateField_AFM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh);
	}

	else {

		MOpticalCUDA_UpdateField_FMDM <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh);
	}
}

#endif

#endif