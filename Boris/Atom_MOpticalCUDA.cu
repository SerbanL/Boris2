#include "Atom_MOpticalCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_MOPTICAL) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"
#include "Atom_MeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

__global__ void MOpticalCUDA_UpdateField_Cubic(ManagedAtom_MeshCUDA& cuaMesh)
{
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Heff1.linear_size()) {

		cuBReal cHmo = *cuaMesh.pcHmo;
		cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pcHmo, cHmo);

		//magneto-optical field along z direction only : spatial and time dependence set through the usual material parameter mechanism
		Heff1[idx] += cuReal3(0, 0, cHmo);
	}
}

//----------------------- UpdateField LAUNCHER

void Atom_MOpticalCUDA::UpdateField(void)
{
	MOpticalCUDA_UpdateField_Cubic <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh);
}

#endif

#endif