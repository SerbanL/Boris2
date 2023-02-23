#include "Atom_MOpticalCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_MOPTICAL) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"
#include "Atom_MeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

__global__ void MOpticalCUDA_UpdateField_Cubic(ManagedAtom_MeshCUDA& cuaMesh, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff1.linear_size()) {

		cuBReal cHmo = *cuaMesh.pcHmo;
		cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pcHmo, cHmo);

		//magneto-optical field along z direction only : spatial and time dependence set through the usual material parameter mechanism
		Heff1[idx] += cuReal3(0, 0, cHmo);

		if (do_reduction) {

			//energy density
			int non_empty_cells = M1.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuBReal)MUB * M1[idx] * MU0 * cuReal3(0, 0, cHmo) / (non_empty_cells * M1.h.dim());
		}

		if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = cuReal3(0, 0, cHmo);
		if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = -(cuBReal)MUB * M1[idx] * MU0 * cuReal3(0, 0, cHmo) / M1.h.dim();
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//----------------------- UpdateField LAUNCHER

void Atom_MOpticalCUDA::UpdateField(void)
{
	if (paMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();

		MOpticalCUDA_UpdateField_Cubic <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, cuModule, true);
	}
	else MOpticalCUDA_UpdateField_Cubic <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, cuModule, false);
}

#endif

#endif