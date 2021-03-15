#include "Atom_Demag_NCUDA.h"
#include "Atom_MeshCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_DEMAG_N) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

__global__ void Demag_NCUDA_Cubic_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;
	MatPCUDA<cuReal2, cuBReal>& Nxy = *cuaMesh.pNxy;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal energy_ = 0.0;

	if (idx < Heff1.linear_size()) {

		//used to convert moment to magnetization in each atomistic unit cell
		cuBReal conversion = (cuBReal)MUB / M1.h.dim();

		cuReal3 Heff_value = cuReal3();

		if (M1.is_not_empty(idx)) {

			Heff_value = cuReal3(-cuReal2(Nxy).x * M1[idx].x, -cuReal2(Nxy).y * M1[idx].y, -(1 - cuReal2(Nxy).x - cuReal2(Nxy).y) * M1[idx].z) * conversion;

			if (do_reduction) {

				int non_empty_cells = M1.get_nonempty_cells();
				if (non_empty_cells) energy_ = -(cuBReal)MU0 * conversion * M1[idx] * Heff_value / (2 * non_empty_cells);
			}

			if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = Heff_value;
			if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[idx] = -(cuBReal)MUB_MU0 * M1[idx] * Heff_value / (2 * M1.h.dim());
		}

		Heff1[idx] += Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//----------------------- UpdateField LAUNCHER

void Atom_Demag_NCUDA::UpdateField(void)
{
	if (paMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();

		Demag_NCUDA_Cubic_UpdateField <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, cuModule, true);
	}
	else {

		Demag_NCUDA_Cubic_UpdateField <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, cuModule, false);
	}
}

#endif

#endif