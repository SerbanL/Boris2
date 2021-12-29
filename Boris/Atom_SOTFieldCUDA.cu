#include "Atom_SOTFieldCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_SOTFIELD) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

#include "Atom_Mesh_CubicCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

__global__ void SOTFieldCUDA_ASC_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;
	cuVEC_VC<cuReal3>& E = *cuaMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Heff1.linear_size()) {

		cuReal3 SOTField;

		if (M1.is_not_empty(idx)) {

			cuBReal grel = *cuaMesh.pgrel;
			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal SHA = *cuaMesh.pSHA;
			cuBReal flSOT = *cuaMesh.pflSOT;
			cuReal3 STp = *cuaMesh.pSTp;
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pgrel, grel, *cuaMesh.pmu_s, mu_s, *cuaMesh.pSHA, SHA, *cuaMesh.pflSOT, flSOT, *cuaMesh.pSTp, STp);

			if (cuIsNZ(grel)) {

				cuBReal conv = M1.h.dim() / MUB;
				cuBReal a_const = -conv * (SHA * MUB_E / (GAMMA * grel)) / (mu_s * mu_s * (M1.rect.e.z - M1.rect.s.z));

				int idx_E = E.position_to_cellidx(M1.cellidx_to_position(idx));
				cuReal3 p_vec = (STp ^ E[idx_E]) * elC[idx_E];

				cuReal3 SOTField = a_const * ((M1[idx] ^ p_vec) + flSOT * mu_s * p_vec);
				Heff1[idx] += SOTField;
			}
		}

		if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = SOTField;
	}
}

//----------------------- UpdateField LAUNCHER

void Atom_SOTFieldCUDA::UpdateField(void)
{
	if (!paMeshCUDA->EComputation_Enabled()) return;

	SOTFieldCUDA_ASC_UpdateField <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, cuModule);
}

#endif

#endif