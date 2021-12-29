#include "Atom_STFieldCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_STFIELD) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

#include "Atom_Mesh_CubicCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

__global__ void Atom_STFieldCUDA_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;
	cuVEC_VC<cuReal3>& E = *cuaMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Heff1.linear_size()) {

		cuReal3 STField;

		if (M1.is_not_empty(idx)) {

			cuBReal grel = *cuaMesh.pgrel;
			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal flSOT = *cuaMesh.pflSOT;
			cuReal2 STq = *cuaMesh.pSTq;
			cuReal2 STa = *cuaMesh.pSTa;
			cuReal3 STp = *cuaMesh.pSTp;
			cuaMesh.update_parameters_mcoarse(idx, *cuaMesh.pgrel, grel, *cuaMesh.pmu_s, mu_s, *cuaMesh.pflSOT, flSOT, *cuaMesh.pSTq, STq, *cuaMesh.pSTa, STa, *cuaMesh.pSTp, STp);

			if (cuIsNZ(grel)) {

				cuBReal dotprod = (M1[idx] * STp) / mu_s;
				cuBReal neta = STq.i / (STa.i + STa.j * dotprod) + STq.j / (STa.i - STa.j * dotprod);

				int idx_E = E.position_to_cellidx(M1.cellidx_to_position(idx));
				//z component of Jc
				cuBReal Jc = (elC[idx_E] * E[idx_E].z);

				cuBReal conv = M1.h.dim() / MUB;
				cuBReal a_const = -conv * (neta * (cuBReal)MUB_E * Jc / ((cuBReal)GAMMA * grel)) / (mu_s * mu_s * (M1.rect.e.z - M1.rect.s.z));

				STField = a_const * ((M1[idx] ^ STp) + flSOT * mu_s * STp);
				Heff1[idx] += STField;
			}
		}

		if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = STField;
	}
}

//----------------------- UpdateField LAUNCHER

void Atom_STFieldCUDA::UpdateField(void)
{
	if (!paMeshCUDA->EComputation_Enabled()) return;

	Atom_STFieldCUDA_UpdateField <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, cuModule);
}

#endif

#endif