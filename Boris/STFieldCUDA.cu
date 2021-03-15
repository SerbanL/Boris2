#include "STFieldCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STFIELD

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

#include "Mesh_FerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

__global__ void STFieldCUDA_UpdateField(ManagedMeshCUDA& cuMesh, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC_VC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Heff.linear_size()) {

		cuReal3 STField;

		if (M.is_not_empty(idx)) {

			cuBReal grel = *cuMesh.pgrel;
			cuBReal Ms = *cuMesh.pMs;
			cuBReal flSOT = *cuMesh.pflSOT;
			cuReal2 STq = *cuMesh.pSTq;
			cuReal2 STa = *cuMesh.pSTa;
			cuReal3 STp = *cuMesh.pSTp;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pgrel, grel, *cuMesh.pMs, Ms, *cuMesh.pflSOT, flSOT, *cuMesh.pSTq, STq, *cuMesh.pSTa, STa, *cuMesh.pSTp, STp);

			if (cuIsNZ(grel)) {

				cuBReal dotprod = (M[idx] * STp) / Ms;
				cuBReal neta = STq.i / (STa.i + STq.j * dotprod) + STq.j / (STa.i - STq.j * dotprod);

				int idx_E = E.position_to_cellidx(M.cellidx_to_position(idx));
				//z component of Jc
				cuBReal Jc = (elC[idx_E] * E[idx_E].z);

				cuBReal a_const = -(neta * (cuBReal)MUB_E * Jc / ((cuBReal)GAMMA * grel)) / (Ms * Ms * (M.rect.e.z - M.rect.s.z));

				STField = a_const * ((M[idx] ^ STp) + flSOT * Ms * STp);
				Heff[idx] += STField;
			}
		}

		if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = STField;
	}
}

//----------------------- UpdateField LAUNCHER

void STFieldCUDA::UpdateField(void)
{
	if (!pMeshCUDA->EComputation_Enabled()) return;

	STFieldCUDA_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule);
}

#endif

#endif