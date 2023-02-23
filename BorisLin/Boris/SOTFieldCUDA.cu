#include "SOTFieldCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SOTFIELD

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

#include "Mesh_FerromagneticCUDA.h"
#include "Mesh_AntiFerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

__global__ void SOTFieldCUDA_FM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC_VC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Heff.linear_size()) {

		cuReal3 SOTField;

		if (M.is_not_empty(idx)) {

			cuBReal grel = *cuMesh.pgrel;
			cuBReal Ms = *cuMesh.pMs;
			cuBReal SHA = *cuMesh.pSHA;
			cuBReal flSOT = *cuMesh.pflSOT;
			cuReal3 STp = *cuMesh.pSTp;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pgrel, grel, *cuMesh.pMs, Ms, *cuMesh.pSHA, SHA, *cuMesh.pflSOT, flSOT, *cuMesh.pSTp, STp);

			if (cuIsNZ(grel)) {

				cuBReal a_const = -(SHA * MUB_E / (GAMMA * grel)) / (Ms * Ms * (M.rect.e.z - M.rect.s.z));

				int idx_E = E.position_to_cellidx(M.cellidx_to_position(idx));
				cuReal3 p_vec = STp ^ (elC[idx_E] * E[idx_E]);

				SOTField = a_const * ((M[idx] ^ p_vec) + flSOT * Ms * p_vec);
				Heff[idx] += SOTField;
			}
		}

		if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = SOTField;
	}
}

__global__ void SOTFieldCUDA_AFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;
	cuVEC_VC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Heff.linear_size()) {

		cuReal3 SOTField_A, SOTField_B;

		if (M.is_not_empty(idx)) {

			cuReal2 grel_AFM = *cuMesh.pgrel_AFM;
			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuBReal SHA = *cuMesh.pSHA;
			cuBReal flSOT = *cuMesh.pflSOT;
			cuReal3 STp = *cuMesh.pSTp;
			cuMesh.update_parameters_mcoarse(idx, *cuMesh.pgrel_AFM, grel_AFM, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pSHA, SHA, *cuMesh.pflSOT, flSOT, *cuMesh.pSTp, STp);

			if (cuIsNZ(grel_AFM.i + grel_AFM.j)) {

				cuBReal a_const_A = -(SHA * (cuBReal)MUB_E / ((cuBReal)GAMMA * grel_AFM.i)) / (Ms_AFM.i * Ms_AFM.i * (M.rect.e.z - M.rect.s.z));
				cuBReal a_const_B = -(SHA * (cuBReal)MUB_E / ((cuBReal)GAMMA * grel_AFM.j)) / (Ms_AFM.j * Ms_AFM.j * (M.rect.e.z - M.rect.s.z));

				int idx_E = E.position_to_cellidx(M.cellidx_to_position(idx));
				cuReal3 p_vec = STp ^ (elC[idx_E] * E[idx_E]);

				SOTField_A = a_const_A * ((M[idx] ^ p_vec) + flSOT * Ms_AFM.i * p_vec);
				SOTField_B = a_const_B * ((M2[idx] ^ p_vec) + flSOT * Ms_AFM.j * p_vec);

				Heff[idx] += SOTField_A;
				Heff2[idx] += SOTField_B;
			}
		}

		if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = SOTField_A;
		if (cuModule.pModule_Heff2->linear_size()) (*cuModule.pModule_Heff2)[idx] = SOTField_B;
	}
}

//----------------------- UpdateField LAUNCHER

void SOTFieldCUDA::UpdateField(void)
{
	if (!pMeshCUDA->EComputation_Enabled()) return;

	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		SOTFieldCUDA_AFM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule);
	}
	else {

		SOTFieldCUDA_FM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule);
	}
}

#endif

#endif