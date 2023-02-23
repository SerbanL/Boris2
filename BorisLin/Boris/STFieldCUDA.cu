#include "STFieldCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STFIELD

#include "BorisCUDALib.cuh"

#include "MeshDefs.h"

#include "Mesh_FerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

#include "ManagedAtom_MeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

////////////////////////////////////////
//      FIXED POLARIZATION VERSION    //
////////////////////////////////////////

__global__ void STFieldCUDA_fixedpolarization_UpdateField(ManagedMeshCUDA& cuMesh, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC_VC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Heff.linear_size()) {

		cuReal3 ST_Field;

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
				cuBReal neta = STq.i / (STa.i + STa.j * dotprod) + STq.j / (STa.i - STa.j * dotprod);

				int idx_E = E.position_to_cellidx(M.cellidx_to_position(idx));
				//z component of Jc
				cuBReal Jc = (elC[idx_E] * E[idx_E].z);

				cuBReal a_const = -(neta * (cuBReal)MUB_E * Jc / ((cuBReal)GAMMA * grel)) / (Ms * Ms * (M.rect.e.z - M.rect.s.z));

				ST_Field = a_const * ((M[idx] ^ STp) + flSOT * Ms * STp);
				Heff[idx] += ST_Field;
			}
		}

		if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = ST_Field;
	}
}

/////////////////////////////////////////////
// POLARIZATION FROM TOP AND BOTTOM MESHES //
/////////////////////////////////////////////

//Top mesh is ferromagnetic
__global__ void STFieldCUDA_TopFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Top, size_t coupled_meshes, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC_VC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

		//skip empty cells
		if (M.is_not_empty(cell_idx)) {

			cuBReal grel = *cuMesh.pgrel;
			cuBReal Ms = *cuMesh.pMs;
			cuBReal flSOT2 = *cuMesh.pflSOT2;
			cuReal2 STq2 = *cuMesh.pSTq2;
			cuReal2 STa2 = *cuMesh.pSTa2;
			cuMesh.update_parameters_mcoarse(cell_idx, *cuMesh.pgrel, grel, *cuMesh.pMs, Ms, *cuMesh.pflSOT2, flSOT2, *cuMesh.pSTq2, STq2, *cuMesh.pSTa2, STa2);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M_Top = *(pMesh_Top[mesh_idx].pM);

				//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M.rect.s.x - M_Top.rect.s.x,
					(j + 0.5) * h.y + M.rect.s.y - M_Top.rect.s.y,
					M_Top.h.z / 2);

				//can't couple to an empty cell
				if (!M_Top.rect.contains(cell_rel_pos + M_Top.rect.s) || M_Top.is_empty(cell_rel_pos)) continue;

				//get magnetization value in top mesh cell to couple with
				cuReal3 p = cu_normalize(M_Top[cell_rel_pos]);
				cuReal3 m = cu_normalize(M[cell_idx]);

				cuBReal dotprod = m * p;
				cuBReal neta = STq2.i / (STa2.i + STa2.j * dotprod) + STq2.j / (STa2.i - STa2.j * dotprod);

				int idx_E = E.position_to_cellidx(M.cellidx_to_position(cell_idx));
				//z component of J
				cuBReal Jc = E[idx_E].z * elC[idx_E];

				cuBReal a_const = -(neta * (cuBReal)MUB_E * Jc / ((cuBReal)GAMMA * grel)) / (Ms * h.z);

				cuReal3 ST_Field = a_const * ((m ^ p) + flSOT2 * p);
				Heff[cell_idx] += ST_Field;

				//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
				if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += ST_Field;

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}
}

//Top mesh is atomistic
__global__ void STFieldCUDA_TopAtom_UpdateField(ManagedMeshCUDA& cuMesh, ManagedAtom_MeshCUDA* pMesh_Top, size_t coupled_meshes, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC_VC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

		//skip empty cells
		if (M.is_not_empty(cell_idx)) {

			cuBReal grel = *cuMesh.pgrel;
			cuBReal Ms = *cuMesh.pMs;
			cuBReal flSOT2 = *cuMesh.pflSOT2;
			cuReal2 STq2 = *cuMesh.pSTq2;
			cuReal2 STa2 = *cuMesh.pSTa2;
			cuMesh.update_parameters_mcoarse(cell_idx, *cuMesh.pgrel, grel, *cuMesh.pMs, Ms, *cuMesh.pflSOT2, flSOT2, *cuMesh.pSTq2, STq2, *cuMesh.pSTa2, STa2);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M1 = *(pMesh_Top[mesh_idx].pM1);

				//coupling rectangle in atomistic mesh in absolute coordinates
				cuRect rect_c = cuRect(
					cuReal3(i * h.x, j * h.y, M.rect.e.z),
					cuReal3((i + 1) * h.x, (j + 1) * h.y, M1.h.z + M.rect.e.z));
				rect_c += cuReal3(M.rect.s.x, M.rect.s.y, 0.0);

				//get magnetization value in top mesh cell (the polarization)
				cuReal3 p = cu_normalize(M1.average(rect_c - M1.rect.get_s()));
				cuReal3 m = cu_normalize(M[cell_idx]);

				cuBReal dotprod = m * p;
				cuBReal neta = STq2.i / (STa2.i + STa2.j * dotprod) + STq2.j / (STa2.i - STa2.j * dotprod);

				int idx_E = E.position_to_cellidx(M.cellidx_to_position(cell_idx));
				//z component of J
				cuBReal Jc = E[idx_E].z * elC[idx_E];

				cuBReal a_const = -(neta * (cuBReal)MUB_E * Jc / ((cuBReal)GAMMA * grel)) / (Ms * h.z);

				cuReal3 ST_Field = a_const * ((m ^ p) + flSOT2 * p);

				Heff[cell_idx] += ST_Field;

				//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
				if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += ST_Field;

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}
}

//Bottom mesh is ferromagnetic
__global__ void STFieldCUDA_BotFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Bot, size_t coupled_meshes, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC_VC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x;

		//skip empty cells
		if (M.is_not_empty(cell_idx)) {

			cuBReal grel = *cuMesh.pgrel;
			cuBReal Ms = *cuMesh.pMs;
			cuBReal flSOT = *cuMesh.pflSOT;
			cuReal2 STq = *cuMesh.pSTq;
			cuReal2 STa = *cuMesh.pSTa;
			cuMesh.update_parameters_mcoarse(cell_idx, *cuMesh.pgrel, grel, *cuMesh.pMs, Ms, *cuMesh.pflSOT, flSOT, *cuMesh.pSTq, STq, *cuMesh.pSTa, STa);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M_Bot = *(pMesh_Bot[mesh_idx].pM);

				//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M.rect.s.x - M_Bot.rect.s.x,
					(j + 0.5) * h.y + M.rect.s.y - M_Bot.rect.s.y,
					M_Bot.rect.height() - M_Bot.h.z / 2);

				//can't couple to an empty cell
				if (!M_Bot.rect.contains(cell_rel_pos + M_Bot.rect.s) || M_Bot.is_empty(cell_rel_pos)) continue;

				//yes, then get value of magnetization used in coupling with current cell at cell_idx
				cuReal3 p = cu_normalize(M_Bot[cell_rel_pos]);
				cuReal3 m = cu_normalize(M[cell_idx]);

				cuBReal dotprod = m * p;
				cuBReal neta = STq.i / (STa.i + STa.j * dotprod) + STq.j / (STa.i - STa.j * dotprod);

				int idx_E = E.position_to_cellidx(M.cellidx_to_position(cell_idx));
				//z component of J
				cuBReal Jc = E[idx_E].z * elC[idx_E];

				cuBReal a_const = -(neta * (cuBReal)MUB_E * Jc / ((cuBReal)GAMMA * grel)) / (Ms * h.z);

				cuReal3 ST_Field = a_const * ((m ^ p) + flSOT * p);

				Heff[cell_idx] += ST_Field;

				if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += ST_Field;

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}
}

//Bottom mesh is atomistic
__global__ void STFieldCUDA_BotAtom_UpdateField(ManagedMeshCUDA& cuMesh, ManagedAtom_MeshCUDA* pMesh_Bot, size_t coupled_meshes, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC_VC<cuReal3>& E = *cuMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x;

		//skip empty cells
		if (M.is_not_empty(cell_idx)) {

			cuBReal grel = *cuMesh.pgrel;
			cuBReal Ms = *cuMesh.pMs;
			cuBReal flSOT = *cuMesh.pflSOT;
			cuReal2 STq = *cuMesh.pSTq;
			cuReal2 STa = *cuMesh.pSTa;
			cuMesh.update_parameters_mcoarse(cell_idx, *cuMesh.pgrel, grel, *cuMesh.pMs, Ms, *cuMesh.pflSOT, flSOT, *cuMesh.pSTq, STq, *cuMesh.pSTa, STa);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M1 = *pMesh_Bot[mesh_idx].pM1;

				//coupling rectangle in atomistic mesh in absolute coordinates
				cuRect rect_c = cuRect(
					cuReal3(i * h.x, j * h.y, M1.rect.e.z - M1.h.z),
					cuReal3((i + 1) * h.x, (j + 1) * h.y, M1.rect.e.z));
				rect_c += cuReal3(M.rect.s.x, M.rect.s.y, 0.0);

				//get magnetization value in top mesh cell (the polarization)
				cuReal3 p = cu_normalize(M1.average(rect_c - M1.rect.get_s()));
				cuReal3 m = cu_normalize(M[cell_idx]);

				cuBReal dotprod = m * p;
				cuBReal neta = STq.i / (STa.i + STa.j * dotprod) + STq.j / (STa.i - STa.j * dotprod);

				int idx_E = E.position_to_cellidx(M.cellidx_to_position(cell_idx));
				//z component of J
				cuBReal Jc = E[idx_E].z * elC[idx_E];

				cuBReal a_const = -(neta * (cuBReal)MUB_E * Jc / ((cuBReal)GAMMA * grel)) / (Ms * h.z);

				cuReal3 ST_Field = a_const * ((m ^ p) + flSOT * p);

				Heff[cell_idx] += ST_Field;

				if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += ST_Field;

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}
}

//----------------------- UpdateField LAUNCHER

void STFieldCUDA::UpdateField(void)
{
	if (!pMeshCUDA->EComputation_Enabled()) return;

	if (fixed_polarization) {

		////////////////////////////////////////
		//      FIXED POLARIZATION VERSION    //
		////////////////////////////////////////

		STFieldCUDA_fixedpolarization_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (pMeshCUDA->cuMesh, cuModule);
	}
	else {

		/////////////////////////////////////////////
		// POLARIZATION FROM TOP AND BOTTOM MESHES //
		/////////////////////////////////////////////

		//Coupling from ferromagnetic meshes

		//Top
		if (pMeshFM_Top.size()) {

			STFieldCUDA_TopFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshFM_Top, pMeshFM_Top.size(), cuModule);
		}

		//Bottom
		if (pMeshFM_Bot.size()) {

			STFieldCUDA_BotFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshFM_Bot, pMeshFM_Bot.size(), cuModule);
		}

		//Coupling from atomistic meshes

		//Top
		if (pMeshAtom_Top.size()) {

			STFieldCUDA_TopAtom_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAtom_Top, pMeshAtom_Top.size(), cuModule);
		}

		//Bottom
		if (pMeshAtom_Bot.size()) {

			STFieldCUDA_BotAtom_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAtom_Bot, pMeshAtom_Bot.size(), cuModule);
		}
	}
}

#endif

#endif