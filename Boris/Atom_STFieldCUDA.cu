#include "Atom_STFieldCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_STFIELD) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "Atom_MeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

#include "ManagedMeshCUDA.h"
#include "MeshParamsControlCUDA.h"

#include "MeshDefs.h"

__global__ void Atom_STFieldCUDA_fixedpolarization_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;
	cuVEC_VC<cuReal3>& E = *cuaMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < Heff1.linear_size()) {

		cuReal3 ST_Field;

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

				ST_Field = a_const * ((M1[idx] ^ STp) + flSOT * mu_s * STp);
				Heff1[idx] += ST_Field;
			}
		}

		if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[idx] = ST_Field;
	}
}

//Top mesh is atomistic
__global__ void Atom_STFieldCUDA_TopAtom_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedAtom_MeshCUDA* paMesh_Top, size_t coupled_meshes, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;
	cuVEC_VC<cuReal3>& E = *cuaMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M1.n;
	cuReal3 h = M1.h;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

		//skip empty cells
		if (M1.is_not_empty(cell_idx)) {

			cuBReal grel = *cuaMesh.pgrel;
			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal flSOT2 = *cuaMesh.pflSOT2;
			cuReal2 STq2 = *cuaMesh.pSTq2;
			cuReal2 STa2 = *cuaMesh.pSTa2;
			cuaMesh.update_parameters_mcoarse(cell_idx, *cuaMesh.pgrel, grel, *cuaMesh.pmu_s, mu_s, *cuaMesh.pflSOT2, flSOT2, *cuaMesh.pSTq2, STq2, *cuaMesh.pSTa2, STa2);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M1_Top = *(paMesh_Top[mesh_idx].pM1);

				//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M1.rect.s.x - M1_Top.rect.s.x,
					(j + 0.5) * h.y + M1.rect.s.y - M1_Top.rect.s.y,
					M1_Top.h.z / 2);

				//can't couple to an empty cell
				if (!M1_Top.rect.contains(cell_rel_pos + M1_Top.rect.s) || M1_Top.is_empty(cell_rel_pos)) continue;

				//get magnetization value in top mesh cell to couple with
				cuReal3 p = M1_Top[cell_rel_pos].normalized();
				cuReal3 m = M1[cell_idx] / mu_s;

				cuBReal dotprod = m * p;
				cuBReal neta = STq2.i / (STa2.i + STa2.j * dotprod) + STq2.j / (STa2.i - STa2.j * dotprod);

				int idx_E = E.position_to_cellidx(M1.cellidx_to_position(cell_idx));
				//z component of J
				cuBReal Jc = E[idx_E].z * elC[idx_E];

				cuBReal conv = M1.h.dim() / MUB;
				cuBReal a_const = -conv * (neta * MUB_E * Jc / (GAMMA * grel)) / (mu_s * h.z);

				cuReal3 ST_Field = a_const * ((m ^ p) + flSOT2 * p);
				Heff1[cell_idx] += ST_Field;

				//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
				if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += ST_Field;

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}
}

//Top mesh is ferromagnetic
__global__ void Atom_STFieldCUDA_TopFM_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedMeshCUDA* pMesh_Top, size_t coupled_meshes, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;
	cuVEC_VC<cuReal3>& E = *cuaMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M1.n;
	cuReal3 h = M1.h;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

		//skip empty cells
		if (M1.is_not_empty(cell_idx)) {

			cuBReal grel = *cuaMesh.pgrel;
			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal flSOT2 = *cuaMesh.pflSOT2;
			cuReal2 STq2 = *cuaMesh.pSTq2;
			cuReal2 STa2 = *cuaMesh.pSTa2;
			cuaMesh.update_parameters_mcoarse(cell_idx, *cuaMesh.pgrel, grel, *cuaMesh.pmu_s, mu_s, *cuaMesh.pflSOT2, flSOT2, *cuaMesh.pSTq2, STq2, *cuaMesh.pSTa2, STa2);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M = *(pMesh_Top[mesh_idx].pM);

				cuRect tmeshRect = M.rect;

				//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M1.rect.s.x - tmeshRect.s.x,
					(j + 0.5) * h.y + M1.rect.s.y - tmeshRect.s.y,
					M.h.z / 2);

				//can't couple to an empty cell
				if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s) || M.is_empty(cell_rel_pos)) continue;

				//get magnetization value in top mesh cell to couple with
				cuReal3 p = M[cell_rel_pos].normalized();
				cuReal3 m = M1[cell_idx] / mu_s;

				cuBReal dotprod = m * p;
				cuBReal neta = STq2.i / (STa2.i + STa2.j * dotprod) + STq2.j / (STa2.i - STa2.j * dotprod);

				int idx_E = E.position_to_cellidx(M1.cellidx_to_position(cell_idx));
				//z component of J
				cuBReal Jc = E[idx_E].z * elC[idx_E];

				cuBReal conv = M1.h.dim() / MUB;
				cuBReal a_const = -conv * (neta * MUB_E * Jc / (GAMMA * grel)) / (mu_s * h.z);

				cuReal3 ST_Field = a_const * ((m ^ p) + flSOT2 * p);
				Heff1[cell_idx] += ST_Field;

				//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
				if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += ST_Field;

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}
}

//Bottom mesh is atomistic
__global__ void Atom_STFieldCUDA_BotAtom_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedAtom_MeshCUDA* paMesh_Bot, size_t coupled_meshes, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;
	cuVEC_VC<cuReal3>& E = *cuaMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M1.n;
	cuReal3 h = M1.h;

	cuBReal energy_ = 0.0;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x;

		//skip empty cells
		if (M1.is_not_empty(cell_idx)) {

			cuBReal grel = *cuaMesh.pgrel;
			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal flSOT = *cuaMesh.pflSOT;
			cuReal2 STq = *cuaMesh.pSTq;
			cuReal2 STa = *cuaMesh.pSTa;
			cuaMesh.update_parameters_mcoarse(cell_idx, *cuaMesh.pgrel, grel, *cuaMesh.pmu_s, mu_s, *cuaMesh.pflSOT, flSOT, *cuaMesh.pSTq, STq, *cuaMesh.pSTa, STa);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M1_Bot = *(paMesh_Bot[mesh_idx].pM1);

				//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M1.rect.s.x - M1_Bot.rect.s.x,
					(j + 0.5) * h.y + M1.rect.s.y - M1_Bot.rect.s.y,
					M1_Bot.rect.height() - M1_Bot.h.z / 2);

				//can't couple to an empty cell
				if (!M1_Bot.rect.contains(cell_rel_pos + M1_Bot.rect.s) || M1_Bot.is_empty(cell_rel_pos)) continue;

				//yes, then get value of magnetization used in coupling with current cell at cell_idx
				cuReal3 p = M1_Bot[cell_rel_pos].normalized();
				cuReal3 m = M1[cell_idx] / mu_s;

				cuBReal dotprod = m * p;
				cuBReal neta = STq.i / (STa.i + STa.j * dotprod) + STq.j / (STa.i - STa.j * dotprod);

				int idx_E = E.position_to_cellidx(M1.cellidx_to_position(cell_idx));
				//z component of J
				cuBReal Jc = E[idx_E].z * elC[idx_E];

				cuBReal conv = M1.h.dim() / MUB;
				cuBReal a_const = -conv * (neta * MUB_E * Jc / (GAMMA * grel)) / (mu_s * h.z);

				cuReal3 ST_Field = a_const * ((m ^ p) + flSOT * p);
				Heff1[cell_idx] += ST_Field;

				//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
				if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += ST_Field;

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}
}

//Bottom mesh is ferromagnetic
__global__ void Atom_STFieldCUDA_BotFM_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedMeshCUDA* pMesh_Bot, size_t coupled_meshes, ManagedModulesCUDA& cuModule)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;
	cuVEC_VC<cuReal3>& E = *cuaMesh.pE;
	cuVEC_VC<cuBReal>& elC = *cuaMesh.pelC;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M1.n;
	cuReal3 h = M1.h;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x;

		//skip empty cells
		if (M1.is_not_empty(cell_idx)) {

			cuBReal grel = *cuaMesh.pgrel;
			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal flSOT = *cuaMesh.pflSOT;
			cuReal2 STq = *cuaMesh.pSTq;
			cuReal2 STa = *cuaMesh.pSTa;
			cuaMesh.update_parameters_mcoarse(cell_idx, *cuaMesh.pgrel, grel, *cuaMesh.pmu_s, mu_s, *cuaMesh.pflSOT, flSOT, *cuaMesh.pSTq, STq, *cuaMesh.pSTa, STa);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M = *pMesh_Bot[mesh_idx].pM;

				cuRect bmeshRect = M.rect;

				//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M1.rect.s.x - bmeshRect.s.x,
					(j + 0.5) * h.y + M1.rect.s.y - bmeshRect.s.y,
					bmeshRect.height() - M.h.z / 2);

				//can't couple to an empty cell
				if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s) || M.is_empty(cell_rel_pos)) continue;

				//get magnetization value in top mesh cell to couple with
				cuReal3 p = M[cell_rel_pos].normalized();
				cuReal3 m = M1[cell_idx] / mu_s;

				cuBReal dotprod = m * p;
				cuBReal neta = STq.i / (STa.i + STa.j * dotprod) + STq.j / (STa.i - STa.j * dotprod);

				int idx_E = E.position_to_cellidx(M1.cellidx_to_position(cell_idx));
				//z component of J
				cuBReal Jc = E[idx_E].z * elC[idx_E];

				cuBReal conv = M1.h.dim() / MUB;
				cuBReal a_const = -conv * (neta * MUB_E * Jc / (GAMMA * grel)) / (mu_s * h.z);

				cuReal3 ST_Field = a_const * ((m ^ p) + flSOT * p);
				Heff1[cell_idx] += ST_Field;

				//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
				if (cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += ST_Field;

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}
}

//----------------------- UpdateField LAUNCHER

void Atom_STFieldCUDA::UpdateField(void)
{
	if (!paMeshCUDA->EComputation_Enabled()) return;

	if (fixed_polarization) {

		////////////////////////////////////////
		//      FIXED POLARIZATION VERSION    //
		////////////////////////////////////////

		Atom_STFieldCUDA_fixedpolarization_UpdateField <<< (paMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (paMeshCUDA->cuaMesh, cuModule);
	}
	else {

		/////////////////////////////////////////////
		// POLARIZATION FROM TOP AND BOTTOM MESHES //
		/////////////////////////////////////////////

		//Coupling from ferromagnetic meshes

		//Top
		if (pMeshFM_Top.size()) {

			Atom_STFieldCUDA_TopFM_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, pMeshFM_Top, pMeshFM_Top.size(), cuModule);
		}

		//Bottom
		if (pMeshFM_Bot.size()) {

			Atom_STFieldCUDA_BotFM_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, pMeshFM_Bot, pMeshFM_Bot.size(), cuModule);
		}

		//Coupling from atomistic meshes

		//Top
		if (paMesh_Top.size()) {

			Atom_STFieldCUDA_TopAtom_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, paMesh_Top, paMesh_Top.size(), cuModule);
		}

		//Bottom
		if (paMesh_Bot.size()) {

			Atom_STFieldCUDA_BotAtom_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, paMesh_Bot, paMesh_Bot.size(), cuModule);
		}
	}
}

#endif

#endif