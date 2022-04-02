#include "Atom_SurfExchangeCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_SURFEXCHANGE) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "Atom_MeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

#include "ManagedMeshCUDA.h"
#include "MeshParamsControlCUDA.h"

#include "MeshDefs.h"

//Top mesh is atomistic
__global__ void SurfExchangeCUDA_Top_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedAtom_MeshCUDA* paMesh_Top, size_t coupled_meshes, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M1.n;
	cuReal3 h = M1.h;

	cuBReal energy_ = 0.0;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

		//skip empty cells
		if (M1.is_not_empty(cell_idx)) {

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuaMesh.update_parameters_mcoarse(cell_idx, *cuaMesh.pmu_s, mu_s);

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

				cuBReal Js = *(paMesh_Top[mesh_idx].pJs);
				paMesh_Top[mesh_idx].update_parameters_atposition(cell_rel_pos, *(paMesh_Top[mesh_idx].pJs), Js);

				//get magnetization value in top mesh cell to couple with
				cuReal3 m_j = M1_Top[cell_rel_pos].normalized();
				cuReal3 m_i = M1[cell_idx] / mu_s;

				cuBReal dot_prod = m_i * m_j;

				cuReal3 Hsurfexch = m_j * Js / ((cuBReal)MUB_MU0 * mu_s);

				if (do_reduction) {

					energy_ = -Js * dot_prod / (M1.get_nonempty_cells() * h.dim());
				}

				Heff1[cell_idx] += Hsurfexch;

				//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
				if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexch;
				if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy_ * M1.get_nonempty_cells();

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Top mesh is ferromagnetic
__global__ void SurfExchangeCUDA_TopFM_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedMeshCUDA* pMesh_Top, size_t coupled_meshes, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M1.n;
	cuReal3 h = M1.h;

	cuBReal energy_ = 0.0;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

		//skip empty cells
		if (M1.is_not_empty(cell_idx)) {

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal Js = *cuaMesh.pJs;
			cuaMesh.update_parameters_mcoarse(cell_idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pJs, Js);

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
				cuReal3 m_j = M[cell_rel_pos].normalized();
				cuReal3 m_i = M1[cell_idx] / mu_s;

				cuBReal dot_prod = m_i * m_j;

				cuReal3 Hsurfexch = m_j * Js / (MUB_MU0 * mu_s);

				if (do_reduction) {

					energy_ = -Js * dot_prod / (M1.get_nonempty_cells() * h.dim());
				}

				Heff1[cell_idx] += Hsurfexch;

				//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
				if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexch;
				if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy_ * M1.get_nonempty_cells();

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Top mesh is antiferromagnetic
__global__ void SurfExchangeCUDA_TopAFM_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedMeshCUDA* pMesh_Top, size_t coupled_meshes, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M1.n;
	cuReal3 h = M1.h;

	cuBReal energy_ = 0.0;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

		//skip empty cells
		if (M1.is_not_empty(cell_idx)) {

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal Js = *cuaMesh.pJs;
			cuBReal Js2 = *cuaMesh.pJs2;
			cuaMesh.update_parameters_mcoarse(cell_idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pJs, Js, *cuaMesh.pJs2, Js2);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M = *(pMesh_Top[mesh_idx].pM);
				cuVEC_VC<cuReal3>& M2 = *(pMesh_Top[mesh_idx].pM2);

				cuRect tmeshRect = M.rect;

				//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M1.rect.s.x - tmeshRect.s.x,
					(j + 0.5) * h.y + M1.rect.s.y - tmeshRect.s.y,
					M.h.z / 2);

				//can't couple to an empty cell
				if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s) || M.is_empty(cell_rel_pos)) continue;

				//get magnetization value in top mesh cell to couple with
				cuReal3 m_j1 = M[cell_rel_pos].normalized();
				cuReal3 m_j2 = M2[cell_rel_pos].normalized();
				cuReal3 m_i = M1[cell_idx] / mu_s;

				cuBReal dot_prod1 = m_i * m_j1;
				cuBReal dot_prod2 = m_i * m_j2;

				cuReal3 Hsurfexch = m_j1 * Js / (MUB_MU0 * mu_s);
				Hsurfexch += m_j2 * Js2 / (MUB_MU0 * mu_s);

				if (do_reduction) {

					energy_ = -Js * dot_prod1 / (M1.get_nonempty_cells() * h.dim());
					energy_ += -Js2 * dot_prod2 / (M1.get_nonempty_cells() * h.dim());
				}

				Heff1[cell_idx] += Hsurfexch;

				//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
				if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexch;
				if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy_ * M1.get_nonempty_cells();

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Bottom mesh is atomistic
__global__ void SurfExchangeCUDA_Bot_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedAtom_MeshCUDA* paMesh_Bot, size_t coupled_meshes, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

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

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal Js = *cuaMesh.pJs;
			cuaMesh.update_parameters_mcoarse(cell_idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pJs, Js);

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
				cuReal3 m_j = M1_Bot[cell_rel_pos].normalized();
				cuReal3 m_i = M1[cell_idx] / mu_s;

				cuBReal dot_prod = m_i * m_j;

				cuReal3 Hsurfexch = m_j * Js / ((cuBReal)MUB_MU0 * mu_s);

				if (do_reduction) {

					energy_ = -Js * dot_prod / (M1.get_nonempty_cells()* h.dim());
				}

				Heff1[cell_idx] += Hsurfexch;

				//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
				if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexch;
				if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy_ * M1.get_nonempty_cells();

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Bottom mesh is ferromagnetic
__global__ void SurfExchangeCUDA_BotFM_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedMeshCUDA* pMesh_Bot, size_t coupled_meshes, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

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

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal Js = *cuaMesh.pJs;
			cuaMesh.update_parameters_mcoarse(cell_idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pJs, Js);

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
				cuReal3 m_j = M[cell_rel_pos].normalized();
				cuReal3 m_i = M1[cell_idx] / mu_s;

				cuBReal dot_prod = m_i * m_j;

				cuReal3 Hsurfexch = m_j * Js / (MUB_MU0 * mu_s);

				if (do_reduction) {

					energy_ = -Js * dot_prod / (M1.get_nonempty_cells()* h.dim());
				}

				Heff1[cell_idx] += Hsurfexch;

				//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
				if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexch;
				if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy_ * M1.get_nonempty_cells();

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Bottom mesh is antiferromagnetic
__global__ void SurfExchangeCUDA_BotAFM_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedMeshCUDA* pMesh_Bot, size_t coupled_meshes, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;
	cuVEC<cuReal3>& Heff1 = *cuaMesh.pHeff1;

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

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal Js = *cuaMesh.pJs;
			cuBReal Js2 = *cuaMesh.pJs2;
			cuaMesh.update_parameters_mcoarse(cell_idx, *cuaMesh.pmu_s, mu_s, *cuaMesh.pJs, Js, *cuaMesh.pJs2, Js2);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M = *pMesh_Bot[mesh_idx].pM;
				cuVEC_VC<cuReal3>& M2 = *pMesh_Bot[mesh_idx].pM2;

				cuRect bmeshRect = M.rect;

				//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M1.rect.s.x - bmeshRect.s.x,
					(j + 0.5) * h.y + M1.rect.s.y - bmeshRect.s.y,
					bmeshRect.height() - M.h.z / 2);

				//can't couple to an empty cell
				if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s) || M.is_empty(cell_rel_pos)) continue;

				//get magnetization value in top mesh cell to couple with
				cuReal3 m_j1 = M[cell_rel_pos].normalized();
				cuReal3 m_j2 = M2[cell_rel_pos].normalized();
				cuReal3 m_i = M1[cell_idx] / mu_s;

				cuBReal dot_prod1 = m_i * m_j1;
				cuBReal dot_prod2 = m_i * m_j2;

				cuReal3 Hsurfexch = m_j1 * Js / (MUB_MU0 * mu_s);
				Hsurfexch += m_j2 * Js2 / (MUB_MU0 * mu_s);

				if (do_reduction) {

					energy_ = -Js * dot_prod1 / (M1.get_nonempty_cells()* h.dim());
					energy_ += -Js2 * dot_prod2 / (M1.get_nonempty_cells()* h.dim());
				}

				Heff1[cell_idx] += Hsurfexch;

				//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
				if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexch;
				if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy_ * M1.get_nonempty_cells();

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//----------------------- UpdateField LAUNCHER

void Atom_SurfExchangeCUDA::UpdateField(void)
{
	if (paMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();
		ZeroModuleVECs();

		//Top - Atomistic
		if (paMesh_Top.size()) {

			SurfExchangeCUDA_Top_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, paMesh_Top, paMesh_Top.size(), cuModule, true);
		}

		//Top - Ferromagnetic
		if (pMeshFM_Top.size()) {

			SurfExchangeCUDA_TopFM_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, pMeshFM_Top, pMeshFM_Top.size(), cuModule, true);
		}

		//Top - AntiFerromagnetic
		if (pMeshAFM_Top.size()) {

			SurfExchangeCUDA_TopAFM_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, pMeshAFM_Top, pMeshAFM_Top.size(), cuModule, true);
		}

		//Bottom - Atomistic
		if (paMesh_Bot.size()) {

			SurfExchangeCUDA_Bot_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, paMesh_Bot, paMesh_Bot.size(), cuModule, true);
		}

		//Bottom - Ferromagnetic
		if (pMeshFM_Bot.size()) {

			SurfExchangeCUDA_BotFM_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, pMeshFM_Bot, pMeshFM_Bot.size(), cuModule, true);
		}

		//Bottom - AntiFerromagnetic
		if (pMeshAFM_Bot.size()) {

			SurfExchangeCUDA_BotAFM_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, pMeshAFM_Bot, pMeshAFM_Bot.size(), cuModule, true);
		}
	}
	else {

		//Top - Atomistic
		if (paMesh_Top.size()) {

			SurfExchangeCUDA_Top_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, paMesh_Top, paMesh_Top.size(), cuModule, false);
		}

		//Top - Ferromagnetic
		if (pMeshFM_Top.size()) {

			SurfExchangeCUDA_TopFM_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, pMeshFM_Top, pMeshFM_Top.size(), cuModule, false);
		}

		//Top - AntiFerromagnetic
		if (pMeshAFM_Top.size()) {

			SurfExchangeCUDA_TopAFM_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, pMeshAFM_Top, pMeshAFM_Top.size(), cuModule, false);
		}

		//Bottom - Atomistic
		if (paMesh_Bot.size()) {

			SurfExchangeCUDA_Bot_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, paMesh_Bot, paMesh_Bot.size(), cuModule, false);
		}

		//Bottom - Ferromagnetic
		if (pMeshFM_Bot.size()) {

			SurfExchangeCUDA_BotFM_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, pMeshFM_Bot, pMeshFM_Bot.size(), cuModule, false);
		}

		//Bottom - AntiFerromagnetic
		if (pMeshAFM_Bot.size()) {

			SurfExchangeCUDA_BotAFM_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, pMeshAFM_Bot, pMeshAFM_Bot.size(), cuModule, false);
		}
	}
}

//----------------------- Initialization

//Top mesh is ferromagnetic
__global__ void set_Atom_SurfExchangeCUDA_pointers_kernel(
	ManagedAtom_MeshCUDA& cuaMesh,
	ManagedAtom_MeshCUDA* paMesh_Bot, size_t coupled_bot_meshes,
	ManagedAtom_MeshCUDA* paMesh_Top, size_t coupled_top_meshes,
	ManagedMeshCUDA* pMeshFM_Bot, size_t coupledFM_bot_meshes,
	ManagedMeshCUDA* pMeshFM_Top, size_t coupledFM_top_meshes,
	ManagedMeshCUDA* pMeshAFM_Bot, size_t coupledAFM_bot_meshes,
	ManagedMeshCUDA* pMeshAFM_Top, size_t coupledAFM_top_meshes)
{
	if (threadIdx.x == 0) cuaMesh.paMesh_Bot = paMesh_Bot;
	if (threadIdx.x == 1) cuaMesh.paMesh_Bot_size = coupled_bot_meshes;
	if (threadIdx.x == 2) cuaMesh.paMesh_Top = paMesh_Top;
	if (threadIdx.x == 3) cuaMesh.paMesh_Top_size = coupled_top_meshes;

	if (threadIdx.x == 4) cuaMesh.pMeshFM_Bot = pMeshFM_Bot;
	if (threadIdx.x == 5) cuaMesh.pMeshFM_Bot_size = coupledFM_bot_meshes;
	if (threadIdx.x == 6) cuaMesh.pMeshFM_Top = pMeshFM_Top;
	if (threadIdx.x == 7) cuaMesh.pMeshFM_Top_size = coupledFM_top_meshes;

	if (threadIdx.x == 8) cuaMesh.pMeshAFM_Bot = pMeshAFM_Bot;
	if (threadIdx.x == 9) cuaMesh.pMeshAFM_Bot_size = coupledAFM_bot_meshes;
	if (threadIdx.x == 10) cuaMesh.pMeshAFM_Top = pMeshAFM_Top;
	if (threadIdx.x == 11) cuaMesh.pMeshAFM_Top_size = coupledAFM_top_meshes;
}

//Called by Atom_SurfExchangeCUDA module
void Atom_SurfExchangeCUDA::set_Atom_SurfExchangeCUDA_pointers(void)
{
	set_Atom_SurfExchangeCUDA_pointers_kernel <<< 1, CUDATHREADS >>>
		(paMeshCUDA->cuaMesh, 
		paMesh_Bot, paMesh_Bot.size(), paMesh_Top, paMesh_Top.size(),
		pMeshFM_Bot, pMeshFM_Bot.size(), pMeshFM_Top, pMeshFM_Top.size(),
		pMeshAFM_Bot, pMeshAFM_Bot.size(), pMeshAFM_Top, pMeshAFM_Top.size());
}

#endif

#endif