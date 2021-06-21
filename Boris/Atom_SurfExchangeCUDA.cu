#include "Atom_SurfExchangeCUDA.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_SURFEXCHANGE) && ATOMISTIC == 1

#include "BorisCUDALib.cuh"

#include "Atom_MeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"
#include "MeshDefs.h"

//Top mesh is ferromagnetic
__global__ void SurfExchangeCUDA_Top_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedAtom_MeshCUDA* paMesh_Top, size_t coupled_meshes, ManagedModulesCUDA& cuModule, int& coupled_cells, bool do_reduction)
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
		int cell_idx_top = i + j * n.x + (n.z - 1) * n.x*n.y;

		//skip empty cells
		if (M1.is_not_empty(cell_idx_top)) {

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuaMesh.update_parameters_mcoarse(cell_idx_top, *cuaMesh.pmu_s, mu_s);

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
				cuReal3 m_i = M1[cell_idx_top] / mu_s;

				cuBReal dot_prod = m_i * m_j;

				cuReal3 Hsurfexh = m_j * Js / ((cuBReal)MUB_MU0 * mu_s);

				if (do_reduction) {

					energy_ = -Js * dot_prod / (coupled_cells * h.dim());
				}

				Heff1[cell_idx_top] += Hsurfexh;

				//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
				if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx_top] += Hsurfexh;
				if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx_top] += energy_ * coupled_cells;

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Bottom mesh is ferromagnetic
__global__ void SurfExchangeCUDA_Bot_UpdateField(ManagedAtom_MeshCUDA& cuaMesh, ManagedAtom_MeshCUDA* paMesh_Bot, size_t coupled_meshes, ManagedModulesCUDA& cuModule, int& coupled_cells, bool do_reduction)
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
		int cell_idx_bot = i + j * n.x;

		//skip empty cells
		if (M1.is_not_empty(cell_idx_bot)) {

			cuBReal mu_s = *cuaMesh.pmu_s;
			cuBReal Js = *cuaMesh.pJs;
			cuaMesh.update_parameters_mcoarse(cell_idx_bot, *cuaMesh.pmu_s, mu_s, *cuaMesh.pJs, Js);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M1_Bot = *(paMesh_Bot[mesh_idx].pM1);

				//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M1.rect.s.x - M1_Bot.rect.s.x,
					(j + 0.5) * h.y + M1.rect.s.y - M1_Bot.rect.s.y,
					M1_Bot.rect.e.z - M1_Bot.rect.s.z - M1_Bot.h.z / 2);

				//can't couple to an empty cell
				if (!M1_Bot.rect.contains(cell_rel_pos + M1_Bot.rect.s) || M1_Bot.is_empty(cell_rel_pos)) continue;

				//yes, then get value of magnetization used in coupling with current cell at cell_idx
				cuReal3 m_j = M1_Bot[cell_rel_pos].normalized();
				cuReal3 m_i = M1[cell_idx_bot] / mu_s;

				cuBReal dot_prod = m_i * m_j;

				cuReal3 Hsurfexh = m_j * Js / ((cuBReal)MUB_MU0 * mu_s);

				if (do_reduction) {

					energy_ = -Js * dot_prod / (coupled_cells * h.dim());
				}

				Heff1[cell_idx_bot] += Hsurfexh;

				//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
				if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx_bot] += Hsurfexh;
				if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx_bot] += energy_ * coupled_cells;

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

		//Top
		if (paMesh_Top.size()) {

			SurfExchangeCUDA_Top_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, paMesh_Top, paMesh_Top.size(), cuModule, coupled_cells, true);
		}

		//Bottom
		if (paMesh_Bot.size()) {

			SurfExchangeCUDA_Bot_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, paMesh_Bot, paMesh_Bot.size(), cuModule, coupled_cells, true);
		}
	}
	else {

		//Top
		if (paMesh_Top.size()) {

			SurfExchangeCUDA_Top_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, paMesh_Top, paMesh_Top.size(), cuModule, coupled_cells, false);
		}

		//Bottom
		if (paMesh_Bot.size()) {

			SurfExchangeCUDA_Bot_UpdateField <<< (paMeshCUDA->n.x * paMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(paMeshCUDA->cuaMesh, paMesh_Bot, paMesh_Bot.size(), cuModule, coupled_cells, false);
		}
	}
}

//----------------------- Initialization

//Top mesh is ferromagnetic
__global__ void set_Atom_SurfExchangeCUDA_pointers_kernel(
	ManagedAtom_MeshCUDA& cuaMesh,
	ManagedAtom_MeshCUDA* paMesh_Bot, size_t coupled_bot_meshes,
	ManagedAtom_MeshCUDA* paMesh_Top, size_t coupled_top_meshes)
{
	if (threadIdx.x == 0) cuaMesh.paMesh_Bot = paMesh_Bot;
	if (threadIdx.x == 1) cuaMesh.paMesh_Bot_size = coupled_bot_meshes;
	if (threadIdx.x == 2) cuaMesh.paMesh_Top = paMesh_Top;
	if (threadIdx.x == 3) cuaMesh.paMesh_Top_size = coupled_top_meshes;
}

//Called by Atom_SurfExchangeCUDA module
void Atom_SurfExchangeCUDA::set_Atom_SurfExchangeCUDA_pointers(void)
{
	set_Atom_SurfExchangeCUDA_pointers_kernel <<< 1, CUDATHREADS >>>
		(paMeshCUDA->cuaMesh, paMesh_Bot, paMesh_Bot.size(), paMesh_Top, paMesh_Top.size());
}

#endif

#endif

