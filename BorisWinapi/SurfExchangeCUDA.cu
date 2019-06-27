#include "SurfExchangeCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_SURFEXCHANGE

#include "BorisCUDALib.cuh"

#include "Mesh_FerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

__global__ void SurfExchangeCUDA_Top_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Top, size_t coupled_meshes, cuReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	cuReal energy_ = 0.0;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

		//skip empty cells
		if (M.is_not_empty(cell_idx)) {

			cuReal Ms = *cuMesh.pMs;
			cuMesh.update_parameters_mcoarse(cell_idx, *cuMesh.pMs, Ms);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				//relative coordinates to read value from top mesh (the one we're coupling to here)
				cuReal3 cell_rel_pos = cuReal3((i + 0.5) * h.x, (j + 0.5) * h.y, pMesh_Top[mesh_idx].pM->h.z / 2);

				cuVEC_VC<cuReal3>& M_Top = *(pMesh_Top[mesh_idx].pM);

				//can't couple to an empty cell
				if (cuIsZ(M_Top[cell_rel_pos].norm())) continue;

				cuReal J1 = *(pMesh_Top[mesh_idx].pJ1);
				cuReal J2 = *(pMesh_Top[mesh_idx].pJ2);
				pMesh_Top[mesh_idx].update_parameters_atposition(cell_rel_pos, *(pMesh_Top[mesh_idx].pJ1), J1, *(pMesh_Top[mesh_idx].pJ2), J2);

				//get magnetization value in top mesh cell to couple with
				cuReal3 m_j = M_Top[cell_rel_pos].normalized();
				cuReal3 m_i = M[cell_idx] / Ms;

				cuReal dot_prod = m_i * m_j;

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexh = (m_j / ((cuReal)MU0 * Ms * h.z)) * (J1 + 2 * J2 * dot_prod);

				Heff[cell_idx] += Hsurfexh;

				if (do_reduction) {

					int non_empty_cells = M.get_nonempty_cells();
					if (non_empty_cells) energy_ = ((-1 * J1 - 2 * J2 * dot_prod) * dot_prod / h.z) / non_empty_cells;
				}
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

__global__ void SurfExchangeCUDA_Bot_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Bot, size_t coupled_meshes, cuReal& energy, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	cuReal energy_ = 0.0;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x;

		//skip empty cells
		if (M.is_not_empty(cell_idx)) {

			cuReal Ms = *cuMesh.pMs;
			cuReal J1 = *cuMesh.pJ1;
			cuReal J2 = *cuMesh.pJ2;
			cuMesh.update_parameters_mcoarse(cell_idx, *cuMesh.pMs, Ms, *cuMesh.pJ1, J1, *cuMesh.pJ2, J2);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				//relative coordinates to read value from bottom mesh (the one we're coupling to here)
				cuReal3 cell_rel_pos = cuReal3((i + 0.5) * h.x, (j + 0.5) * h.y, pMesh_Bot[mesh_idx].pM->rect.e.z - (pMesh_Bot[mesh_idx].pM->h.z / 2));

				cuVEC_VC<cuReal3>& M_Bot = *(pMesh_Bot[mesh_idx].pM);

				//can't couple to an empty cell
				if (cuIsZ(M_Bot[cell_rel_pos].norm())) continue;

				//yes, then get value of magnetization used in coupling with current cell at cell_idx
				cuReal3 m_j = M_Bot[cell_rel_pos].normalized();
				cuReal3 m_i = M[cell_idx] / Ms;

				cuReal dot_prod = m_i * m_j;

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexh = (m_j / ((cuReal)MU0 * Ms * h.z)) * (J1 + 2 * J2 * dot_prod);

				Heff[cell_idx] += Hsurfexh;

				if (do_reduction) {

					int non_empty_cells = M.get_nonempty_cells();
					if (non_empty_cells) energy_ = ((-1 * J1 - 2 * J2 * dot_prod) * dot_prod / h.z) / non_empty_cells;
				}
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void SurfExchangeCUDA::UpdateField(void)
{
	if (pMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();

		//Top
		if (pMesh_Top.size()) {

			SurfExchangeCUDA_Top_UpdateField << < (pMeshCUDA->n.x*pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >
				(pMeshCUDA->cuMesh, pMesh_Top, pMesh_Top.size(), energy, true);
		}

		//Bottom
		if (pMesh_Bot.size()) {

			SurfExchangeCUDA_Bot_UpdateField << < (pMeshCUDA->n.x*pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >
				(pMeshCUDA->cuMesh, pMesh_Bot, pMesh_Bot.size(), energy, true);
		}
	}
	else {

		//Top
		if (pMesh_Top.size()) {

			SurfExchangeCUDA_Top_UpdateField << < (pMeshCUDA->n.x*pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >
				(pMeshCUDA->cuMesh, pMesh_Top, pMesh_Top.size(), energy, false);
		}

		//Bottom
		if (pMesh_Bot.size()) {

			SurfExchangeCUDA_Bot_UpdateField << < (pMeshCUDA->n.x*pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >
				(pMeshCUDA->cuMesh, pMesh_Bot, pMesh_Bot.size(), energy, false);
		}
	}
}

#endif

#endif

