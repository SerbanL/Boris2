#include "SurfExchangeCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_SURFEXCHANGE

#include "BorisCUDALib.cuh"

#include "Mesh_FerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"
#include "MeshDefs.h"

__global__ void SurfExchangeCUDA_Top_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Top, size_t coupled_meshes, cuBReal& energy, int& coupled_cells, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	//thickness of layer - SurfExchange applies for layers in the xy plane
	cuBReal thickness = M.rect.e.z - M.rect.s.z;

	cuBReal energy_ = 0.0;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int k = idx / (n.x * n.y);
		int cell_idx_top = i + j * n.x + (n.z - 1) * n.x*n.y;

		//skip empty cells
		if (M.is_not_empty(cell_idx_top)) {

			cuBReal Ms = *cuMesh.pMs;
			cuMesh.update_parameters_mcoarse(cell_idx_top, *cuMesh.pMs, Ms);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				//relative coordinates to read value from top mesh (the one we're coupling to here)
				cuReal3 cell_rel_pos = cuReal3((i + 0.5) * h.x, (j + 0.5) * h.y, pMesh_Top[mesh_idx].pM->h.z / 2);

				cuVEC_VC<cuReal3>& M_Top = *(pMesh_Top[mesh_idx].pM);

				//can't couple to an empty cell
				if (cuIsZ(M_Top[cell_rel_pos].norm())) continue;

				cuBReal J1 = *(pMesh_Top[mesh_idx].pJ1);
				cuBReal J2 = *(pMesh_Top[mesh_idx].pJ2);
				pMesh_Top[mesh_idx].update_parameters_atposition(cell_rel_pos, *(pMesh_Top[mesh_idx].pJ1), J1, *(pMesh_Top[mesh_idx].pJ2), J2);

				//get magnetization value in top mesh cell to couple with
				cuReal3 m_j = M_Top[cell_rel_pos].normalized();
				cuReal3 m_i = M[cell_idx_top] / Ms;

				cuBReal dot_prod = m_i * m_j;

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexh = (m_j / ((cuBReal)MU0 * Ms * thickness)) * (J1 + 2 * J2 * dot_prod);

				//couple all cells through the layer thickness : the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
				Heff[i + j * n.x + k * n.x*n.y] += Hsurfexh;

				if (do_reduction && k == n.z - 1) {

					energy_ = (-1 * J1 - 2 * J2 * dot_prod) * dot_prod / (thickness * coupled_cells);
				}
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

__global__ void SurfExchangeCUDA_Bot_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Bot, size_t coupled_meshes, cuBReal& energy, int& coupled_cells, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	//thickness of layer - SurfExchange applies for layers in the xy plane
	cuBReal thickness = M.rect.e.z - M.rect.s.z;

	cuBReal energy_ = 0.0;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int k = idx / (n.x * n.y);
		int cell_idx_bot = i + j * n.x;

		//skip empty cells
		if (M.is_not_empty(cell_idx_bot)) {

			cuBReal Ms = *cuMesh.pMs;
			cuBReal J1 = *cuMesh.pJ1;
			cuBReal J2 = *cuMesh.pJ2;
			cuMesh.update_parameters_mcoarse(cell_idx_bot, *cuMesh.pMs, Ms, *cuMesh.pJ1, J1, *cuMesh.pJ2, J2);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				//relative coordinates to read value from bottom mesh (the one we're coupling to here)
				cuReal3 cell_rel_pos = cuReal3((i + 0.5) * h.x, (j + 0.5) * h.y, pMesh_Bot[mesh_idx].pM->rect.e.z - (pMesh_Bot[mesh_idx].pM->h.z / 2));

				cuVEC_VC<cuReal3>& M_Bot = *(pMesh_Bot[mesh_idx].pM);

				//can't couple to an empty cell
				if (cuIsZ(M_Bot[cell_rel_pos].norm())) continue;

				//yes, then get value of magnetization used in coupling with current cell at cell_idx
				cuReal3 m_j = M_Bot[cell_rel_pos].normalized();
				cuReal3 m_i = M[cell_idx_bot] / Ms;

				cuBReal dot_prod = m_i * m_j;

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexh = (m_j / ((cuBReal)MU0 * Ms * thickness)) * (J1 + 2 * J2 * dot_prod);

				Heff[i + j * n.x + k * n.x*n.y] += Hsurfexh;

				if (do_reduction && k == 0) {

					energy_ = (-1 * J1 - 2 * J2 * dot_prod) * dot_prod / (thickness * coupled_cells);
				}
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void SurfExchangeCUDA::UpdateField(void)
{
	//antiferromagnetic meshes can have a SurfExchange module, but in this case it's just a dummy module to enable exchange biasing of a neighboring ferromagnetic mesh
	//we don't include a surface exchange field in the antiferromagnetic mesh
	if (pMeshCUDA->GetMeshType() == MESH_ANTIFERROMAGNETIC) return;

	if (pMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();

		//Top
		if (pMesh_Top.size()) {

			SurfExchangeCUDA_Top_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMesh_Top, pMesh_Top.size(), energy, coupled_cells, true);
		}

		//Bottom
		if (pMesh_Bot.size()) {

			SurfExchangeCUDA_Bot_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMesh_Bot, pMesh_Bot.size(), energy, coupled_cells, true);
		}
	}
	else {

		//Top
		if (pMesh_Top.size()) {

			SurfExchangeCUDA_Top_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMesh_Top, pMesh_Top.size(), energy, coupled_cells, false);
		}

		//Bottom
		if (pMesh_Bot.size()) {

			SurfExchangeCUDA_Bot_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMesh_Bot, pMesh_Bot.size(), energy, coupled_cells, false);
		}
	}
}

#endif

#endif

