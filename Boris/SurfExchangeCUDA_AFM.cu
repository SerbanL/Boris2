#include "SurfExchangeCUDA_AFM.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SURFEXCHANGE

#include "BorisCUDALib.cuh"

#include "Mesh_AntiFerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"
#include "MeshDefs.h"

//Top mesh is ferromagnetic
__global__ void SurfExchangeCUDA_AFM_TopFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Top, size_t coupled_meshes, cuBReal& energy, int& coupled_cells, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

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

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuMesh.update_parameters_mcoarse(cell_idx_top, *cuMesh.pMs_AFM, Ms_AFM);

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

				cuBReal J1 = *(pMesh_Top[mesh_idx].pJ1);
				cuBReal J2 = *(pMesh_Top[mesh_idx].pJ2);
				pMesh_Top[mesh_idx].update_parameters_atposition(cell_rel_pos, *(pMesh_Top[mesh_idx].pJ1), J1, *(pMesh_Top[mesh_idx].pJ2), J2);

				//get magnetization value in top mesh cell to couple with
				cuReal3 m_j = M_Top[cell_rel_pos].normalized();
				cuReal3 m_i1 = M[cell_idx_top] / Ms_AFM.i;
				cuReal3 m_i2 = M2[cell_idx_top] / Ms_AFM.j;

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexh = (m_j / ((cuBReal)MU0 * Ms_AFM.i * thickness)) * J1;
				cuReal3 Hsurfexh2 = (m_j / ((cuBReal)MU0 * Ms_AFM.j * thickness)) * J2;

				//couple all cells through the layer thickness : the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
				Heff[i + j * n.x + k * n.x*n.y] += Hsurfexh;
				Heff2[i + j * n.x + k * n.x*n.y] += Hsurfexh2;

				if (do_reduction && k == n.z - 1) {

					energy_ = (-J1 * (m_i1 * m_j) - J2 * (m_i2 * m_j)) / (thickness * coupled_cells);
				}
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//Top mesh is antiferromagnetic
__global__ void SurfExchangeCUDA_AFM_TopAFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Top, size_t coupled_meshes, cuBReal& energy, int& coupled_cells, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

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

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuMesh.update_parameters_mcoarse(cell_idx_top, *cuMesh.pMs_AFM, Ms_AFM);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M_Top = *(pMesh_Top[mesh_idx].pM);
				cuVEC_VC<cuReal3>& M2_Top = *(pMesh_Top[mesh_idx].pM2);

				//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M.rect.s.x - M_Top.rect.s.x,
					(j + 0.5) * h.y + M.rect.s.y - M_Top.rect.s.y,
					M_Top.h.z / 2);

				//can't couple to an empty cell
				if (!M_Top.rect.contains(cell_rel_pos + M_Top.rect.s) || M_Top.is_empty(cell_rel_pos)) continue;

				cuBReal J1 = *(pMesh_Top[mesh_idx].pJ1);
				cuBReal J2 = *(pMesh_Top[mesh_idx].pJ2);
				pMesh_Top[mesh_idx].update_parameters_atposition(cell_rel_pos, *(pMesh_Top[mesh_idx].pJ1), J1, *(pMesh_Top[mesh_idx].pJ2), J2);

				//get magnetization value in top mesh cell to couple with
				cuReal3 m_j1 = M_Top[cell_rel_pos].normalized();
				cuReal3 m_j2 = M2_Top[cell_rel_pos].normalized();
				cuReal3 m_i1 = M[cell_idx_top] / Ms_AFM.i;
				cuReal3 m_i2 = M2[cell_idx_top] / Ms_AFM.j;

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexh = (m_j1 / ((cuBReal)MU0 * Ms_AFM.i * thickness)) * J1;
				cuReal3 Hsurfexh2 = (m_j2 / ((cuBReal)MU0 * Ms_AFM.j * thickness)) * J2;

				//couple all cells through the layer thickness : the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
				Heff[i + j * n.x + k * n.x*n.y] += Hsurfexh;
				Heff2[i + j * n.x + k * n.x*n.y] += Hsurfexh2;

				if (do_reduction && k == n.z - 1) {

					energy_ = (-J1 * (m_i1 * m_j1) - J2 * (m_i2 * m_j2)) / (thickness * coupled_cells);
				}
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//Bottom mesh is ferromagnetic
__global__ void SurfExchangeCUDA_AFM_BotFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Bot, size_t coupled_meshes, cuBReal& energy, int& coupled_cells, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

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

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuBReal J1 = *cuMesh.pJ1;
			cuBReal J2 = *cuMesh.pJ2;
			cuMesh.update_parameters_mcoarse(cell_idx_bot, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pJ1, J1, *cuMesh.pJ2, J2);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M_Bot = *(pMesh_Bot[mesh_idx].pM);

				//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M.rect.s.x - M_Bot.rect.s.x,
					(j + 0.5) * h.y + M.rect.s.y - M_Bot.rect.s.y,
					M_Bot.rect.e.z - M_Bot.rect.s.z - M_Bot.h.z / 2);

				//can't couple to an empty cell
				if (!M_Bot.rect.contains(cell_rel_pos + M_Bot.rect.s) || M_Bot.is_empty(cell_rel_pos)) continue;

				//yes, then get value of magnetization used in coupling with current cell at cell_idx
				cuReal3 m_j = M_Bot[cell_rel_pos].normalized();
				cuReal3 m_i1 = M[cell_idx_bot] / Ms_AFM.i;
				cuReal3 m_i2 = M2[cell_idx_bot] / Ms_AFM.j;

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexh = (m_j / ((cuBReal)MU0 * Ms_AFM.i * thickness)) * J1;
				cuReal3 Hsurfexh2 = (m_j / ((cuBReal)MU0 * Ms_AFM.j * thickness)) * J2;

				Heff[i + j * n.x + k * n.x*n.y] += Hsurfexh;
				Heff2[i + j * n.x + k * n.x*n.y] += Hsurfexh2;

				if (do_reduction && k == 0) {

					energy_ = (-J1 * (m_i1 * m_j) - J2 * (m_i2 * m_j)) / (thickness * coupled_cells);
				}
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//Bottom mesh is antiferromagnetic
__global__ void SurfExchangeCUDA_AFM_BotAFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Bot, size_t coupled_meshes, cuBReal& energy, int& coupled_cells, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

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

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuBReal J1 = *cuMesh.pJ1;
			cuBReal J2 = *cuMesh.pJ2;
			cuMesh.update_parameters_mcoarse(cell_idx_bot, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pJ1, J1, *cuMesh.pJ2, J2);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M_Bot = *(pMesh_Bot[mesh_idx].pM);
				cuVEC_VC<cuReal3>& M2_Bot = *(pMesh_Bot[mesh_idx].pM2);

				//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M.rect.s.x - M_Bot.rect.s.x,
					(j + 0.5) * h.y + M.rect.s.y - M_Bot.rect.s.y,
					M_Bot.rect.e.z - M_Bot.rect.s.z - M_Bot.h.z / 2);

				//can't couple to an empty cell
				if (!M_Bot.rect.contains(cell_rel_pos + M_Bot.rect.s) || M_Bot.is_empty(cell_rel_pos)) continue;

				//yes, then get value of magnetization used in coupling with current cell at cell_idx
				cuReal3 m_j1 = M_Bot[cell_rel_pos].normalized();
				cuReal3 m_j2 = M2_Bot[cell_rel_pos].normalized();
				cuReal3 m_i1 = M[cell_idx_bot] / Ms_AFM.i;
				cuReal3 m_i2 = M2[cell_idx_bot] / Ms_AFM.j;

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexh = (m_j1 / ((cuBReal)MU0 * Ms_AFM.i * thickness)) * J1;
				cuReal3 Hsurfexh2 = (m_j2 / ((cuBReal)MU0 * Ms_AFM.j * thickness)) * J2;

				Heff[i + j * n.x + k * n.x*n.y] += Hsurfexh;
				Heff2[i + j * n.x + k * n.x*n.y] += Hsurfexh2;

				if (do_reduction && k == 0) {

					energy_ = (-J1 * (m_i1 * m_j1) - J2 * (m_i2 * m_j2)) / (thickness * coupled_cells);
				}
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//----------------------- UpdateField LAUNCHER

void SurfExchangeCUDA_AFM::UpdateField(void)
{
	if (pMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();

		//Coupling from ferromagnetic meshes

		//Top
		if (pMeshFM_Top.size()) {

			SurfExchangeCUDA_AFM_TopFM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshFM_Top, pMeshFM_Top.size(), energy, coupled_cells, true);
		}

		//Bottom
		if (pMeshFM_Bot.size()) {

			SurfExchangeCUDA_AFM_BotFM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshFM_Bot, pMeshFM_Bot.size(), energy, coupled_cells, true);
		}

		//Coupling from antiferromagnetic meshes

		//Top
		if (pMeshAFM_Top.size()) {

			SurfExchangeCUDA_AFM_TopAFM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAFM_Top, pMeshAFM_Top.size(), energy, coupled_cells, true);
		}

		//Bottom
		if (pMeshAFM_Bot.size()) {

			SurfExchangeCUDA_AFM_BotAFM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAFM_Bot, pMeshAFM_Bot.size(), energy, coupled_cells, true);
		}
	}
	else {

		//Coupling from ferromagnetic meshes

		//Top
		if (pMeshFM_Top.size()) {

			SurfExchangeCUDA_AFM_TopFM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshFM_Top, pMeshFM_Top.size(), energy, coupled_cells, false);
		}

		//Bottom
		if (pMeshFM_Bot.size()) {

			SurfExchangeCUDA_AFM_BotFM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshFM_Bot, pMeshFM_Bot.size(), energy, coupled_cells, false);
		}

		//Coupling from antiferromagnetic meshes

		//Top
		if (pMeshAFM_Top.size()) {

			SurfExchangeCUDA_AFM_TopAFM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAFM_Top, pMeshAFM_Top.size(), energy, coupled_cells, false);
		}

		//Bottom
		if (pMeshAFM_Bot.size()) {

			SurfExchangeCUDA_AFM_BotAFM_UpdateField <<< (pMeshCUDA->n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAFM_Bot, pMeshAFM_Bot.size(), energy, coupled_cells, false);
		}
	}
}

#endif

#endif

