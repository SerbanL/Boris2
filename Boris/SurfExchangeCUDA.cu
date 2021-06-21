#include "SurfExchangeCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SURFEXCHANGE

#include "BorisCUDALib.cuh"

#include "Mesh_FerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"
#include "MeshDefs.h"

//Top mesh is ferromagnetic
__global__ void SurfExchangeCUDA_TopFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Top, size_t coupled_meshes, ManagedModulesCUDA& cuModule, int& coupled_cells, bool do_reduction)
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
		int cell_idx_top = i + j * n.x + (n.z - 1) * n.x*n.y;

		//skip empty cells
		if (M.is_not_empty(cell_idx_top)) {

			cuBReal Ms = *cuMesh.pMs;
			cuMesh.update_parameters_mcoarse(cell_idx_top, *cuMesh.pMs, Ms);

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
				cuReal3 m_i = M[cell_idx_top] / Ms;

				cuBReal dot_prod = m_i * m_j;

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexh = (m_j / ((cuBReal)MU0 * Ms * thickness)) * (J1 + 2 * J2 * dot_prod);

				if (do_reduction) {

					energy_ = (-1 * J1 - 2 * J2 * dot_prod) * dot_prod / (thickness * coupled_cells);
				}

				//couple all cells through the layer thickness : the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
				for (int k = 0; k < n.z; k++) {

					int cell_idx = i + j * n.x + k * n.x*n.y;
					Heff[cell_idx] += Hsurfexh;

					//NOTE : we must add into the module display VECs, since there could be 2 contributions for some cells (top and bottom). This is why we had to zero the VECs before calling this kernel.
					if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexh;
					if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy_ * coupled_cells;
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Top mesh is antiferromagnetic
__global__ void SurfExchangeCUDA_TopAFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Top, size_t coupled_meshes, ManagedModulesCUDA& cuModule, int& coupled_cells, bool do_reduction)
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
				cuReal3 m_i = M[cell_idx_top] / Ms;

				//total surface exchange field in coupling cells, including contributions from both sub-lattices
				cuReal3 Hsurfexh = (m_j1 / ((cuBReal)MU0 * Ms * thickness)) * J1 + (m_j2 / ((cuBReal)MU0 * Ms * thickness)) * J2;

				if (do_reduction) {

					energy_ = (-J1 * (m_i * m_j1) - J2 * (m_i * m_j2)) / (thickness * coupled_cells);
				}

				//couple all cells through the layer thickness : the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
				for (int k = 0; k < n.z; k++) {

					int cell_idx = i + j * n.x + k * n.x*n.y;
					Heff[cell_idx] += Hsurfexh;

					if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexh;
					if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy_ * coupled_cells;
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//The top mesh is diamagnetic
__global__ void SurfExchangeCUDA_TopDiamagnetic_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Top, size_t coupled_meshes, ManagedModulesCUDA& cuModule, int& coupled_cells, bool do_reduction)
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

				cuVEC_VC<cuReal3>& M_Top = *(pMesh_Top[mesh_idx].pM);

				//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M.rect.s.x - M_Top.rect.s.x,
					(j + 0.5) * h.y + M.rect.s.y - M_Top.rect.s.y,
					M_Top.h.z / 2);

				//can't couple to an empty cell
				if (!M_Top.rect.contains(cell_rel_pos + M_Top.rect.s) || M_Top.is_empty(cell_rel_pos)) continue;

				cuBReal neta_dia = *(pMesh_Top[mesh_idx].pneta_dia);
				pMesh_Top[mesh_idx].update_parameters_atposition(cell_rel_pos, *(pMesh_Top[mesh_idx].pneta_dia), neta_dia);

				//Hse = neta_dia * sus * Hext / (mu0 * Ms * tF)
				//energy_density  = -neta_dia * sus * Hext.mF / tF

				cuReal3 Mdia = M_Top[cell_rel_pos];
				cuReal3 m_i = M[cell_idx_top] / Ms;

				cuReal3 Hsurfexh = neta_dia * Mdia / (MU0 * Ms * thickness);

				if (do_reduction) {

					energy_ = -neta_dia * (Mdia*m_i) / (thickness * coupled_cells);
				}

				//couple all cells through the layer thickness : the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
				for (int k = 0; k < n.z; k++) {

					int cell_idx = i + j * n.x + k * n.x*n.y;
					Heff[cell_idx] += Hsurfexh;

					if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexh;
					if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy_ * coupled_cells;
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Bottom mesh is ferromagnetic
__global__ void SurfExchangeCUDA_BotFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Bot, size_t coupled_meshes, ManagedModulesCUDA& cuModule, int& coupled_cells, bool do_reduction)
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
		int cell_idx_bot = i + j * n.x;

		//skip empty cells
		if (M.is_not_empty(cell_idx_bot)) {

			cuBReal Ms = *cuMesh.pMs;
			cuBReal J1 = *cuMesh.pJ1;
			cuBReal J2 = *cuMesh.pJ2;
			cuMesh.update_parameters_mcoarse(cell_idx_bot, *cuMesh.pMs, Ms, *cuMesh.pJ1, J1, *cuMesh.pJ2, J2);

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
				cuReal3 m_i = M[cell_idx_bot] / Ms;

				cuBReal dot_prod = m_i * m_j;

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexh = (m_j / ((cuBReal)MU0 * Ms * thickness)) * (J1 + 2 * J2 * dot_prod);

				if (do_reduction) {

					energy_ = (-1 * J1 - 2 * J2 * dot_prod) * dot_prod / (thickness * coupled_cells);
				}

				//couple all cells through the layer thickness : the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
				for (int k = 0; k < n.z; k++) {

					int cell_idx = i + j * n.x + k * n.x*n.y;
					Heff[cell_idx] += Hsurfexh;

					if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexh;
					if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy_ * coupled_cells;
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Bottom mesh is antiferromagnetic
__global__ void SurfExchangeCUDA_BotAFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Bot, size_t coupled_meshes, ManagedModulesCUDA& cuModule, int& coupled_cells, bool do_reduction)
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
				cuReal3 m_i = M[cell_idx_bot] / Ms;

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexh = (m_j1 / ((cuBReal)MU0 * Ms * thickness)) * J1 + (m_j2 / ((cuBReal)MU0 * Ms * thickness)) * J2;

				if (do_reduction) {

					energy_ = (-J1 * (m_i * m_j1) - J2 * (m_i * m_j2)) / (thickness * coupled_cells);
				}

				//couple all cells through the layer thickness : the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
				for (int k = 0; k < n.z; k++) {

					int cell_idx = i + j * n.x + k * n.x*n.y;
					Heff[cell_idx] += Hsurfexh;

					if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexh;
					if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy_ * coupled_cells;
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Bottom mesh is diamagnetic
__global__ void SurfExchangeCUDA_BotDiamagnetic_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Bot, size_t coupled_meshes, ManagedModulesCUDA& cuModule, int& coupled_cells, bool do_reduction)
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
			cuMesh.update_parameters_mcoarse(cell_idx_bot, *cuMesh.pMs, Ms);

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

				cuBReal neta_dia = *(pMesh_Bot[mesh_idx].pneta_dia);
				pMesh_Bot[mesh_idx].update_parameters_atposition(cell_rel_pos, *(pMesh_Bot[mesh_idx].pneta_dia), neta_dia);

				//Hse = neta_dia * sus * Hext / (mu0 * Ms * tF)
				//energy_density  = -neta_dia * sus * Hext.mF / tF

				cuReal3 Mdia = M_Bot[cell_rel_pos];
				cuReal3 m_i = M[cell_idx_bot] / Ms;

				cuReal3 Hsurfexh = neta_dia * Mdia / (MU0 * Ms * thickness);

				if (do_reduction) {

					energy_ = -neta_dia * (Mdia*m_i) / (thickness * coupled_cells);
				}

				//couple all cells through the layer thickness : the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
				for (int k = 0; k < n.z; k++) {

					int cell_idx = i + j * n.x + k * n.x*n.y;
					Heff[cell_idx] += Hsurfexh;

					if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexh;
					if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy_ * coupled_cells;
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//----------------------- UpdateField LAUNCHER

void SurfExchangeCUDA::UpdateField(void)
{
	if (pMeshCUDA->GetMeshType() == MESH_DIAMAGNETIC) return;

	if (pMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();
		ZeroModuleVECs();

		//Coupling from ferromagnetic meshes

		//Top
		if (pMeshFM_Top.size()) {

			SurfExchangeCUDA_TopFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshFM_Top, pMeshFM_Top.size(), cuModule, coupled_cells, true);
		}

		//Bottom
		if (pMeshFM_Bot.size()) {

			SurfExchangeCUDA_BotFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshFM_Bot, pMeshFM_Bot.size(), cuModule, coupled_cells, true);
		}

		//Coupling from antiferromagnetic meshes

		//Top
		if (pMeshAFM_Top.size()) {

			SurfExchangeCUDA_TopAFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAFM_Top, pMeshAFM_Top.size(), cuModule, coupled_cells, true);
		}

		//Bottom
		if (pMeshAFM_Bot.size()) {

			SurfExchangeCUDA_BotAFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAFM_Bot, pMeshAFM_Bot.size(), cuModule, coupled_cells, true);
		}

		//Coupling from diamagnetic meshes

		//Top
		if (pMeshDia_Top.size()) {

			SurfExchangeCUDA_TopDiamagnetic_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshDia_Top, pMeshDia_Top.size(), cuModule, coupled_cells, true);
		}

		//Bottom
		if (pMeshDia_Bot.size()) {

			SurfExchangeCUDA_BotDiamagnetic_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >
				(pMeshCUDA->cuMesh, pMeshDia_Bot, pMeshDia_Bot.size(), cuModule, coupled_cells, true);
		}
	}
	else {

		//Coupling from ferromagnetic

		//Top
		if (pMeshFM_Top.size()) {

			SurfExchangeCUDA_TopFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshFM_Top, pMeshFM_Top.size(), cuModule, coupled_cells, false);
		}

		//Bottom
		if (pMeshFM_Bot.size()) {

			SurfExchangeCUDA_BotFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshFM_Bot, pMeshFM_Bot.size(), cuModule, coupled_cells, false);
		}

		//Coupling from antiferromagnetic

		//Top
		if (pMeshAFM_Top.size()) {

			SurfExchangeCUDA_TopAFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAFM_Top, pMeshAFM_Top.size(), cuModule, coupled_cells, false);
		}

		//Bottom
		if (pMeshAFM_Bot.size()) {

			SurfExchangeCUDA_BotAFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAFM_Bot, pMeshAFM_Bot.size(), cuModule, coupled_cells, false);
		}

		//Coupling from diamagnetic meshes

		//Top
		if (pMeshDia_Top.size()) {

			SurfExchangeCUDA_TopDiamagnetic_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshDia_Top, pMeshDia_Top.size(), cuModule, coupled_cells, false);
		}

		//Bottom
		if (pMeshDia_Bot.size()) {

			SurfExchangeCUDA_BotDiamagnetic_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshDia_Bot, pMeshDia_Bot.size(), cuModule, coupled_cells, false);
		}
	}
}

//----------------------- Initialization

//Top mesh is ferromagnetic
__global__ void set_SurfExchangeCUDA_pointers_kernel(
	ManagedMeshCUDA& cuMesh,
	ManagedMeshCUDA* pMeshFM_Bot, size_t coupledFM_bot_meshes,
	ManagedMeshCUDA* pMeshFM_Top, size_t coupledFM_top_meshes,
	ManagedMeshCUDA* pMeshAFM_Bot, size_t coupledAFM_bot_meshes,
	ManagedMeshCUDA* pMeshAFM_Top, size_t coupledAFM_top_meshes)
{
	if (threadIdx.x == 0) cuMesh.pMeshFM_Bot = pMeshFM_Bot;
	if (threadIdx.x == 1) cuMesh.pMeshFM_Bot_size = coupledFM_bot_meshes;
	if (threadIdx.x == 2) cuMesh.pMeshFM_Top = pMeshFM_Top;
	if (threadIdx.x == 3) cuMesh.pMeshFM_Top_size = coupledFM_top_meshes;

	if (threadIdx.x == 4) cuMesh.pMeshAFM_Bot = pMeshAFM_Bot;
	if (threadIdx.x == 5) cuMesh.pMeshAFM_Bot_size = coupledAFM_bot_meshes;
	if (threadIdx.x == 6) cuMesh.pMeshAFM_Top = pMeshAFM_Top;
	if (threadIdx.x == 7) cuMesh.pMeshAFM_Top_size = coupledAFM_top_meshes;
}

//Called by SurfExchangeCUDA module
void SurfExchangeCUDA::set_SurfExchangeCUDA_pointers(void)
{
	set_SurfExchangeCUDA_pointers_kernel <<< 1, CUDATHREADS >>> 
		(pMeshCUDA->cuMesh,
		pMeshFM_Bot, pMeshFM_Bot.size(), pMeshFM_Top, pMeshFM_Top.size(),
		pMeshAFM_Bot, pMeshAFM_Bot.size(), pMeshAFM_Top, pMeshAFM_Top.size());
}

#endif

#endif

