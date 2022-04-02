#include "SurfExchangeCUDA_AFM.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SURFEXCHANGE

#include "BorisCUDALib.cuh"

#include "Mesh_AntiFerromagneticCUDA.h"
#include "MeshParamsControlCUDA.h"

#include "ManagedAtom_MeshCUDA.h"
#include "Atom_MeshParamsControlCUDA.h"

#include "MeshDefs.h"

//Top mesh is ferromagnetic
__global__ void SurfExchangeCUDA_AFM_TopFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Top, size_t coupled_meshes, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	cuBReal energy_ = 0.0;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

		//skip empty cells
		if (M.is_not_empty(cell_idx)) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuMesh.update_parameters_mcoarse(cell_idx, *cuMesh.pMs_AFM, Ms_AFM);

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
				cuReal3 m_j = cu_normalize(M_Top[cell_rel_pos]);
				cuReal3 m_i1 = cu_normalize(M[cell_idx]);
				cuReal3 m_i2 = cu_normalize(M2[cell_idx]);

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexch = (m_j / ((cuBReal)MU0 * Ms_AFM.i * h.z)) * J1;
				cuReal3 Hsurfexch2 = (m_j / ((cuBReal)MU0 * Ms_AFM.j * h.z)) * J2;

				cuBReal energy1_ = 0.0, energy2_ = 0.0;

				if (do_reduction) {

					energy1_ = (-J1 * (m_i1 * m_j)) / (h.z * M.get_nonempty_cells());
					energy2_ = (-J2 * (m_i2 * m_j)) / (h.z * M.get_nonempty_cells());
					energy_ = (energy1_ + energy2_) / 2;
				}

				Heff[cell_idx] += Hsurfexch;
				Heff2[cell_idx] += Hsurfexch2;

				if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexch;
				if (do_reduction && cuModule.pModule_Heff2->linear_size()) (*cuModule.pModule_Heff2)[cell_idx] += Hsurfexch2;
				if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy1_ * M.get_nonempty_cells();
				if (do_reduction && cuModule.pModule_energy2->linear_size()) (*cuModule.pModule_energy2)[cell_idx] += energy2_ * M.get_nonempty_cells();

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Top mesh is antiferromagnetic
__global__ void SurfExchangeCUDA_AFM_TopAFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Top, size_t coupled_meshes, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	cuBReal energy_ = 0.0;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

		//skip empty cells
		if (M.is_not_empty(cell_idx)) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuMesh.update_parameters_mcoarse(cell_idx, *cuMesh.pMs_AFM, Ms_AFM);

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
				cuReal3 m_j1 = cu_normalize(M_Top[cell_rel_pos]);
				cuReal3 m_j2 = cu_normalize(M2_Top[cell_rel_pos]);
				cuReal3 m_i1 = cu_normalize(M[cell_idx]);
				cuReal3 m_i2 = cu_normalize(M2[cell_idx]);

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexch = (m_j1 / ((cuBReal)MU0 * Ms_AFM.i * h.z)) * J1;
				cuReal3 Hsurfexch2 = (m_j2 / ((cuBReal)MU0 * Ms_AFM.j * h.z)) * J2;

				cuBReal energy1_ = 0.0, energy2_ = 0.0;

				if (do_reduction) {

					energy1_ = (-J1 * (m_i1 * m_j1)) / (h.z * M.get_nonempty_cells());
					energy2_ = (-J2 * (m_i2 * m_j2)) / (h.z * M.get_nonempty_cells());
					energy_ = (energy1_ + energy2_) / 2;
				}

				Heff[cell_idx] += Hsurfexch;
				Heff2[cell_idx] += Hsurfexch2;

				if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexch;
				if (do_reduction && cuModule.pModule_Heff2->linear_size()) (*cuModule.pModule_Heff2)[cell_idx] += Hsurfexch2;
				if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy1_ * M.get_nonempty_cells();
				if (do_reduction && cuModule.pModule_energy2->linear_size()) (*cuModule.pModule_energy2)[cell_idx] += energy2_ * M.get_nonempty_cells();

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Top mesh is atomistic
__global__ void SurfExchangeCUDA_AFM_TopAtom_UpdateField(ManagedMeshCUDA& cuMesh, ManagedAtom_MeshCUDA* pMesh_Top, size_t coupled_meshes, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	cuBReal energy_ = 0.0;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

		//skip empty cells
		if (M.is_not_empty(cell_idx)) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuMesh.update_parameters_mcoarse(cell_idx, *cuMesh.pMs_AFM, Ms_AFM);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M1 = *(pMesh_Top[mesh_idx].pM1);

				//coupling rectangle in atomistic mesh in absolute coordinates
				cuRect rect_c = cuRect(
					cuReal3(i * h.x, j * h.y, M.rect.e.z),
					cuReal3((i + 1) * h.x, (j + 1) * h.y, M1.h.z + M.rect.e.z));
				rect_c += cuReal3(M.rect.s.x, M.rect.s.y, 0.0);

				//cells box in atomistic mesh
				cuBox acells = M1.box_from_rect_min(rect_c);

				//find total "directed energy" contribution from atomistic mesh : i.e. sum all mj * Js contributions from atomistic moments in the coupling area at the interface
				cuReal3 total_directed_coupling_energy1 = cuReal3();
				cuReal3 total_directed_coupling_energy2 = cuReal3();
				for (int ai = acells.s.i; ai < acells.e.i; ai++) {
					for (int aj = acells.s.j; aj < acells.e.j; aj++) {

						int acell_idx = ai + aj * M1.n.x;

						if (M1.is_empty(acell_idx)) continue;

						//Js value from atomistic mesh
						cuBReal Js = *pMesh_Top[mesh_idx].pJs;
						cuBReal Js2 = *pMesh_Top[mesh_idx].pJs2;
						cuBReal mu_s = *pMesh_Top[mesh_idx].pmu_s;
						pMesh_Top[mesh_idx].update_parameters_mcoarse(acell_idx, *pMesh_Top[mesh_idx].pJs, Js, *pMesh_Top[mesh_idx].pJs2, Js2, *pMesh_Top[mesh_idx].pmu_s, mu_s);

						total_directed_coupling_energy1 += M1[acell_idx] * Js / mu_s;
						total_directed_coupling_energy2 += M1[acell_idx] * Js2 / mu_s;
					}
				}

				//now obtain coupling field from atomistic mesh at micromagnetic cell
				cuReal3 Hsurfexch1 = (total_directed_coupling_energy1 / (h.x * h.y)) / (MU0 * Ms_AFM.i * h.z);
				cuReal3 Hsurfexch2 = (total_directed_coupling_energy2 / (h.x * h.y)) / (MU0 * Ms_AFM.j * h.z);

				cuBReal energy1_ = 0.0, energy2_ = 0.0;

				if (do_reduction) {

					energy1_ = -MU0 * M[cell_idx] * Hsurfexch1 / M.get_nonempty_cells();
					energy2_ = -MU0 * M2[cell_idx] * Hsurfexch2 / M.get_nonempty_cells();
					energy_ = (energy1_ + energy2_) / 2;
				}

				Heff[cell_idx] += Hsurfexch1;
				Heff2[cell_idx] += Hsurfexch2;

				if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexch1;
				if (do_reduction && cuModule.pModule_Heff2->linear_size()) (*cuModule.pModule_Heff2)[cell_idx] += Hsurfexch2;
				if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy1_ * M.get_nonempty_cells();
				if (do_reduction && cuModule.pModule_energy2->linear_size()) (*cuModule.pModule_energy2)[cell_idx] += energy2_ * M.get_nonempty_cells();

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Bottom mesh is ferromagnetic
__global__ void SurfExchangeCUDA_AFM_BotFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Bot, size_t coupled_meshes, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	cuBReal energy_ = 0.0;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x;

		//skip empty cells
		if (M.is_not_empty(cell_idx)) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuBReal J1 = *cuMesh.pJ1;
			cuBReal J2 = *cuMesh.pJ2;
			cuMesh.update_parameters_mcoarse(cell_idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pJ1, J1, *cuMesh.pJ2, J2);

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
				cuReal3 m_j = cu_normalize(M_Bot[cell_rel_pos]);
				cuReal3 m_i1 = cu_normalize(M[cell_idx]);
				cuReal3 m_i2 = cu_normalize(M2[cell_idx]);

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexch = (m_j / ((cuBReal)MU0 * Ms_AFM.i * h.z)) * J1;
				cuReal3 Hsurfexch2 = (m_j / ((cuBReal)MU0 * Ms_AFM.j * h.z)) * J2;

				cuBReal energy1_ = 0.0, energy2_ = 0.0;

				if (do_reduction) {

					energy1_ = (-J1 * (m_i1 * m_j)) / (h.z * M.get_nonempty_cells());
					energy2_ = (-J2 * (m_i2 * m_j)) / (h.z * M.get_nonempty_cells());
					energy_ = (energy1_ + energy2_) / 2;
				}

				Heff[cell_idx] += Hsurfexch;
				Heff2[cell_idx] += Hsurfexch2;

				if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexch;
				if (do_reduction && cuModule.pModule_Heff2->linear_size()) (*cuModule.pModule_Heff2)[cell_idx] += Hsurfexch2;
				if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy1_ * M.get_nonempty_cells();
				if (do_reduction && cuModule.pModule_energy2->linear_size()) (*cuModule.pModule_energy2)[cell_idx] += energy2_ * M.get_nonempty_cells();

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Bottom mesh is antiferromagnetic
__global__ void SurfExchangeCUDA_AFM_BotAFM_UpdateField(ManagedMeshCUDA& cuMesh, ManagedMeshCUDA* pMesh_Bot, size_t coupled_meshes, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	cuBReal energy_ = 0.0;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x;

		//skip empty cells
		if (M.is_not_empty(cell_idx)) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuBReal J1 = *cuMesh.pJ1;
			cuBReal J2 = *cuMesh.pJ2;
			cuMesh.update_parameters_mcoarse(cell_idx, *cuMesh.pMs_AFM, Ms_AFM, *cuMesh.pJ1, J1, *cuMesh.pJ2, J2);

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
				cuReal3 m_j1 = cu_normalize(M_Bot[cell_rel_pos]);
				cuReal3 m_j2 = cu_normalize(M2_Bot[cell_rel_pos]);
				cuReal3 m_i1 = cu_normalize(M[cell_idx]);
				cuReal3 m_i2 = cu_normalize(M2[cell_idx]);

				//total surface exchange field in coupling cells, including bilinear and biquadratic terms
				cuReal3 Hsurfexch = (m_j1 / ((cuBReal)MU0 * Ms_AFM.i * h.z)) * J1;
				cuReal3 Hsurfexch2 = (m_j2 / ((cuBReal)MU0 * Ms_AFM.j * h.z)) * J2;

				cuBReal energy1_ = 0.0, energy2_ = 0.0;

				if (do_reduction) {

					energy1_ = (-J1 * (m_i1 * m_j1)) / (h.z * M.get_nonempty_cells());
					energy2_ = (-J2 * (m_i2 * m_j2)) / (h.z * M.get_nonempty_cells());
					energy_ = (energy1_ + energy2_) / 2;
				}

				Heff[cell_idx] += Hsurfexch;
				Heff2[cell_idx] += Hsurfexch2;

				if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexch;
				if (do_reduction && cuModule.pModule_Heff2->linear_size()) (*cuModule.pModule_Heff2)[cell_idx] += Hsurfexch2;
				if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy1_ * M.get_nonempty_cells();
				if (do_reduction && cuModule.pModule_energy2->linear_size()) (*cuModule.pModule_energy2)[cell_idx] += energy2_ * M.get_nonempty_cells();

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//Bottom mesh is atomistic
__global__ void SurfExchangeCUDA_AFM_BotAtom_UpdateField(ManagedMeshCUDA& cuMesh, ManagedAtom_MeshCUDA* pMesh_Bot, size_t coupled_meshes, ManagedModulesCUDA& cuModule, bool do_reduction)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;
	cuVEC<cuReal3>& Heff = *cuMesh.pHeff;
	cuVEC<cuReal3>& Heff2 = *cuMesh.pHeff2;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = M.n;
	cuReal3 h = M.h;

	cuBReal energy_ = 0.0;

	if (idx < n.x * n.y) {

		int i = idx % n.x;
		int j = idx / n.x;
		int cell_idx = i + j * n.x;

		//skip empty cells
		if (M.is_not_empty(cell_idx)) {

			cuReal2 Ms_AFM = *cuMesh.pMs_AFM;
			cuMesh.update_parameters_mcoarse(cell_idx, *cuMesh.pMs_AFM, Ms_AFM);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < coupled_meshes; mesh_idx++) {

				cuVEC_VC<cuReal3>& M1 = *pMesh_Bot[mesh_idx].pM1;

				//coupling rectangle in atomistic mesh in absolute coordinates
				cuRect rect_c = cuRect(
					cuReal3(i * h.x, j * h.y, M1.rect.e.z - M1.h.z),
					cuReal3((i + 1) * h.x, (j + 1) * h.y, M1.rect.e.z));
				rect_c += cuReal3(M.rect.s.x, M.rect.s.y, 0.0);

				//cells box in atomistic mesh
				cuBox acells = M1.box_from_rect_min(rect_c);

				//find total "directed energy" contribution from atomistic mesh : i.e. sum all mj * Js contributions from atomistic moments in the coupling area at the interface
				//NOTE : at atomistic/micromagnetic coupling, it's the atomistic mesh which sets coupling constant, not the top mesh
				cuReal3 total_directed_coupling_energy1 = cuReal3();
				cuReal3 total_directed_coupling_energy2 = cuReal3();
				for (int ai = acells.s.i; ai < acells.e.i; ai++) {
					for (int aj = acells.s.j; aj < acells.e.j; aj++) {

						int acell_idx = ai + aj * M1.n.x + (M1.n.z - 1) * M1.n.x * M1.n.y;

						if (M1.is_empty(acell_idx)) continue;

						//Js value from atomistic mesh
						cuBReal Js = *pMesh_Bot[mesh_idx].pJs;
						cuBReal Js2 = *pMesh_Bot[mesh_idx].pJs2;
						cuBReal mu_s = *pMesh_Bot[mesh_idx].pmu_s;
						pMesh_Bot[mesh_idx].update_parameters_mcoarse(acell_idx, *pMesh_Bot[mesh_idx].pJs, Js, *pMesh_Bot[mesh_idx].pJs2, Js2, *pMesh_Bot[mesh_idx].pmu_s, mu_s);

						total_directed_coupling_energy1 += M1[acell_idx] * Js / mu_s;
						total_directed_coupling_energy2 += M1[acell_idx] * Js2 / mu_s;
					}
				}

				//now obtain coupling field from atomistic mesh at micromagnetic cell
				cuReal3 Hsurfexch1 = (total_directed_coupling_energy1 / (h.x * h.y)) / (MU0 * Ms_AFM.i * h.z);
				cuReal3 Hsurfexch2 = (total_directed_coupling_energy2 / (h.x * h.y)) / (MU0 * Ms_AFM.j * h.z);

				cuBReal energy1_ = 0.0, energy2_ = 0.0;

				if (do_reduction) {

					energy1_ = -MU0 * M[cell_idx] * Hsurfexch1 / M.get_nonempty_cells();
					energy2_ = -MU0 * M2[cell_idx] * Hsurfexch2 / M.get_nonempty_cells();
					energy_ = (energy1_ + energy2_) / 2;
				}

				Heff[cell_idx] += Hsurfexch1;
				Heff2[cell_idx] += Hsurfexch2;

				if (do_reduction && cuModule.pModule_Heff->linear_size()) (*cuModule.pModule_Heff)[cell_idx] += Hsurfexch1;
				if (do_reduction && cuModule.pModule_Heff2->linear_size()) (*cuModule.pModule_Heff2)[cell_idx] += Hsurfexch2;
				if (do_reduction && cuModule.pModule_energy->linear_size()) (*cuModule.pModule_energy)[cell_idx] += energy1_ * M.get_nonempty_cells();
				if (do_reduction && cuModule.pModule_energy2->linear_size()) (*cuModule.pModule_energy2)[cell_idx] += energy2_ * M.get_nonempty_cells();

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, *cuModule.penergy);
}

//----------------------- UpdateField LAUNCHER

void SurfExchangeCUDA_AFM::UpdateField(void)
{
	if (pMeshCUDA->CurrentTimeStepSolved()) {

		ZeroEnergy();
		ZeroModuleVECs();

		//Coupling from ferromagnetic meshes

		//Top
		if (pMeshFM_Top.size()) {

			SurfExchangeCUDA_AFM_TopFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshFM_Top, pMeshFM_Top.size(), cuModule, true);
		}

		//Bottom
		if (pMeshFM_Bot.size()) {

			SurfExchangeCUDA_AFM_BotFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshFM_Bot, pMeshFM_Bot.size(), cuModule, true);
		}

		//Coupling from antiferromagnetic meshes

		//Top
		if (pMeshAFM_Top.size()) {

			SurfExchangeCUDA_AFM_TopAFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAFM_Top, pMeshAFM_Top.size(), cuModule, true);
		}

		//Bottom
		if (pMeshAFM_Bot.size()) {

			SurfExchangeCUDA_AFM_BotAFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAFM_Bot, pMeshAFM_Bot.size(), cuModule, true);
		}

		//Coupling from atomistic meshes

		//Top
		if (pMeshAtom_Top.size()) {

			SurfExchangeCUDA_AFM_TopAtom_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAtom_Top, pMeshAtom_Top.size(), cuModule, true);
		}

		//Bottom
		if (pMeshAtom_Bot.size()) {

			SurfExchangeCUDA_AFM_BotAtom_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAtom_Bot, pMeshAtom_Bot.size(), cuModule, true);
		}
	}
	else {

		//Coupling from ferromagnetic meshes

		//Top
		if (pMeshFM_Top.size()) {

			SurfExchangeCUDA_AFM_TopFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshFM_Top, pMeshFM_Top.size(), cuModule, false);
		}

		//Bottom
		if (pMeshFM_Bot.size()) {

			SurfExchangeCUDA_AFM_BotFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshFM_Bot, pMeshFM_Bot.size(), cuModule, false);
		}

		//Coupling from antiferromagnetic meshes

		//Top
		if (pMeshAFM_Top.size()) {

			SurfExchangeCUDA_AFM_TopAFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAFM_Top, pMeshAFM_Top.size(), cuModule, false);
		}

		//Bottom
		if (pMeshAFM_Bot.size()) {

			SurfExchangeCUDA_AFM_BotAFM_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAFM_Bot, pMeshAFM_Bot.size(), cuModule, false);
		}

		//Coupling from atomistic meshes

		//Top
		if (pMeshAtom_Top.size()) {

			SurfExchangeCUDA_AFM_TopAtom_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAtom_Top, pMeshAtom_Top.size(), cuModule, false);
		}

		//Bottom
		if (pMeshAtom_Bot.size()) {

			SurfExchangeCUDA_AFM_BotAtom_UpdateField <<< (pMeshCUDA->n.x * pMeshCUDA->n.y + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(pMeshCUDA->cuMesh, pMeshAtom_Bot, pMeshAtom_Bot.size(), cuModule, false);
		}
	}
}

//----------------------- Initialization

//Current mesh is anti-ferromagnetic
__global__ void set_SurfExchangeCUDA_AFM_pointers_kernel(
	ManagedMeshCUDA& cuMesh,
	ManagedMeshCUDA* pMeshFM_Bot, size_t coupledFM_bot_meshes,
	ManagedMeshCUDA* pMeshFM_Top, size_t coupledFM_top_meshes,
	ManagedMeshCUDA* pMeshAFM_Bot, size_t coupledAFM_bot_meshes,
	ManagedMeshCUDA* pMeshAFM_Top, size_t coupledAFM_top_meshes,
	ManagedAtom_MeshCUDA* pMeshAtom_Bot, size_t coupledAtom_bot_meshes,
	ManagedAtom_MeshCUDA* pMeshAtom_Top, size_t coupledAtom_top_meshes)
{
	if (threadIdx.x == 0) cuMesh.pMeshFM_Bot = pMeshFM_Bot;
	if (threadIdx.x == 1) cuMesh.pMeshFM_Bot_size = coupledFM_bot_meshes;
	if (threadIdx.x == 2) cuMesh.pMeshFM_Top = pMeshFM_Top;
	if (threadIdx.x == 3) cuMesh.pMeshFM_Top_size = coupledFM_top_meshes;

	if (threadIdx.x == 4) cuMesh.pMeshAFM_Bot = pMeshAFM_Bot;
	if (threadIdx.x == 5) cuMesh.pMeshAFM_Bot_size = coupledAFM_bot_meshes;
	if (threadIdx.x == 6) cuMesh.pMeshAFM_Top = pMeshAFM_Top;
	if (threadIdx.x == 7) cuMesh.pMeshAFM_Top_size = coupledAFM_top_meshes;

	if (threadIdx.x == 8) cuMesh.pMeshAtom_Bot = pMeshAtom_Bot;
	if (threadIdx.x == 9) cuMesh.pMeshAtom_Bot_size = coupledAtom_bot_meshes;
	if (threadIdx.x == 10) cuMesh.pMeshAtom_Top = pMeshAtom_Top;
	if (threadIdx.x == 11) cuMesh.pMeshAtom_Top_size = coupledAtom_top_meshes;
}

//Called by SurfExchangeCUDA module
void SurfExchangeCUDA_AFM::set_SurfExchangeCUDA_AFM_pointers(void)
{
	set_SurfExchangeCUDA_AFM_pointers_kernel <<< 1, CUDATHREADS >>>
		(pMeshCUDA->cuMesh,
		pMeshFM_Bot, pMeshFM_Bot.size(), pMeshFM_Top, pMeshFM_Top.size(),
		pMeshAFM_Bot, pMeshAFM_Bot.size(), pMeshAFM_Top, pMeshAFM_Top.size(),
		pMeshAtom_Bot, pMeshAtom_Bot.size(), pMeshAtom_Top, pMeshAtom_Top.size());
}

#endif

#endif