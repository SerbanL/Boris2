#include "ManagedAtom_MeshCUDA.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"
#include "Atom_MeshParamsControlCUDA.h"

#include "ManagedMeshCUDA.h"
#include "MeshParamsControlCUDA.h"

#include "ModulesDefs.h"

//----------------------------------- MONTE-CARLO METHODS FOR ENERGY COMPUTATION

//Energy Deltas

//SIMPLE CUBIC

//switch function which adds all assigned energy contributions in this mesh
//If Mnew is passed in as cuReal3(), then this function returns the current spin energy only - all functions in the switch statement below implement this eventuality.
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC(int spin_index, cuReal3 Mnew, int*& cuaModules, int numModules, cuReal3& Ha)
{
	cuBReal energy = 0.0;

	for (int idx = 0; idx < numModules; idx++) {

		switch (cuaModules[idx]) {

		case MOD_DEMAG_N:
			energy += Get_Atomistic_EnergyChange_SC_DemagNCUDA(spin_index, Mnew);
			break;

		case MOD_DEMAG:
			energy += Get_Atomistic_EnergyChange_SC_DemagCUDA(spin_index, Mnew);
			break;

		case MOD_ATOM_DIPOLEDIPOLE:
			energy += Get_Atomistic_EnergyChange_SC_DipoleDipoleCUDA(spin_index, Mnew);
			break;

		case MOD_STRAYFIELD_MESH:
			energy += Get_Atomistic_EnergyChange_SC_StrayField_AtomMeshCUDA(spin_index, Mnew);
			break;

		case MOD_EXCHANGE:
			energy += Get_Atomistic_EnergyChange_SC_ExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_DMEXCHANGE:
			energy += Get_Atomistic_EnergyChange_SC_DMExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_IDMEXCHANGE:
			energy += Get_Atomistic_EnergyChange_SC_iDMExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_VIDMEXCHANGE:
			energy += Get_Atomistic_EnergyChange_SC_viDMExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_SURFEXCHANGE:
			energy += Get_Atomistic_EnergyChange_SC_SurfExchangeCUDA(spin_index, Mnew);
			break;

		case MOD_ZEEMAN:
			energy += Get_Atomistic_EnergyChange_SC_ZeemanCUDA(spin_index, Mnew, Ha);
			break;

		case MOD_MOPTICAL:
			energy += Get_Atomistic_EnergyChange_SC_MOpticalCUDA(spin_index, Mnew);
			break;

		case MOD_ANIUNI:
			energy += Get_Atomistic_EnergyChange_SC_AnisotropyCUDA(spin_index, Mnew);
			break;

		case MOD_ANICUBI:
			energy += Get_Atomistic_EnergyChange_SC_AnisotropyCubiCUDA(spin_index, Mnew);
			break;

		case MOD_ANIBI:
			energy += Get_Atomistic_EnergyChange_SC_AnisotropyBiaxialCUDA(spin_index, Mnew);
			break;

		case MOD_ANITENS:
			energy += Get_Atomistic_EnergyChange_SC_AnisotropyTensorialCUDA(spin_index, Mnew);
			break;

		default:
			break;
		}
	}

	return energy;
}

//Atom_Demag_N
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_DemagNCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	//Nxy shouldn't have a temperature (or spatial) dependence so not using update_parameters_mcoarse here
	cuReal2 Nxy = *pNxy;

	cuBReal Nz = (1 - Nxy.x - Nxy.y);

	cuBReal r = (cuBReal)MUB / M1.h.dim();
	cuReal3 S = M1[spin_index];

	if (Mnew != cuReal3()) return ((cuBReal)MUB_MU0 / 2) * r * (Nxy.x * (Mnew.x*Mnew.x - S.x*S.x) + Nxy.y * (Mnew.y*Mnew.y - S.y*S.y) + Nz * (Mnew.z*Mnew.z - S.z*S.z));
	else return ((cuBReal)MUB_MU0 / 2) * r * (Nxy.x * S.x*S.x + Nxy.y * S.y*S.y + Nz * S.z*S.z);
}

//Atom_Demag
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_DemagCUDA(int spin_index, cuReal3 Mnew)
{
	if (pAtom_Demag_Heff && pAtom_Demag_Heff->linear_size()) {

		cuVEC_VC<cuReal3>& M1 = *pM1;

		if (Mnew != cuReal3()) return -(cuBReal)MUB_MU0 * (*pAtom_Demag_Heff)[M1.cellidx_to_position(spin_index)] * (Mnew - M1[spin_index]);
		else return -(cuBReal)MUB_MU0 * (*pAtom_Demag_Heff)[M1.cellidx_to_position(spin_index)] * M1[spin_index];
	}
	else return 0.0;
}

//Atom_DipoleDipole
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_DipoleDipoleCUDA(int spin_index, cuReal3 Mnew)
{
	if (pAtom_Demag_Heff && pAtom_Demag_Heff->linear_size()) {

		cuVEC_VC<cuReal3>& M1 = *pM1;

		if (Mnew != cuReal3()) return -(cuBReal)MUB_MU0 * (*pAtom_Demag_Heff)[M1.cellidx_to_position(spin_index)] * (Mnew - M1[spin_index]);
		else return -(cuBReal)MUB_MU0 * (*pAtom_Demag_Heff)[M1.cellidx_to_position(spin_index)] * M1[spin_index];
	}
	else return 0.0;
}

//StrayField_AtomMesh
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_StrayField_AtomMeshCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuReal3 Hstray = cuReal3();

	if (pstrayField && pstrayField->linear_size()) {

		Hstray = (*pstrayField)[spin_index];
	}

	if (Mnew != cuReal3()) return -(cuBReal)MUB_MU0 * (Mnew - M1[spin_index]) * Hstray;
	else return -(cuBReal)MUB_MU0 * M1[spin_index] * Hstray;
}

//Atom_ExchangeCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_ExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal J = *pJ;
	update_parameters_mcoarse(spin_index, *pJ, J);

	if (Mnew != cuReal3()) return -J * ((Mnew.normalized() - M1[spin_index].normalized()) * M1.ngbr_dirsum(spin_index));
	else return -J * (M1[spin_index].normalized() * M1.ngbr_dirsum(spin_index));
}

//Atom_DMExchangeCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_DMExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal J = *pJ;
	cuBReal D = *pD;
	update_parameters_mcoarse(spin_index, *pJ, J, *pD, D);

	//local spin energy given:
	//1) exchange: -J * Sum_over_neighbors_j (Si . Sj), where Si, Sj are unit vectors
	//2) DM exchange: D * Sum_over_neighbors_j(rij . (Si x Sj))

	//Now anisotropic_ngbr_dirsum returns rij x Sj, and Si . (rij x Sj) = -Si. (Sj x rij) = -rij . (Si x Sj)
	//Thus compute: -D * (paMesh->M1[spin_index].normalized() * paMesh->M1.anisotropic_ngbr_dirsum(spin_index));

	if (Mnew != cuReal3()) return (Mnew.normalized() - M1[spin_index].normalized()) * (-J * M1.ngbr_dirsum(spin_index) - D * M1.anisotropic_ngbr_dirsum(spin_index));
	else return M1[spin_index].normalized() * (-J * M1.ngbr_dirsum(spin_index) - D * M1.anisotropic_ngbr_dirsum(spin_index));
}

//Atom_iDMExchangeCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_iDMExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal J = *pJ;
	cuBReal D = *pD;
	update_parameters_mcoarse(spin_index, *pJ, J, *pD, D);

	//local spin energy given:
	//1) exchange: -J * Sum_over_neighbors_j (Si . Sj), where Si, Sj are unit vectors
	//2) iDM exchange: D * Sum_over_neighbors_j((rij x z) . (Si x Sj))

	//Now zanisotropic_ngbr_dirsum returns (rij x z) x Sj, and Si . ((rij x z) x Sj) = -Si. (Sj x (rij x z)) = -(rij x z) . (Si x Sj)
	//Thus compute: -D * (paMesh->M1[spin_index].normalized() * paMesh->M1.zanisotropic_ngbr_dirsum(spin_index));

	if (Mnew != cuReal3()) return (Mnew.normalized() - M1[spin_index].normalized()) * (-J * M1.ngbr_dirsum(spin_index) - D * M1.zanisotropic_ngbr_dirsum(spin_index));
	else return M1[spin_index].normalized() * (-J * M1.ngbr_dirsum(spin_index) - D * M1.zanisotropic_ngbr_dirsum(spin_index));
}

//Atom_viDMExchangeCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_viDMExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal J = *pJ;
	cuBReal D = *pD;
	cuReal3 D_dir = *pD_dir;
	update_parameters_mcoarse(spin_index, *pJ, J, *pD, D, *pD_dir, D_dir);

	cuReal3 hexch_D_x = M1.xanisotropic_ngbr_dirsum(spin_index);
	cuReal3 hexch_D_y = M1.yanisotropic_ngbr_dirsum(spin_index);
	cuReal3 hexch_D_z = M1.zanisotropic_ngbr_dirsum(spin_index);
	cuReal3 hexch_D = D_dir.x * hexch_D_x + D_dir.y * hexch_D_y + D_dir.z * hexch_D_z;

	if (Mnew != cuReal3()) return (Mnew.normalized() - M1[spin_index].normalized()) * (-J * M1.ngbr_dirsum(spin_index) - D * hexch_D);
	else return M1[spin_index].normalized() * (-J * M1.ngbr_dirsum(spin_index) - D * hexch_D);
}

//Atom_SurfExchangeCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_SurfExchangeCUDA(int spin_index, cuReal3 Mnew)
{
	cuBReal energy_new = 0, energy_old = 0;

	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuSZ3 n = M1.n;
	cuReal3 h = M1.h;

	//if spin is on top surface then look at paMesh_Top
	if (spin_index / (n.x * n.y) == n.z - 1 && (paMesh_Top_size || pMeshFM_Top_size)) {

		if (!M1.is_empty(spin_index)) {

			int i = spin_index % n.x;
			int j = (spin_index / n.x) % n.y;
			bool cell_coupled = false;

			//check all meshes for coupling
			//1. coupling from other atomistic meshes
			for (int mesh_idx = 0; mesh_idx < paMesh_Top_size; mesh_idx++) {

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
				cuReal3 m_i = M1[spin_index].normalized();

				cuBReal dot_prod = m_i * m_j;
				energy_old = -Js * dot_prod;

				if (Mnew != cuReal3()) {

					cuReal3 mnew_i = Mnew.normalized();
					cuBReal dot_prod_new = mnew_i * m_j;
					energy_new = -Js * dot_prod_new;
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				cell_coupled = true;
				break;
			}

			if (!cell_coupled) {

				//2. coupling from micromagnetic meshes - ferromagnetic
				for (int mesh_idx = 0; mesh_idx < pMeshFM_Top_size; mesh_idx++) {

					cuVEC_VC<cuReal3>& M = *(pMeshFM_Top[mesh_idx].pM);

					cuRect tmeshRect = M.rect;

					//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
					cuReal3 cell_rel_pos = cuReal3(
						(i + 0.5) * h.x + M1.rect.s.x - tmeshRect.s.x,
						(j + 0.5) * h.y + M1.rect.s.y - tmeshRect.s.y,
						M.h.z / 2);

					//can't couple to an empty cell
					if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s) || M.is_empty(cell_rel_pos)) continue;

					cuBReal Js = *pJs;
					update_parameters_mcoarse(spin_index, *pJs, Js);

					//get magnetization value in top mesh cell to couple with
					cuReal3 m_j = M[cell_rel_pos].normalized();
					cuReal3 m_i = M1[spin_index].normalized();

					cuBReal dot_prod = m_i * m_j;
					energy_old = -Js * dot_prod;

					if (Mnew != cuReal3()) {

						cuReal3 mnew_i = Mnew.normalized();
						cuBReal dot_prod_new = mnew_i * m_j;
						energy_new = -Js * dot_prod_new;
					}

					//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
					cell_coupled = true;
					break;
				}
			}

			if (!cell_coupled) {

				//2. coupling from micromagnetic meshes - antiferromagnetic
				for (int mesh_idx = 0; mesh_idx < pMeshAFM_Top_size; mesh_idx++) {

					cuVEC_VC<cuReal3>& M = *(pMeshAFM_Top[mesh_idx].pM);
					cuVEC_VC<cuReal3>& M2 = *(pMeshAFM_Top[mesh_idx].pM2);

					cuRect tmeshRect = M.rect;

					//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
					cuReal3 cell_rel_pos = cuReal3(
						(i + 0.5) * h.x + M1.rect.s.x - tmeshRect.s.x,
						(j + 0.5) * h.y + M1.rect.s.y - tmeshRect.s.y,
						M.h.z / 2);

					//can't couple to an empty cell
					if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s) || M.is_empty(cell_rel_pos)) continue;

					cuBReal Js = *pJs;
					cuBReal Js2 = *pJs2;
					update_parameters_mcoarse(spin_index, *pJs, Js, *pJs2, Js2);

					//get magnetization value in top mesh cell to couple with
					cuReal3 m_j1 = M[cell_rel_pos].normalized();
					cuReal3 m_j2 = M2[cell_rel_pos].normalized();
					cuReal3 m_i = M1[spin_index].normalized();

					cuBReal dot_prod1 = m_i * m_j1;
					cuBReal dot_prod2 = m_i * m_j2;
					energy_old = -Js * dot_prod1;
					energy_old += -Js2 * dot_prod2;

					if (Mnew != cuReal3()) {

						cuReal3 mnew_i = Mnew.normalized();
						cuBReal dot_prod_new1 = mnew_i * m_j1;
						cuBReal dot_prod_new2 = mnew_i * m_j2;
						energy_new = -Js * dot_prod_new1;
						energy_new += -Js2 * dot_prod_new2;
					}

					//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
					break;
				}
			}
		}
	}

	//if spin is on bottom surface then look at paMesh_Top
	if (spin_index / (n.x * n.y) == 0 && (paMesh_Bot_size || pMeshFM_Bot_size)) {

		if (!M1.is_empty(spin_index)) {

			int i = spin_index % n.x;
			int j = (spin_index / n.x) % n.y;
			bool cell_coupled = false;

			cuBReal Js = *pJs;
			update_parameters_mcoarse(spin_index, *pJs, Js);

			//check all meshes for coupling
			//1. coupling from other atomistic meshes
			for (int mesh_idx = 0; mesh_idx < paMesh_Bot_size; mesh_idx++) {

				cuVEC_VC<cuReal3>& M1_Bot = *(paMesh_Bot[mesh_idx].pM1);

				//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
				cuReal3 cell_rel_pos = cuReal3(
					(i + 0.5) * h.x + M1.rect.s.x - M1_Bot.rect.s.x,
					(j + 0.5) * h.y + M1.rect.s.y - M1_Bot.rect.s.y,
					M1_Bot.rect.e.z - M1_Bot.rect.s.z - M1_Bot.h.z / 2);

				//can't couple to an empty cell
				if (!M1_Bot.rect.contains(cell_rel_pos + M1_Bot.rect.s) || M1_Bot.is_empty(cell_rel_pos)) continue;

				//get magnetization value in top mesh cell to couple with
				cuReal3 m_j = M1_Bot[cell_rel_pos].normalized();
				cuReal3 m_i = M1[spin_index].normalized();

				cuBReal dot_prod = m_i * m_j;

				energy_old += -Js * dot_prod;

				if (Mnew != cuReal3()) {

					cuReal3 mnew_i = Mnew.normalized();
					cuBReal dot_prod_new = mnew_i * m_j;
					energy_new += -Js * dot_prod_new;
				}
				
				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				cell_coupled = true;
				break;
			}

			if (!cell_coupled) {

				//2. coupling from micromagnetic meshes - ferromagnetic
				for (int mesh_idx = 0; mesh_idx < pMeshFM_Bot_size; mesh_idx++) {

					cuVEC_VC<cuReal3>& M = *pMeshFM_Bot[mesh_idx].pM;

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
					cuReal3 m_i = M1[spin_index].normalized();

					cuBReal dot_prod = m_i * m_j;

					energy_old += -Js * dot_prod;

					if (Mnew != cuReal3()) {

						cuReal3 mnew_i = Mnew.normalized();
						cuBReal dot_prod_new = mnew_i * m_j;
						energy_new += -Js * dot_prod_new;
					}

					//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
					cell_coupled = true;
					break;
				}
			}

			if (!cell_coupled) {

				//2. coupling from micromagnetic meshes - antiferromagnetic
				for (int mesh_idx = 0; mesh_idx < pMeshAFM_Bot_size; mesh_idx++) {

					cuVEC_VC<cuReal3>& M = *pMeshAFM_Bot[mesh_idx].pM;
					cuVEC_VC<cuReal3>& M2 = *pMeshAFM_Bot[mesh_idx].pM2;

					cuRect bmeshRect = M.rect;

					//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
					cuReal3 cell_rel_pos = cuReal3(
						(i + 0.5) * h.x + M1.rect.s.x - bmeshRect.s.x,
						(j + 0.5) * h.y + M1.rect.s.y - bmeshRect.s.y,
						bmeshRect.height() - M.h.z / 2);

					//can't couple to an empty cell
					if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s) || M.is_empty(cell_rel_pos)) continue;

					cuBReal Js2 = *pJs2;
					update_parameters_mcoarse(spin_index, *pJs2, Js2);

					//get magnetization value in top mesh cell to couple with
					cuReal3 m_j1 = M[cell_rel_pos].normalized();
					cuReal3 m_j2 = M2[cell_rel_pos].normalized();
					cuReal3 m_i = M1[spin_index].normalized();

					cuBReal dot_prod1 = m_i * m_j1;
					cuBReal dot_prod2 = m_i * m_j2;

					energy_old += -Js * dot_prod1;
					energy_old += -Js2 * dot_prod2;

					if (Mnew != cuReal3()) {

						cuReal3 mnew_i = Mnew.normalized();
						cuBReal dot_prod_new1 = mnew_i * m_j1;
						cuBReal dot_prod_new2 = mnew_i * m_j2;
						energy_new += -Js * dot_prod_new1;
						energy_new += -Js2 * dot_prod_new2;
					}

					//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
					break;
				}
			}
		}
	}

	if (Mnew != cuReal3()) return energy_new - energy_old;
	else return energy_old;
}

//Atom_ZeemanCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_ZeemanCUDA(int spin_index, cuReal3 Mnew, cuReal3& Ha)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuReal3 Hext = cuReal3();

	if (pHavec && pHavec->linear_size()) {

		Hext = (*pHavec)[spin_index];
	}
	else {

		cuBReal cHA = *pcHA;
		update_parameters_mcoarse(spin_index, *pcHA, cHA);

		Hext = cHA * Ha;
	}

	if (Mnew != cuReal3()) return -(cuBReal)MUB_MU0 * (Mnew - M1[spin_index]) * Hext;
	else return -(cuBReal)MUB_MU0 * M1[spin_index] * Hext;
}

//Atom_MOpticalCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_MOpticalCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal cHmo = *pcHmo;
	update_parameters_mcoarse(spin_index, *pcHmo, cHmo);

	if (Mnew != cuReal3()) return -(cuBReal)MUB * (Mnew - M1[spin_index]) * (cuBReal)MU0 * cuReal3(0, 0, cHmo);
	else return -(cuBReal)MUB * M1[spin_index] * (cuBReal)MU0 * cuReal3(0, 0, cHmo);
}

//Atom_AnisotropyCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_AnisotropyCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal K1 = *pK1;
	cuBReal K2 = *pK2;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	update_parameters_mcoarse(spin_index, *pK2, K2, *pmcanis_ea1, mcanis_ea1);

	//calculate m.ea dot product
	cuBReal dotprod = M1[spin_index].normalized() * mcanis_ea1;
	cuBReal dpsq = dotprod * dotprod;
	
	if (Mnew != cuReal3()) {

		cuBReal dotprod_new = Mnew.normalized() * mcanis_ea1;
		cuBReal dpsq_new = dotprod_new * dotprod_new;

		//Hamiltonian contribution as -Ku * (S * ea)^2, where S is the local spin direction
		return -K1 * (dotprod_new * dotprod_new - dotprod * dotprod) - K2 * (dpsq_new * dpsq_new - dpsq * dpsq);
	}
	else return -K1 * dotprod * dotprod - K2 * dpsq * dpsq;
}

//Atom_AnisotropyCubiCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_AnisotropyCubiCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal K1 = *pK1;
	cuBReal K2 = *pK2;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	cuReal3 mcanis_ea2 = *pmcanis_ea2;
	cuReal3 mcanis_ea3 = *pmcanis_ea3;
	update_parameters_mcoarse(spin_index, *pK1, K1, *pK2, K2, *pmcanis_ea1, mcanis_ea1, *pmcanis_ea2, mcanis_ea2, *pmcanis_ea3, mcanis_ea3);

	cuReal3 S = M1[spin_index].normalized();
	cuReal3 Snew = Mnew.normalized();

	//calculate m.ea1, m.ea2 and m.ea3 dot products
	cuBReal d1 = S * mcanis_ea1;
	cuBReal d2 = S * mcanis_ea2;
	cuBReal d3 = S * mcanis_ea3;

	if (Mnew != cuReal3()) {

		cuBReal d1_new = Snew * mcanis_ea1;
		cuBReal d2_new = Snew * mcanis_ea2;
		cuBReal d3_new = Snew * mcanis_ea3;

		//Hamiltonian contribution as K * (Sx^2*Sy^2 + Sx^2*Sz^2 + Sy^2*Sz^2), where S is the local spin direction (for easy axes coinciding with the xyz system)
		//This is equivalent to the form -K/2 * (Sx^4 + Sy^4 + Sz^4) - energy zero point differs but that's immaterial.
		//Also note the correct signs here for given easy axes (need to be careful, some publications have this wrong).
		return K1 * (d1_new*d1_new*d2_new*d2_new + d1_new*d1_new*d3_new*d3_new + d2_new*d2_new*d3_new*d3_new - d1*d1*d2*d2 - d1*d1*d3*d3 - d2*d2*d3*d3)
			+ K2 * (d1_new*d2_new*d3_new*d1_new*d2_new*d3_new - d1*d2*d3*d1*d2*d3);
	}
	else return K1 * (d1*d1*d2*d2 + d1*d1*d3*d3 + d2*d2*d3*d3) + K2 * d1*d2*d3*d1*d2*d3;
}

//Atom_AnisotropyBiaxialCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_AnisotropyBiaxialCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal K1 = *pK1;
	cuBReal K2 = *pK2;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	cuReal3 mcanis_ea2 = *pmcanis_ea2;
	cuReal3 mcanis_ea3 = *pmcanis_ea3;
	update_parameters_mcoarse(spin_index, *pK1, K1, *pK2, K2, *pmcanis_ea1, mcanis_ea1, *pmcanis_ea2, mcanis_ea2, *pmcanis_ea3, mcanis_ea3);

	//calculate m.ea1 dot product (uniaxial contribution)
	cuBReal u1 = M1[spin_index].normalized() * mcanis_ea1;

	//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
	cuBReal b1 = M1[spin_index].normalized() * mcanis_ea2;
	cuBReal b2 = M1[spin_index].normalized() * mcanis_ea3;

	if (Mnew != cuReal3()) {

		//calculate m.ea1 dot product (uniaxial contribution)
		cuBReal u1new = Mnew.normalized() * mcanis_ea1;

		//calculate m.ea2 and m.ea3 dot products (biaxial contribution)
		cuBReal b1new = Mnew.normalized() * mcanis_ea2;
		cuBReal b2new = Mnew.normalized() * mcanis_ea3;

		return (K1 * (1 - u1new*u1new) + K2 * b1new*b1new*b2new*b2new) - (K1 * (1 - u1 * u1) + K2 * b1*b1*b2*b2);
	}
	else return K1 * (1 - u1*u1) + K2 * b1*b1*b2*b2;
}

//Atom_AnisotropyTensorialCUDA
__device__ cuBReal ManagedAtom_MeshCUDA::Get_Atomistic_EnergyChange_SC_AnisotropyTensorialCUDA(int spin_index, cuReal3 Mnew)
{
	cuVEC_VC<cuReal3>& M1 = *pM1;

	cuBReal K1 = *pK1;
	cuBReal K2 = *pK2;
	cuBReal K3 = *pK3;
	cuReal3 mcanis_ea1 = *pmcanis_ea1;
	cuReal3 mcanis_ea2 = *pmcanis_ea2;
	cuReal3 mcanis_ea3 = *pmcanis_ea3;
	update_parameters_mcoarse(spin_index, *pK1, K1, *pK2, K2, *pK3, K3, *pmcanis_ea1, mcanis_ea1, *pmcanis_ea2, mcanis_ea2, *pmcanis_ea3, mcanis_ea3);

	auto Get_Energy = [&](cuBReal a, cuBReal b, cuBReal c) -> cuBReal {

		cuBReal energy_ = 0.0;

		for (int tidx = 0; tidx < pKt->linear_size(); tidx++) {

			//for each energy density term d*a^n1 b^n2 c^n3 we have an effective field contribution as:
			//(-d / mu0 mu_s) * [n1 * a^(n1-1) * b^n2 * c^n3 * mcanis_ea1 + n2 * a^n1) * b^(n2-1) * c^n3 * mcanis_ea2 + n3 * a^n1 * b^n2 * c^(n3-1) * mcanis_ea3] - for each n1, n2, n3 > 0

			cuBReal coeff;
			int order = (*pKt)[tidx].j + (*pKt)[tidx].k + (*pKt)[tidx].l;
			if (order == 2) coeff = K1 * (*pKt)[tidx].i;
			else if (order == 4) coeff = K2 * (*pKt)[tidx].i;
			else if (order == 6) coeff = K3 * (*pKt)[tidx].i;
			else coeff = (*pKt)[tidx].i;

			energy_ += coeff * pow(a, (*pKt)[tidx].j)*pow(b, (*pKt)[tidx].k)*pow(c, (*pKt)[tidx].l);
		}

		return energy_;
	};

	//calculate dot products
	cuBReal a = M1[spin_index].normalized() * mcanis_ea1;
	cuBReal b = M1[spin_index].normalized() * mcanis_ea2;
	cuBReal c = M1[spin_index].normalized() * mcanis_ea3;

	cuBReal energy_ = Get_Energy(a, b, c);

	if (Mnew != cuReal3()) {

		cuBReal anew = Mnew.normalized() * mcanis_ea1;
		cuBReal bnew = Mnew.normalized() * mcanis_ea2;
		cuBReal cnew = Mnew.normalized() * mcanis_ea3;

		cuBReal energynew_ = Get_Energy(anew, bnew, cnew);

		return energynew_ - energy_;
	}
	else return energy_;
}

#endif