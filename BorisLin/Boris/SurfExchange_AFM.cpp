#include "stdafx.h"
#include "SurfExchange_AFM.h"

#ifdef MODULE_COMPILATION_SURFEXCHANGE

#include "SuperMesh.h"

#include "Mesh_AntiFerromagnetic.h"
#include "MeshParamsControl.h"

#include "Atom_Mesh.h"
#include "Atom_MeshParamsControl.h"

/////////////////////////////////////////////////////////////////
//SurfExchange_AFM
//

SurfExchange_AFM::SurfExchange_AFM(Mesh *pMesh_) :
	Modules(),
	ProgramStateNames(this, {}, {})
{
	pMesh = pMesh_;

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

SurfExchange_AFM::~SurfExchange_AFM()
{
}

BError SurfExchange_AFM::Initialize(void)
{
	BError error(CLASS_STR(SurfExchange_AFM));

	//Need to identify all magnetic meshes participating in surface exchange coupling with this module:
	//1. Must be ferromagnetic or anti-ferromagnetic and have the SurfExchange module set -> this results in surface exchange coupling (or equivalently exchange bias - same formula, except in this case you should have J2 set to zero)
	//2. Must overlap in the x-y plane only with the mesh holding this module (either top or bottom) but not along z (i.e. mesh rectangles must not intersect)
	//3. No other magnetic meshes can be sandwiched in between - there could be other types of non-magnetic meshes in between of course (e.g. insulator, conductive layers etc).

	SuperMesh* pSMesh = pMesh->pSMesh;

	pMesh_Bot.clear();
	pMesh_Top.clear();
	paMesh_Bot.clear();
	paMesh_Top.clear();

	//---

	//lambda used to check condition 3
	auto check_candidate = [&](Rect xy_intersection, double z1, double z2) -> bool {

		//check all meshes to find a magnetic mesh with SurfExchange modules set, which intersects in the xy plane with xy_intersection, and has z coordinates between z1 and z2.
		//if found then current candidate not good
		for (int idx = 0; idx < pSMesh->size(); idx++) {

			//consider all meshes in turn - condition 1 first
			if ((*pSMesh)[idx]->MComputation_Enabled() && (*pSMesh)[idx]->IsModuleSet(MOD_SURFEXCHANGE)) {

				//get xy_meshrect (with z coordinates set to zero)
				Rect xy_meshrect = (*pSMesh)[idx]->GetMeshRect();
				xy_meshrect.s.z = 0; xy_meshrect.e.z = 0;

				if (xy_meshrect.intersects(xy_intersection)) {

					//intersection found. are the z coordinates in range also?
					if (IsGE((*pSMesh)[idx]->GetMeshRect().s.z, z1) && IsSE((*pSMesh)[idx]->GetMeshRect().e.z, z2)) {

						//new candidate found - note, new candidate cannot be the mesh being checked or the current candidate, so this is guranteed to be a better candidate
						return false;
					}
				}
			}
		}

		//no new candidate found - current candidate has been validated as the best one of its type (with given intersection)
		return true;
	};

	//---

	Rect meshRect = pMesh->GetMeshRect();

	Rect xy_meshRect = meshRect;
	xy_meshRect.s.z = 0; xy_meshRect.e.z = 0;

	for (int idx = 0; idx < pSMesh->size(); idx++) {

		//consider all meshes in turn - condition 1 first
		if ((*pSMesh)[idx]->MComputation_Enabled() && (*pSMesh)[idx]->IsModuleSet(MOD_SURFEXCHANGE)) {

			Rect candidate_meshRect = (*pSMesh)[idx]->GetMeshRect();

			//candidate mesh at the top
			if (IsGE(candidate_meshRect.s.z, meshRect.e.z)) {

				double z1 = meshRect.e.z;				//start z
				double z2 = candidate_meshRect.s.z;		//end z
				candidate_meshRect.s.z = 0; candidate_meshRect.e.z = 0;		//leave just the xy plane rect
				
				if (candidate_meshRect.intersects(xy_meshRect)) {

					//passes condition 2 - identified a candidate mesh at idx index. Does it pass condition 3?
					if (check_candidate(candidate_meshRect.get_intersection(xy_meshRect), z1, z2)) {

						if (!(*pSMesh)[idx]->is_atomistic()) pMesh_Top.push_back(dynamic_cast<Mesh*>((*pSMesh)[idx]));
						else paMesh_Top.push_back(dynamic_cast<Atom_Mesh*>((*pSMesh)[idx]));
					}
				}
			}

			//candidate mesh at the botttom
			else if (IsSE(candidate_meshRect.e.z, meshRect.s.z)) {

				double z1 = candidate_meshRect.e.z;		//start z
				double z2 = meshRect.s.z;				//end z
				candidate_meshRect.s.z = 0; candidate_meshRect.e.z = 0;		//leave just the xy plane rect

				if (candidate_meshRect.intersects(xy_meshRect)) {

					//passes condition 2 - identified a candidate mesh at idx index. Does it pass condition 3?
					if (check_candidate(candidate_meshRect.get_intersection(xy_meshRect), z1, z2)) {

						if (!(*pSMesh)[idx]->is_atomistic()) pMesh_Bot.push_back(dynamic_cast<Mesh*>((*pSMesh)[idx]));
						else paMesh_Bot.push_back(dynamic_cast<Atom_Mesh*>((*pSMesh)[idx]));
					}
				}
			}
		}
	}
	
	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect,
		(MOD_)pMesh->Get_Module_Heff_Display() == MOD_SURFEXCHANGE || pMesh->IsOutputDataSet_withRect(DATA_E_SURFEXCH) || pMesh->IsOutputDataSet(DATA_T_SURFEXCH),
		(MOD_)pMesh->Get_Module_Energy_Display() == MOD_SURFEXCHANGE || pMesh->IsOutputDataSet_withRect(DATA_E_SURFEXCH),
		true);
	if (!error) initialized = true;

	return error;
}

BError SurfExchange_AFM::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SurfExchange_AFM));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError SurfExchange_AFM::MakeCUDAModule(void)
{
	BError error(CLASS_STR(SurfExchange_AFM));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new SurfExchangeCUDA_AFM(pMesh->pMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double SurfExchange_AFM::UpdateField(void)
{
	double energy = 0;

	SZ3 n = pMesh->n;

	//zero module display VECs if needed, since contributions must be added into them to account for possiblility of 2 contributions (top and bottom)
	ZeroModuleVECs();

	if (pMesh_Top.size()) {

		//surface exchange coupling at the top
#pragma omp parallel for reduction(+:energy)
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

				//empty cell here ... next
				if (pMesh->M.is_empty(cell_idx)) continue;

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				pMesh->update_parameters_mcoarse(cell_idx, pMesh->Ms_AFM, Ms_AFM);

				DBL3 Hsurfexch1 = DBL3(), Hsurfexch2 = DBL3();
				double cell_energy1 = 0.0, cell_energy2 = 0.0;
				bool cell_coupled = false;

				//check all meshes for coupling
				//1 : coupling into this micromagnetic mesh from other micromagnetic meshes (FM or AFM)
				for (int mesh_idx = 0; mesh_idx < (int)pMesh_Top.size(); mesh_idx++) {

					if (!check_cell_coupling(pMesh->M, pMesh_Top[mesh_idx]->M,
						(i + 0.5) * pMesh->h.x, (j + 0.5) * pMesh->h.y, pMesh_Top[mesh_idx]->h.z / 2)) continue;

					DBL3 cell_rel_pos = DBL3(
						(i + 0.5) * pMesh->h.x + pMesh->M.rect.s.x - pMesh_Top[mesh_idx]->M.rect.s.x,
						(j + 0.5) * pMesh->h.y + pMesh->M.rect.s.y - pMesh_Top[mesh_idx]->M.rect.s.y,
						pMesh_Top[mesh_idx]->h.z / 2);

					if (pMesh_Top[mesh_idx]->GetMeshType() == MESH_FERROMAGNETIC) {

						//Surface exchange field from a ferromagnetic mesh

						//Top mesh sets J1 and J2 values
						double J1 = pMesh_Top[mesh_idx]->J1;
						double J2 = pMesh_Top[mesh_idx]->J2;
						pMesh_Top[mesh_idx]->update_parameters_atposition(cell_rel_pos, pMesh_Top[mesh_idx]->J1, J1, pMesh_Top[mesh_idx]->J2, J2);

						//get magnetization value in top mesh cell to couple with
						DBL3 m_j = normalize(pMesh_Top[mesh_idx]->M[cell_rel_pos]);
						DBL3 m_i1 = normalize(pMesh->M[cell_idx]);
						DBL3 m_i2 = normalize(pMesh->M2[cell_idx]);

						//total surface exchange field in coupling cells
						Hsurfexch1 = (m_j / (MU0 * Ms_AFM.i * pMesh->h.z)) * J1;
						Hsurfexch2 = (m_j / (MU0 * Ms_AFM.j * pMesh->h.z)) * J2;
						cell_energy1 = (-J1 * (m_i1 * m_j)) / pMesh->h.z;
						cell_energy2 = (-J2 * (m_i2 * m_j)) / pMesh->h.z;
					}
					else if (pMesh_Top[mesh_idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

						//Surface exchange field from an antiferromagnetic mesh

						//Top mesh sets J1 and J2 values
						double J1 = pMesh_Top[mesh_idx]->J1;
						double J2 = pMesh_Top[mesh_idx]->J2;
						pMesh_Top[mesh_idx]->update_parameters_atposition(cell_rel_pos, pMesh_Top[mesh_idx]->J1, J1, pMesh_Top[mesh_idx]->J2, J2);

						//get magnetization value in top mesh cell to couple with
						DBL3 m_j1 = normalize(pMesh_Top[mesh_idx]->M[cell_rel_pos]);
						DBL3 m_j2 = normalize(pMesh_Top[mesh_idx]->M2[cell_rel_pos]);
						DBL3 m_i1 = normalize(pMesh->M[cell_idx]);
						DBL3 m_i2 = normalize(pMesh->M2[cell_idx]);

						//total surface exchange field in coupling cells
						Hsurfexch1 = (m_j1 / (MU0 * Ms_AFM.i * pMesh->h.z)) * J1;
						Hsurfexch2 = (m_j2 / (MU0 * Ms_AFM.j * pMesh->h.z)) * J2;
						cell_energy1 = (-J1 * (m_i1 * m_j1)) / pMesh->h.z;
						cell_energy2 = (-J2 * (m_i2 * m_j2)) / pMesh->h.z;
					}
					
					//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
					cell_coupled = true;
					break;
				}

				if (!cell_coupled) {

					//2 : coupling into this micromagnetic mesh from atomistic meshes
					for (int mesh_idx = 0; mesh_idx < (int)paMesh_Top.size(); mesh_idx++) {

						VEC_VC<DBL3>& M1 = paMesh_Top[mesh_idx]->M1;

						//coupling rectangle in atomistic mesh in absolute coordinates
						Rect rect_c = Rect(
							DBL3(i * pMesh->h.x, j * pMesh->h.y, pMesh->M.rect.e.z),
							DBL3((i + 1) * pMesh->h.x, (j + 1) * pMesh->h.y, M1.h.z + pMesh->M.rect.e.z));
						rect_c += DBL3(pMesh->M.rect.s.x, pMesh->M.rect.s.y, 0.0);

						//cells box in atomistic mesh
						Box acells = M1.box_from_rect_min(rect_c);

						//find total "directed energy" contribution from atomistic mesh : i.e. sum all mj * Js contributions from atomistic moments in the coupling area at the interface
						DBL3 total_directed_coupling_energy1 = DBL3();
						DBL3 total_directed_coupling_energy2 = DBL3();
						for (int ai = acells.s.i; ai < acells.e.i; ai++) {
							for (int aj = acells.s.j; aj < acells.e.j; aj++) {

								int acell_idx = ai + aj * M1.n.x;

								if (M1.is_empty(acell_idx)) continue;

								//Js value from atomistic mesh
								double Js = paMesh_Top[mesh_idx]->Js;
								double Js2 = paMesh_Top[mesh_idx]->Js2;
								double mu_s = paMesh_Top[mesh_idx]->mu_s;
								paMesh_Top[mesh_idx]->update_parameters_mcoarse(acell_idx, paMesh_Top[mesh_idx]->Js, Js, paMesh_Top[mesh_idx]->Js2, Js2, paMesh_Top[mesh_idx]->mu_s, mu_s);

								total_directed_coupling_energy1 += M1[acell_idx] * Js / mu_s;
								total_directed_coupling_energy2 += M1[acell_idx] * Js2 / mu_s;
							}
						}

						//now obtain coupling field from atomistic mesh at micromagnetic cell
						Hsurfexch1 = (total_directed_coupling_energy1 / (pMesh->h.x * pMesh->h.y)) / (MU0 * Ms_AFM.i * pMesh->h.z);
						Hsurfexch2 = (total_directed_coupling_energy2 / (pMesh->h.x * pMesh->h.y)) / (MU0 * Ms_AFM.j * pMesh->h.z);
						cell_energy1 = -MU0 * pMesh->M[cell_idx] * Hsurfexch1;
						cell_energy2 = -MU0 * pMesh->M2[cell_idx] * Hsurfexch2;

						//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
						break;
					}
				}

				pMesh->Heff[cell_idx] += Hsurfexch1;
				pMesh->Heff2[cell_idx] += Hsurfexch2;

				if (Module_Heff.linear_size()) Module_Heff[cell_idx] += Hsurfexch1;
				if (Module_Heff2.linear_size()) Module_Heff2[cell_idx] += Hsurfexch2;
				if (Module_energy.linear_size()) Module_energy[cell_idx] += cell_energy1;
				if (Module_energy2.linear_size()) Module_energy2[cell_idx] += cell_energy2;

				energy += (cell_energy1 + cell_energy2) / 2;
			}
		}
	}

	if (pMesh_Bot.size()) {

		//surface exchange coupling at the bottom
#pragma omp parallel for reduction(+:energy)
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				int cell_idx = i + j * n.x;

				//empty cell here ... next
				if (pMesh->M.is_empty(cell_idx)) continue;

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				double J1 = pMesh->J1;
				double J2 = pMesh->J2;
				pMesh->update_parameters_mcoarse(cell_idx, pMesh->Ms_AFM, Ms_AFM, pMesh->J1, J1, pMesh->J2, J2);

				DBL3 Hsurfexch1 = DBL3(), Hsurfexch2 = DBL3();
				double cell_energy1 = 0.0, cell_energy2 = 0.0;
				bool cell_coupled = false;

				//check all meshes for coupling
				for (int mesh_idx = 0; mesh_idx < (int)pMesh_Bot.size(); mesh_idx++) {

					if (!check_cell_coupling(pMesh->M, pMesh_Bot[mesh_idx]->M,
						(i + 0.5) * pMesh->h.x, (j + 0.5) * pMesh->h.y, pMesh_Bot[mesh_idx]->meshRect.height() - pMesh_Bot[mesh_idx]->h.z / 2)) continue;

					DBL3 cell_rel_pos = DBL3(
						(i + 0.5) * pMesh->h.x + pMesh->M.rect.s.x - pMesh_Bot[mesh_idx]->M.rect.s.x,
						(j + 0.5) * pMesh->h.y + pMesh->M.rect.s.y - pMesh_Bot[mesh_idx]->M.rect.s.y,
						pMesh_Bot[mesh_idx]->meshRect.height() - pMesh_Bot[mesh_idx]->h.z / 2);

					if (pMesh_Bot[mesh_idx]->GetMeshType() == MESH_FERROMAGNETIC) {

						//Surface exchange field from a ferromagnetic mesh

						//get magnetization value in top mesh cell to couple with
						DBL3 m_j = normalize(pMesh_Bot[mesh_idx]->M[cell_rel_pos]);
						DBL3 m_i1 = normalize(pMesh->M[cell_idx]);
						DBL3 m_i2 = normalize(pMesh->M2[cell_idx]);

						//total surface exchange field in coupling cells
						Hsurfexch1 = (m_j / (MU0 * Ms_AFM.i * pMesh->h.z)) * J1;
						Hsurfexch2 = (m_j / (MU0 * Ms_AFM.j * pMesh->h.z)) * J2;
						cell_energy1 = (-J1 * (m_i1 * m_j)) / pMesh->h.z;
						cell_energy2 = (-J2 * (m_i2 * m_j)) / pMesh->h.z;
					}
					else if (pMesh_Bot[mesh_idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

						//Surface exchange field from a antiferromagnetic mesh

						//get magnetization value in top mesh cell to couple with
						DBL3 m_j1 = normalize(pMesh_Bot[mesh_idx]->M[cell_rel_pos]);
						DBL3 m_j2 = normalize(pMesh_Bot[mesh_idx]->M2[cell_rel_pos]);
						DBL3 m_i1 = normalize(pMesh->M[cell_idx]);
						DBL3 m_i2 = normalize(pMesh->M2[cell_idx]);

						//total surface exchange field in coupling cells
						Hsurfexch1 = (m_j1 / (MU0 * Ms_AFM.i * pMesh->h.z)) * J1;
						Hsurfexch2 = (m_j2 / (MU0 * Ms_AFM.j * pMesh->h.z)) * J2;
						cell_energy1 = (-J1 * (m_i1 * m_j1)) / pMesh->h.z;
						cell_energy2 = (-J2 * (m_i2 * m_j2)) / pMesh->h.z;
					}

					//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
					cell_coupled = true;
					break;
				}

				if (!cell_coupled) {

					//2 : coupling into this micromagnetic mesh from atomistic meshes
					for (int mesh_idx = 0; mesh_idx < (int)paMesh_Bot.size(); mesh_idx++) {

						VEC_VC<DBL3>& M1 = paMesh_Bot[mesh_idx]->M1;

						//coupling rectangle in atomistic mesh in absolute coordinates
						Rect rect_c = Rect(
							DBL3(i * pMesh->h.x, j * pMesh->h.y, M1.rect.e.z - M1.h.z),
							DBL3((i + 1) * pMesh->h.x, (j + 1) * pMesh->h.y, M1.rect.e.z));
						rect_c += DBL3(pMesh->M.rect.s.x, pMesh->M.rect.s.y, 0.0);

						//cells box in atomistic mesh
						Box acells = M1.box_from_rect_min(rect_c);

						//find total "directed energy" contribution from atomistic mesh : i.e. sum all mj * Js contributions from atomistic moments in the coupling area at the interface
						//NOTE : at atomistic/micromagnetic coupling, it's the atomistic mesh which sets coupling constant, not the top mesh
						DBL3 total_directed_coupling_energy1 = DBL3();
						DBL3 total_directed_coupling_energy2 = DBL3();
						for (int ai = acells.s.i; ai < acells.e.i; ai++) {
							for (int aj = acells.s.j; aj < acells.e.j; aj++) {

								int acell_idx = ai + aj * M1.n.x + (M1.n.z - 1) * M1.n.x * M1.n.y;

								if (M1.is_empty(acell_idx)) continue;

								//Js value from atomistic mesh
								double Js = paMesh_Bot[mesh_idx]->Js;
								double Js2 = paMesh_Bot[mesh_idx]->Js2;
								double mu_s = paMesh_Bot[mesh_idx]->mu_s;
								paMesh_Bot[mesh_idx]->update_parameters_mcoarse(acell_idx, paMesh_Bot[mesh_idx]->Js, Js, paMesh_Bot[mesh_idx]->Js2, Js2, paMesh_Bot[mesh_idx]->mu_s, mu_s);

								total_directed_coupling_energy1 += M1[acell_idx] * Js / mu_s;
								total_directed_coupling_energy2 += M1[acell_idx] * Js2 / mu_s;
							}
						}

						//now obtain coupling field from atomistic mesh at micromagnetic cell
						Hsurfexch1 = (total_directed_coupling_energy1 / (pMesh->h.x * pMesh->h.y)) / (MU0 * Ms_AFM.i * pMesh->h.z);
						Hsurfexch2 = (total_directed_coupling_energy2 / (pMesh->h.x * pMesh->h.y)) / (MU0 * Ms_AFM.j * pMesh->h.z);
						cell_energy1 = -MU0 * pMesh->M[cell_idx] * Hsurfexch1;
						cell_energy2 = -MU0 * pMesh->M2[cell_idx] * Hsurfexch2;

						//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
						break;
					}
				}

				pMesh->Heff[cell_idx] += Hsurfexch1;
				pMesh->Heff2[cell_idx] += Hsurfexch2;

				if (Module_Heff.linear_size()) Module_Heff[cell_idx] += Hsurfexch1;
				if (Module_Heff2.linear_size()) Module_Heff2[cell_idx] += Hsurfexch2;
				if (Module_energy.linear_size()) Module_energy[cell_idx] += cell_energy1;
				if (Module_energy2.linear_size()) Module_energy2[cell_idx] += cell_energy2;

				energy += (cell_energy1 + cell_energy2) / 2;
			}
		}
	}

	energy /= pMesh->M.get_nonempty_cells();
	this->energy = energy;

	return this->energy;
}

//-------------------Energy methods

//AFM
DBL2 SurfExchange_AFM::Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B)
{
	DBL2 energy_new = DBL2(), energy_old = DBL2();

	SZ3 n = pMesh->n;

	//if spin is on top surface then look at paMesh_Top
	if (spin_index / (n.x * n.y) == n.z - 1 && pMesh_Top.size()) {

		if (!pMesh->M.is_empty(spin_index)) {

			int i = spin_index % n.x;
			int j = (spin_index / n.x) % n.y;

			DBL2 Ms_AFM = pMesh->Ms_AFM;
			pMesh->update_parameters_mcoarse(spin_index, pMesh->Ms_AFM, Ms_AFM);

			bool cell_coupled = false;

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < (int)pMesh_Top.size(); mesh_idx++) {

				if (!check_cell_coupling(pMesh->M, pMesh_Top[mesh_idx]->M,
					(i + 0.5) * pMesh->h.x, (j + 0.5) * pMesh->h.y, pMesh_Top[mesh_idx]->h.z / 2)) continue;

				DBL3 cell_rel_pos = DBL3(
					(i + 0.5) * pMesh->h.x + pMesh->M.rect.s.x - pMesh_Top[mesh_idx]->M.rect.s.x,
					(j + 0.5) * pMesh->h.y + pMesh->M.rect.s.y - pMesh_Top[mesh_idx]->M.rect.s.y,
					pMesh_Top[mesh_idx]->h.z / 2);

				if (pMesh_Top[mesh_idx]->GetMeshType() == MESH_FERROMAGNETIC) {

					//Surface exchange field from a ferromagnetic mesh

					//Top mesh sets J1 and J2 values
					double J1 = pMesh_Top[mesh_idx]->J1;
					double J2 = pMesh_Top[mesh_idx]->J2;
					pMesh_Top[mesh_idx]->update_parameters_atposition(cell_rel_pos, pMesh_Top[mesh_idx]->J1, J1, pMesh_Top[mesh_idx]->J2, J2);

					//get magnetization value in top mesh cell to couple with
					DBL3 m_j = normalize(pMesh_Top[mesh_idx]->M[cell_rel_pos]);
					DBL3 m_i1 = normalize(pMesh->M[spin_index]);
					DBL3 m_i2 = normalize(pMesh->M2[spin_index]);

					energy_old.i = (-J1 * (m_i1 * m_j)) / pMesh->h.z;
					energy_old.j = (-J2 * (m_i2 * m_j)) / pMesh->h.z;

					if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

						DBL3 mnew_i1 = normalize(Mnew_A);
						DBL3 mnew_i2 = normalize(Mnew_B);

						energy_new.i = (-J1 * (mnew_i1 * m_j)) / pMesh->h.z;
						energy_new.j = (-J2 * (mnew_i2 * m_j)) / pMesh->h.z;
					}
				}
				else if (pMesh_Top[mesh_idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					//Surface exchange field from an antiferromagnetic mesh

					//Top mesh sets J1 and J2 values
					double J1 = pMesh_Top[mesh_idx]->J1;
					double J2 = pMesh_Top[mesh_idx]->J2;
					pMesh_Top[mesh_idx]->update_parameters_atposition(cell_rel_pos, pMesh_Top[mesh_idx]->J1, J1, pMesh_Top[mesh_idx]->J2, J2);

					//get magnetization value in top mesh cell to couple with
					DBL3 m_j1 = normalize(pMesh_Top[mesh_idx]->M[cell_rel_pos]);
					DBL3 m_j2 = normalize(pMesh_Top[mesh_idx]->M2[cell_rel_pos]);
					DBL3 m_i1 = normalize(pMesh->M[spin_index]);
					DBL3 m_i2 = normalize(pMesh->M2[spin_index]);

					energy_old.i = (-J1 * (m_i1 * m_j1)) / pMesh->h.z;
					energy_old.j = (-J2 * (m_i2 * m_j2)) / pMesh->h.z;

					if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

						DBL3 mnew_i1 = normalize(Mnew_A);
						DBL3 mnew_i2 = normalize(Mnew_B);

						energy_new.i = (-J1 * (mnew_i1 * m_j1)) / pMesh->h.z;
						energy_new.j = (-J2 * (mnew_i2 * m_j2)) / pMesh->h.z;
					}
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				cell_coupled = true;
				break;
			}

			if (!cell_coupled) {

				//2 : coupling into this micromagnetic mesh from atomistic meshes
				for (int mesh_idx = 0; mesh_idx < (int)paMesh_Top.size(); mesh_idx++) {

					VEC_VC<DBL3>& M1 = paMesh_Top[mesh_idx]->M1;

					//coupling rectangle in atomistic mesh in absolute coordinates
					Rect rect_c = Rect(
						DBL3(i * pMesh->h.x, j * pMesh->h.y, pMesh->M.rect.e.z),
						DBL3((i + 1) * pMesh->h.x, (j + 1) * pMesh->h.y, M1.h.z + pMesh->M.rect.e.z));
					rect_c += DBL3(pMesh->M.rect.s.x, pMesh->M.rect.s.y, 0.0);

					//cells box in atomistic mesh
					Box acells = M1.box_from_rect_min(rect_c);

					//find total "directed energy" contribution from atomistic mesh : i.e. sum all mj * Js contributions from atomistic moments in the coupling area at the interface
					DBL3 total_directed_coupling_energy1 = DBL3();
					DBL3 total_directed_coupling_energy2 = DBL3();
					for (int ai = acells.s.i; ai < acells.e.i; ai++) {
						for (int aj = acells.s.j; aj < acells.e.j; aj++) {

							int acell_idx = ai + aj * M1.n.x;

							if (M1.is_empty(acell_idx)) continue;

							//Js value from atomistic mesh
							double Js = paMesh_Top[mesh_idx]->Js;
							double Js2 = paMesh_Top[mesh_idx]->Js2;
							double mu_s = paMesh_Top[mesh_idx]->mu_s;
							paMesh_Top[mesh_idx]->update_parameters_mcoarse(acell_idx, paMesh_Top[mesh_idx]->Js, Js, paMesh_Top[mesh_idx]->Js2, Js2, paMesh_Top[mesh_idx]->mu_s, mu_s);

							total_directed_coupling_energy1 += M1[acell_idx] * Js / mu_s;
							total_directed_coupling_energy2 += M1[acell_idx] * Js2 / mu_s;
						}
					}

					//now obtain coupling field from atomistic mesh at micromagnetic cell
					DBL3 Hsurfexch1 = (total_directed_coupling_energy1 / (pMesh->h.x * pMesh->h.y)) / (MU0 * Ms_AFM.i * pMesh->h.z);
					DBL3 Hsurfexch2 = (total_directed_coupling_energy2 / (pMesh->h.x * pMesh->h.y)) / (MU0 * Ms_AFM.j * pMesh->h.z);
					energy_old.i = -MU0 * pMesh->M[spin_index] * Hsurfexch1;
					energy_old.j = -MU0 * pMesh->M2[spin_index] * Hsurfexch2;

					if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

						energy_new.i = -MU0 * Mnew_A * Hsurfexch1;
						energy_new.j = -MU0 * Mnew_B * Hsurfexch2;
					}

					//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
					break;
				}
			}
		}
	}

	if (spin_index / (n.x * n.y) == 0 && pMesh_Bot.size()) {

		//surface exchange coupling at the bottom
		
		if (!pMesh->M.is_empty(spin_index)) {

			int i = spin_index % n.x;
			int j = (spin_index / n.x) % n.y;

			DBL2 Ms_AFM = pMesh->Ms_AFM;
			double J1 = pMesh->J1;
			double J2 = pMesh->J2;
			pMesh->update_parameters_mcoarse(spin_index, pMesh->Ms_AFM, Ms_AFM, pMesh->J1, J1, pMesh->J2, J2);

			bool cell_coupled = false;

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < (int)pMesh_Bot.size(); mesh_idx++) {

				if (!check_cell_coupling(pMesh->M, pMesh_Bot[mesh_idx]->M,
					(i + 0.5) * pMesh->h.x, (j + 0.5) * pMesh->h.y, pMesh_Bot[mesh_idx]->meshRect.height() - pMesh_Bot[mesh_idx]->h.z / 2)) continue;

				DBL3 cell_rel_pos = DBL3(
					(i + 0.5) * pMesh->h.x + pMesh->M.rect.s.x - pMesh_Bot[mesh_idx]->M.rect.s.x,
					(j + 0.5) * pMesh->h.y + pMesh->M.rect.s.y - pMesh_Bot[mesh_idx]->M.rect.s.y,
					pMesh_Bot[mesh_idx]->meshRect.height() - pMesh_Bot[mesh_idx]->h.z / 2);

				if (pMesh_Bot[mesh_idx]->GetMeshType() == MESH_FERROMAGNETIC) {

					//Surface exchange field from a ferromagnetic mesh

					//get magnetization value in top mesh cell to couple with
					DBL3 m_j = normalize(pMesh_Bot[mesh_idx]->M[cell_rel_pos]);
					DBL3 m_i1 = normalize(pMesh->M[spin_index]);
					DBL3 m_i2 = normalize(pMesh->M2[spin_index]);

					energy_old.i += (-J1 * (m_i1 * m_j)) / pMesh->h.z;
					energy_old.j += (-J2 * (m_i2 * m_j)) / pMesh->h.z;

					if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

						DBL3 mnew_i1 = normalize(Mnew_A);
						DBL3 mnew_i2 = normalize(Mnew_B);

						energy_new.i += (-J1 * (mnew_i1 * m_j)) / pMesh->h.z;
						energy_new.j += (-J2 * (mnew_i2 * m_j)) / pMesh->h.z;
					}
				}
				else if (pMesh_Bot[mesh_idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					//Surface exchange field from a antiferromagnetic mesh

					//get magnetization value in top mesh cell to couple with
					DBL3 m_j1 = normalize(pMesh_Bot[mesh_idx]->M[cell_rel_pos]);
					DBL3 m_j2 = normalize(pMesh_Bot[mesh_idx]->M2[cell_rel_pos]);
					DBL3 m_i1 = normalize(pMesh->M[spin_index]);
					DBL3 m_i2 = normalize(pMesh->M2[spin_index]);

					energy_old.i += (-J1 * (m_i1 * m_j1)) / pMesh->h.z;
					energy_old.j += (-J2 * (m_i2 * m_j2)) / pMesh->h.z;

					if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

						DBL3 mnew_i1 = normalize(Mnew_A);
						DBL3 mnew_i2 = normalize(Mnew_B);

						energy_new.i += (-J1 * (mnew_i1 * m_j1)) / pMesh->h.z;
						energy_new.j += (-J2 * (mnew_i2 * m_j2)) / pMesh->h.z;
					}
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				cell_coupled = true;
				break;
			}

			if (!cell_coupled) {

				//2 : coupling into this micromagnetic mesh from atomistic meshes
				for (int mesh_idx = 0; mesh_idx < (int)paMesh_Bot.size(); mesh_idx++) {

					VEC_VC<DBL3>& M1 = paMesh_Bot[mesh_idx]->M1;

					//coupling rectangle in atomistic mesh in absolute coordinates
					Rect rect_c = Rect(
						DBL3(i * pMesh->h.x, j * pMesh->h.y, M1.rect.e.z - M1.h.z),
						DBL3((i + 1) * pMesh->h.x, (j + 1) * pMesh->h.y, M1.rect.e.z));
					rect_c += DBL3(pMesh->M.rect.s.x, pMesh->M.rect.s.y, 0.0);

					//cells box in atomistic mesh
					Box acells = M1.box_from_rect_min(rect_c);

					//find total "directed energy" contribution from atomistic mesh : i.e. sum all mj * Js contributions from atomistic moments in the coupling area at the interface
					//NOTE : at atomistic/micromagnetic coupling, it's the atomistic mesh which sets coupling constant, not the top mesh
					DBL3 total_directed_coupling_energy1 = DBL3();
					DBL3 total_directed_coupling_energy2 = DBL3();
					for (int ai = acells.s.i; ai < acells.e.i; ai++) {
						for (int aj = acells.s.j; aj < acells.e.j; aj++) {

							int acell_idx = ai + aj * M1.n.x + (M1.n.z - 1) * M1.n.x * M1.n.y;

							if (M1.is_empty(acell_idx)) continue;

							//Js value from atomistic mesh
							double Js = paMesh_Bot[mesh_idx]->Js;
							double Js2 = paMesh_Bot[mesh_idx]->Js2;
							double mu_s = paMesh_Bot[mesh_idx]->mu_s;
							paMesh_Bot[mesh_idx]->update_parameters_mcoarse(acell_idx, paMesh_Bot[mesh_idx]->Js, Js, paMesh_Bot[mesh_idx]->Js2, Js2, paMesh_Bot[mesh_idx]->mu_s, mu_s);

							total_directed_coupling_energy1 += M1[acell_idx] * Js / mu_s;
							total_directed_coupling_energy2 += M1[acell_idx] * Js2 / mu_s;
						}
					}

					//now obtain coupling field from atomistic mesh at micromagnetic cell
					DBL3 Hsurfexch1 = (total_directed_coupling_energy1 / (pMesh->h.x * pMesh->h.y)) / (MU0 * Ms_AFM.i * pMesh->h.z);
					DBL3 Hsurfexch2 = (total_directed_coupling_energy2 / (pMesh->h.x * pMesh->h.y)) / (MU0 * Ms_AFM.j * pMesh->h.z);
					energy_old.i += -MU0 * pMesh->M[spin_index] * Hsurfexch1;
					energy_old.j += -MU0 * pMesh->M2[spin_index] * Hsurfexch2;

					if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

						energy_new.i += -MU0 * Mnew_A * Hsurfexch1;
						energy_new.j += -MU0 * Mnew_B * Hsurfexch2;
					}

					//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
					break;
				}
			}
		}
	}

	if (Mnew_A != DBL3() && Mnew_B != DBL3()) return pMesh->h.dim() * (energy_new - energy_old);
	else return pMesh->h.dim() * energy_old;
}

//-------------------Torque methods

DBL3 SurfExchange_AFM::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return pModuleCUDA->GetTorque(avRect);
#endif

	return CalculateTorque(pMesh->M, avRect);
}

#endif