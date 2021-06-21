#include "stdafx.h"
#include "SurfExchange_AFM.h"

#ifdef MODULE_COMPILATION_SURFEXCHANGE

#include "Mesh_AntiFerromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

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

	coupled_cells = 0;

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
		if ((*pSMesh)[idx]->MComputation_Enabled() && (*pSMesh)[idx]->IsModuleSet(MOD_SURFEXCHANGE) && !(*pSMesh)[idx]->is_atomistic()) {

			Rect candidate_meshRect = (*pSMesh)[idx]->GetMeshRect();

			//candidate mesh at the top
			if (IsGE(candidate_meshRect.s.z, meshRect.e.z)) {

				double z1 = meshRect.e.z;				//start z
				double z2 = candidate_meshRect.s.z;		//end z
				candidate_meshRect.s.z = 0; candidate_meshRect.e.z = 0;		//leave just the xy plane rect
				
				if (candidate_meshRect.intersects(xy_meshRect)) {

					//passes condition 2 - identified a candidate mesh at idx index. Does it pass condition 3?
					if (check_candidate(candidate_meshRect.get_intersection(xy_meshRect), z1, z2)) {

						pMesh_Top.push_back(dynamic_cast<Mesh*>((*pSMesh)[idx]));
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

						pMesh_Bot.push_back(dynamic_cast<Mesh*>((*pSMesh)[idx]));
					}
				}
			}
		}
	}
	
	//count number of coupled cells (either top or bottom) in this mesh

	SZ3 n = pMesh->n;
	
	if (pMesh_Top.size()) {

		//surface exchange coupling at the top
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

				//empty cell here ... next
				if (pMesh->M.is_empty(cell_idx)) continue;

				//check all meshes for coupling
				for (int mesh_idx = 0; mesh_idx < (int)pMesh_Top.size(); mesh_idx++) {

					Rect tmeshRect = pMesh_Top[mesh_idx]->GetMeshRect();

					//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
					DBL3 cell_rel_pos = DBL3(
						(i + 0.5) * pMesh->h.x + meshRect.s.x - tmeshRect.s.x, 
						(j + 0.5) * pMesh->h.y + meshRect.s.y - tmeshRect.s.y, 
						pMesh_Top[mesh_idx]->h.z / 2);

					//can't couple to an empty cell
					if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s) || pMesh_Top[mesh_idx]->M.is_empty(cell_rel_pos)) continue;

					//if we are here then the cell in this mesh at cell_idx has something to couple to so count it : it will contribute to the surface exchange energy density
					coupled_cells++;
				}
			}
		}
	}

	if (pMesh_Bot.size()) {

		//surface exchange coupling at the bottom
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				int cell_idx = i + j * n.x;

				//empty cell here ... next
				if (pMesh->M.is_empty(cell_idx)) continue;

				//check all meshes for coupling
				for (int mesh_idx = 0; mesh_idx < (int)pMesh_Bot.size(); mesh_idx++) {

					Rect bmeshRect = pMesh_Bot[mesh_idx]->GetMeshRect();

					//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
					DBL3 cell_rel_pos = DBL3(
						(i + 0.5) * pMesh->h.x + meshRect.s.x - bmeshRect.s.x,
						(j + 0.5) * pMesh->h.y + meshRect.s.y - bmeshRect.s.y,
						pMesh_Bot[mesh_idx]->meshRect.e.z - pMesh_Bot[mesh_idx]->meshRect.s.z - (pMesh_Bot[mesh_idx]->h.z / 2));

					//can't couple to an empty cell
					if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s) || pMesh_Bot[mesh_idx]->M.is_empty(cell_rel_pos)) continue;

					//if we are here then the cell in this mesh at cell_idx has something to couple to so count it : it will contribute to the surface exchange energy density
					coupled_cells++;
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

	//thickness of layer - SurfExchange applies for layers in the xy plane
	double thickness = pMesh->meshRect.e.z - pMesh->meshRect.s.z;

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

				//check all meshes for coupling
				for (int mesh_idx = 0; mesh_idx < (int)pMesh_Top.size(); mesh_idx++) {

					Rect tmeshRect = pMesh_Top[mesh_idx]->GetMeshRect();

					//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
					DBL3 cell_rel_pos = DBL3(
						(i + 0.5) * pMesh->h.x + pMesh->meshRect.s.x - tmeshRect.s.x,
						(j + 0.5) * pMesh->h.y + pMesh->meshRect.s.y - tmeshRect.s.y,
						pMesh_Top[mesh_idx]->h.z / 2);

					//can't couple to an empty cell
					if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s) || pMesh_Top[mesh_idx]->M.is_empty(cell_rel_pos)) continue;

					DBL3 Hsurfexh1, Hsurfexh2;
					double cell_energy1 = 0.0, cell_energy2 = 0.0;

					if (pMesh_Top[mesh_idx]->GetMeshType() == MESH_FERROMAGNETIC) {

						//Surface exchange field from a ferromagnetic mesh

						//Top mesh sets J1 and J2 values
						double J1 = pMesh_Top[mesh_idx]->J1;
						double J2 = pMesh_Top[mesh_idx]->J2;
						pMesh_Top[mesh_idx]->update_parameters_atposition(cell_rel_pos, pMesh_Top[mesh_idx]->J1, J1, pMesh_Top[mesh_idx]->J2, J2);

						//get magnetization value in top mesh cell to couple with
						DBL3 m_j = pMesh_Top[mesh_idx]->M[cell_rel_pos].normalized();
						DBL3 m_i1 = pMesh->M[cell_idx] / Ms_AFM.i;
						DBL3 m_i2 = pMesh->M2[cell_idx] / Ms_AFM.j;

						//total surface exchange field in coupling cells
						Hsurfexh1 = (m_j / (MU0 * Ms_AFM.i * thickness)) * J1;
						Hsurfexh2 = (m_j / (MU0 * Ms_AFM.j * thickness)) * J2;
						cell_energy1 = (-J1 * (m_i1 * m_j)) / thickness;
						cell_energy2 = (-J2 * (m_i2 * m_j)) / thickness;
					}
					else if (pMesh_Top[mesh_idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

						//Surface exchange field from an antiferromagnetic mesh

						//Top mesh sets J1 and J2 values
						double J1 = pMesh_Top[mesh_idx]->J1;
						double J2 = pMesh_Top[mesh_idx]->J2;
						pMesh_Top[mesh_idx]->update_parameters_atposition(cell_rel_pos, pMesh_Top[mesh_idx]->J1, J1, pMesh_Top[mesh_idx]->J2, J2);

						//get magnetization value in top mesh cell to couple with
						DBL3 m_j1 = pMesh_Top[mesh_idx]->M[cell_rel_pos].normalized();
						DBL3 m_j2 = pMesh_Top[mesh_idx]->M2[cell_rel_pos].normalized();
						DBL3 m_i1 = pMesh->M[cell_idx] / Ms_AFM.i;
						DBL3 m_i2 = pMesh->M2[cell_idx] / Ms_AFM.j;

						//total surface exchange field in coupling cells
						Hsurfexh1 = (m_j1 / (MU0 * Ms_AFM.i * thickness)) * J1;
						Hsurfexh2 = (m_j2 / (MU0 * Ms_AFM.j * thickness)) * J2;
						cell_energy1 = (-J1 * (m_i1 * m_j1)) / thickness;
						cell_energy2 = (-J2 * (m_i2 * m_j2)) / thickness;
					}

					//couple all cells through the layer thickness : the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
					for (int k = 0; k < n.z; k++) {

						int idx = i + j * n.x + k * n.x*n.y;
						pMesh->Heff[idx] += Hsurfexh1;
						pMesh->Heff2[idx] += Hsurfexh2;

						if (Module_Heff.linear_size()) Module_Heff[idx] += Hsurfexh1;
						if (Module_Heff2.linear_size()) Module_Heff2[idx] += Hsurfexh2;
						if (Module_energy.linear_size()) Module_energy[idx] += cell_energy1;
						if (Module_energy2.linear_size()) Module_energy2[idx] += cell_energy2;
					}

					energy += (cell_energy1 + cell_energy2) / 2;
					
					//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
					break;
				}
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

				//check all meshes for coupling
				for (int mesh_idx = 0; mesh_idx < (int)pMesh_Bot.size(); mesh_idx++) {

					Rect bmeshRect = pMesh_Bot[mesh_idx]->GetMeshRect();

					//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
					DBL3 cell_rel_pos = DBL3(
						(i + 0.5) * pMesh->h.x + pMesh->meshRect.s.x - bmeshRect.s.x,
						(j + 0.5) * pMesh->h.y + pMesh->meshRect.s.y - bmeshRect.s.y,
						pMesh_Bot[mesh_idx]->meshRect.e.z - pMesh_Bot[mesh_idx]->meshRect.s.z - (pMesh_Bot[mesh_idx]->h.z / 2));

					//can't couple to an empty cell
					if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s) || pMesh_Bot[mesh_idx]->M.is_empty(cell_rel_pos)) continue;

					DBL3 Hsurfexh1, Hsurfexh2;
					double cell_energy1 = 0.0, cell_energy2 = 0.0;

					if (pMesh_Bot[mesh_idx]->GetMeshType() == MESH_FERROMAGNETIC) {

						//Surface exchange field from a ferromagnetic mesh

						//get magnetization value in top mesh cell to couple with
						DBL3 m_j = pMesh_Bot[mesh_idx]->M[cell_rel_pos].normalized();
						DBL3 m_i1 = pMesh->M[cell_idx] / Ms_AFM.i;
						DBL3 m_i2 = pMesh->M2[cell_idx] / Ms_AFM.j;

						//total surface exchange field in coupling cells
						Hsurfexh1 = (m_j / (MU0 * Ms_AFM.i * thickness)) * J1;
						Hsurfexh2 = (m_j / (MU0 * Ms_AFM.j * thickness)) * J2;
						cell_energy1 = (-J1 * (m_i1 * m_j)) / thickness;
						cell_energy2 = (-J2 * (m_i2 * m_j)) / thickness;
					}
					else if (pMesh_Bot[mesh_idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

						//Surface exchange field from a antiferromagnetic mesh

						//get magnetization value in top mesh cell to couple with
						DBL3 m_j1 = pMesh_Bot[mesh_idx]->M[cell_rel_pos].normalized();
						DBL3 m_j2 = pMesh_Bot[mesh_idx]->M2[cell_rel_pos].normalized();
						DBL3 m_i1 = pMesh->M[cell_idx] / Ms_AFM.i;
						DBL3 m_i2 = pMesh->M2[cell_idx] / Ms_AFM.j;

						//total surface exchange field in coupling cells
						Hsurfexh1 = (m_j1 / (MU0 * Ms_AFM.i * thickness)) * J1;
						Hsurfexh2 = (m_j2 / (MU0 * Ms_AFM.j * thickness)) * J2;
						cell_energy1 = (-J1 * (m_i1 * m_j1)) / thickness;
						cell_energy2 = (-J2 * (m_i2 * m_j2)) / thickness;
					}

					//couple all cells through the layer thickness : the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
					for (int k = 0; k < n.z; k++) {

						int idx = i + j * n.x + k * n.x*n.y;
						pMesh->Heff[idx] += Hsurfexh1;
						pMesh->Heff2[idx] += Hsurfexh2;

						if (Module_Heff.linear_size()) Module_Heff[idx] += Hsurfexh1;
						if (Module_Heff2.linear_size()) Module_Heff2[idx] += Hsurfexh2;
						if (Module_energy.linear_size()) Module_energy[idx] += cell_energy1;
						if (Module_energy2.linear_size()) Module_energy2[idx] += cell_energy2;
					}

					energy += (cell_energy1 + cell_energy2) / 2;

					//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
					break;
				}
			}
		}
	}

	if (coupled_cells) energy /= coupled_cells;
	else energy = 0.0;

	this->energy = energy;

	return this->energy;
}

//-------------------Torque methods

DBL3 SurfExchange_AFM::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return reinterpret_cast<SurfExchangeCUDA_AFM*>(pModuleCUDA)->GetTorque(avRect);
#endif

	return CalculateTorque(pMesh->M, avRect);
}

#endif