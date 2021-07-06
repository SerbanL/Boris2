#include "stdafx.h"
#include "Atom_SurfExchange.h"

#if defined(MODULE_COMPILATION_SURFEXCHANGE) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "SuperMesh.h"
#include "Atom_MeshParamsControl.h"

/////////////////////////////////////////////////////////////////
//SurfExchange
//

Atom_SurfExchange::Atom_SurfExchange(Atom_Mesh *paMesh_) :
	Modules(),
	ProgramStateNames(this, {}, {})
{
	paMesh = paMesh_;

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (paMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

Atom_SurfExchange::~Atom_SurfExchange()
{
}

BError Atom_SurfExchange::Initialize(void)
{
	BError error(CLASS_STR(Atom_SurfExchange));

	//Need to identify all magnetic meshes participating in surface exchange coupling with this module:
	//1. Must be ferromagnetic or anti-ferromagnetic and have the SurfExchange module set -> this results in surface exchange coupling (or equivalently exchange bias - same formula, except in this case you should have J2 set to zero)
	//2. Must overlap in the x-y plane only with the mesh holding this module (either top or bottom) but not along z (i.e. mesh rectangles must not intersect)
	//3. No other magnetic meshes can be sandwiched in between - there could be other types of non-magnetic meshes in between of course (e.g. insulator, conductive layers etc).

	SuperMesh* pSMesh = paMesh->pSMesh;

	paMesh_Bot.clear();
	paMesh_Top.clear();

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

	Rect meshRect = paMesh->GetMeshRect();

	Rect xy_meshRect = meshRect;
	xy_meshRect.s.z = 0; xy_meshRect.e.z = 0;

	for (int idx = 0; idx < pSMesh->size(); idx++) {

		//consider all meshes in turn - condition 1 first
		if ((*pSMesh)[idx]->MComputation_Enabled() && (*pSMesh)[idx]->IsModuleSet(MOD_SURFEXCHANGE) && (*pSMesh)[idx]->is_atomistic()) {

			Rect candidate_meshRect = (*pSMesh)[idx]->GetMeshRect();

			//candidate mesh at the top
			if (IsGE(candidate_meshRect.s.z, meshRect.e.z)) {

				double z1 = meshRect.e.z;				//start z
				double z2 = candidate_meshRect.s.z;		//end z
				candidate_meshRect.s.z = 0; candidate_meshRect.e.z = 0;		//leave just the xy plane rect

				if (candidate_meshRect.intersects(xy_meshRect)) {

					//passes condition 2 - identified a candidate mesh at idx index. Does it pass condition 3?
					if (check_candidate(candidate_meshRect.get_intersection(xy_meshRect), z1, z2)) {

						paMesh_Top.push_back(dynamic_cast<Atom_Mesh*>((*pSMesh)[idx]));
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

						paMesh_Bot.push_back(dynamic_cast<Atom_Mesh*>((*pSMesh)[idx]));
					}
				}
			}
		}
	}

	//count number of coupled cells (either top or bottom) in this mesh

	SZ3 n = paMesh->n;

	if (paMesh_Top.size()) {

		//surface exchange coupling at the top
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

				//empty cell here ... next
				if (paMesh->M1.is_empty(cell_idx)) continue;

				//check all meshes for coupling
				for (int mesh_idx = 0; mesh_idx < (int)paMesh_Top.size(); mesh_idx++) {

					Rect tmeshRect = paMesh_Top[mesh_idx]->GetMeshRect();

					//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
					DBL3 cell_rel_pos = DBL3(
						(i + 0.5) * paMesh->h.x + meshRect.s.x - tmeshRect.s.x,
						(j + 0.5) * paMesh->h.y + meshRect.s.y - tmeshRect.s.y,
						paMesh_Top[mesh_idx]->h.z / 2);

					//can't couple to an empty cell
					if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s) || paMesh_Top[mesh_idx]->M1.is_empty(cell_rel_pos)) continue;

					//if we are here then the cell in this mesh at cell_idx has something to couple to so count it : it will contribute to the surface exchange energy density
					coupled_cells++;
				}
			}
		}
	}

	if (paMesh_Bot.size()) {

		//surface exchange coupling at the bottom
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				int cell_idx = i + j * n.x;

				//empty cell here ... next
				if (paMesh->M1.is_empty(cell_idx)) continue;

				//check all meshes for coupling
				for (int mesh_idx = 0; mesh_idx < (int)paMesh_Bot.size(); mesh_idx++) {

					Rect bmeshRect = paMesh_Bot[mesh_idx]->GetMeshRect();

					//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
					DBL3 cell_rel_pos = DBL3(
						(i + 0.5) * paMesh->h.x + meshRect.s.x - bmeshRect.s.x,
						(j + 0.5) * paMesh->h.y + meshRect.s.y - bmeshRect.s.y,
						paMesh_Bot[mesh_idx]->meshRect.e.z - paMesh_Bot[mesh_idx]->meshRect.s.z - (paMesh_Bot[mesh_idx]->h.z / 2));

					//can't couple to an empty cell
					if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s) || paMesh_Bot[mesh_idx]->M1.is_empty(cell_rel_pos)) continue;

					//if we are here then the cell in this mesh at cell_idx has something to couple to so count it : it will contribute to the surface exchange energy density
					coupled_cells++;
				}
			}
		}
	}

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		paMesh->h, paMesh->meshRect,
		(MOD_)paMesh->Get_Module_Heff_Display() == MOD_SURFEXCHANGE || paMesh->IsOutputDataSet_withRect(DATA_E_SURFEXCH) || paMesh->IsOutputDataSet(DATA_T_SURFEXCH),
		(MOD_)paMesh->Get_Module_Energy_Display() == MOD_SURFEXCHANGE || paMesh->IsOutputDataSet_withRect(DATA_E_SURFEXCH));
	if (!error) initialized = true;

	return error;
}

BError Atom_SurfExchange::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_SurfExchange));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError Atom_SurfExchange::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Atom_SurfExchange));

#if COMPILECUDA == 1

	if (paMesh->paMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Atom_SurfExchangeCUDA(paMesh->paMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_SurfExchange::UpdateField(void)
{
	double energy = 0;
	
	SZ3 n = paMesh->n;

	//zero module display VECs if needed, since contributions must be added into them to account for possiblility of 2 contributions (top and bottom)
	ZeroModuleVECs();

	if (paMesh_Top.size()) {

		//surface exchange coupling at the top
#pragma omp parallel for reduction(+:energy)
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

				//empty cell here ... next
				if (paMesh->M1.is_empty(cell_idx)) continue;

				double mu_s = paMesh->mu_s;
				paMesh->update_parameters_mcoarse(cell_idx, paMesh->mu_s, mu_s);

				//check all meshes for coupling
				for (int mesh_idx = 0; mesh_idx < (int)paMesh_Top.size(); mesh_idx++) {

					Rect tmeshRect = paMesh_Top[mesh_idx]->GetMeshRect();

					//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
					DBL3 cell_rel_pos = DBL3(
						(i + 0.5) * paMesh->h.x + paMesh->meshRect.s.x - tmeshRect.s.x,
						(j + 0.5) * paMesh->h.y + paMesh->meshRect.s.y - tmeshRect.s.y,
						paMesh_Top[mesh_idx]->h.z / 2);

					//can't couple to an empty cell
					if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s) || paMesh_Top[mesh_idx]->M1.is_empty(cell_rel_pos)) continue;

					//Top mesh sets Js
					double Js = paMesh_Top[mesh_idx]->Js;
					paMesh_Top[mesh_idx]->update_parameters_atposition(cell_rel_pos, paMesh_Top[mesh_idx]->Js, Js);

					//get magnetization value in top mesh cell to couple with
					DBL3 m_j = paMesh_Top[mesh_idx]->M1[cell_rel_pos].normalized();
					DBL3 m_i = paMesh->M1[cell_idx] / mu_s;

					double dot_prod = m_i * m_j;

					DBL3 Hsurfexh = m_j * Js / (MUB_MU0 * mu_s);
					double cell_energy = -Js * dot_prod / paMesh->M1.h.dim();

					paMesh->Heff1[cell_idx] += Hsurfexh;

					if (Module_Heff.linear_size()) Module_Heff[cell_idx] += Hsurfexh;
					if (Module_energy.linear_size()) Module_energy[cell_idx] += cell_energy;

					energy += cell_energy;

					//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
					break;
				}
			}
		}
	}

	if (paMesh_Bot.size()) {

		//surface exchange coupling at the bottom
#pragma omp parallel for reduction(+:energy)
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				int cell_idx = i + j * n.x;

				//empty cell here ... next
				if (paMesh->M1.is_empty(cell_idx)) continue;

				double mu_s = paMesh->mu_s;
				double Js = paMesh->Js;
				paMesh->update_parameters_mcoarse(cell_idx, paMesh->mu_s, mu_s, paMesh->Js, Js);

				//check all meshes for coupling
				for (int mesh_idx = 0; mesh_idx < (int)paMesh_Bot.size(); mesh_idx++) {

					Rect bmeshRect = paMesh_Bot[mesh_idx]->GetMeshRect();

					//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
					DBL3 cell_rel_pos = DBL3(
						(i + 0.5) * paMesh->h.x + paMesh->meshRect.s.x - bmeshRect.s.x,
						(j + 0.5) * paMesh->h.y + paMesh->meshRect.s.y - bmeshRect.s.y,
						paMesh_Bot[mesh_idx]->meshRect.e.z - paMesh_Bot[mesh_idx]->meshRect.s.z - (paMesh_Bot[mesh_idx]->h.z / 2));

					//can't couple to an empty cell
					if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s) || paMesh_Bot[mesh_idx]->M1.is_empty(cell_rel_pos)) continue;

					//get magnetization value in top mesh cell to couple with
					DBL3 m_j = paMesh_Bot[mesh_idx]->M1[cell_rel_pos].normalized();
					DBL3 m_i = paMesh->M1[cell_idx] / mu_s;

					double dot_prod = m_i * m_j;

					DBL3 Hsurfexh = m_j * Js / (MUB_MU0 * mu_s);
					double cell_energy = -Js * dot_prod / paMesh->M1.h.dim();

					paMesh->Heff1[cell_idx] += Hsurfexh;

					energy += cell_energy;

					if (Module_Heff.linear_size()) Module_Heff[cell_idx] = Hsurfexh;
					if (Module_energy.linear_size()) Module_energy[cell_idx] = cell_energy;

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

DBL3 Atom_SurfExchange::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return reinterpret_cast<Atom_SurfExchangeCUDA*>(pModuleCUDA)->GetTorque(avRect);
#endif

	return CalculateTorque(paMesh->M1, avRect);
}

//-------------------Energy methods

//For simple cubic mesh spin_index coincides with index in M1
double Atom_SurfExchange::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	double energy_new = 0, energy_old = 0;

	SZ3 n = paMesh->n;

	//if spin is on top surface then look at paMesh_Top
	if (spin_index / (n.x * n.y) == n.z - 1 && paMesh_Top.size()) {

		if (!paMesh->M1.is_empty(spin_index)) {

			int i = spin_index % n.x;
			int j = (spin_index / n.x) % n.y;

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < (int)paMesh_Top.size(); mesh_idx++) {

				Rect tmeshRect = paMesh_Top[mesh_idx]->GetMeshRect();

				//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
				DBL3 cell_rel_pos = DBL3(
					(i + 0.5) * paMesh->h.x + paMesh->meshRect.s.x - tmeshRect.s.x,
					(j + 0.5) * paMesh->h.y + paMesh->meshRect.s.y - tmeshRect.s.y,
					paMesh_Top[mesh_idx]->h.z / 2);

				//can't couple to an empty cell
				if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s) || paMesh_Top[mesh_idx]->M1.is_empty(cell_rel_pos)) continue;

				//Top mesh sets Js
				double Js = paMesh_Top[mesh_idx]->Js;
				paMesh_Top[mesh_idx]->update_parameters_atposition(cell_rel_pos, paMesh_Top[mesh_idx]->Js, Js);

				//get magnetization value in top mesh cell to couple with
				DBL3 m_j = paMesh_Top[mesh_idx]->M1[cell_rel_pos].normalized();
				DBL3 m_i = paMesh->M1[spin_index].normalized();
				double dot_prod = m_i * m_j;
				energy_old = -Js * dot_prod;

				if (Mnew != DBL3()) {

					DBL3 mnew_i = Mnew.normalized();
					double dot_prod_new = mnew_i * m_j;
					energy_new = -Js * dot_prod_new;
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	//if spin is on bottom surface then look at paMesh_Top
	if (spin_index / (n.x * n.y) == 0 && paMesh_Bot.size()) {

		if (!paMesh->M1.is_empty(spin_index)) {

			int i = spin_index % n.x;
			int j = (spin_index / n.x) % n.y;

			double Js = paMesh->Js;
			paMesh->update_parameters_mcoarse(spin_index, paMesh->Js, Js);

			//check all meshes for coupling
			for (int mesh_idx = 0; mesh_idx < (int)paMesh_Bot.size(); mesh_idx++) {

				Rect bmeshRect = paMesh_Bot[mesh_idx]->GetMeshRect();

				//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
				DBL3 cell_rel_pos = DBL3(
					(i + 0.5) * paMesh->h.x + paMesh->meshRect.s.x - bmeshRect.s.x,
					(j + 0.5) * paMesh->h.y + paMesh->meshRect.s.y - bmeshRect.s.y,
					paMesh_Bot[mesh_idx]->meshRect.e.z - paMesh_Bot[mesh_idx]->meshRect.s.z - (paMesh_Bot[mesh_idx]->h.z / 2));

				//can't couple to an empty cell
				if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s) || paMesh_Bot[mesh_idx]->M1.is_empty(cell_rel_pos)) continue;

				//get magnetization value in top mesh cell to couple with
				DBL3 m_j = paMesh_Bot[mesh_idx]->M1[cell_rel_pos].normalized();
				DBL3 m_i = paMesh->M1[spin_index].normalized();
				double dot_prod = m_i * m_j;
				energy_old = -Js * dot_prod;

				if (Mnew != DBL3()) {

					DBL3 mnew_i = Mnew.normalized();
					double dot_prod_new = mnew_i * m_j;
					energy_new = -Js * dot_prod_new;
				}

				//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
				break;
			}
		}
	}

	if (Mnew != DBL3()) return energy_new - energy_old;
	else return energy_old;
}

#endif