#include "stdafx.h"
#include "SurfExchange.h"

#ifdef MODULE_COMPILATION_SURFEXCHANGE

#include "Mesh_Ferromagnetic.h"
#include "Mesh_AntiFerromagnetic.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

/////////////////////////////////////////////////////////////////
//SurfExchange
//

SurfExchange::SurfExchange(Mesh *pMesh_) :
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

SurfExchange::~SurfExchange()
{
}

BError SurfExchange::Initialize(void)
{
	BError error(CLASS_STR(SurfExchange));

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
	auto check_candidate = [&](Rect xy_intersection, double z1, double z2) -> int {

		//xy_intersection must have z coordinates set to zero, but set them here to be sure
		xy_intersection.s.z = 0; xy_intersection.e.z = 0;

		//check all meshes to find a magnetic mesh with SurfExchange modules set, which intersects in the xy plane with xy_intersection, and has z coordinates between z1 and z2.
		//return index of found mesh, else return -1 (-1 return means candidate being checked - which generated this call - is valid, else a new, better candidate has been found.
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
						return idx;
					}
				}
			}
		}

		//no new candidate found - current candidate has been validated as the best one of its type (with given intersection)
		return -1;
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
					int coupledmesh_idx = idx;
					while (true) {

						int check_idx = check_candidate(candidate_meshRect.get_intersection(xy_meshRect), z1, z2);

						if (check_idx < 0) {

							pMesh_Top.push_back(dynamic_cast<Mesh*>((*pSMesh)[coupledmesh_idx]));
							
							break;
						}
						else {

							coupledmesh_idx = check_idx;
							candidate_meshRect = (*pSMesh)[check_idx]->GetMeshRect();
							z2 = candidate_meshRect.s.z;
							candidate_meshRect.s.z = 0; candidate_meshRect.e.z = 0;
						}
					};
				}
			}

			//candidate mesh at the botttom
			else if (IsSE(candidate_meshRect.e.z, meshRect.s.z)) {

				double z1 = candidate_meshRect.e.z;		//start z
				double z2 = meshRect.s.z;				//end z
				candidate_meshRect.s.z = 0; candidate_meshRect.e.z = 0;		//leave just the xy plane rect

				if (candidate_meshRect.intersects(xy_meshRect)) {

					//passes condition 2 - identified a candidate mesh at idx index. Does it pass condition 3?
					int coupledmesh_idx = idx;
					while (true) {

						int check_idx = check_candidate(candidate_meshRect.get_intersection(xy_meshRect), z1, z2);

						if (check_idx < 0) {

							pMesh_Bot.push_back(dynamic_cast<Mesh*>((*pSMesh)[coupledmesh_idx]));
							
							break;
						}
						else {

							coupledmesh_idx = check_idx;
							candidate_meshRect = (*pSMesh)[check_idx]->GetMeshRect();
							z1 = candidate_meshRect.e.z;
							candidate_meshRect.s.z = 0; candidate_meshRect.e.z = 0;
						}
					};
				}
			}
		}
	}

	//count number of coupled cells (either top or bottom) in this mesh

	INT3 n = pMesh->n;

	if (pMesh_Top.size()) {

		//surface exchange coupling at the top
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

				//empty cell here ... next
				if (pMesh->M.is_empty(cell_idx)) continue;

				//check all meshes for coupling
				for (int mesh_idx = 0; mesh_idx < (int)pMesh_Top.size(); mesh_idx++) {

					//relative coordinates to read value from top mesh (the one we're coupling to here)
					DBL3 cell_rel_pos = DBL3((i + 0.5) * pMesh->h.x, (j + 0.5) * pMesh->h.y, pMesh_Top[mesh_idx]->h.z / 2);

					//can't couple to an empty cell
					if (IsZ(pMesh_Top[mesh_idx]->M[cell_rel_pos].norm())) continue;

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

					//relative coordinates to read value from bottom mesh (the one we're coupling to here)
					DBL3 cell_rel_pos = DBL3((i + 0.5) * pMesh->h.x, (j + 0.5) * pMesh->h.y, pMesh_Bot[mesh_idx]->meshRect.e.z - (pMesh_Bot[mesh_idx]->h.z / 2));

					//can't couple to an empty cell
					if (IsZ(pMesh_Bot[mesh_idx]->M[cell_rel_pos].norm())) continue;

					//if we are here then the cell in this mesh at cell_idx has something to couple to so count it : it will contribute to the surface exchange energy density
					coupled_cells++;
				}
			}
		}
	}

	initialized = true;

	return error;
}

BError SurfExchange::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SurfExchange));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError SurfExchange::MakeCUDAModule(void)
{
	BError error(CLASS_STR(SurfExchange));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new SurfExchangeCUDA(pMesh->pMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double SurfExchange::UpdateField(void)
{
	if (pMesh->GetMeshType() == MESH_DIAMAGNETIC) return 0.0;

	double energy = 0;

	INT3 n = pMesh->n;

	//thickness of layer - SurfExchange applies for layers in the xy plane
	double thickness = pMesh->meshRect.e.z - pMesh->meshRect.s.z;

	//FM to FM coupling at the top
	if (pMesh_Top.size()) {

		//surface exchange coupling at the top
		#pragma omp parallel for reduction(+:energy)
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

				//empty cell here ... next
				if (pMesh->M.is_empty(cell_idx)) continue;

				double Ms = pMesh->Ms;
				pMesh->update_parameters_mcoarse(cell_idx, pMesh->Ms, Ms);

				//check all meshes for coupling
				for (int mesh_idx = 0; mesh_idx < (int)pMesh_Top.size(); mesh_idx++) {

					//relative coordinates to read value from top mesh (the one we're coupling to here)
					DBL3 cell_rel_pos = DBL3((i + 0.5) * pMesh->h.x, (j + 0.5) * pMesh->h.y, pMesh_Top[mesh_idx]->h.z / 2);

					//can't couple to an empty cell
					if (IsZ(pMesh_Top[mesh_idx]->M[cell_rel_pos].norm())) continue;

					DBL3 Hsurfexh;

					if (pMesh_Top[mesh_idx]->GetMeshType() == MESH_DIAMAGNETIC) {

						//Surface exchange field from a diamagnetic mesh - EXPERIMENTAL

						//diamagnetic mesh sets coupling constant
						double neta_dia = pMesh_Top[mesh_idx]->neta_dia;
						pMesh_Top[mesh_idx]->update_parameters_atposition(cell_rel_pos, pMesh_Top[mesh_idx]->neta_dia, neta_dia);

						//Hse = neta_dia * sus * Hext / (mu0 * Ms * tF)
						//energy_density  = -neta_dia * sus * Hext.mF / tF

						DBL3 Mdia = pMesh_Top[mesh_idx]->M[cell_rel_pos];
						DBL3 m_i = pMesh->M[cell_idx] / Ms;

						Hsurfexh = neta_dia * Mdia / (MU0 * Ms * thickness);
						energy += -neta_dia * (Mdia*m_i) / thickness;
					}
					else {

						//Surface exchange field from a ferromagnetic (RKKY) or an an antiferromagnetic mesh (exchange bias)

						//Top mesh sets J1 and J2 values
						double J1 = pMesh_Top[mesh_idx]->J1;
						double J2 = pMesh_Top[mesh_idx]->J2;
						pMesh_Top[mesh_idx]->update_parameters_atposition(cell_rel_pos, pMesh_Top[mesh_idx]->J1, J1, pMesh_Top[mesh_idx]->J2, J2);

						//get magnetization value in top mesh cell to couple with
						DBL3 m_j = pMesh_Top[mesh_idx]->M[cell_rel_pos].normalized();
						DBL3 m_i = pMesh->M[cell_idx] / Ms;

						double dot_prod = m_i * m_j;

						//total surface exchange field in coupling cells, including bilinear and biquadratic terms
						Hsurfexh = (m_j / (MU0 * Ms * thickness)) * (J1 + 2 * J2 * dot_prod);
						energy += (-1 * J1 - 2 * J2 * dot_prod) * dot_prod / thickness;
					}

					//couple all cells through the layer thickness : the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
					for (int k = 0; k < n.z; k++) pMesh->Heff[i + j * n.x + k * n.x*n.y] += Hsurfexh;
				}
			}
		}
	}

	//FM to FM coupling at the bottom
	if (pMesh_Bot.size()) {

		//surface exchange coupling at the bottom
		#pragma omp parallel for reduction(+:energy)
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				int cell_idx = i + j * n.x;

				//empty cell here ... next
				if (pMesh->M.is_empty(cell_idx)) continue;

				double Ms = pMesh->Ms;
				double J1 = pMesh->J1;
				double J2 = pMesh->J2;
				pMesh->update_parameters_mcoarse(cell_idx, pMesh->Ms, Ms, pMesh->J1, J1, pMesh->J2, J2);

				//check all meshes for coupling
				for (int mesh_idx = 0; mesh_idx < (int)pMesh_Bot.size(); mesh_idx++) {

					//relative coordinates to read value from bottom mesh (the one we're coupling to here)
					DBL3 cell_rel_pos = DBL3((i + 0.5) * pMesh->h.x, (j + 0.5) * pMesh->h.y, pMesh_Bot[mesh_idx]->meshRect.e.z - (pMesh_Bot[mesh_idx]->h.z / 2));

					//can't couple to an empty cell
					if (IsZ(pMesh_Bot[mesh_idx]->M[cell_rel_pos].norm())) continue;

					DBL3 Hsurfexh;

					if (pMesh_Bot[mesh_idx]->GetMeshType() == MESH_DIAMAGNETIC) {

						//Surface exchange field from a diamagnetic mesh - EXPERIMENTAL

						//diamagnetic mesh sets coupling constant
						double neta_dia = pMesh_Bot[mesh_idx]->neta_dia;
						pMesh_Bot[mesh_idx]->update_parameters_atposition(cell_rel_pos, pMesh_Bot[mesh_idx]->neta_dia, neta_dia);

						//Hse = neta_dia * sus * Hext / (mu0 * Ms * tF)
						//energy_density  = -neta_dia * sus * Hext.mF / tF

						DBL3 Mdia = pMesh_Bot[mesh_idx]->M[cell_rel_pos];
						DBL3 m_i = pMesh->M[cell_idx] / Ms;

						Hsurfexh = neta_dia * Mdia / (MU0 * Ms * thickness);
						energy += -neta_dia * (Mdia*m_i) / thickness;
					}
					else {

						//Surface exchange field from a ferromagnetic (RKKY) or an an antiferromagnetic mesh (exchange bias)

						//get magnetization value in top mesh cell to couple with
						DBL3 m_j = pMesh_Bot[mesh_idx]->M[cell_rel_pos].normalized();
						DBL3 m_i = pMesh->M[cell_idx] / Ms;

						double dot_prod = m_i * m_j;

						//total surface exchange field in coupling cells, including bilinear and biquadratic terms
						Hsurfexh = (m_j / (MU0 * Ms * thickness)) * (J1 + 2 * J2 * dot_prod);
						energy += (-1 * J1 - 2 * J2 * dot_prod) * dot_prod / thickness;
					}

					//couple all cells through the layer thickness : the surface exchange field is applicable for effectively 2D layers, but the simulation allows 3D meshes.
					for (int k = 0; k < n.z; k++) pMesh->Heff[i + j * n.x + k * n.x*n.y] += Hsurfexh;
				}
			}
		}
	}

	if (coupled_cells) energy /= coupled_cells;
	else energy = 0.0;

	this->energy = energy;

	return this->energy;
}

#endif