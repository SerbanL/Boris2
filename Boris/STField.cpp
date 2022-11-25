#include "stdafx.h"
#include "STField.h"

#ifdef MODULE_COMPILATION_STFIELD

#include "SuperMesh.h"
#include "Mesh_Ferromagnetic.h"
#include "MeshParamsControl.h"

#include "Atom_Mesh.h"
#include "Atom_MeshParamsControl.h"

#if COMPILECUDA == 1
#include "STFieldCUDA.h"
#endif

STField::STField(Mesh *pMesh_) :
	Modules(),
	ProgramStateNames(this, {}, {})
{
	pMesh = pMesh_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

STField::~STField()
{
}

BError STField::Initialize(void)
{
	BError error(CLASS_STR(STField));

	SuperMesh* pSMesh = pMesh->pSMesh;

	pMesh_Bot.clear();
	pMesh_Top.clear();
	paMesh_Bot.clear();
	paMesh_Top.clear();

	if (pMesh->STp.get0() == DBL3()) {

		//not fixed polarization : this is instead provided directly by other magnetic meshes top and bottom along z direction

		//Need to identify all magnetic meshes which provide polarization vectors:
		//1. Must be ferromagnetic (micromagnetic or atomistic)
		//2. Must overlap in the x-y plane only with the mesh holding this module (either top or bottom) but not along z (i.e. mesh rectangles must not intersect)
		//3. No other magnetic meshes can be sandwiched in between - there could be other types of non-magnetic meshes in between of course (e.g. insulator, conductive layers etc).

		//---

		//lambda used to check condition 3
		auto check_candidate = [&](Rect xy_intersection, double z1, double z2) -> bool {

			//check all meshes to find a magnetic mesh, except for 2-sublattice, which intersects in the xy plane with xy_intersection, and has z coordinates between z1 and z2.
			//if found then current candidate not good
			for (int idx = 0; idx < pSMesh->size(); idx++) {

				//consider all meshes in turn - condition 1 first
				if ((*pSMesh)[idx]->MComputation_Enabled() && (*pSMesh)[idx]->GetMeshType() != MESH_ANTIFERROMAGNETIC) {

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
			if ((*pSMesh)[idx]->MComputation_Enabled() && (*pSMesh)[idx]->GetMeshType() != MESH_ANTIFERROMAGNETIC) {

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
	}

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect, 
		(MOD_)pMesh->Get_Module_Heff_Display() == MOD_STFIELD, 
		(MOD_)pMesh->Get_Module_Energy_Display() == MOD_STFIELD, 
		pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	return error;
}

BError STField::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(STField));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError STField::MakeCUDAModule(void)
{
	BError error(CLASS_STR(STField));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new STFieldCUDA(pMesh->pMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double STField::UpdateField(void)
{
	if (!pMesh->E.linear_size()) return 0.0;

	if (pMesh->STp.get0() != DBL3()) {

		////////////////////////////////////////
		//      FIXED POLARIZATION VERSION    //
		////////////////////////////////////////

#pragma omp parallel for
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			double grel = pMesh->grel;
			double Ms = pMesh->Ms;
			double flSOT = pMesh->flSOT;
			DBL2 STq = pMesh->STq;
			DBL2 STa = pMesh->STa;
			DBL3 STp = pMesh->STp;
			pMesh->update_parameters_mcoarse(idx, pMesh->grel, grel, pMesh->Ms, Ms, pMesh->flSOT, flSOT, pMesh->STq, STq, pMesh->STa, STa, pMesh->STp, STp);

			if (IsNZ(grel)) {

				double dotprod = (pMesh->M[idx] * STp) / Ms;
				double neta = STq.i / (STa.i + STa.j * dotprod) + STq.j / (STa.i - STa.j * dotprod);

				int idx_E = pMesh->E.position_to_cellidx(pMesh->M.cellidx_to_position(idx));
				//z component of J
				double Jc = pMesh->E[idx_E].z * pMesh->elC[idx_E];

				double a_const = -(neta * MUB_E * Jc / (GAMMA * grel)) / (Ms * Ms * pMesh->GetMeshDimensions().z);

				DBL3 ST_Field = a_const * ((pMesh->M[idx] ^ STp) + flSOT * Ms * STp);
				pMesh->Heff[idx] += ST_Field;

				if (Module_Heff.linear_size()) Module_Heff[idx] = ST_Field;
			}
		}
	}
	else {

		/////////////////////////////////////////////
		// POLARIZATION FROM TOP AND BOTTOM MESHES //
		/////////////////////////////////////////////

		SZ3 n = pMesh->n;

		//zero module display VECs if needed, since contributions must be added into them to account for possiblility of 2 contributions (top and bottom)
		ZeroModuleVECs();

		if (pMesh_Top.size() || paMesh_Top.size()) {

			//polarization from top mesh
#pragma omp parallel for
			for (int j = 0; j < n.y; j++) {
				for (int i = 0; i < n.x; i++) {

					int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

					//empty cell here ... next
					if (pMesh->M.is_empty(cell_idx)) continue;

					double grel = pMesh->grel;
					double Ms = pMesh->Ms;
					double flSOT2 = pMesh->flSOT2;
					DBL2 STq2 = pMesh->STq2;
					DBL2 STa2 = pMesh->STa2;
					pMesh->update_parameters_mcoarse(cell_idx, 
						pMesh->grel, grel, pMesh->Ms, Ms, 
						pMesh->flSOT2, flSOT2, pMesh->STq2, STq2, pMesh->STa2, STa2);

					//effective field for this cell
					DBL3 ST_Field = DBL3();
					bool cell_coupled = false;

					//check all meshes for coupling
					//1 : coupling into this micromagnetic mesh from other micromagnetic meshes (FM)
					for (int mesh_idx = 0; mesh_idx < (int)pMesh_Top.size(); mesh_idx++) {

						if (!check_cell_coupling(pMesh->M, pMesh_Top[mesh_idx]->M,
							(i + 0.5) * pMesh->h.x, (j + 0.5) * pMesh->h.y, pMesh_Top[mesh_idx]->h.z / 2)) continue;

						DBL3 cell_rel_pos = DBL3(
							(i + 0.5) * pMesh->h.x + pMesh->M.rect.s.x - pMesh_Top[mesh_idx]->M.rect.s.x,
							(j + 0.5) * pMesh->h.y + pMesh->M.rect.s.y - pMesh_Top[mesh_idx]->M.rect.s.y,
							pMesh_Top[mesh_idx]->h.z / 2);

						if (pMesh_Top[mesh_idx]->GetMeshType() == MESH_FERROMAGNETIC) {

							//get magnetization value in top mesh cell (the polarization)
							DBL3 p = normalize(pMesh_Top[mesh_idx]->M[cell_rel_pos]);
							DBL3 m = normalize(pMesh->M[cell_idx]);

							double dotprod = m * p;
							double neta = STq2.i / (STa2.i + STa2.j * dotprod) + STq2.j / (STa2.i - STa2.j * dotprod);

							int idx_E = pMesh->E.position_to_cellidx(pMesh->M.cellidx_to_position(cell_idx));
							//z component of J
							double Jc = pMesh->E[idx_E].z * pMesh->elC[idx_E];

							double a_const = -(neta * MUB_E * Jc / (GAMMA * grel)) / (Ms * pMesh->h.z);

							ST_Field = a_const * ((m ^ p) + flSOT2 * p);
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

							//get magnetization value in top mesh cell (the polarization)
							DBL3 p = normalize(M1.average(rect_c - M1.rect.get_s()));
							DBL3 m = normalize(pMesh->M[cell_idx]);

							double dotprod = m * p;
							double neta = STq2.i / (STa2.i + STa2.j * dotprod) + STq2.j / (STa2.i - STa2.j * dotprod);

							int idx_E = pMesh->E.position_to_cellidx(pMesh->M.cellidx_to_position(cell_idx));
							//z component of J
							double Jc = pMesh->E[idx_E].z * pMesh->elC[idx_E];

							double a_const = -(neta * MUB_E * Jc / (GAMMA * grel)) / (Ms * pMesh->h.z);

							ST_Field = a_const * ((m ^ p) + flSOT2 * p);

							//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
							break;
						}
					}

					pMesh->Heff[cell_idx] += ST_Field;
					if (Module_Heff.linear_size()) Module_Heff[cell_idx] += ST_Field;
				}
			}
		}

		if (pMesh_Bot.size() || paMesh_Bot.size()) {

			//polarization from bottom mesh
#pragma omp parallel for
			for (int j = 0; j < n.y; j++) {
				for (int i = 0; i < n.x; i++) {

					int cell_idx = i + j * n.x;

					//empty cell here ... next
					if (pMesh->M.is_empty(cell_idx)) continue;

					double grel = pMesh->grel;
					double Ms = pMesh->Ms;
					double flSOT = pMesh->flSOT;
					DBL2 STq = pMesh->STq;
					DBL2 STa = pMesh->STa;
					pMesh->update_parameters_mcoarse(cell_idx,
						pMesh->grel, grel, pMesh->Ms, Ms,
						pMesh->flSOT, flSOT, pMesh->STq, STq, pMesh->STa, STa);

					//effective field for this cell
					DBL3 ST_Field = DBL3();
					bool cell_coupled = false;

					//check all meshes for coupling
					//1 : coupling into this micromagnetic mesh from other micromagnetic meshes (FM)
					for (int mesh_idx = 0; mesh_idx < (int)pMesh_Bot.size(); mesh_idx++) {

						if (!check_cell_coupling(pMesh->M, pMesh_Bot[mesh_idx]->M,
							(i + 0.5) * pMesh->h.x, (j + 0.5) * pMesh->h.y, pMesh_Bot[mesh_idx]->meshRect.height() - pMesh_Bot[mesh_idx]->h.z / 2)) continue;

						DBL3 cell_rel_pos = DBL3(
							(i + 0.5) * pMesh->h.x + pMesh->M.rect.s.x - pMesh_Bot[mesh_idx]->M.rect.s.x,
							(j + 0.5) * pMesh->h.y + pMesh->M.rect.s.y - pMesh_Bot[mesh_idx]->M.rect.s.y,
							pMesh_Bot[mesh_idx]->meshRect.height() - pMesh_Bot[mesh_idx]->h.z / 2);

						if (pMesh_Bot[mesh_idx]->GetMeshType() == MESH_FERROMAGNETIC) {

							//get magnetization value in bottom mesh cell (the polarization)
							DBL3 p = normalize(pMesh_Bot[mesh_idx]->M[cell_rel_pos]);
							DBL3 m = normalize(pMesh->M[cell_idx]);

							double dotprod = m * p;
							double neta = STq.i / (STa.i + STa.j * dotprod) + STq.j / (STa.i - STa.j * dotprod);

							int idx_E = pMesh->E.position_to_cellidx(pMesh->M.cellidx_to_position(cell_idx));
							//z component of J
							double Jc = pMesh->E[idx_E].z * pMesh->elC[idx_E];

							double a_const = -(neta * MUB_E * Jc / (GAMMA * grel)) / (Ms * pMesh->h.z);

							ST_Field = a_const * ((m ^ p) + flSOT * p);
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

							//get magnetization value in top mesh cell (the polarization)
							DBL3 p = normalize(M1.average(rect_c - M1.rect.get_s()));
							DBL3 m = normalize(pMesh->M[cell_idx]);

							double dotprod = m * p;
							double neta = STq.i / (STa.i + STa.j * dotprod) + STq.j / (STa.i - STa.j * dotprod);

							int idx_E = pMesh->E.position_to_cellidx(pMesh->M.cellidx_to_position(cell_idx));
							//z component of J
							double Jc = pMesh->E[idx_E].z * pMesh->elC[idx_E];

							double a_const = -(neta * MUB_E * Jc / (GAMMA * grel)) / (Ms * pMesh->h.z);

							ST_Field = a_const * ((m ^ p) + flSOT * p);

							//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
							break;
						}
					}

					pMesh->Heff[cell_idx] += ST_Field;
					if (Module_Heff.linear_size()) Module_Heff[cell_idx] += ST_Field;
				}
			}
		}
	}

	//don't count this as a contribution to the total energy density
	return 0.0;
}

#endif
