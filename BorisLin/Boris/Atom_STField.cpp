#include "stdafx.h"
#include "Atom_STField.h"

#if defined(MODULE_COMPILATION_STFIELD) && ATOMISTIC == 1

#include "SuperMesh.h"
#include "Atom_Mesh.h"
#include "Atom_MeshParamsControl.h"
#include "Mesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "Atom_STFieldCUDA.h"
#endif

Atom_STField::Atom_STField(Atom_Mesh *paMesh_) :
	Modules(),
	ProgramStateNames(this, {}, {})
{
	paMesh = paMesh_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (paMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

Atom_STField::~Atom_STField()
{
}

BError Atom_STField::Initialize(void)
{
	BError error(CLASS_STR(Atom_STField));

	SuperMesh* pSMesh = paMesh->pSMesh;

	paMesh_Bot.clear();
	paMesh_Top.clear();
	pMesh_Bot.clear();
	pMesh_Top.clear();

	if (paMesh->STp.get0() == DBL3()) {

		//not fixed polarization : this is instead provided directly by other magnetic meshes top and bottom along z direction

		//Need to identify all magnetic meshes which provide polarization vectors:
		//1. Must be ferromagnetic (micromagnetic or atomistic)
		//2. Must overlap in the x-y plane only with the mesh holding this module (either top or bottom) but not along z (i.e. mesh rectangles must not intersect)
		//3. No other magnetic meshes can be sandwiched in between - there could be other types of non-magnetic meshes in between of course (e.g. insulator, conductive layers etc).

		//---

		//lambda used to check condition 3
		auto check_candidate = [&](Rect xy_intersection, double z1, double z2) -> bool {

			//check all meshes to find a magnetic mesh with SurfExchange modules set, which intersects in the xy plane with xy_intersection, and has z coordinates between z1 and z2.
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

		Rect meshRect = paMesh->GetMeshRect();

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

							if ((*pSMesh)[idx]->is_atomistic()) paMesh_Top.push_back(dynamic_cast<Atom_Mesh*>((*pSMesh)[idx]));
							else pMesh_Top.push_back(dynamic_cast<Mesh*>((*pSMesh)[idx]));
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

							if ((*pSMesh)[idx]->is_atomistic()) paMesh_Bot.push_back(dynamic_cast<Atom_Mesh*>((*pSMesh)[idx]));
							else pMesh_Bot.push_back(dynamic_cast<Mesh*>((*pSMesh)[idx]));
						}
					}
				}
			}
		}
	}

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		paMesh->h, paMesh->meshRect, 
		(MOD_)paMesh->Get_Module_Heff_Display() == MOD_STFIELD, 
		(MOD_)paMesh->Get_Module_Energy_Display() == MOD_STFIELD);
	if (!error)	initialized = true;

	return error;
}

BError Atom_STField::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Atom_STField));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError Atom_STField::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Atom_STField));

#if COMPILECUDA == 1

	if (paMesh->paMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new Atom_STFieldCUDA(paMesh->paMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Atom_STField::UpdateField(void)
{
	if (!paMesh->E.linear_size()) return 0.0;

	double conv = paMesh->M1.h.dim() / MUB;
	
	if (paMesh->STp.get0() != DBL3()) {

		////////////////////////////////////////
		//      FIXED POLARIZATION VERSION    //
		////////////////////////////////////////

#pragma omp parallel for
		for (int idx = 0; idx < paMesh->n.dim(); idx++) {

			double grel = paMesh->grel;
			double mu_s = paMesh->mu_s;
			double flSOT = paMesh->flSOT;
			DBL2 STq = paMesh->STq;
			DBL2 STa = paMesh->STa;
			DBL3 STp = paMesh->STp;
			paMesh->update_parameters_mcoarse(idx, paMesh->grel, grel, paMesh->mu_s, mu_s, paMesh->flSOT, flSOT, paMesh->STq, STq, paMesh->STa, STa, paMesh->STp, STp);

			if (IsNZ(grel)) {

				double dotprod = (paMesh->M1[idx] * STp) / mu_s;
				double neta = STq.i / (STa.i + STa.j * dotprod) + STq.j / (STa.i - STa.j * dotprod);

				int idx_E = paMesh->E.position_to_cellidx(paMesh->M1.cellidx_to_position(idx));
				//z component of J
				double Jc = paMesh->E[idx_E].z * paMesh->elC[idx_E];

				double a_const = -conv * (neta * MUB_E * Jc / (GAMMA * grel)) / (mu_s * mu_s * paMesh->GetMeshDimensions().z);

				DBL3 ST_Field = a_const * ((paMesh->M1[idx] ^ STp) + flSOT * mu_s * STp);
				paMesh->Heff1[idx] += ST_Field;

				if (Module_Heff.linear_size()) Module_Heff[idx] = ST_Field;
			}
		}
	}
	else {

		/////////////////////////////////////////////
		// POLARIZATION FROM TOP AND BOTTOM MESHES //
		/////////////////////////////////////////////

		SZ3 n = paMesh->n;

		//zero module display VECs if needed, since contributions must be added into them to account for possiblility of 2 contributions (top and bottom)
		ZeroModuleVECs();

		if (paMesh_Top.size() || pMesh_Top.size()) {

			//polarization from top mesh
#pragma omp parallel for
			for (int j = 0; j < n.y; j++) {
				for (int i = 0; i < n.x; i++) {

					int cell_idx = i + j * n.x + (n.z - 1) * n.x*n.y;

					//empty cell here ... next
					if (paMesh->M1.is_empty(cell_idx)) continue;

					double grel = paMesh->grel;
					double mu_s = paMesh->mu_s;
					double flSOT2 = paMesh->flSOT2;
					DBL2 STq2 = paMesh->STq2;
					DBL2 STa2 = paMesh->STa2;
					paMesh->update_parameters_mcoarse(cell_idx, paMesh->grel, grel, paMesh->mu_s, mu_s, paMesh->flSOT2, flSOT2, paMesh->STq2, STq2, paMesh->STa2, STa2);

					//effective field for this cell
					DBL3 ST_Field = DBL3();
					bool cell_coupled = false;

					//check all meshes for coupling
					//1. coupling from other atomistic meshes
					for (int mesh_idx = 0; mesh_idx < (int)paMesh_Top.size(); mesh_idx++) {

						Rect tmeshRect = paMesh_Top[mesh_idx]->GetMeshRect();

						//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
						DBL3 cell_rel_pos = DBL3(
							(i + 0.5) * paMesh->h.x + paMesh->meshRect.s.x - tmeshRect.s.x,
							(j + 0.5) * paMesh->h.y + paMesh->meshRect.s.y - tmeshRect.s.y,
							paMesh_Top[mesh_idx]->h.z / 2);

						//can't couple to an empty cell
						if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s) || paMesh_Top[mesh_idx]->M1.is_empty(cell_rel_pos)) continue;

						//get magnetization value in top mesh cell (the polarization)
						DBL3 p = paMesh_Top[mesh_idx]->M1[cell_rel_pos].normalized();
						DBL3 m = paMesh->M1[cell_idx] / mu_s;

						double dotprod = m * p;
						double neta = STq2.i / (STa2.i + STa2.j * dotprod) + STq2.j / (STa2.i - STa2.j * dotprod);

						int idx_E = paMesh->E.position_to_cellidx(paMesh->M1.cellidx_to_position(cell_idx));
						//z component of J
						double Jc = paMesh->E[idx_E].z * paMesh->elC[idx_E];

						double a_const = -conv * (neta * MUB_E * Jc / (GAMMA * grel)) / (mu_s * paMesh->h.z);

						ST_Field = a_const * ((m ^ p) + flSOT2 * p);

						//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
						cell_coupled = true;
						break;
					}

					if (!cell_coupled) {

						//2. coupling from micromagnetic meshes
						for (int mesh_idx = 0; mesh_idx < (int)pMesh_Top.size(); mesh_idx++) {

							Rect tmeshRect = pMesh_Top[mesh_idx]->GetMeshRect();

							//relative coordinates to read value from top mesh (the one we're coupling to here) - relative to top mesh
							DBL3 cell_rel_pos = DBL3(
								(i + 0.5) * paMesh->h.x + paMesh->meshRect.s.x - tmeshRect.s.x,
								(j + 0.5) * paMesh->h.y + paMesh->meshRect.s.y - tmeshRect.s.y,
								pMesh_Top[mesh_idx]->h.z / 2);

							//can't couple to an empty cell
							if (!tmeshRect.contains(cell_rel_pos + tmeshRect.s) || pMesh_Top[mesh_idx]->M.is_empty(cell_rel_pos)) continue;

							//get magnetization value in top mesh cell (the polarization)
							DBL3 p = pMesh_Top[mesh_idx]->M[cell_rel_pos].normalized();
							DBL3 m = paMesh->M1[cell_idx] / mu_s;

							double dotprod = m * p;
							double neta = STq2.i / (STa2.i + STa2.j * dotprod) + STq2.j / (STa2.i - STa2.j * dotprod);

							int idx_E = paMesh->E.position_to_cellidx(paMesh->M1.cellidx_to_position(cell_idx));
							//z component of J
							double Jc = paMesh->E[idx_E].z * paMesh->elC[idx_E];

							double a_const = -conv * (neta * MUB_E * Jc / (GAMMA * grel)) / (mu_s * paMesh->h.z);

							ST_Field = a_const * ((m ^ p) + flSOT2 * p);

							//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
							break;
						}
					}

					paMesh->Heff1[cell_idx] += ST_Field;
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
					if (paMesh->M1.is_empty(cell_idx)) continue;

					double grel = paMesh->grel;
					double mu_s = paMesh->mu_s;
					double flSOT = paMesh->flSOT;
					DBL2 STq = paMesh->STq;
					DBL2 STa = paMesh->STa;
					paMesh->update_parameters_mcoarse(cell_idx, paMesh->grel, grel, paMesh->mu_s, mu_s, paMesh->flSOT, flSOT, paMesh->STq, STq, paMesh->STa, STa);

					//effective field for this cell
					DBL3 ST_Field = DBL3();
					bool cell_coupled = false;

					//check all meshes for coupling
					//1. coupling from other atomistic meshes
					for (int mesh_idx = 0; mesh_idx < (int)paMesh_Bot.size(); mesh_idx++) {

						Rect bmeshRect = paMesh_Bot[mesh_idx]->GetMeshRect();

						//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
						DBL3 cell_rel_pos = DBL3(
							(i + 0.5) * paMesh->h.x + paMesh->meshRect.s.x - bmeshRect.s.x,
							(j + 0.5) * paMesh->h.y + paMesh->meshRect.s.y - bmeshRect.s.y,
							bmeshRect.height() - (paMesh_Bot[mesh_idx]->h.z / 2));

						//can't couple to an empty cell
						if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s) || paMesh_Bot[mesh_idx]->M1.is_empty(cell_rel_pos)) continue;

						//get magnetization value in top mesh cell (the polarization)
						DBL3 p = paMesh_Bot[mesh_idx]->M1[cell_rel_pos].normalized();
						DBL3 m = paMesh->M1[cell_idx] / mu_s;

						double dotprod = m * p;
						double neta = STq.i / (STa.i + STa.j * dotprod) + STq.j / (STa.i - STa.j * dotprod);

						int idx_E = paMesh->E.position_to_cellidx(paMesh->M1.cellidx_to_position(cell_idx));
						//z component of J
						double Jc = paMesh->E[idx_E].z * paMesh->elC[idx_E];

						double a_const = -conv * (neta * MUB_E * Jc / (GAMMA * grel)) / (mu_s * paMesh->h.z);

						ST_Field = a_const * ((m ^ p) + flSOT * p);

						//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
						cell_coupled = true;
						break;
					}

					if (!cell_coupled) {

						//2. coupling from micromagnetic meshes
						for (int mesh_idx = 0; mesh_idx < (int)pMesh_Bot.size(); mesh_idx++) {

							Rect bmeshRect = pMesh_Bot[mesh_idx]->GetMeshRect();

							//relative coordinates to read value from bottom mesh (the one we're coupling to here) - relative to bottom mesh
							DBL3 cell_rel_pos = DBL3(
								(i + 0.5) * paMesh->h.x + paMesh->meshRect.s.x - bmeshRect.s.x,
								(j + 0.5) * paMesh->h.y + paMesh->meshRect.s.y - bmeshRect.s.y,
								bmeshRect.height() - pMesh_Bot[mesh_idx]->h.z / 2);

							//can't couple to an empty cell
							if (!bmeshRect.contains(cell_rel_pos + bmeshRect.s) || pMesh_Bot[mesh_idx]->M.is_empty(cell_rel_pos)) continue;

							//get magnetization value in top mesh cell (the polarization)
							DBL3 p = pMesh_Bot[mesh_idx]->M[cell_rel_pos].normalized();
							DBL3 m = paMesh->M1[cell_idx] / mu_s;

							double dotprod = m * p;
							double neta = STq.i / (STa.i + STa.j * dotprod) + STq.j / (STa.i - STa.j * dotprod);

							int idx_E = paMesh->E.position_to_cellidx(paMesh->M1.cellidx_to_position(cell_idx));
							//z component of J
							double Jc = paMesh->E[idx_E].z * paMesh->elC[idx_E];

							double a_const = -conv * (neta * MUB_E * Jc / (GAMMA * grel)) / (mu_s * paMesh->h.z);

							ST_Field = a_const * ((m ^ p) + flSOT * p);

							//for each cell, either it's not coupled to any other mesh cell (so we never get here), or else it's coupled to exactly one cell on this surface (thus can stop looping over meshes now)
							break;
						}
					}

					paMesh->Heff1[cell_idx] += ST_Field;
					if (Module_Heff.linear_size()) Module_Heff[cell_idx] = ST_Field;
				}
			}
		}
	}

	//don't count this as a contribution to the total energy density
	return 0.0;
}

#endif