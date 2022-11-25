#include "stdafx.h"
#include "SMElasticCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_MELASTIC

#include "SMElastic.h"
#include "SuperMesh.h"

#include "MElasticCUDA.h"

SMElasticCUDA::SMElasticCUDA(SuperMesh* pSMesh_, SMElastic* pSMElastic_) :
	ModulesCUDA()
{
	pSMesh = pSMesh_;
	pSMElastic = pSMElastic_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
}

SMElasticCUDA::~SMElasticCUDA()
{

}

//-------------------Abstract base class method implementations

BError SMElasticCUDA::Initialize(void)
{
	BError error(CLASS_STR(SHeaSMElasticCUDAtCUDA));

	ZeroEnergy();

	error = pSMElastic->Initialize();
	if (error) return error;

	//clear everything then rebuild
	pMElastic.clear();
	CMBNDcontactsCUDA.clear();
	CMBNDcontacts.clear();
	pu_disp.clear();

	//now build pMElastic (and pu_disp)
	for (int idx = 0; idx < pSMesh->size(); idx++) {

		if ((*pSMesh)[idx]->IsModuleSet(MOD_MELASTIC)) {

			pMElastic.push_back(dynamic_cast<MElasticCUDA*>((*pSMesh)[idx]->GetCUDAModule(MOD_MELASTIC)));
			pu_disp.push_back(&(*pSMesh)[idx]->pMeshBaseCUDA->u_disp);
		}
	}

	//set cmbnd flags
	for (int idx = 0; idx < (int)pMElastic.size(); idx++) {

		//it's easier to just copy the flags entirely from the cpu versions.
		//Notes :
		//1. By design the cpu versions are required to keep size and flags up to date (but not mesh values)
		//2. pMElastic in SMelastic has exactly the same size and order
		//3. SMelastic UpdateConfiguration was called just before, which called this CUDA version at the end.

		if (!(*pu_disp[idx])()->copyflags_from_cpuvec(*pSMElastic->pu_disp[idx])) error(BERROR_GPUERROR_CRIT);
	}

	for (int idx = 0; idx < pSMElastic->CMBNDcontacts.size(); idx++) {

		std::vector<cu_obj<CMBNDInfoCUDA>> mesh_contacts;
		std::vector<CMBNDInfoCUDA> mesh_contacts_cpu;

		for (int idx_contact = 0; idx_contact < pSMElastic->CMBNDcontacts[idx].size(); idx_contact++) {

			cu_obj<CMBNDInfoCUDA> contact;

			contact()->copy_from_CMBNDInfo<CMBNDInfo>(pSMElastic->CMBNDcontacts[idx][idx_contact]);

			mesh_contacts.push_back(contact);

			mesh_contacts_cpu.push_back(pSMElastic->CMBNDcontacts[idx][idx_contact]);
		}

		CMBNDcontactsCUDA.push_back(mesh_contacts);
		CMBNDcontacts.push_back(mesh_contacts_cpu);
	}

	initialized = true;

	return error;
}

BError SMElasticCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SMElasticCUDA));

	Uninitialize();

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED, UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED, UPDATECONFIG_SWITCHCUDASTATE)) {

		//clear everything then rebuild
		pMElastic.clear();
		pu_disp.clear();

		//now build pMElastic (and pu_disp)
		for (int idx = 0; idx < pSMesh->size(); idx++) {

			if ((*pSMesh)[idx]->IsModuleSet(MOD_MELASTIC)) {

				pMElastic.push_back(dynamic_cast<MElasticCUDA*>((*pSMesh)[idx]->GetCUDAModule(MOD_MELASTIC)));
				pu_disp.push_back(&(*pSMesh)[idx]->pMeshBaseCUDA->u_disp);
			}
		}
	}

	return error;
}

void SMElasticCUDA::UpdateField(void)
{
	//only need to update this after an entire magnetization equation time step is solved
	//also if el_dT is set to zero skip the elastic solver : this will maintain a fixed strain
	if (!pSMesh->CurrentTimeStepSolved() || pSMElastic->el_dT < MINTIMESTEP) {

		for (int idx = 0; idx < pMElastic.size(); idx++) {

			//calculate effective field in each mesh
			pMElastic[idx]->Calculate_MElastic_Field();
		}

		return;
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	//
	//If magneto-elastic solver is enabled (el_dT not zero) then iterate solver

	if (pSMElastic->el_dT >= MINTIMESTEP) {

		//using el_dT time steps until the full magnetic dT is covered:
		//1a. evaluate v_vel for next el_dT time step
		//1b. evaluate sig_d, sig_od for next el_dT time step using updated v_vel (v_vel and sig are staggered on the time axis, so immediately use result from 1a)
		//1c. calculate strain using stress, by inverting the stifness coefficients matrix (cubic crystal symmetry)

		//2. calculate magnetoelastic effective field using calculated strain

		////////////////////////////////
		//1.
		double dT = pSMElastic->el_dT;

		//number of sub_steps to cover magnetic_dT required when advancing in smaller el_dT steps
		int sub_steps = (int)floor_epsilon(pSMElastic->magnetic_dT / pSMElastic->el_dT);

		//any left-over epsilon_dT < el_dT
		double epsilon_dT = pSMElastic->magnetic_dT - pSMElastic->el_dT * sub_steps;

		for (int step_idx = 0; step_idx < sub_steps + 1; step_idx++) {

			//the last step may have a different time step - take this epsilon_dT step (if not zero)
			if (step_idx == sub_steps) {

				if (epsilon_dT) dT = epsilon_dT;
				else continue;
			}

			//1a. Update velocity (and displacement for display)
			for (int idx = 0; idx < pMElastic.size(); idx++) {

				pMElastic[idx]->Iterate_Elastic_Solver_Velocity(dT);
			}
			
			//CMBND for continuous velocity
			for (int idx1 = 0; idx1 < (int)CMBNDcontacts.size(); idx1++) {
				for (int idx2 = 0; idx2 < (int)CMBNDcontacts[idx1].size(); idx2++) {

					//contact descriptor
					CMBNDInfoCUDA& contact = CMBNDcontacts[idx1][idx2];

					//mexh indexes
					int idx_sec = contact.mesh_idx.i;
					int idx_pri = contact.mesh_idx.j;

					size_t size = contact.cells_box.size().dim();

					int axis;
					if (contact.cell_shift.x) axis = 1;
					else if (contact.cell_shift.y) axis = 2;
					else axis = 3;

					pMElastic[idx_pri]->make_velocity_continuous(
						size, axis, CMBNDcontactsCUDA[idx1][idx2],
						pMElastic[idx_sec]->vx, pMElastic[idx_sec]->vy, pMElastic[idx_sec]->vz, pMElastic[idx_sec]->pMeshCUDA->u_disp);
				}
			}

			//1b. Update stress
			for (int idx = 0; idx < pMElastic.size(); idx++) {

				pMElastic[idx]->Iterate_Elastic_Solver_Stress(dT);
			}
			
			//CMBND for continuous stress
			for (int idx1 = 0; idx1 < (int)CMBNDcontacts.size(); idx1++) {
				for (int idx2 = 0; idx2 < (int)CMBNDcontacts[idx1].size(); idx2++) {

					//contact descriptor
					CMBNDInfoCUDA& contact = CMBNDcontacts[idx1][idx2];

					//mexh indexes
					int idx_sec = contact.mesh_idx.i;
					int idx_pri = contact.mesh_idx.j;

					size_t size = contact.cells_box.size().dim();

					int axis;
					if (contact.cell_shift.x) axis = 1;
					else if (contact.cell_shift.y) axis = 2;
					else axis = 3;

					pMElastic[idx_pri]->make_stress_continuous(
						size, axis, CMBNDcontactsCUDA[idx1][idx2],
						pMElastic[idx_sec]->sdd, pMElastic[idx_sec]->sxy, pMElastic[idx_sec]->sxz, pMElastic[idx_sec]->syz,
						pMElastic[idx_sec]->pMeshCUDA->u_disp);
				}
			}
		}

		//1c. Update strain from stress
		for (int idx = 0; idx < pMElastic.size(); idx++) {

			pMElastic[idx]->Calculate_Strain_From_Stress();
		}

		////////////////////////////////
		//2.

		for (int idx = 0; idx < pMElastic.size(); idx++) {

			pMElastic[idx]->Calculate_MElastic_Field();
		}

		//update the magnetic dT that will be used next time around to increment the elastic solver by
		pSMElastic->magnetic_dT = pSMesh->GetTimeStep();
		if (pSMElastic->linked_el_dT) pSMElastic->el_dT = pSMElastic->magnetic_dT;
	}
}

#endif

#endif