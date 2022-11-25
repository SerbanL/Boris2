#include "stdafx.h"
#include "SMElastic.h"

#ifdef MODULE_COMPILATION_MELASTIC

#include "SuperMesh.h"
#include "MElastic.h"
#include "MElastic_Boundaries.h"

SMElastic::SMElastic(SuperMesh *pSMesh_) :
	Modules(),
	ProgramStateNames(this, { VINFO(el_dT), VINFO(linked_el_dT), VINFO(fixed_u_surfaces), VINFO(stress_surfaces_rect), VINFO(stress_surfaces_equations) }, {})
{
	pSMesh = pSMesh_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pSMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

//-------------------Abstract base class method implementations

BError SMElastic::Initialize(void)
{
	BError error(CLASS_STR(SMElastic));

	energy = 0.0;

	//Must have at least one fixed surface defined.
	if (!fixed_u_surfaces.size()) return error(BERROR_INCORRECTCONFIG);

	//el_dT must be set correctly using the magnetic time step
	magnetic_dT = pSMesh->GetTimeStep();
	if (linked_el_dT) el_dT = magnetic_dT;

	//clear everything then rebuild
	pMElastic.clear();
	pu_disp.clear();
	CMBNDcontacts.clear();

	//now build pMElastic (and pu_disp)
	for (int idx = 0; idx < pSMesh->size(); idx++) {

		if ((*pSMesh)[idx]->IsModuleSet(MOD_MELASTIC)) {

			pMElastic.push_back(dynamic_cast<MElastic*>((*pSMesh)[idx]->GetModule(MOD_MELASTIC)));
			pu_disp.push_back(&(*pSMesh)[idx]->u_disp);
		}
	}
	
	//set cmbnd flags (also building contacts)
	for (int idx = 0; idx < (int)pMElastic.size(); idx++) {

		//build CMBND contacts and set flags for u_disp (we don't actually need to use cmbnd flags in u_disp as boundary values are just over-written after to ensure continuity)
		CMBNDcontacts.push_back(pu_disp[idx]->set_cmbnd_flags(idx, pu_disp));
	}
	
	initialized = true;

	return error;
}

BError SMElastic::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SMElastic));

	Uninitialize();

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHADDED, UPDATECONFIG_MESHDELETED, UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED, UPDATECONFIG_SWITCHCUDASTATE)) {

		//clear everything then rebuild
		pMElastic.clear();
		pu_disp.clear();

		//now build pMElastic (and pu_disp)
		for (int idx = 0; idx < pSMesh->size(); idx++) {

			if ((*pSMesh)[idx]->IsModuleSet(MOD_MELASTIC)) {

				pMElastic.push_back(dynamic_cast<MElastic*>((*pSMesh)[idx]->GetModule(MOD_MELASTIC)));
				pu_disp.push_back(&(*pSMesh)[idx]->u_disp);
			}
		}

		for (int idx = 0; idx < pMElastic.size(); idx++) {

			//setup SMElastic pointer in each MElastic module
			pMElastic[idx]->pSMEl = this;
		}
	}

	//el_dT must be set correctly using the magnetic time step
	magnetic_dT = pSMesh->GetTimeStep();
	if (linked_el_dT) el_dT = magnetic_dT;

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError SMElastic::MakeCUDAModule(void)
{
	BError error(CLASS_STR(SMElastic));

#if COMPILECUDA == 1

	pModuleCUDA = new SMElasticCUDA(pSMesh, this);
	error = pModuleCUDA->Error_On_Create();

#endif

	return error;
}

double SMElastic::UpdateField(void)
{
	//only need to update this after an entire magnetization equation time step is solved
	//also if el_dT is set to zero skip the elastic solver : this will maintain a fixed strain
	if (!pSMesh->CurrentTimeStepSolved() || el_dT < MINTIMESTEP) {

		for (int idx = 0; idx < pMElastic.size(); idx++) {

			//calculate effective field in each mesh
			pMElastic[idx]->Calculate_MElastic_Field();
		}

		return 0.0;
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	//
	//If magneto-elastic solver is enabled (el_dT not zero) then iterate solver

	if (el_dT >= MINTIMESTEP) {

		//using el_dT time steps until the full magnetic dT is covered:
		//1a. evaluate v_vel for next el_dT time step
		//1b. evaluate sig_d, sig_od for next el_dT time step using updated v_vel (v_vel and sig are staggered on the time axis, so immediately use result from 1a)
		//1c. calculate strain using stress, by inverting the stifness coefficients matrix (cubic crystal symmetry)

		//2. calculate magnetoelastic effective field using calculated strain

		////////////////////////////////
		//1.
		double dT = el_dT;

		//number of sub_steps to cover magnetic_dT required when advancing in smaller el_dT steps
		int sub_steps = (int)floor_epsilon(magnetic_dT / el_dT);

		//any left-over epsilon_dT < el_dT
		double epsilon_dT = magnetic_dT - el_dT * sub_steps;

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
					CMBNDInfo& contact = CMBNDcontacts[idx1][idx2];

					//mexh indexes
					int idx_sec = contact.mesh_idx.i;
					int idx_pri = contact.mesh_idx.j;

					pMElastic[idx_pri]->make_velocity_continuous(
						contact, 
						pMElastic[idx_sec]->vx, pMElastic[idx_sec]->vy, pMElastic[idx_sec]->vz, pMElastic[idx_sec]->pMesh->u_disp);
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
					CMBNDInfo& contact = CMBNDcontacts[idx1][idx2];

					//mexh indexes
					int idx_sec = contact.mesh_idx.i;
					int idx_pri = contact.mesh_idx.j;

					pMElastic[idx_pri]->make_stress_continuous(
						contact, 
						pMElastic[idx_sec]->sdd,
						pMElastic[idx_sec]->sxy, pMElastic[idx_sec]->sxz, pMElastic[idx_sec]->syz,
						pMElastic[idx_sec]->pMesh->u_disp);
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
		magnetic_dT = pSMesh->GetTimeStep();
		if (linked_el_dT) el_dT = magnetic_dT;
	}

	return 0.0;
}

void SMElastic::set_linked_el_dT(bool flag)
{
	linked_el_dT = flag;
	if (linked_el_dT) el_dT = pSMesh->GetTimeStep();
}

//------------------- Fixed and Stress Surfaces

BError SMElastic::Add_Fixed_Surface(Rect surface_rect)
{
	BError error(CLASS_STR(SMElastic));

	if (!surface_rect.IsPlane()) return error(BERROR_INCORRECTVALUE);

	fixed_u_surfaces.push_back(surface_rect);

	return error;
}

BError SMElastic::Add_Stress_Surface(Rect surface_rect, std::string equation)
{
	BError error(CLASS_STR(SMElastic));

	if (!surface_rect.IsPlane()) return error(BERROR_INCORRECTVALUE);

	stress_surfaces_rect.push_back(surface_rect);
	stress_surfaces_equations.push_back(equation);

	return error;
}

#endif