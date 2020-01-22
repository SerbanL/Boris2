#include "stdafx.h"
#include "HeatCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_HEAT

#include "Heat.h"
#include "MeshCUDA.h"
#include "Mesh.h"
#include "SuperMesh.h"

HeatCUDA::HeatCUDA(Mesh* pMesh_, SuperMesh* pSMesh_, Heat* pHeat_)
	: ModulesCUDA()
{
	pMesh = pMesh_;
	pMeshCUDA = pMesh->pMeshCUDA;
	pSMesh = pSMesh_;
	pHeat = pHeat_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	if (!error_on_create) error_on_create = temp_cmbnd_funcs()->set_pointers(pMeshCUDA);
}

HeatCUDA::~HeatCUDA()
{
	//Transport managed quantities must be transferred back to their cpu versions - values only as sizes must match
	if (Holder_Module_Available()) {

		pMeshCUDA->Temp()->copy_to_cpuvec(pMesh->Temp);
	}

	//clear mesh quantities as they're not used any more
	pMeshCUDA->Temp()->clear();
}


//-------------------Abstract base class method implementations

BError HeatCUDA::Initialize(void)
{
	BError error(CLASS_STR(HeatCUDA));

	//no energy density contribution here
	ZeroEnergy();

	pHeat->SetRobinBoundaryConditions();

	initialized = true;

	return error;
}

BError HeatCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(HeatCUDA));

	Uninitialize();
	
	bool success = true;
	
	//need this when we switch cuda mode
	if (!Q_equation.is_set() && pHeat->Q_equation.is_set()) error = SetQEquation(pHeat->Q_equation.get_scalar_fspec());

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE, UPDATECONFIG_MODULEADDED, UPDATECONFIG_MODULEDELETED)) {

		//update mesh dimensions in equation constants
		if (Q_equation.is_set()) {

			error = SetQEquation(pHeat->Q_equation.get_scalar_fspec());
		}

		if (pMeshCUDA->M()->size_cpu().dim()) {

			//in a ferromagnet the magnetization sets the shape only on initialization. If already initialized then shape already set.
			if (pMeshCUDA->Temp()->size_cpu().dim()) {

				success = pMeshCUDA->Temp()->resize(pMesh->h_t, pMesh->meshRect);
			}
			else {

				success = pMeshCUDA->Temp()->set_from_cpuvec(pMesh->Temp);
			}
		}
		else if (pMeshCUDA->elC()->size_cpu().dim()) {

			//in a normal metal the electrical conductivity sets the shape

			//before doing this must make sure elC was itself set in the Transport module (it could be this module is being updated before the Transport module)
			//Transport module is guaranteed to be set otherwise elC would have zero size - it does mean Transport has UpdateConfiguration called twice but it doesn't matter.
			if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

				error = (*pMesh)(MOD_TRANSPORT)->UpdateConfiguration(cfgMessage);
				if (error) return error;
			}

			if (success) {

				if (pMeshCUDA->Temp()->size_cpu().dim()) {

					success = pMeshCUDA->Temp()->resize(pMesh->h_t, pMesh->meshRect, (cuVEC_VC<cuBReal>&)pMeshCUDA->elC);
				}
				else {

					success = pMeshCUDA->Temp()->set_from_cpuvec(pMesh->Temp);
				}
			}
		}
		else {

			//in an insulator (or conductor with Transport module not set) the shape is set directly in Temp
			if (pMeshCUDA->Temp()->size_cpu().dim()) {

				success = pMeshCUDA->Temp()->resize(pMesh->h_t, pMesh->meshRect);
			}
			else {

				success = pMeshCUDA->Temp()->set_from_cpuvec(pMesh->Temp);
			}
		}

		//allocate memory for the heatEq_RHS auxiliary vector
		if (success) {

			if (heatEq_RHS.size() != pMesh->n_t.dim()) {

				success = heatEq_RHS.resize(pMesh->n_t.dim());
				if (success) heatEq_RHS.set(0.0);
			}
		}
	}

	if (!success) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	
	return error;
}

void HeatCUDA::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		if (Q_equation.is_set()) {

			SetQEquation(pHeat->Q_equation.get_scalar_fspec());
		}
	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

		if (Q_equation.is_set()) Q_equation.clear();
	}
}

void HeatCUDA::UpdateField(void)
{
}

#endif

#endif