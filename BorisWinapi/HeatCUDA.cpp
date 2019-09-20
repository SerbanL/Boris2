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

	error_on_create = UpdateConfiguration();

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

	initialized = true;

	//no energy density contribution here
	ZeroEnergy();

	return error;
}

BError HeatCUDA::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(HeatCUDA));

	Uninitialize();
	
	bool success = true;
	
	if (pMeshCUDA->M()->size_cpu().dim()) {

		//in a ferromagnet the magnetization sets the shape only on initialization. If already initialized then shape already set.
		if (pMeshCUDA->Temp()->size_cpu().dim())
			success = pMeshCUDA->Temp()->resize(pMesh->h_t, pMesh->meshRect);
		else
			success = pMeshCUDA->Temp()->set_from_cpuvec(pMesh->Temp);
	}
	else if (pMeshCUDA->elC()->size_cpu().dim()) {

		//in a normal metal the electrical conductivity sets the shape

		//before doing this must make sure elC was itself set in the Transport module (it could be this module is being updated before the Transport module)
		//Transport module is guaranteed to be set otherwise elC would have zero size - it does mean Transport has UpdateConfiguration called twice but it doesn't matter.
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			error = pMesh->pMod(MOD_TRANSPORT)->UpdateConfiguration();
			if (error) return error;
		}

		if (success) {

			if (pMeshCUDA->Temp()->size_cpu().dim())
				success = pMeshCUDA->Temp()->resize(pMesh->h_t, pMesh->meshRect, (cuVEC_VC<cuBReal>&)pMeshCUDA->elC);
			else
				success = pMeshCUDA->Temp()->set_from_cpuvec(pMesh->Temp);
		}
	}
	else {

		//in an insulator (or conductor with Transport module not set) the shape is set directly in Temp
		if (pMeshCUDA->Temp()->size_cpu().dim())
			success = pMeshCUDA->Temp()->resize(pMesh->h_t, pMesh->meshRect);
		else
			success = pMeshCUDA->Temp()->set_from_cpuvec(pMesh->Temp);
	}
	
	//allocate memory for the heatEq_RHS auxiliary vector
	if (success) {

		if (heatEq_RHS.size() != pMesh->n_t.dim()) {

			success = heatEq_RHS.resize(pMesh->n_t.dim());
			if (success) heatEq_RHS.set(0.0);
		}
	}
	
	pHeat->SetRobinBoundaryConditions();
	
	if (!success) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	
	return error;
}

void HeatCUDA::UpdateField(void)
{
}

#endif

#endif