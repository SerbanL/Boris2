#include "stdafx.h"

#include "SuperMeshCUDA.h"
#include "SuperMesh.h"
#include "PhysQRep.h"

#if COMPILECUDA == 1

SuperMeshCUDA::SuperMeshCUDA(SuperMesh* pSMesh_) :
	MeshDisplayCUDA(),
	sMeshRect(pSMesh_->sMeshRect),
	sMeshRect_fm(pSMesh_->sMeshRect_fm), n_fm(pSMesh_->n_fm), h_fm(pSMesh_->h_fm),
	sMeshRect_e(pSMesh_->sMeshRect_e), n_e(pSMesh_->n_e), h_e(pSMesh_->h_e)
{
	pSMesh = pSMesh_;
}

SuperMeshCUDA::~SuperMeshCUDA()
{
}

//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

vector<PhysQ> SuperMeshCUDA::FetchOnScreenPhysicalQuantity(double detail_level)
{
	vector<PhysQ> physQ;
	
	//get anything displayed on super-mesh
	switch (pSMesh->displayedPhysicalQuantity) {

	case MESHDISPLAY_SM_DEMAG:

		if (pSMesh->IsSuperMeshModuleSet(MODS_SDEMAG)) {
			
			if (prepare_display(n_fm, sMeshRect_fm, detail_level, reinterpret_cast<SDemag*>(pSMesh->pSMod(MODS_SDEMAG))->GetDemagFieldCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				physQ.push_back(PhysQ(pdisplay_vec_vec, pSMesh->displayedPhysicalQuantity).set_focus(true, pSMesh->superMeshHandle));
			}
		}
		break;
		
	case MESHDISPLAY_SM_OERSTED:

		if (pSMesh->IsSuperMeshModuleSet(MODS_OERSTED)) {

			if (prepare_display(n_e, sMeshRect_e, detail_level, reinterpret_cast<Oersted*>(pSMesh->pSMod(MODS_OERSTED))->GetOerstedFieldCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				physQ.push_back(PhysQ(pdisplay_vec_vec, pSMesh->displayedPhysicalQuantity).set_focus(true, pSMesh->superMeshHandle));
			}
		}
		break;
		
	case MESHDISPLAY_SM_STRAYH:

		if (pSMesh->IsSuperMeshModuleSet(MODS_STRAYFIELD)) {

			if (prepare_display(n_fm, sMeshRect_fm, detail_level, reinterpret_cast<StrayField*>(pSMesh->pSMod(MODS_STRAYFIELD))->GetStrayFieldCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				physQ.push_back(PhysQ(pdisplay_vec_vec, pSMesh->displayedPhysicalQuantity).set_focus(true, pSMesh->superMeshHandle));
			}
		}
		break;
		
	}

	return physQ;
}

//----------------------------------- Getters

//check if ODE solver needs spin accumulation solved
bool SuperMeshCUDA::SolveSpinCurrent(void)
{
	return pSMesh->SolveSpinCurrent();
}

#endif