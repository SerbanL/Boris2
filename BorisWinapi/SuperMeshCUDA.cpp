#include "stdafx.h"

#include "SuperMeshCUDA.h"
#include "SuperMesh.h"
#include "PhysQRep.h"
#include "BorisLib.h"

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
				physQ.push_back(PhysQ(pdisplay_vec_vec, pSMesh->displayedPhysicalQuantity, (VEC3REP_)pSMesh->vec3rep).set_focus(true, pSMesh->superMeshHandle));
			}
		}
		break;
		
	case MESHDISPLAY_SM_OERSTED:

		if (pSMesh->IsSuperMeshModuleSet(MODS_OERSTED)) {

			if (prepare_display(n_e, sMeshRect_e, detail_level, reinterpret_cast<Oersted*>(pSMesh->pSMod(MODS_OERSTED))->GetOerstedFieldCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				physQ.push_back(PhysQ(pdisplay_vec_vec, pSMesh->displayedPhysicalQuantity, (VEC3REP_)pSMesh->vec3rep).set_focus(true, pSMesh->superMeshHandle));
			}
		}
		break;
		
	case MESHDISPLAY_SM_STRAYH:

		if (pSMesh->IsSuperMeshModuleSet(MODS_STRAYFIELD)) {

			if (prepare_display(n_fm, sMeshRect_fm, detail_level, reinterpret_cast<StrayField*>(pSMesh->pSMod(MODS_STRAYFIELD))->GetStrayFieldCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				physQ.push_back(PhysQ(pdisplay_vec_vec, pSMesh->displayedPhysicalQuantity, (VEC3REP_)pSMesh->vec3rep).set_focus(true, pSMesh->superMeshHandle));
			}
		}
		break;
		
	}

	return physQ;
}

//save the quantity currently displayed on screen in an ovf2 file using the specified format
BError SuperMeshCUDA::SaveOnScreenPhysicalQuantity(string fileName, string ovf2_dataType)
{
	BError error(__FUNCTION__);

	OVF2 ovf2;

	switch (pSMesh->displayedPhysicalQuantity) {

	case MESHDISPLAY_NONE:
		error = pSMesh->active_mesh()->SaveOnScreenPhysicalQuantity(fileName, ovf2_dataType);
		break;

	case MESHDISPLAY_SM_DEMAG:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_SDEMAG)) {

			prepare_display(n_fm, sMeshRect_fm, h_fm.mindim(), reinterpret_cast<SDemag*>(pSMesh->pSMod(MODS_SDEMAG))->GetDemagFieldCUDA());
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		break;

	case MESHDISPLAY_SM_OERSTED:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_OERSTED)) {

			prepare_display(n_e, sMeshRect_e, h_e.mindim(), reinterpret_cast<Oersted*>(pSMesh->pSMod(MODS_OERSTED))->GetOerstedFieldCUDA());
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		break;

	case MESHDISPLAY_SM_STRAYH:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRAYFIELD)) {

			prepare_display(n_fm, sMeshRect_fm, h_fm.mindim(), reinterpret_cast<StrayField*>(pSMesh->pSMod(MODS_STRAYFIELD))->GetStrayFieldCUDA());
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		break;
	}

	return error;
}

//Before calling a run of GetDisplayedMeshValue, make sure to call PrepareDisplayedMeshValue : this calculates and stores in displayVEC storage and quantities which don't have memory allocated directly, but require computation and temporary storage.
void SuperMeshCUDA::PrepareDisplayedMeshValue(void)
{
	switch (pSMesh->displayedPhysicalQuantity) {

	case MESHDISPLAY_SM_DEMAG:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_SDEMAG)) {

			prepare_display(n_fm, sMeshRect_fm, h_fm.mindim(), reinterpret_cast<SDemag*>(pSMesh->pSMod(MODS_SDEMAG))->GetDemagFieldCUDA());
		}
		break;

	case MESHDISPLAY_SM_OERSTED:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_OERSTED)) {

			prepare_display(n_e, sMeshRect_e, h_e.mindim(), reinterpret_cast<Oersted*>(pSMesh->pSMod(MODS_OERSTED))->GetOerstedFieldCUDA());
		}
		break;

	case MESHDISPLAY_SM_STRAYH:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRAYFIELD)) {

			prepare_display(n_fm, sMeshRect_fm, h_fm.mindim(), reinterpret_cast<StrayField*>(pSMesh->pSMod(MODS_STRAYFIELD))->GetStrayFieldCUDA());
		}
		break;
	}
}

//return value of currently displayed mesh quantity at the given absolute position; the value is read directly from the storage VEC, not from the displayed PhysQ.
//Return an Any as the displayed quantity could be either a scalar or a vector.
Any SuperMeshCUDA::GetDisplayedMeshValue(DBL3 abs_pos)
{
	//get anything displayed on super-mesh
	switch (pSMesh->displayedPhysicalQuantity) {

		////////////////
		//no quantity displayed on the supermesh, so use individual mesh displayed quantities
		////////////////

	case MESHDISPLAY_NONE:
	{
		//find which mesh holds abs_pos, if any, and return value displayed for that mesh
		for (int idx = 0; idx < pSMesh->pMesh.size(); idx++) {

			if (pSMesh->pMesh[idx]->meshRect.contains(abs_pos)) return pSMesh->pMesh[idx]->GetDisplayedMeshValue(abs_pos);
		}
	}
	break;

	////////////////
	//use a quantity displayed on the supermesh
	////////////////

	case MESHDISPLAY_SM_DEMAG:
		return (*pdisplay_vec_vec)[abs_pos - sMeshRect_fm.s];
		break;

	case MESHDISPLAY_SM_OERSTED:
		return (*pdisplay_vec_vec)[abs_pos - sMeshRect_e.s];
		break;

	case MESHDISPLAY_SM_STRAYH:
		return (*pdisplay_vec_vec)[abs_pos - sMeshRect_fm.s];
		break;
	}

	return Any();
}

//return average value for currently displayed mesh quantity in the given relative rectangle
Any SuperMeshCUDA::GetAverageDisplayedMeshValue(Rect rel_rect)
{
	//get anything displayed on super-mesh
	switch (pSMesh->displayedPhysicalQuantity) {

		////////////////
		//no quantity displayed on the supermesh, so use individual mesh displayed quantities
		////////////////

	case MESHDISPLAY_NONE:
	{
		//get average value from currently focused mesh instead
		return pSMesh->active_mesh()->GetAverageDisplayedMeshValue(rel_rect);
	}
	break;

	////////////////
	//use a quantity displayed on the supermesh
	////////////////

	case MESHDISPLAY_SM_DEMAG:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_SDEMAG)) {

			prepare_display(n_fm, sMeshRect_fm, h_fm.mindim(), reinterpret_cast<SDemag*>(pSMesh->pSMod(MODS_SDEMAG))->GetDemagFieldCUDA());
			return pdisplay_vec_vec->average_nonempty_omp(rel_rect);
		}
		break;

	case MESHDISPLAY_SM_OERSTED:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_OERSTED)) {

			prepare_display(n_e, sMeshRect_e, h_e.mindim(), reinterpret_cast<Oersted*>(pSMesh->pSMod(MODS_OERSTED))->GetOerstedFieldCUDA());
			return pdisplay_vec_vec->average_nonempty_omp(rel_rect);
		}
		break;

	case MESHDISPLAY_SM_STRAYH:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRAYFIELD)) {

			prepare_display(n_fm, sMeshRect_fm, h_fm.mindim(), reinterpret_cast<StrayField*>(pSMesh->pSMod(MODS_STRAYFIELD))->GetStrayFieldCUDA());
			return pdisplay_vec_vec->average_nonempty_omp(rel_rect);
		}
		break;
	}

	return Any();
}

//----------------------------------- Getters

//check if ODE solver needs spin accumulation solved
bool SuperMeshCUDA::SolveSpinCurrent(void)
{
	return pSMesh->SolveSpinCurrent();
}

#endif