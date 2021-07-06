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

std::vector<PhysQ> SuperMeshCUDA::FetchOnScreenPhysicalQuantity(double detail_level)
{
	std::vector<PhysQ> physQ;
	
	//get anything displayed on super-mesh
	switch (pSMesh->displayedPhysicalQuantity) {

	case MESHDISPLAY_SM_DEMAG:

		if (pSMesh->IsSuperMeshModuleSet(MODS_SDEMAG)) {
			
			if (prepare_display(n_fm, sMeshRect_fm, detail_level, dynamic_cast<SDemag*>(pSMesh->pSMod(MODS_SDEMAG))->GetDemagFieldCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				physQ.push_back(PhysQ(pdisplay_vec_vec, pSMesh->displayedPhysicalQuantity, (VEC3REP_)pSMesh->vec3rep).set_focus(true, pSMesh->superMeshHandle));
			}
		}
		break;
		
	case MESHDISPLAY_SM_OERSTED:

		if (pSMesh->IsSuperMeshModuleSet(MODS_OERSTED)) {

			if (prepare_display(n_e, sMeshRect_e, detail_level, dynamic_cast<Oersted*>(pSMesh->pSMod(MODS_OERSTED))->GetOerstedFieldCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				physQ.push_back(PhysQ(pdisplay_vec_vec, pSMesh->displayedPhysicalQuantity, (VEC3REP_)pSMesh->vec3rep).set_focus(true, pSMesh->superMeshHandle));
			}
		}
		break;
		
	case MESHDISPLAY_SM_STRAYH:

		if (pSMesh->IsSuperMeshModuleSet(MODS_STRAYFIELD)) {

			if (prepare_display(n_fm, sMeshRect_fm, detail_level, dynamic_cast<StrayField*>(pSMesh->pSMod(MODS_STRAYFIELD))->GetStrayFieldCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				physQ.push_back(PhysQ(pdisplay_vec_vec, pSMesh->displayedPhysicalQuantity, (VEC3REP_)pSMesh->vec3rep).set_focus(true, pSMesh->superMeshHandle));
			}
		}
		break;
		
	}

	return physQ;
}

//save the quantity currently displayed on screen for named mesh in an ovf2 file using the specified format
BError SuperMeshCUDA::SaveOnScreenPhysicalQuantity(std::string meshName, std::string fileName, std::string ovf2_dataType)
{
	BError error(__FUNCTION__);

	if (!pSMesh->contains(meshName) && meshName != pSMesh->superMeshHandle) return error(BERROR_INCORRECTNAME);

	OVF2 ovf2;

	switch (pSMesh->displayedPhysicalQuantity) {

	case MESHDISPLAY_NONE:
		if (meshName != pSMesh->superMeshHandle) error = pSMesh->pMesh[meshName]->SaveOnScreenPhysicalQuantity(fileName, ovf2_dataType);
		break;

	case MESHDISPLAY_SM_DEMAG:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_SDEMAG)) {

			prepare_display(n_fm, sMeshRect_fm, h_fm.mindim(), dynamic_cast<SDemag*>(pSMesh->pSMod(MODS_SDEMAG))->GetDemagFieldCUDA());
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		break;

	case MESHDISPLAY_SM_OERSTED:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_OERSTED)) {

			prepare_display(n_e, sMeshRect_e, h_e.mindim(), dynamic_cast<Oersted*>(pSMesh->pSMod(MODS_OERSTED))->GetOerstedFieldCUDA());
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		break;

	case MESHDISPLAY_SM_STRAYH:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRAYFIELD)) {

			prepare_display(n_fm, sMeshRect_fm, h_fm.mindim(), dynamic_cast<StrayField*>(pSMesh->pSMod(MODS_STRAYFIELD))->GetStrayFieldCUDA());
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		break;
	}

	return error;
}

//extract profile from named mesh, from currently display mesh quantity, but reading directly from the quantity
//Displayed mesh quantity can be scalar or a vector; pass in std::vector pointers, then check for nullptr to determine what type is displayed
//if do_average = true then build average and don't return anything, else return just a single-shot profile. If read_average = true then simply read out the internally stored averaged profile by assigning to pointer.
void SuperMeshCUDA::GetPhysicalQuantityProfile(DBL3 start, DBL3 end, double step, DBL3 stencil, std::vector<DBL3>*& pprofile_dbl3, std::vector<double>*& pprofile_dbl, std::string meshName, bool do_average, bool read_average)
{
	size_t size = round((end - start).norm() / step) + 1;

	//get anything displayed on super-mesh
	switch (pSMesh->displayedPhysicalQuantity) {

		////////////////
		//no quantity displayed on the supermesh, so use individual mesh displayed quantities
		////////////////

	case MESHDISPLAY_NONE:
		if (pSMesh->contains(meshName)) pSMesh->pMesh[meshName]->GetPhysicalQuantityProfile(start, end, step, stencil, pprofile_dbl3, pprofile_dbl, do_average, read_average);
		break;

	////////////////
	//use a quantity displayed on the supermesh
	////////////////

	case MESHDISPLAY_SM_DEMAG:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_SDEMAG)) {

			if (profile_storage_vec.size() != size) { if (!profile_storage_vec.resize(size)) return; }
			dynamic_cast<SDemag*>(pSMesh->pSMod(MODS_SDEMAG))->GetDemagFieldCUDA()()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil, profile_storage_vec);

			if (pSMesh->profile_storage_dbl3.size() != size) { pSMesh->profile_storage_dbl3.resize(size); }
			profile_storage_vec.copy_to_vector(pSMesh->profile_storage_dbl3);
			pprofile_dbl3 = &pSMesh->profile_storage_dbl3;
		}
		break;

	case MESHDISPLAY_SM_OERSTED:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_OERSTED)) {

			if (profile_storage_vec.size() != size) { if (!profile_storage_vec.resize(size)) return; }
			dynamic_cast<Oersted*>(pSMesh->pSMod(MODS_OERSTED))->GetOerstedFieldCUDA()()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil, profile_storage_vec);

			if (pSMesh->profile_storage_dbl3.size() != size) { pSMesh->profile_storage_dbl3.resize(size); }
			profile_storage_vec.copy_to_vector(pSMesh->profile_storage_dbl3);
			pprofile_dbl3 = &pSMesh->profile_storage_dbl3;
		}
		break;

	case MESHDISPLAY_SM_STRAYH:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRAYFIELD)) {

			if (profile_storage_vec.size() != size) { if (!profile_storage_vec.resize(size)) return; }
			dynamic_cast<StrayField*>(pSMesh->pSMod(MODS_STRAYFIELD))->GetStrayFieldCUDA()()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil, profile_storage_vec);

			if (pSMesh->profile_storage_dbl3.size() != size) { pSMesh->profile_storage_dbl3.resize(size); }
			profile_storage_vec.copy_to_vector(pSMesh->profile_storage_dbl3);
			pprofile_dbl3 = &pSMesh->profile_storage_dbl3;
		}
		break;
	}
}

//return average value for currently displayed mesh quantity for named mesh in the given relative rectangle
Any SuperMeshCUDA::GetAverageDisplayedMeshValue(std::string meshName, Rect rel_rect, std::vector<MeshShape> shapes)
{
	if (!pSMesh->contains(meshName) && meshName != pSMesh->superMeshHandle) return Any();

	//get anything displayed on super-mesh
	switch (pSMesh->displayedPhysicalQuantity) {

		////////////////
		//no quantity displayed on the supermesh, so use individual mesh displayed quantities
		////////////////

	case MESHDISPLAY_NONE:
	{
		if (pSMesh->contains(meshName)) return pSMesh->pMesh[meshName]->GetAverageDisplayedMeshValue(rel_rect, shapes);
	}
	break;

	////////////////
	//use a quantity displayed on the supermesh
	////////////////

	case MESHDISPLAY_SM_DEMAG:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_SDEMAG)) {

			return (DBL3)dynamic_cast<SDemag*>(pSMesh->pSMod(MODS_SDEMAG))->GetDemagFieldCUDA()()->average_nonempty(n_fm.dim(), (cuRect)rel_rect);
		}
		break;

	case MESHDISPLAY_SM_OERSTED:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_OERSTED)) {

			return (DBL3)dynamic_cast<Oersted*>(pSMesh->pSMod(MODS_OERSTED))->GetOerstedFieldCUDA()()->average_nonempty(n_e.dim(), (cuRect)rel_rect);
		}
		break;

	case MESHDISPLAY_SM_STRAYH:

		//pdisplay_vec_vec at maximum resolution
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRAYFIELD)) {

			return (DBL3)dynamic_cast<StrayField*>(pSMesh->pSMod(MODS_STRAYFIELD))->GetStrayFieldCUDA()()->average_nonempty(n_fm.dim(), (cuRect)rel_rect);
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

//check disabled_transport_solver flag
bool SuperMeshCUDA::DisabledTransportSolver(void)
{
	return pSMesh->disabled_transport_solver;
}

#endif