#include "stdafx.h"
#include "Atom_MeshCUDA.h"
#include "Atom_Mesh.h"
#include "SuperMesh.h"
#include "PhysQRep.h"
#include "BorisLib.h"
#include "OVF2_Handlers.h"

#if COMPILECUDA == 1

Atom_MeshCUDA::Atom_MeshCUDA(Atom_Mesh* paMesh) :
	MeshBaseCUDA(paMesh),
	Atom_MeshParamsCUDA(dynamic_cast<Atom_MeshParams*>(paMesh)),
	MeshDisplayCUDA(),
	n_dm(paMesh->n_dm), h_dm(paMesh->h_dm)
{
	this->paMesh = paMesh;

	//make cuda objects in gpu memory from their cpu memory equivalents

	//-----Magnetic properties

	//Moment
	if (!M1()->set_from_cpuvec(paMesh->M1)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//effective field - sum total field of all the added modules
	if (!Heff1()->set_from_cpuvec(paMesh->Heff1)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//-----Electric conduction properties (Electron charge and spin Transport)

	//electrical potential - on n_e, h_e mesh
	if (!V()->set_from_cpuvec(paMesh->V)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//electrical conductivity - on n_e, h_e mesh
	if (!elC()->set_from_cpuvec(paMesh->elC)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//electrical field - on n_e, h_e mesh
	if (!E()->set_from_cpuvec(paMesh->E)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//electrical current density - on n_e, h_e mesh
	if (!S()->set_from_cpuvec(paMesh->S)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//-----Thermal conduction properties

	//temperature calculated by Heat module
	if (!Temp()->set_from_cpuvec(paMesh->Temp)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!Temp_l()->set_from_cpuvec(paMesh->Temp_l)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//-----Mechanical properties

	//mechanical displacement and strain calculated by MElastic module
	if (!u_disp()->set_from_cpuvec(paMesh->u_disp)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!strain_diag()->set_from_cpuvec(paMesh->strain_diag)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!strain_odiag()->set_from_cpuvec(paMesh->strain_odiag)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);
}

Atom_MeshCUDA::~Atom_MeshCUDA()
{
	if (Holder_Mesh_Available()) {

		//if this mesh is being deleted, then cuda could be switching off. We need to copy over data in cuVECs back to cpu memory

		//-----Magnetic properties

		//Moment
		M1()->copy_to_cpuvec(paMesh->M1);

		//effective field - sum total field of all the added modules
		Heff1()->copy_to_cpuvec(paMesh->Heff1);

		//-----Electric conduction properties (Electron charge and spin Transport)

		//Transport module not currently used for atomistic meshes : TO DO

		//-----Thermal conduction properties

		Temp()->copy_to_cpuvec(paMesh->Temp);
		Temp_l()->copy_to_cpuvec(paMesh->Temp_l);

		//-----Mechanical properties

		//MElastic module not currently used for atomistic meshes : TO DO
	}
}

//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

PhysQ Atom_MeshCUDA::FetchOnScreenPhysicalQuantity(double detail_level, bool getBackground)
{
	int physicalQuantity = paMesh->displayedPhysicalQuantity;
	if (getBackground) physicalQuantity = paMesh->displayedBackgroundPhysicalQuantity;

	switch (physicalQuantity) {

	case MESHDISPLAY_NONE:
		return PhysQ(meshRect, h, physicalQuantity);

	case MESHDISPLAY_MOMENT:

		if (prepare_display(n, meshRect, detail_level, M1)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_vec, physicalQuantity, (VEC3REP_)paMesh->vec3rep);
		}
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:
		if ((MOD_)paMesh->Get_Module_Heff_Display() == MOD_ALL || (MOD_)paMesh->Get_Module_Heff_Display() == MOD_ERROR) {

			if (prepare_display(n, meshRect, detail_level, Heff1)) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, physicalQuantity, (VEC3REP_)paMesh->vec3rep);
			}
		}
		else {

			MOD_ Module_Heff = (MOD_)paMesh->Get_ActualModule_Heff_Display();
			if (paMesh->IsModuleSet(Module_Heff)) {

				if (prepare_display(n, meshRect, detail_level, paMesh->pMod(Module_Heff)->Get_Module_HeffCUDA())) {

					//return PhysQ made from the cpu version of coarse mesh display.
					return PhysQ(pdisplay_vec_vec, physicalQuantity, (VEC3REP_)paMesh->vec3rep);
				}
			}
		}
		break;

	case MESHDISPLAY_ENERGY:
	{
		MOD_ Module_Energy = (MOD_)paMesh->Get_ActualModule_Energy_Display();
		if (paMesh->IsModuleSet(Module_Energy)) {

			if (prepare_display(n, meshRect, detail_level, paMesh->pMod(Module_Energy)->Get_Module_EnergyCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_sca, physicalQuantity);
			}
		}
	}
		break;

	case MESHDISPLAY_TEMPERATURE:

		if (prepare_display(n_t, meshRect, detail_level, Temp)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_sca, physicalQuantity);
		}
		break;

	case MESHDISPLAY_PARAMVAR:
	{
		void* s_scaling;

		if (paMesh->is_paramvarequation_set((PARAM_)paMesh->displayedParamVar)) {

			//if text equation is set, then we need to calculate the output into a display VEC
			//We could of course calculate this inside the MatP object in its s_scaling VEC, then get it here through reference
			//This is wasteful however as without some nasty book-keeping we could end up with many s_scaling VECs allocated when they are not needed
			//better to just use the single VEC in Mesh intended for display purposes - just means a little bit more work here.

			if (paMesh->is_paramvar_scalar((PARAM_)paMesh->displayedParamVar)) {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				paMesh->displayVEC_SCA.resize(paMesh->get_paramtype_cellsize((PARAM_)paMesh->displayedParamVar), paMesh->meshRect);
				//now calculate it based on the set text equation
				paMesh->calculate_meshparam_s_scaling((PARAM_)paMesh->displayedParamVar, paMesh->displayVEC_SCA, paMesh->pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &paMesh->displayVEC_SCA;
			}
			else {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				paMesh->displayVEC_VEC.resize(paMesh->get_paramtype_cellsize((PARAM_)paMesh->displayedParamVar), paMesh->meshRect);
				//now calculate it based on the set text equation
				paMesh->calculate_meshparam_s_scaling((PARAM_)paMesh->displayedParamVar, paMesh->displayVEC_VEC, paMesh->pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &paMesh->displayVEC_VEC;
			}
		}
		else {

			//..otherwise we can just get the s_scaling VEC from the MatP object directly.
			s_scaling = paMesh->get_meshparam_s_scaling((PARAM_)paMesh->displayedParamVar);
		}

		if (paMesh->is_paramvar_scalar((PARAM_)paMesh->displayedParamVar)) {

			return PhysQ(reinterpret_cast<VEC<double>*>(s_scaling), physicalQuantity);
		}
		else {

			return PhysQ(reinterpret_cast<VEC<DBL3>*>(s_scaling), physicalQuantity, (VEC3REP_)paMesh->vec3rep);
		}
	}
	break;

	case MESHDISPLAY_CUSTOM_VEC:
		return PhysQ(&paMesh->displayVEC_VEC, physicalQuantity, (VEC3REP_)paMesh->vec3rep);
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		return PhysQ(&paMesh->displayVEC_SCA, physicalQuantity);
		break;
	}

	return PhysQ(meshRect, h, physicalQuantity);
}

//save the quantity currently displayed on screen in an ovf2 file using the specified format
BError Atom_MeshCUDA::SaveOnScreenPhysicalQuantity(std::string fileName, std::string ovf2_dataType)
{
	BError error(__FUNCTION__);

	OVF2 ovf2;

	switch (paMesh->displayedPhysicalQuantity) {

	case MESHDISPLAY_NONE:
		return error(BERROR_COULDNOTSAVEFILE);
		break;

	case MESHDISPLAY_MOMENT:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), M1);
		error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vc_vec, ovf2_dataType);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:
		if ((MOD_)paMesh->Get_Module_Heff_Display() == MOD_ALL || (MOD_)paMesh->Get_Module_Heff_Display() == MOD_ERROR) {

			prepare_display(n, meshRect, h.mindim(), Heff1);
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		else {

			MOD_ Module_Heff = (MOD_)paMesh->Get_ActualModule_Heff_Display();
			if (paMesh->IsModuleSet(Module_Heff)) {

				prepare_display(n, meshRect, h.mindim(), paMesh->pMod(Module_Heff)->Get_Module_HeffCUDA());
				error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
			}
		}
		break;

	case MESHDISPLAY_ENERGY:
	{
		MOD_ Module_Energy = (MOD_)paMesh->Get_ActualModule_Energy_Display();
		if (paMesh->IsModuleSet(Module_Energy)) {

			prepare_display(n, meshRect, h.mindim(), paMesh->pMod(Module_Energy)->Get_Module_EnergyCUDA());
			error = ovf2.Write_OVF2_SCA(fileName, *pdisplay_vec_sca, ovf2_dataType);
		}
	}
		break;

	case MESHDISPLAY_TEMPERATURE:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_t, meshRect, h_t.mindim(), Temp);
		error = ovf2.Write_OVF2_SCA(fileName, *pdisplay_vec_vc_sca, ovf2_dataType);
		break;

	case MESHDISPLAY_CUSTOM_VEC:
		error = ovf2.Write_OVF2_VEC(fileName, paMesh->displayVEC_VEC, ovf2_dataType);
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		error = ovf2.Write_OVF2_SCA(fileName, paMesh->displayVEC_SCA, ovf2_dataType);
		break;
	}

	return error;
}

//extract profile from focused mesh, from currently display mesh quantity, but reading directly from the quantity
//Displayed	mesh quantity can be scalar or a vector; pass in std::vector pointers, then check for nullptr to determine what type is displayed
//if do_average = true then build average and don't return anything, else return just a single-shot profile. If read_average = true then simply read out the internally stored averaged profile by assigning to pointer.
void Atom_MeshCUDA::GetPhysicalQuantityProfile(DBL3 start, DBL3 end, double step, DBL3 stencil, std::vector<DBL3>*& pprofile_dbl3, std::vector<double>*& pprofile_dbl, bool do_average, bool read_average)
{
	size_t size = round((end - start).norm() / step) + 1;

	auto read_profile_vec = [&](void) -> void
	{
		if (paMesh->profile_storage_dbl3.size() != profile_storage_vec.size()) { paMesh->profile_storage_dbl3.resize(profile_storage_vec.size()); }
		profile_storage_vec.copy_to_vector(paMesh->profile_storage_dbl3);
		pprofile_dbl3 = &paMesh->profile_storage_dbl3;
	};

	auto read_profile_sca = [&](void) -> void
	{
		if (paMesh->profile_storage_dbl.size() != profile_storage_sca.size()) { paMesh->profile_storage_dbl.resize(profile_storage_sca.size()); }
		profile_storage_sca.copy_to_vector(paMesh->profile_storage_dbl);
		pprofile_dbl = &paMesh->profile_storage_dbl;
	};

	auto setup_profile_cuvecvc_vec = [&](cu_obj<cuVEC_VC<cuReal3>>& vec) -> void
	{
		if (profile_storage_vec.size() != size) { if (!profile_storage_vec.resize(size)) return; }
		if (do_average) {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil);
			average_mesh_profile(vec, paMesh->num_profile_averages);
		}
		else {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil, profile_storage_vec);
			if (paMesh->profile_storage_dbl3.size() != size) { paMesh->profile_storage_dbl3.resize(size); }
			profile_storage_vec.copy_to_vector(paMesh->profile_storage_dbl3);
			pprofile_dbl3 = &paMesh->profile_storage_dbl3;
		}
	};

	auto setup_profile_cuvec_vec = [&](cu_obj<cuVEC<cuReal3>>& vec) -> void
	{
		if (profile_storage_vec.size() != size) { if (!profile_storage_vec.resize(size)) return; }
		if (do_average) {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil);
			average_mesh_profile(vec, paMesh->num_profile_averages);
		}
		else {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil, profile_storage_vec);
			if (paMesh->profile_storage_dbl3.size() != size) { paMesh->profile_storage_dbl3.resize(size); }
			profile_storage_vec.copy_to_vector(paMesh->profile_storage_dbl3);
			pprofile_dbl3 = &paMesh->profile_storage_dbl3;
		}
	};

	auto setup_profile_cuvecvc_sca = [&](cu_obj<cuVEC_VC<cuBReal>>& vec) -> void
	{
		if (profile_storage_sca.size() != size) { if (!profile_storage_sca.resize(size)) return; }
		if (do_average) {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil);
			average_mesh_profile(vec, paMesh->num_profile_averages);
		}
		else {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil, profile_storage_sca);
			if (paMesh->profile_storage_dbl.size() != size) { paMesh->profile_storage_dbl.resize(size); }
			profile_storage_sca.copy_to_vector(paMesh->profile_storage_dbl);
			pprofile_dbl = &paMesh->profile_storage_dbl;
		}
	};

	auto setup_profile_cuvec_sca = [&](cu_obj<cuVEC<cuBReal>>& vec) -> void
	{
		if (profile_storage_sca.size() != size) { if (!profile_storage_sca.resize(size)) return; }
		if (do_average) {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil);
			average_mesh_profile(vec, paMesh->num_profile_averages);
		}
		else {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil, profile_storage_sca);
			if (paMesh->profile_storage_dbl.size() != size) { paMesh->profile_storage_dbl.resize(size); }
			profile_storage_sca.copy_to_vector(paMesh->profile_storage_dbl);
			pprofile_dbl = &paMesh->profile_storage_dbl;
		}
	};

	if (read_average) paMesh->num_profile_averages = 0;

	switch (paMesh->displayedPhysicalQuantity) {

	default:
	case MESHDISPLAY_NONE:
		break;

	case MESHDISPLAY_MOMENT:
		if (read_average) { read_profile_vec(); return; }
		setup_profile_cuvecvc_vec(M1);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:
		if (read_average) { read_profile_vec(); return; }
		if ((MOD_)paMesh->Get_Module_Heff_Display() == MOD_ALL || (MOD_)paMesh->Get_Module_Heff_Display() == MOD_ERROR) {

			setup_profile_cuvec_vec(Heff1);
		}
		else {

			MOD_ Module_Heff = (MOD_)paMesh->Get_ActualModule_Heff_Display();
			if (paMesh->IsModuleSet(Module_Heff)) {

				setup_profile_cuvec_vec(paMesh->pMod(Module_Heff)->Get_Module_HeffCUDA());
			}
		}
		break;

	case MESHDISPLAY_ENERGY:
	{
		if (read_average) { read_profile_sca(); return; }
		MOD_ Module_Energy = (MOD_)paMesh->Get_ActualModule_Energy_Display();
		if (paMesh->IsModuleSet(Module_Energy)) {

			setup_profile_cuvec_sca(paMesh->pMod(Module_Energy)->Get_Module_EnergyCUDA());
		}
	}
	break;

	case MESHDISPLAY_TEMPERATURE:
		if (read_average) { read_profile_sca(); return; }
		setup_profile_cuvecvc_sca(Temp);
		break;

	case MESHDISPLAY_PARAMVAR:
	{
		void* s_scaling;

		if (paMesh->is_paramvarequation_set((PARAM_)paMesh->displayedParamVar)) {

			//if text equation is set, then we need to calculate the output into a display VEC
			//We could of course calculate this inside the MatP object in its s_scaling VEC, then get it here through reference
			//This is wasteful however as without some nasty book-keeping we could end up with many s_scaling VECs allocated when they are not needed
			//better to just use the single VEC in Mesh intended for display purposes - just means a little bit more work here.

			if (paMesh->is_paramvar_scalar((PARAM_)paMesh->displayedParamVar)) {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				paMesh->displayVEC_SCA.resize(paMesh->get_paramtype_cellsize((PARAM_)paMesh->displayedParamVar), paMesh->meshRect);
				//now calculate it based on the set text equation
				paMesh->calculate_meshparam_s_scaling((PARAM_)paMesh->displayedParamVar, paMesh->displayVEC_SCA, paMesh->pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &paMesh->displayVEC_SCA;
			}
			else {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				paMesh->displayVEC_VEC.resize(paMesh->get_paramtype_cellsize((PARAM_)paMesh->displayedParamVar), paMesh->meshRect);
				//now calculate it based on the set text equation
				paMesh->calculate_meshparam_s_scaling((PARAM_)paMesh->displayedParamVar, paMesh->displayVEC_VEC, paMesh->pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &paMesh->displayVEC_VEC;
			}
		}
		else {

			//..otherwise we can just get the s_scaling VEC from the MatP object directly.
			s_scaling = paMesh->get_meshparam_s_scaling((PARAM_)paMesh->displayedParamVar);
		}

		if (paMesh->is_paramvar_scalar((PARAM_)paMesh->displayedParamVar)) {

			paMesh->profile_storage_dbl = reinterpret_cast<VEC<double>*>(s_scaling)->extract_profile(start, end, step, stencil);
			pprofile_dbl = &paMesh->profile_storage_dbl;
		}
		else {

			paMesh->profile_storage_dbl3 = reinterpret_cast<VEC<DBL3>*>(s_scaling)->extract_profile(start, end, step, stencil);
			pprofile_dbl3 = &paMesh->profile_storage_dbl3;
		}
	}
	break;

	case MESHDISPLAY_CUSTOM_VEC:
		paMesh->profile_storage_dbl3 = paMesh->displayVEC_VEC.extract_profile(start, end, step, stencil);
		pprofile_dbl3 = &paMesh->profile_storage_dbl3;
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		paMesh->profile_storage_dbl = paMesh->displayVEC_SCA.extract_profile(start, end, step, stencil);
		pprofile_dbl = &paMesh->profile_storage_dbl;
		break;
	}
}

//return average value for currently displayed mesh quantity in the given relative rectangle
Any Atom_MeshCUDA::GetAverageDisplayedMeshValue(Rect rel_rect)
{
	switch (paMesh->displayedPhysicalQuantity) {

	default:
	case MESHDISPLAY_NONE:
		break;

	case MESHDISPLAY_MOMENT:
		return (DBL3)M1()->average_nonempty(paMesh->M1.linear_size(), (cuRect)rel_rect);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:
		return (DBL3)Heff1()->average_nonempty(paMesh->Heff1.linear_size(), (cuRect)rel_rect);
		break;

	case MESHDISPLAY_ENERGY:
	{
		MOD_ Module_Energy = (MOD_)paMesh->Get_ActualModule_Energy_Display();
		if (paMesh->IsModuleSet(Module_Energy)) {

			return (double)paMesh->pMod(Module_Energy)->Get_Module_EnergyCUDA()()->average_nonempty(paMesh->pMod(Module_Energy)->Get_Module_Energy().linear_size(), (cuRect)rel_rect);
		}
	}
	break;

	case MESHDISPLAY_TEMPERATURE:
		return (double)Temp()->average_nonempty(paMesh->Temp.linear_size(), (cuRect)rel_rect);
		break;

	case MESHDISPLAY_PARAMVAR:
	{
		void* s_scaling;

		if (paMesh->is_paramvarequation_set((PARAM_)paMesh->displayedParamVar)) {

			//if text equation is set, then we need to calculate the output into a display VEC
			//We could of course calculate this inside the MatP object in its s_scaling VEC, then get it here through reference
			//This is wasteful however as without some nasty book-keeping we could end up with many s_scaling VECs allocated when they are not needed
			//better to just use the single VEC in Mesh intended for display purposes - just means a little bit more work here.

			if (paMesh->is_paramvar_scalar((PARAM_)paMesh->displayedParamVar)) {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				paMesh->displayVEC_SCA.resize(paMesh->get_paramtype_cellsize((PARAM_)paMesh->displayedParamVar), paMesh->meshRect);
				//now calculate it based on the set text equation
				paMesh->calculate_meshparam_s_scaling((PARAM_)paMesh->displayedParamVar, paMesh->displayVEC_SCA, paMesh->pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &paMesh->displayVEC_SCA;
			}
			else {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				paMesh->displayVEC_VEC.resize(paMesh->get_paramtype_cellsize((PARAM_)paMesh->displayedParamVar), paMesh->meshRect);
				//now calculate it based on the set text equation
				paMesh->calculate_meshparam_s_scaling((PARAM_)paMesh->displayedParamVar, paMesh->displayVEC_VEC, paMesh->pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &paMesh->displayVEC_VEC;
			}
		}
		else {

			//..otherwise we can just get the s_scaling VEC from the MatP object directly.
			s_scaling = paMesh->get_meshparam_s_scaling((PARAM_)paMesh->displayedParamVar);
		}

		if (paMesh->is_paramvar_scalar((PARAM_)paMesh->displayedParamVar)) {

			return reinterpret_cast<VEC<double>*>(s_scaling)->average_nonempty_omp(rel_rect);
		}
		else {

			return reinterpret_cast<VEC<DBL3>*>(s_scaling)->average_nonempty_omp(rel_rect);
		}
	}
	break;

	case MESHDISPLAY_CUSTOM_VEC:
		return paMesh->displayVEC_VEC.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		return paMesh->displayVEC_SCA.average_nonempty_omp(rel_rect);
		break;
	}

	return (double)0.0;
}

//copy aux_vec_sca in GPU memory to displayVEC in CPU memory
void Atom_MeshCUDA::copy_aux_vec_sca(VEC<double>& displayVEC)
{
	aux_vec_sca()->copy_to_cpuvec(displayVEC);
}

//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

//magnetization dynamics computation enabled
bool Atom_MeshCUDA::MComputation_Enabled(void)
{
	return paMesh->Heff1.linear_size();
}

bool Atom_MeshCUDA::Magnetism_Enabled(void)
{
	return paMesh->M1.linear_size();
}

//electrical conduction computation enabled
bool Atom_MeshCUDA::EComputation_Enabled(void)
{
	return paMesh->V.linear_size();
}

//thermal conduction computation enabled
bool Atom_MeshCUDA::TComputation_Enabled(void)
{
	return paMesh->Temp.linear_size();
}

//mechanical computation enabled
bool Atom_MeshCUDA::MechComputation_Enabled(void)
{
	return paMesh->u_disp.linear_size();
}

bool Atom_MeshCUDA::GInterface_Enabled(void)
{
	//TO DO
	return false;
	//return (DBL2(paMesh->Gmix.get0()).norm() > 0);
}

//----------------------------------- OTHER MESH SHAPE CONTROL

//copy all meshes controlled using change_mesh_shape from cpu to gpu versions
BError Atom_MeshCUDA::copy_shapes_from_cpu(void)
{
	//Primary quantities are : M, elC, Temp, u_disp

	BError error(__FUNCTION__);

	bool success = true;

	//1. shape moments
	if (M1()->size_cpu().dim()) success &= M1()->set_from_cpuvec(paMesh->M1);

	//2. shape electrical conductivity
	if (elC()->size_cpu().dim()) success &= elC()->set_from_cpuvec(paMesh->elC);

	//3. shape temperature
	if (Temp()->size_cpu().dim()) success &= Temp()->set_from_cpuvec(paMesh->Temp);

	//4. shape mechanical properties
	if (u_disp()->size_cpu().dim()) success &= u_disp()->set_from_cpuvec(paMesh->u_disp);

	//if adding any more here also remember to edit change_mesh_shape

	if (!success) error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

//copy all meshes controlled using change_mesh_shape from gpu to cpu versions
BError Atom_MeshCUDA::copy_shapes_to_cpu(void)
{
	//Primary quantities are : M, elC, Temp, u_disp

	BError error(__FUNCTION__);

	bool success = true;

	//1. shape moments
	if (M1()->size_cpu().dim()) success &= M1()->set_cpuvec(paMesh->M1);

	//2. shape electrical conductivity
	if (elC()->size_cpu().dim()) success &= elC()->set_cpuvec(paMesh->elC);

	//3. shape temperature
	if (Temp()->size_cpu().dim()) success &= Temp()->set_cpuvec(paMesh->Temp);

	//4. shape mechanical properties
	if (u_disp()->size_cpu().dim()) success &= u_disp()->set_cpuvec(paMesh->u_disp);

	//if adding any more here also remember to edit change_mesh_shape

	if (!success) error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

cu_obj<ManagedAtom_DiffEq_CommonCUDA>& Atom_MeshCUDA::Get_ManagedAtom_DiffEq_CommonCUDA(void)
{
	return paMesh->pSMesh->Get_ManagedAtom_DiffEq_CommonCUDA();
}

std::vector<DBL4>& Atom_MeshCUDA::get_tensorial_anisotropy(void)
{
	return paMesh->Kt;
}

#endif