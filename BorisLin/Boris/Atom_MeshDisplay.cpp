#include "stdafx.h"
#include "Atom_Mesh.h"

#include "SuperMesh.h"

//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

//Return physical quantity currently set to display on screen
PhysQ Atom_Mesh::FetchOnScreenPhysicalQuantity(double detail_level, bool getBackground)
{
#if COMPILECUDA == 1
	if (paMeshCUDA) { return paMeshCUDA->FetchOnScreenPhysicalQuantity(detail_level, getBackground); }
#endif

	int physicalQuantity = displayedPhysicalQuantity;
	if (getBackground) physicalQuantity = displayedBackgroundPhysicalQuantity;

	switch (physicalQuantity) {

	case MESHDISPLAY_NONE:
		return PhysQ(meshRect, h, physicalQuantity);
		break;

	case MESHDISPLAY_MOMENT:
		return PhysQ(&M1, physicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:
		if (Module_Heff_Display == MOD_ALL || Module_Heff_Display == MOD_ERROR) {

			return PhysQ(&Heff1, physicalQuantity, (VEC3REP_)vec3rep);
		}
		else {

			MOD_ Module_Heff = (MOD_)Get_ActualModule_Heff_Display();
			if (IsModuleSet(Module_Heff)) return PhysQ(&pMod(Module_Heff)->Get_Module_Heff(), physicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_ENERGY:
	{
		MOD_ Module_Energy = (MOD_)Get_ActualModule_Energy_Display();
		if (IsModuleSet(Module_Energy)) return PhysQ(&pMod(Module_Energy)->Get_Module_Energy(), physicalQuantity);
	}
		break;

	case MESHDISPLAY_CURRDENSITY:
#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1
		if (IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(&dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetChargeCurrent(), physicalQuantity, (VEC3REP_)vec3rep);
		}
#endif
		break;

	case MESHDISPLAY_VOLTAGE:
		return PhysQ(&V, physicalQuantity);
		break;

	case MESHDISPLAY_ELCOND:
		return PhysQ(&elC, physicalQuantity);
		break;

	case MESHDISPLAY_SACCUM:
		return PhysQ(&S, physicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_JSX:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(&dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(0), physicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_JSY:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(&dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(1), physicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_JSZ:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(&dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(2), physicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_TS:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(&dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetSpinTorque(), physicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_TSI:
#ifdef MODULE_COMPILATION_TRANSPORT
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(
				&dynamic_cast<STransport*>(pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorque(dynamic_cast<Transport*>(pMod(MOD_TRANSPORT))),
				physicalQuantity, (VEC3REP_)vec3rep);
		}
#endif
		break;

	case MESHDISPLAY_TEMPERATURE:
		return PhysQ(&Temp, physicalQuantity);
		break;

	case MESHDISPLAY_PARAMVAR:
	{
		void* s_scaling;

		if (is_paramvarequation_set((PARAM_)displayedParamVar)) {

			//if text equation is set, then we need to calculate the output into a display VEC
			//We could of course calculate this inside the MatP object in its s_scaling VEC, then get it here through reference
			//This is wasteful however as without some nasty book-keeping we could end up with many s_scaling VECs allocated when they are not needed
			//better to just use the single VEC in Mesh intended for display purposes - just means a little bit more work here.

			if (is_paramvar_scalar((PARAM_)displayedParamVar)) {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				displayVEC_SCA.resize(get_paramtype_cellsize((PARAM_)displayedParamVar), meshRect);
				//now calculate it based on the set text equation
				calculate_meshparam_s_scaling((PARAM_)displayedParamVar, displayVEC_SCA, pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &displayVEC_SCA;
			}
			else {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				displayVEC_VEC.resize(get_paramtype_cellsize((PARAM_)displayedParamVar), meshRect);
				//now calculate it based on the set text equation
				calculate_meshparam_s_scaling((PARAM_)displayedParamVar, displayVEC_VEC, pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &displayVEC_VEC;
			}
		}
		else {

			//..otherwise we can just get the s_scaling VEC from the MatP object directly.
			s_scaling = get_meshparam_s_scaling((PARAM_)displayedParamVar);
		}

		if (is_paramvar_scalar((PARAM_)displayedParamVar)) {

			return PhysQ(reinterpret_cast<VEC<double>*>(s_scaling), physicalQuantity);
		}
		else {

			return PhysQ(reinterpret_cast<VEC<DBL3>*>(s_scaling), physicalQuantity, (VEC3REP_)vec3rep);
		}
	}
		break;

	case MESHDISPLAY_CUSTOM_VEC:
		return PhysQ(&displayVEC_VEC, physicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		return PhysQ(&displayVEC_SCA, physicalQuantity);
		break;
	}

	return PhysQ(meshRect, h, physicalQuantity);
}

//save the quantity currently displayed on screen in an ovf2 file using the specified format
BError Atom_Mesh::SaveOnScreenPhysicalQuantity(std::string fileName, std::string ovf2_dataType, MESHDISPLAY_ quantity)
{
#if COMPILECUDA == 1
	if (paMeshCUDA) { return paMeshCUDA->SaveOnScreenPhysicalQuantity(fileName, ovf2_dataType, quantity); }
#endif

	BError error(__FUNCTION__);

	OVF2 ovf2;

	switch ((quantity == MESHDISPLAY_NONE ? displayedPhysicalQuantity : quantity)) {

	case MESHDISPLAY_NONE:
		return error(BERROR_COULDNOTSAVEFILE);
		break;

	case MESHDISPLAY_MOMENT:
		error = ovf2.Write_OVF2_VEC(fileName, M1, ovf2_dataType);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:
		if (Module_Heff_Display == MOD_ALL || Module_Heff_Display == MOD_ERROR) error = ovf2.Write_OVF2_VEC(fileName, Heff1, ovf2_dataType);
		else {

			MOD_ Module_Heff = (MOD_)Get_ActualModule_Heff_Display();
			if (IsModuleSet(Module_Heff)) error = ovf2.Write_OVF2_VEC(fileName, pMod(Module_Heff)->Get_Module_Heff(), ovf2_dataType);
		}
		break;

	case MESHDISPLAY_ENERGY:
	{
		MOD_ Module_Energy = (MOD_)Get_ActualModule_Energy_Display();
		if (IsModuleSet(Module_Energy)) error = ovf2.Write_OVF2_SCA(fileName, pMod(Module_Energy)->Get_Module_Energy(), ovf2_dataType);
	}
		break;

	case MESHDISPLAY_CURRDENSITY:
#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1
		if (IsModuleSet(MOD_TRANSPORT)) {

			error = ovf2.Write_OVF2_VEC(fileName, dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetChargeCurrent(), ovf2_dataType);
		}
#endif
		break;

	case MESHDISPLAY_VOLTAGE:
		error = ovf2.Write_OVF2_SCA(fileName, V, ovf2_dataType);
		break;

	case MESHDISPLAY_ELCOND:
		error = ovf2.Write_OVF2_SCA(fileName, elC, ovf2_dataType);
		break;

	case MESHDISPLAY_SACCUM:
		error = ovf2.Write_OVF2_VEC(fileName, S, ovf2_dataType);
		break;

	case MESHDISPLAY_JSX:
		if (IsModuleSet(MOD_TRANSPORT)) {

			error = ovf2.Write_OVF2_VEC(fileName, dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(0), ovf2_dataType);
		}
		break;

	case MESHDISPLAY_JSY:
		if (IsModuleSet(MOD_TRANSPORT)) {

			error = ovf2.Write_OVF2_VEC(fileName, dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(1), ovf2_dataType);
		}
		break;

	case MESHDISPLAY_JSZ:
		if (IsModuleSet(MOD_TRANSPORT)) {

			error = ovf2.Write_OVF2_VEC(fileName, dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(2), ovf2_dataType);
		}
		break;

	case MESHDISPLAY_TS:
		if (IsModuleSet(MOD_TRANSPORT)) {

			error = ovf2.Write_OVF2_VEC(fileName, dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetSpinTorque(), ovf2_dataType);
		}
		break;

	case MESHDISPLAY_TSI:
#ifdef MODULE_COMPILATION_TRANSPORT
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && IsModuleSet(MOD_TRANSPORT)) {

			error = ovf2.Write_OVF2_VEC(fileName, dynamic_cast<STransport*>(pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorque(dynamic_cast<Transport*>(pMod(MOD_TRANSPORT))), ovf2_dataType);
		}
#endif
		break;

	case MESHDISPLAY_TEMPERATURE:
		error = ovf2.Write_OVF2_SCA(fileName, Temp, ovf2_dataType);
		break;

	case MESHDISPLAY_PARAMVAR:
	{
		void* s_scaling;

		if (is_paramvarequation_set((PARAM_)displayedParamVar)) {

			//if text equation is set, then we need to calculate the output into a display VEC
			//We could of course calculate this inside the MatP object in its s_scaling VEC, then get it here through reference
			//This is wasteful however as without some nasty book-keeping we could end up with many s_scaling VECs allocated when they are not needed
			//better to just use the single VEC in Mesh intended for display purposes - just means a little bit more work here.

			if (is_paramvar_scalar((PARAM_)displayedParamVar)) {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				displayVEC_SCA.resize(get_paramtype_cellsize((PARAM_)displayedParamVar), meshRect);
				//now calculate it based on the set text equation
				calculate_meshparam_s_scaling((PARAM_)displayedParamVar, displayVEC_SCA, pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &displayVEC_SCA;
			}
			else {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				displayVEC_VEC.resize(get_paramtype_cellsize((PARAM_)displayedParamVar), meshRect);
				//now calculate it based on the set text equation
				calculate_meshparam_s_scaling((PARAM_)displayedParamVar, displayVEC_VEC, pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &displayVEC_VEC;
			}
		}
		else {

			//..otherwise we can just get the s_scaling VEC from the MatP object directly.
			s_scaling = get_meshparam_s_scaling((PARAM_)displayedParamVar);
		}

		if (is_paramvar_scalar((PARAM_)displayedParamVar)) {

			error = ovf2.Write_OVF2_SCA(fileName, *reinterpret_cast<VEC<double>*>(s_scaling), ovf2_dataType);
		}
		else {

			error = ovf2.Write_OVF2_VEC(fileName, *reinterpret_cast<VEC<DBL3>*>(s_scaling), ovf2_dataType);
		}
	}
	break;

	case MESHDISPLAY_CUSTOM_VEC:
		error = ovf2.Write_OVF2_VEC(fileName, displayVEC_VEC, ovf2_dataType);
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		error = ovf2.Write_OVF2_SCA(fileName, displayVEC_SCA, ovf2_dataType);
		break;
	}

	return error;
}

//extract profile from focused mesh, from currently display mesh quantity, but reading directly from the quantity
//Displayed	mesh quantity can be scalar or a vector; pass in std::vector pointers, then check for nullptr to determine what type is displayed
//if do_average = true then build average and don't return anything, else return just a single-shot profile. If read_average = true then simply read out the internally stored averaged profile by assigning to pointer.
void Atom_Mesh::GetPhysicalQuantityProfile(
	DBL3 start, DBL3 end, double step, DBL3 stencil, 
	std::vector<DBL3>*& pprofile_dbl3, std::vector<double>*& pprofile_dbl, 
	bool do_average, bool read_average, MESHDISPLAY_ quantity)
{
	auto setup_profile_vec_dbl3 = [&](VEC<DBL3>& vec) -> void
	{
		vec.extract_profile(start, end, step, stencil);
		if (do_average) average_mesh_profile(vec.get_last_profile());
		else pprofile_dbl3 = &vec.get_last_profile();
	};

	auto setup_profile_vecvc_dbl3 = [&](VEC_VC<DBL3>& vec) -> void
	{
		vec.extract_profile(start, end, step, stencil);
		if (do_average) average_mesh_profile(vec.get_last_profile());
		else pprofile_dbl3 = &vec.get_last_profile();
	};

	auto setup_profile_vec_dbl = [&](VEC<double>& vec) -> void
	{
		vec.extract_profile(start, end, step, stencil);
		if (do_average) average_mesh_profile(vec.get_last_profile());
		else pprofile_dbl = &vec.get_last_profile();
	};

	auto setup_profile_vecvc_dbl = [&](VEC_VC<double>& vec) -> void
	{
		vec.extract_profile(start, end, step, stencil);
		if (do_average) average_mesh_profile(vec.get_last_profile());
		else pprofile_dbl = &vec.get_last_profile();
	};

	if (stencil.IsNull()) stencil = h;

	if (read_average) num_profile_averages = 0;

#if COMPILECUDA == 1
	if (paMeshCUDA) { paMeshCUDA->GetPhysicalQuantityProfile(
		start, end, step, stencil, 
		pprofile_dbl3, pprofile_dbl, 
		do_average, read_average, quantity); return; }
#endif

	switch ((quantity == MESHDISPLAY_NONE ? displayedPhysicalQuantity : quantity)) {

	default:
	case MESHDISPLAY_NONE:
		break;

	case MESHDISPLAY_MOMENT:
		if (read_average) { pprofile_dbl3 = &profile_storage_dbl3; return; }
		setup_profile_vecvc_dbl3(M1);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:
		if (read_average) { pprofile_dbl3 = &profile_storage_dbl3; return; }
		if (Module_Heff_Display == MOD_ALL || Module_Heff_Display == MOD_ERROR) {
			
			setup_profile_vec_dbl3(Heff1);
		}
		else {

			MOD_ Module_Heff = (MOD_)Get_ActualModule_Heff_Display();
			if (IsModuleSet(Module_Heff)) {

				setup_profile_vec_dbl3(pMod(Module_Heff)->Get_Module_Heff());
			}
		}
		break;

	case MESHDISPLAY_ENERGY:
	{
		if (read_average) { pprofile_dbl = &profile_storage_dbl; return; }
		MOD_ Module_Energy = (MOD_)Get_ActualModule_Energy_Display();
		if (IsModuleSet(Module_Energy)) {

			setup_profile_vec_dbl(pMod(Module_Energy)->Get_Module_Energy());
		}
	}
	break;

	case MESHDISPLAY_CURRDENSITY:
#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1
		if (read_average) { pprofile_dbl3 = &profile_storage_dbl3; return; }
		if (IsModuleSet(MOD_TRANSPORT)) {

			setup_profile_vecvc_dbl3(dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetChargeCurrent());
		}
#endif
		break;

	case MESHDISPLAY_VOLTAGE:
		if (read_average) { pprofile_dbl = &profile_storage_dbl; return; }
		setup_profile_vecvc_dbl(V);
		break;

	case MESHDISPLAY_ELCOND:
		if (read_average) { pprofile_dbl = &profile_storage_dbl; return; }
		setup_profile_vecvc_dbl(elC);
		break;

	case MESHDISPLAY_SACCUM:
		if (read_average) { pprofile_dbl3 = &profile_storage_dbl3; return; }
		setup_profile_vecvc_dbl3(S);
		break;

	case MESHDISPLAY_JSX:
		if (read_average) { pprofile_dbl3 = &profile_storage_dbl3; return; }
		if (IsModuleSet(MOD_TRANSPORT)) {

			setup_profile_vec_dbl3(dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(0));
		}
		break;

	case MESHDISPLAY_JSY:
		if (read_average) { pprofile_dbl3 = &profile_storage_dbl3; return; }
		if (IsModuleSet(MOD_TRANSPORT)) {

			setup_profile_vec_dbl3(dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(1));
		}
		break;

	case MESHDISPLAY_JSZ:
		if (read_average) { pprofile_dbl3 = &profile_storage_dbl3; return; }
		if (IsModuleSet(MOD_TRANSPORT)) {

			setup_profile_vec_dbl3(dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(2));
		}
		break;

	case MESHDISPLAY_TS:
		if (read_average) { pprofile_dbl3 = &profile_storage_dbl3; return; }
		if (IsModuleSet(MOD_TRANSPORT)) {

			setup_profile_vec_dbl3(dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetSpinTorque());
		}
		break;

	case MESHDISPLAY_TSI:
#ifdef MODULE_COMPILATION_TRANSPORT
		if (read_average) { pprofile_dbl3 = &profile_storage_dbl3; return; }
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && IsModuleSet(MOD_TRANSPORT)) {

			setup_profile_vec_dbl3(dynamic_cast<STransport*>(pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorque(dynamic_cast<Transport*>(pMod(MOD_TRANSPORT))));
		}
#endif
		break;

	case MESHDISPLAY_TEMPERATURE:
		if (read_average) { pprofile_dbl = &profile_storage_dbl; return; }
		setup_profile_vecvc_dbl(Temp);
		break;

	case MESHDISPLAY_PARAMVAR:
	{
		void* s_scaling;

		if (is_paramvarequation_set((PARAM_)displayedParamVar)) {

			//if text equation is set, then we need to calculate the output into a display VEC
			//We could of course calculate this inside the MatP object in its s_scaling VEC, then get it here through reference
			//This is wasteful however as without some nasty book-keeping we could end up with many s_scaling VECs allocated when they are not needed
			//better to just use the single VEC in Mesh intended for display purposes - just means a little bit more work here.

			if (is_paramvar_scalar((PARAM_)displayedParamVar)) {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				displayVEC_SCA.resize(get_paramtype_cellsize((PARAM_)displayedParamVar), meshRect);
				//now calculate it based on the set text equation
				calculate_meshparam_s_scaling((PARAM_)displayedParamVar, displayVEC_SCA, pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &displayVEC_SCA;
			}
			else {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				displayVEC_VEC.resize(get_paramtype_cellsize((PARAM_)displayedParamVar), meshRect);
				//now calculate it based on the set text equation
				calculate_meshparam_s_scaling((PARAM_)displayedParamVar, displayVEC_VEC, pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &displayVEC_VEC;
			}
		}
		else {

			//..otherwise we can just get the s_scaling VEC from the MatP object directly.
			s_scaling = get_meshparam_s_scaling((PARAM_)displayedParamVar);
		}

		if (is_paramvar_scalar((PARAM_)displayedParamVar)) {

			profile_storage_dbl = reinterpret_cast<VEC<double>*>(s_scaling)->extract_profile(start, end, step, stencil);
			pprofile_dbl = &profile_storage_dbl;
		}
		else {

			profile_storage_dbl3 = reinterpret_cast<VEC<DBL3>*>(s_scaling)->extract_profile(start, end, step, stencil);
			pprofile_dbl3 = &profile_storage_dbl3;
		}
	}
	break;

	case MESHDISPLAY_CUSTOM_VEC:
		profile_storage_dbl3 = displayVEC_VEC.extract_profile(start, end, step, stencil);
		pprofile_dbl3 = &profile_storage_dbl3;
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		profile_storage_dbl = displayVEC_SCA.extract_profile(start, end, step, stencil);
		pprofile_dbl = &profile_storage_dbl;
		break;
	}
}

//return average value for currently displayed mesh quantity in the given relative rectangle
Any Atom_Mesh::GetAverageDisplayedMeshValue(Rect rel_rect, std::vector<MeshShape> shapes, MESHDISPLAY_ quantity)
{
#if COMPILECUDA == 1
	if (paMeshCUDA) { return paMeshCUDA->GetAverageDisplayedMeshValue(rel_rect, quantity); }
#endif

	switch ((quantity == MESHDISPLAY_NONE ? displayedPhysicalQuantity : quantity)) {
		
	default:
	case MESHDISPLAY_NONE:
		break;

	case MESHDISPLAY_MOMENT:
		if (M1.linear_size()) {

			if (!shapes.size()) return M1.average_nonempty_omp(rel_rect);
			else return M1.shape_getaverage(shapes);
		}
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:
		if (Heff1.linear_size()) {

			if (Module_Heff_Display == MOD_ALL || Module_Heff_Display == MOD_ERROR) {

				if (!shapes.size()) return Heff1.average_nonempty_omp(rel_rect);
				else return Heff1.shape_getaverage(shapes);
			}
			else {

				MOD_ Module_Heff = (MOD_)Get_ActualModule_Heff_Display();
				if (IsModuleSet(Module_Heff)) {

					if (!shapes.size()) return pMod(Module_Heff)->Get_Module_Heff().average_nonempty_omp(rel_rect);
					else return pMod(Module_Heff)->Get_Module_Heff().shape_getaverage(shapes);
				}
			}
		}
		break;

	case MESHDISPLAY_ENERGY:
	{
		MOD_ Module_Energy = (MOD_)Get_ActualModule_Energy_Display();
		if (IsModuleSet(Module_Energy)) {

			if (!shapes.size()) return pMod(Module_Energy)->Get_Module_Energy().average_nonempty_omp(rel_rect);
			else return pMod(Module_Energy)->Get_Module_Energy().shape_getaverage(shapes);
		}
	}
	break;

	case MESHDISPLAY_CURRDENSITY:
#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1
		if (IsModuleSet(MOD_TRANSPORT)) {

			return dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetAverageChargeCurrent(rel_rect, shapes);
		}
#endif
		break;

	case MESHDISPLAY_VOLTAGE:
		if (V.linear_size()) {

			if (!shapes.size()) return V.average_nonempty_omp(rel_rect);
			else return V.shape_getaverage(shapes);
		}
		break;

	case MESHDISPLAY_ELCOND:
		if (elC.linear_size()) {

			if (!shapes.size()) return elC.average_nonempty_omp(rel_rect);
			else return elC.shape_getaverage(shapes);
		}
		break;

	case MESHDISPLAY_SACCUM:
		if (S.linear_size()) {

			if (!shapes.size()) return S.average_nonempty_omp(rel_rect);
			else return S.shape_getaverage(shapes);
		}
		break;

	case MESHDISPLAY_JSX:
#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1
		if (IsModuleSet(MOD_TRANSPORT)) {

			return dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetAverageSpinCurrent(0, rel_rect, shapes);
		}
#endif
		break;

	case MESHDISPLAY_JSY:
#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1
		if (IsModuleSet(MOD_TRANSPORT)) {

			return dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetAverageSpinCurrent(1, rel_rect, shapes);
		}
#endif
		break;

	case MESHDISPLAY_JSZ:
#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1
		if (IsModuleSet(MOD_TRANSPORT)) {

			return dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetAverageSpinCurrent(2, rel_rect, shapes);
		}
#endif
		break;

	case MESHDISPLAY_TS:
#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1
		if (IsModuleSet(MOD_TRANSPORT)) {

			return dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetAverageSpinTorque(rel_rect, shapes);
		}
#endif
		break;

	case MESHDISPLAY_TSI:
#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && IsModuleSet(MOD_TRANSPORT)) {

			//spin torque calculated internally in the Transport module, ready to be read out when needed
			dynamic_cast<STransport*>(pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorque(dynamic_cast<Transport*>(pMod(MOD_TRANSPORT)));
			return dynamic_cast<Atom_Transport*>(pMod(MOD_TRANSPORT))->GetAverageInterfacialSpinTorque(rel_rect, shapes);
		}
#endif
		break;

	case MESHDISPLAY_TEMPERATURE:
		if (Temp.linear_size()) {

			if (!shapes.size()) return Temp.average_nonempty_omp(rel_rect);
			else return Temp.shape_getaverage(shapes);
		}
		break;

	case MESHDISPLAY_PARAMVAR:
	{
		void* s_scaling;

		if (is_paramvarequation_set((PARAM_)displayedParamVar)) {

			//if text equation is set, then we need to calculate the output into a display VEC
			//We could of course calculate this inside the MatP object in its s_scaling VEC, then get it here through reference
			//This is wasteful however as without some nasty book-keeping we could end up with many s_scaling VECs allocated when they are not needed
			//better to just use the single VEC in Mesh intended for display purposes - just means a little bit more work here.

			if (is_paramvar_scalar((PARAM_)displayedParamVar)) {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				displayVEC_SCA.resize(get_paramtype_cellsize((PARAM_)displayedParamVar), meshRect);
				//now calculate it based on the set text equation
				calculate_meshparam_s_scaling((PARAM_)displayedParamVar, displayVEC_SCA, pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &displayVEC_SCA;
			}
			else {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				displayVEC_VEC.resize(get_paramtype_cellsize((PARAM_)displayedParamVar), meshRect);
				//now calculate it based on the set text equation
				calculate_meshparam_s_scaling((PARAM_)displayedParamVar, displayVEC_VEC, pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &displayVEC_VEC;
			}
		}
		else {

			//..otherwise we can just get the s_scaling VEC from the MatP object directly.
			s_scaling = get_meshparam_s_scaling((PARAM_)displayedParamVar);
		}

		if (is_paramvar_scalar((PARAM_)displayedParamVar)) {

			return reinterpret_cast<VEC<double>*>(s_scaling)->average_nonempty_omp(rel_rect);
		}
		else {

			return reinterpret_cast<VEC<DBL3>*>(s_scaling)->average_nonempty_omp(rel_rect);
		}
	}
	break;

	case MESHDISPLAY_CUSTOM_VEC:
		if (displayVEC_VEC.linear_size()) {

			if (!shapes.size()) return displayVEC_VEC.average_nonempty_omp(rel_rect);
			else return displayVEC_VEC.shape_getaverage(shapes);
		}
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		if (displayVEC_SCA.linear_size()) {

			if (!shapes.size()) return displayVEC_SCA.average_nonempty_omp(rel_rect);
			else return displayVEC_SCA.shape_getaverage(shapes);
		}
		break;
	}

	return (double)0.0;
}