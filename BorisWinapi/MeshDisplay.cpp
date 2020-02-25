#include "stdafx.h"
#include "Mesh.h"
#include "SuperMesh.h"

//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

//Return physical quantity currently set to display on screen
PhysQ Mesh::FetchOnScreenPhysicalQuantity(double detail_level, bool getBackground)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) { return pMeshCUDA->FetchOnScreenPhysicalQuantity(detail_level, getBackground); }
#endif

	int physicalQuantity = displayedPhysicalQuantity;
	if (getBackground) physicalQuantity = displayedBackgroundPhysicalQuantity;

	switch (physicalQuantity) {

	case MESHDISPLAY_NONE:
		return PhysQ(meshRect, h, physicalQuantity);
		break;

	case MESHDISPLAY_MAGNETIZATION:
		return PhysQ(&M, physicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_MAGNETIZATION2:
		return PhysQ(&M2, physicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_MAGNETIZATION12:
		return PhysQ(&M, &M2, physicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:
		return PhysQ(&Heff, physicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:
		return PhysQ(&Heff2, physicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:
		return PhysQ(&Heff, &Heff2, physicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_CURRDENSITY:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(&reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetChargeCurrent(), physicalQuantity, (VEC3REP_)vec3rep);
		}
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

			return PhysQ(&reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(0), physicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_JSY:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(&reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(1), physicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_JSZ:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(&reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(2), physicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_TS:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(&reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinTorque(), physicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_TSI:
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(
				&reinterpret_cast<STransport*>(pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorque(reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))),
				physicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_TEMPERATURE:
		return PhysQ(&Temp, physicalQuantity);
		break;

	case MESHDISPLAY_UDISP:
		return PhysQ(&u_disp, physicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_STRAINDIAG:
		return PhysQ(&strain_diag, physicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_STRAINODIAG:
		return PhysQ(&strain_odiag, physicalQuantity, (VEC3REP_)vec3rep);
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

	case MESHDISPLAY_ROUGHNESS:
		if (IsModuleSet(MOD_ROUGHNESS)) {

			return PhysQ(&reinterpret_cast<Roughness*>(pMod(MOD_ROUGHNESS))->GetRoughness(), physicalQuantity);
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
BError Mesh::SaveOnScreenPhysicalQuantity(string fileName, string ovf2_dataType)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) { return pMeshCUDA->SaveOnScreenPhysicalQuantity(fileName, ovf2_dataType); }
#endif

	BError error(__FUNCTION__);

	OVF2 ovf2;

	switch (displayedPhysicalQuantity) {

	case MESHDISPLAY_NONE:
		return error(BERROR_COULDNOTSAVEFILE);
		break;

	case MESHDISPLAY_MAGNETIZATION:
		error = ovf2.Write_OVF2_VEC(fileName, M, ovf2_dataType);
		break;

	case MESHDISPLAY_MAGNETIZATION2:
		error = ovf2.Write_OVF2_VEC(fileName, M2, ovf2_dataType);
		break;

	case MESHDISPLAY_MAGNETIZATION12:
		error = ovf2.Write_OVF2_VEC(fileName, M, ovf2_dataType);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:
		error = ovf2.Write_OVF2_VEC(fileName, Heff, ovf2_dataType);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:
		error = ovf2.Write_OVF2_VEC(fileName, Heff2, ovf2_dataType);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:
		error = ovf2.Write_OVF2_VEC(fileName, Heff, ovf2_dataType);
		break;

	case MESHDISPLAY_CURRDENSITY:
		if (IsModuleSet(MOD_TRANSPORT)) {

			error = ovf2.Write_OVF2_VEC(fileName, reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetChargeCurrent(), ovf2_dataType);
		}
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

			error = ovf2.Write_OVF2_VEC(fileName, reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(0), ovf2_dataType);
		}
		break;

	case MESHDISPLAY_JSY:
		if (IsModuleSet(MOD_TRANSPORT)) {

			error = ovf2.Write_OVF2_VEC(fileName, reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(1), ovf2_dataType);
		}
		break;

	case MESHDISPLAY_JSZ:
		if (IsModuleSet(MOD_TRANSPORT)) {

			error = ovf2.Write_OVF2_VEC(fileName, reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(2), ovf2_dataType);
		}
		break;

	case MESHDISPLAY_TS:
		if (IsModuleSet(MOD_TRANSPORT)) {

			error = ovf2.Write_OVF2_VEC(fileName, reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinTorque(), ovf2_dataType);
		}
		break;

	case MESHDISPLAY_TSI:
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && IsModuleSet(MOD_TRANSPORT)) {

			error = ovf2.Write_OVF2_VEC(fileName, reinterpret_cast<STransport*>(pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorque(reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))), ovf2_dataType);
		}
		break;

	case MESHDISPLAY_TEMPERATURE:
		error = ovf2.Write_OVF2_SCA(fileName, Temp, ovf2_dataType);
		break;

	case MESHDISPLAY_UDISP:
		error = ovf2.Write_OVF2_VEC(fileName, u_disp, ovf2_dataType);
		break;

	case MESHDISPLAY_STRAINDIAG:
		error = ovf2.Write_OVF2_VEC(fileName, strain_diag, ovf2_dataType);
		break;

	case MESHDISPLAY_STRAINODIAG:
		error = ovf2.Write_OVF2_VEC(fileName, strain_odiag, ovf2_dataType);
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

	case MESHDISPLAY_ROUGHNESS:
		if (IsModuleSet(MOD_ROUGHNESS)) {

			error = ovf2.Write_OVF2_SCA(fileName, reinterpret_cast<Roughness*>(pMod(MOD_ROUGHNESS))->GetRoughness(), ovf2_dataType);
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

//Before calling a run of GetDisplayedMeshValue, make sure to call PrepareDisplayedMeshValue : this calculates and stores in displayVEC storage and quantities which don't have memory allocated directly, but require computation and temporary storage.
void Mesh::PrepareDisplayedMeshValue(void)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) { 

		pMeshCUDA->PrepareDisplayedMeshValue(); 
		return;
	}
#endif

	switch (displayedPhysicalQuantity) {

	case MESHDISPLAY_CURRDENSITY:
		if (IsModuleSet(MOD_TRANSPORT)) {

			//charge current calculated internally in the Transport module, ready to be read out when needed
			reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetChargeCurrent();
		}
		break;

	case MESHDISPLAY_JSX:
		if (IsModuleSet(MOD_TRANSPORT)) {

			//spin current calculated internally in the Transport module, ready to be read out when needed
			reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(0);
		}
		break;

	case MESHDISPLAY_JSY:
		if (IsModuleSet(MOD_TRANSPORT)) {

			//spin current calculated internally in the Transport module, ready to be read out when needed
			reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(1);
		}
		break;

	case MESHDISPLAY_JSZ:
		if (IsModuleSet(MOD_TRANSPORT)) {

			//spin current calculated internally in the Transport module, ready to be read out when needed
			reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(2);
		}
		break;

	case MESHDISPLAY_TS:
		if (IsModuleSet(MOD_TRANSPORT)) {

			//spin torque calculated internally in the Transport module, ready to be read out when needed
			reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinTorque();
		}
		break;

	case MESHDISPLAY_TSI:
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && IsModuleSet(MOD_TRANSPORT)) {

			//spin torque calculated internally in the Transport module, ready to be read out when needed
			reinterpret_cast<STransport*>(pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorque(reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT)));
		}
		break;

	default:
		break;
	}
}

//return value of currently displayed mesh quantity at the given absolute position; the value is read directly from the storage VEC, not from the displayed PhysQ.
//Return an Any as the displayed quantity could be either a scalar or a vector.
Any Mesh::GetDisplayedMeshValue(DBL3 abs_pos)
{
	if (!meshRect.contains(abs_pos)) return Any();

#if COMPILECUDA == 1
	if (pMeshCUDA) { return pMeshCUDA->GetDisplayedMeshValue(abs_pos); }
#endif

	DBL3 rel_pos = abs_pos - meshRect.s;

	switch (displayedPhysicalQuantity) {

	case MESHDISPLAY_NONE:
		return (double)0.0;
		break;

	case MESHDISPLAY_MAGNETIZATION:
		if (M.linear_size()) return M[rel_pos];
		break;

	case MESHDISPLAY_MAGNETIZATION2:
		if (M2.linear_size()) return M2[rel_pos];
		break;

	case MESHDISPLAY_MAGNETIZATION12:
		if (M.linear_size()) return M[rel_pos];
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:
		if (Heff.linear_size()) return Heff[rel_pos];
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:
		if (Heff2.linear_size()) return Heff2[rel_pos];
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:
		if (Heff.linear_size()) return Heff[rel_pos];
		break;

	case MESHDISPLAY_CURRDENSITY:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetCalculatedChargeCurrentValue(rel_pos);
		}
		break;

	case MESHDISPLAY_VOLTAGE:
		if (V.linear_size()) return V[rel_pos];
		break;

	case MESHDISPLAY_ELCOND:
		if (elC.linear_size()) return elC[rel_pos];
		break;

	case MESHDISPLAY_SACCUM:
		if (S.linear_size()) return S[rel_pos];
		break;

	case MESHDISPLAY_JSX:
	case MESHDISPLAY_JSY:
	case MESHDISPLAY_JSZ:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetCalculatedSpinCurrentValue(rel_pos);
		}
		break;

	case MESHDISPLAY_TS:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetCalculatedSpinTorqueValue(rel_pos);
		}
		break;

	case MESHDISPLAY_TSI:
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && IsModuleSet(MOD_TRANSPORT)) {

			return reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetCalculatedSpinTorqueValue(rel_pos);
		}
		break;

	case MESHDISPLAY_TEMPERATURE:
		if (Temp.linear_size()) return Temp[rel_pos];
		break;

	case MESHDISPLAY_UDISP:
		if (u_disp.linear_size()) return u_disp[rel_pos];
		break;

	case MESHDISPLAY_STRAINDIAG:
		if (strain_diag.linear_size()) return strain_diag[rel_pos];
		break;

	case MESHDISPLAY_STRAINODIAG:
		if (strain_odiag.linear_size()) return strain_odiag[rel_pos];
		break;

	case MESHDISPLAY_PARAMVAR:
	{
		return get_meshparam_s_scaling_value((PARAM_)displayedParamVar, rel_pos, pSMesh->GetStageTime());
	}
	break;

	case MESHDISPLAY_ROUGHNESS:
		return reinterpret_cast<Roughness*>(pMod(MOD_ROUGHNESS))->GetRoughness()[rel_pos];
		break;

	case MESHDISPLAY_CUSTOM_VEC:
		if (displayVEC_VEC.linear_size()) return displayVEC_VEC[rel_pos];
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		if (displayVEC_SCA.linear_size()) return displayVEC_SCA[rel_pos];
		break;
	}

	return (double)0.0;
}

//return average value for currently displayed mesh quantity in the given relative rectangle
Any Mesh::GetAverageDisplayedMeshValue(Rect rel_rect)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) { return pMeshCUDA->GetAverageDisplayedMeshValue(rel_rect); }
#endif

	switch (displayedPhysicalQuantity) {

	case MESHDISPLAY_NONE:
		return (double)0.0;
		break;

	case MESHDISPLAY_MAGNETIZATION:
		if (M.linear_size()) return M.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_MAGNETIZATION2:
		if (M2.linear_size()) return M2.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_MAGNETIZATION12:
		if (M.linear_size()) return M.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:
		if (Heff.linear_size()) return Heff.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:
		if (Heff2.linear_size()) return Heff2.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:
		if (Heff.linear_size()) return Heff.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_CURRDENSITY:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetAverageChargeCurrent(rel_rect);
		}
		break;

	case MESHDISPLAY_VOLTAGE:
		if (V.linear_size()) return V.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_ELCOND:
		if (elC.linear_size()) return elC.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_SACCUM:
		if (S.linear_size()) return S.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_JSX:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetAverageSpinCurrent(0, rel_rect);
		}
		break;

	case MESHDISPLAY_JSY:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetAverageSpinCurrent(1, rel_rect);
		}
		break;

	case MESHDISPLAY_JSZ:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetAverageSpinCurrent(2, rel_rect);
		}
		break;

	case MESHDISPLAY_TS:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetAverageSpinTorque(rel_rect);
		}
		break;

	case MESHDISPLAY_TSI:
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && IsModuleSet(MOD_TRANSPORT)) {

			//spin torque calculated internally in the Transport module, ready to be read out when needed
			reinterpret_cast<STransport*>(pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorque(reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT)));
			return reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetAverageInterfacialSpinTorque(rel_rect);
		}
		break;

	case MESHDISPLAY_TEMPERATURE:
		if (Temp.linear_size()) return Temp.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_UDISP:
		if (u_disp.linear_size()) return u_disp.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_STRAINDIAG:
		if (strain_diag.linear_size()) return strain_diag.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_STRAINODIAG:
		if (strain_odiag.linear_size()) return strain_odiag.average_nonempty_omp(rel_rect);
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

	case MESHDISPLAY_ROUGHNESS:
		return reinterpret_cast<Roughness*>(pMod(MOD_ROUGHNESS))->GetRoughness().average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_CUSTOM_VEC:
		if (displayVEC_VEC.linear_size()) return displayVEC_VEC.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		if (displayVEC_SCA.linear_size()) return displayVEC_SCA.average_nonempty_omp(rel_rect);
		break;
	}

	return (double)0.0;
}