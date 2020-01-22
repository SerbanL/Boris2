#include "stdafx.h"
#include "Mesh.h"
#include "SuperMesh.h"

//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

//Return physical quantity currently set to display on screen
PhysQ Mesh::FetchOnScreenPhysicalQuantity(double detail_level)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) { return pMeshCUDA->FetchOnScreenPhysicalQuantity(detail_level); }
#endif

	switch (displayedPhysicalQuantity) {

	case MESHDISPLAY_NONE:
		return PhysQ(meshRect, h, displayedPhysicalQuantity);
		break;

	case MESHDISPLAY_MAGNETIZATION:
		return PhysQ(&M, displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_MAGNETIZATION2:
		return PhysQ(&M2, displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_MAGNETIZATION12:
		return PhysQ(&M, &M2, displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:
		return PhysQ(&Heff, displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:
		return PhysQ(&Heff2, displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:
		return PhysQ(&Heff, &Heff2, displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_CURRDENSITY:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(&reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetChargeCurrent(), displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_VOLTAGE:
		return PhysQ(&V, displayedPhysicalQuantity);
		break;

	case MESHDISPLAY_ELCOND:
		return PhysQ(&elC, displayedPhysicalQuantity);
		break;

	case MESHDISPLAY_SACCUM:
		return PhysQ(&S, displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_JSX:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(&reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(0), displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_JSY:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(&reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(1), displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_JSZ:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(&reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinCurrent(2), displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_TS:
		if (IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(&reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))->GetSpinTorque(), displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_TSI:
		if (pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && IsModuleSet(MOD_TRANSPORT)) {

			return PhysQ(
				&reinterpret_cast<STransport*>(pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorque(reinterpret_cast<Transport*>(pMod(MOD_TRANSPORT))),
				displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		}
		break;

	case MESHDISPLAY_TEMPERATURE:
		return PhysQ(&Temp, displayedPhysicalQuantity);
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

			return PhysQ(reinterpret_cast<VEC<double>*>(s_scaling), displayedPhysicalQuantity);
		}
		else {

			return PhysQ(reinterpret_cast<VEC<DBL3>*>(s_scaling), displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		}
	}
		break;

	case MESHDISPLAY_ROUGHNESS:
		if (IsModuleSet(MOD_ROUGHNESS)) {

			return PhysQ(&reinterpret_cast<Roughness*>(pMod(MOD_ROUGHNESS))->GetRoughness(), displayedPhysicalQuantity);
		}
		break;

	case MESHDISPLAY_CUSTOM_VEC:
		return PhysQ(&displayVEC_VEC, displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		return PhysQ(&displayVEC_SCA, displayedPhysicalQuantity);
		break;
	}

	return PhysQ(meshRect, h, displayedPhysicalQuantity);
}