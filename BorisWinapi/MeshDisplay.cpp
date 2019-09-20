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

	case MESHDISPLAY_EFFECTIVEFIELD:
		return PhysQ(&Heff, displayedPhysicalQuantity, (VEC3REP_)vec3rep);
		break;

	case MESHDISPLAY_CURRDENSITY:
		return PhysQ(&Jc, displayedPhysicalQuantity, (VEC3REP_)vec3rep);
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
		void* s_scaling = get_meshparam_s_scaling((PARAM_)displayedParamVar);

		if (is_s_scaling_scalar((PARAM_)displayedParamVar))
			return PhysQ(reinterpret_cast<VEC<double>*>(s_scaling), displayedPhysicalQuantity);
		else
			return PhysQ(reinterpret_cast<VEC<DBL3>*>(s_scaling), displayedPhysicalQuantity, (VEC3REP_)vec3rep);
	}
		break;

	case MESHDISPLAY_ROUGHNESS:
		if (IsModuleSet(MOD_ROUGHNESS)) {

			return PhysQ(&reinterpret_cast<Roughness*>(pMod(MOD_ROUGHNESS))->GetRoughness(), displayedPhysicalQuantity);
		}
		break;
	}

	return PhysQ(meshRect, h, displayedPhysicalQuantity);
}