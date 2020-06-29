#include "stdafx.h"
#include "MeshCUDA.h"
#include "Mesh.h"
#include "SuperMesh.h"
#include "PhysQRep.h"
#include "BorisLib.h"
#include "OVF2_Handlers.h"

#if COMPILECUDA == 1

MeshCUDA::MeshCUDA(Mesh* pMesh) :
	MeshBaseCUDA(pMesh),
	MeshParamsCUDA(dynamic_cast<MeshParams*>(pMesh)),
	MeshDisplayCUDA(),
	n_s(pMesh->n_s), h_s(pMesh->h_s),
	link_stochastic(pMesh->link_stochastic)
{
	this->pMesh = pMesh;

	//make cuda objects in gpu memory from their cpu memory equivalents

	//-----Ferromagnetic properties

	//Magnetization
	if (!M()->set_from_cpuvec(pMesh->M)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!M2()->set_from_cpuvec(pMesh->M2)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);
	
	//effective field - sum total field of all the added modules
	if(!Heff()->set_from_cpuvec(pMesh->Heff)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!Heff2()->set_from_cpuvec(pMesh->Heff2)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//-----Electric conduction properties (Electron charge and spin Transport)

	//electrical potential - on n_e, h_e mesh
	if(!V()->set_from_cpuvec(pMesh->V)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//electrical conductivity - on n_e, h_e mesh
	if(!elC()->set_from_cpuvec(pMesh->elC)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//electrical field - on n_e, h_e mesh
	if (!E()->set_from_cpuvec(pMesh->E)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//electrical current density - on n_e, h_e mesh
	if (!S()->set_from_cpuvec(pMesh->S)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//-----Thermal conduction properties

	//temperature calculated by Heat module
	if(!Temp()->set_from_cpuvec(pMesh->Temp)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!Temp_l()->set_from_cpuvec(pMesh->Temp_l)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//-----Mechanical properties

	//mechanical displacement and strain calculated by MElastic module
	if (!u_disp()->set_from_cpuvec(pMesh->u_disp)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!strain_diag()->set_from_cpuvec(pMesh->strain_diag)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!strain_odiag()->set_from_cpuvec(pMesh->strain_odiag)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);
}

MeshCUDA::~MeshCUDA()
{
	if (Holder_Mesh_Available()) {

		//if this mesh is being deleted, then cuda could be switching off. We need to copy over data in cuVECs back to cpu memory

		//-----Ferromagnetic properties

		//Magnetization
		M()->copy_to_cpuvec(pMesh->M);
		M2()->copy_to_cpuvec(pMesh->M2);

		//effective field - sum total field of all the added modules
		Heff()->copy_to_cpuvec(pMesh->Heff);
		Heff2()->copy_to_cpuvec(pMesh->Heff2);

		//-----Electric conduction properties (Electron charge and spin Transport)

		//copying to cpuvec done in Transport module (if enabled; if not, copying not needed)

		//-----Thermal conduction properties

		//copying to cpuvec done in Heat module (if enabled; if not, copying not needed)

		//-----Mechanical properties

		//copying to cpuvec done in MElastic module (if enabled; if not, copying not needed)
	}
}

//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

PhysQ MeshCUDA::FetchOnScreenPhysicalQuantity(double detail_level, bool getBackground)
{
	int physicalQuantity = pMesh->displayedPhysicalQuantity;
	if (getBackground) physicalQuantity = pMesh->displayedBackgroundPhysicalQuantity;

	switch (physicalQuantity) {

	case MESHDISPLAY_NONE:
		return PhysQ(meshRect, h, physicalQuantity);

	case MESHDISPLAY_MAGNETIZATION:
		
		if (prepare_display(n, meshRect, detail_level, M)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;

	case MESHDISPLAY_MAGNETIZATION2:

		if (prepare_display(n, meshRect, detail_level, M2)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;

	case MESHDISPLAY_MAGNETIZATION12:

		if (prepare_display(n, meshRect, detail_level, M, M2)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_vec, pdisplay2_vec_vc_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;
		
	case MESHDISPLAY_EFFECTIVEFIELD:
		
		if (prepare_display(n, meshRect, detail_level, Heff)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:

		if (prepare_display(n, meshRect, detail_level, Heff2)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:

		if (prepare_display(n, meshRect, detail_level, Heff, Heff2)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vec, pdisplay2_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;
		
	case MESHDISPLAY_CURRDENSITY:

		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n_e, meshRect, detail_level, dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetChargeCurrentCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vc_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		break;

	case MESHDISPLAY_VOLTAGE:

		if (prepare_display(n_e, meshRect, detail_level, V)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_sca, physicalQuantity);
		}
		break;

	case MESHDISPLAY_ELCOND:

		if (prepare_display(n_e, meshRect, detail_level, elC)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_sca, physicalQuantity);
		}
		break;

	case MESHDISPLAY_SACCUM:

		if (prepare_display(n_e, meshRect, detail_level, S)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;

	case MESHDISPLAY_JSX:
		
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n_e, meshRect, detail_level, dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(0))) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		break;

	case MESHDISPLAY_JSY:

		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n_e, meshRect, detail_level, dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(1))) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		break;

	case MESHDISPLAY_JSZ:

		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n_e, meshRect, detail_level, dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(2))) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		break;

	case MESHDISPLAY_TS:

		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n, meshRect, detail_level, dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinTorqueCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		break;

	case MESHDISPLAY_TSI:
		
		if (pMesh->pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n, meshRect, detail_level, 
				dynamic_cast<STransport*>(pMesh->pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorqueCUDA(dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))))) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		break;

	case MESHDISPLAY_TEMPERATURE:
		
		if (prepare_display(n_t, meshRect, detail_level, Temp)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_sca, physicalQuantity);
		}
		break;

	case MESHDISPLAY_UDISP:

		if (prepare_display(n_m, meshRect, detail_level, u_disp)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;

	case MESHDISPLAY_STRAINDIAG:

		if (prepare_display(n_m, meshRect, detail_level, strain_diag)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;

	case MESHDISPLAY_STRAINODIAG:

		if (prepare_display(n_m, meshRect, detail_level, strain_odiag)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;

	case MESHDISPLAY_PARAMVAR:
	{
		void* s_scaling;

		if (pMesh->is_paramvarequation_set((PARAM_)pMesh->displayedParamVar)) {

			//if text equation is set, then we need to calculate the output into a display VEC
			//We could of course calculate this inside the MatP object in its s_scaling VEC, then get it here through reference
			//This is wasteful however as without some nasty book-keeping we could end up with many s_scaling VECs allocated when they are not needed
			//better to just use the single VEC in Mesh intended for display purposes - just means a little bit more work here.

			if (pMesh->is_paramvar_scalar((PARAM_)pMesh->displayedParamVar)) {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				pMesh->displayVEC_SCA.resize(pMesh->get_paramtype_cellsize((PARAM_)pMesh->displayedParamVar), pMesh->meshRect);
				//now calculate it based on the set text equation
				pMesh->calculate_meshparam_s_scaling((PARAM_)pMesh->displayedParamVar, pMesh->displayVEC_SCA, pMesh->pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &pMesh->displayVEC_SCA;
			}
			else {

				//first make sure the display VEC has the right rectangle and cellsize (cellsize appropriate to the type of mesh parameter being displayed - e.g. magnetic, electric, etc..)
				pMesh->displayVEC_VEC.resize(pMesh->get_paramtype_cellsize((PARAM_)pMesh->displayedParamVar), pMesh->meshRect);
				//now calculate it based on the set text equation
				pMesh->calculate_meshparam_s_scaling((PARAM_)pMesh->displayedParamVar, pMesh->displayVEC_VEC, pMesh->pSMesh->GetStageTime());
				//finally set it in s_scaling - come code for setting the PhysQ below
				s_scaling = &pMesh->displayVEC_VEC;
			}
		}
		else {

			//..otherwise we can just get the s_scaling VEC from the MatP object directly.
			s_scaling = pMesh->get_meshparam_s_scaling((PARAM_)pMesh->displayedParamVar);
		}

		if (pMesh->is_paramvar_scalar((PARAM_)pMesh->displayedParamVar)) {

			return PhysQ(reinterpret_cast<VEC<double>*>(s_scaling), physicalQuantity);
		}
		else {

			return PhysQ(reinterpret_cast<VEC<DBL3>*>(s_scaling), physicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
	}
		break;

	case MESHDISPLAY_ROUGHNESS:
		if (pMesh->IsModuleSet(MOD_ROUGHNESS)) {

			return PhysQ(&dynamic_cast<Roughness*>(pMesh->pMod(MOD_ROUGHNESS))->GetRoughness(), physicalQuantity);
		}
		break;

	case MESHDISPLAY_CUSTOM_VEC:
		return PhysQ(&pMesh->displayVEC_VEC, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		return PhysQ(&pMesh->displayVEC_SCA, physicalQuantity);
		break;
	}
	
	return PhysQ(meshRect, h, physicalQuantity);
}

//save the quantity currently displayed on screen in an ovf2 file using the specified format
BError MeshCUDA::SaveOnScreenPhysicalQuantity(string fileName, string ovf2_dataType)
{
	BError error(__FUNCTION__);

	OVF2 ovf2;

	switch (pMesh->displayedPhysicalQuantity) {

	case MESHDISPLAY_NONE:
		return error(BERROR_COULDNOTSAVEFILE);
		break;

	case MESHDISPLAY_MAGNETIZATION:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), M);
		error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vc_vec, ovf2_dataType);
		break;

	case MESHDISPLAY_MAGNETIZATION2:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), M2);
		error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vc_vec, ovf2_dataType);
		break;

	case MESHDISPLAY_MAGNETIZATION12:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), M);
		error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vc_vec, ovf2_dataType);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:

		//pdisplay_vec_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), Heff);
		error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:

		//pdisplay_vec_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), Heff2);
		error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:

		//pdisplay_vec_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), Heff);
		error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		break;

	case MESHDISPLAY_CURRDENSITY:

		//pdisplay_vec_vc_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n_e, meshRect, h_e.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetChargeCurrentCUDA());
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vc_vec, ovf2_dataType);
		}
		break;

	case MESHDISPLAY_VOLTAGE:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_e, meshRect, h_e.mindim(), V);
		error = ovf2.Write_OVF2_SCA(fileName, *pdisplay_vec_vc_sca, ovf2_dataType);
		break;

	case MESHDISPLAY_ELCOND:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_e, meshRect, h_e.mindim(), elC);
		error = ovf2.Write_OVF2_SCA(fileName, *pdisplay_vec_vc_sca, ovf2_dataType);
		break;

	case MESHDISPLAY_SACCUM:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n_e, meshRect, h_e.mindim(), S);
		error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vc_vec, ovf2_dataType);
		break;

	case MESHDISPLAY_JSX:

		//pdisplay_vec_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n_e, meshRect, h_e.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(0));
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		break;

	case MESHDISPLAY_JSY:

		//pdisplay_vec_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n_e, meshRect, h_e.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(1));
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		break;

	case MESHDISPLAY_JSZ:

		//pdisplay_vec_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n_e, meshRect, h_e.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(2));
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		break;

	case MESHDISPLAY_TS:

		//pdisplay_vec_vec at maximumresolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n, meshRect, h.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinTorqueCUDA());
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		break;

	case MESHDISPLAY_TSI:

		//pdisplay_vec_vec at maximumresolution
		if (pMesh->pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n, meshRect, h.mindim(), dynamic_cast<STransport*>(pMesh->pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorqueCUDA(dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))));
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		break;

	case MESHDISPLAY_TEMPERATURE:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_t, meshRect, h_t.mindim(), Temp);
		error = ovf2.Write_OVF2_SCA(fileName, *pdisplay_vec_vc_sca, ovf2_dataType);
		break;

	case MESHDISPLAY_UDISP:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_m, meshRect, h_m.mindim(), u_disp);
		error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vc_vec, ovf2_dataType);
		break;

	case MESHDISPLAY_STRAINDIAG:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_m, meshRect, h_m.mindim(), strain_diag);
		error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vc_vec, ovf2_dataType);
		break;

	case MESHDISPLAY_STRAINODIAG:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_m, meshRect, h_m.mindim(), strain_odiag);
		error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vc_vec, ovf2_dataType);
		break;

	case MESHDISPLAY_ROUGHNESS:
		if (pMesh->IsModuleSet(MOD_ROUGHNESS)) {

			error = ovf2.Write_OVF2_SCA(fileName, dynamic_cast<Roughness*>(pMesh->pMod(MOD_ROUGHNESS))->GetRoughness(), ovf2_dataType);
		}
		break;

	case MESHDISPLAY_CUSTOM_VEC:
		error = ovf2.Write_OVF2_VEC(fileName, pMesh->displayVEC_VEC, ovf2_dataType);
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		error = ovf2.Write_OVF2_SCA(fileName, pMesh->displayVEC_SCA, ovf2_dataType);
		break;
	}

	return error;
}

//Before calling a run of GetDisplayedMeshValue, make sure to call PrepareDisplayedMeshValue : this calculates and stores in displayVEC storage and quantities which don't have memory allocated directly, but require computation and temporary storage.
void MeshCUDA::PrepareDisplayedMeshValue(void)
{
	switch (pMesh->displayedPhysicalQuantity) {

	case MESHDISPLAY_MAGNETIZATION:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), M);
		break;

	case MESHDISPLAY_MAGNETIZATION2:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), M2);
		break;

	case MESHDISPLAY_MAGNETIZATION12:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), M);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:

		//pdisplay_vec_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), Heff);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:

		//pdisplay_vec_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), Heff2);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:

		//pdisplay_vec_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), Heff);
		break;

	case MESHDISPLAY_CURRDENSITY:

		//pdisplay_vec_vc_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n_e, meshRect, h_e.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetChargeCurrentCUDA());
		}
		break;

	case MESHDISPLAY_VOLTAGE:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_e, meshRect, h_e.mindim(), V);
		break;

	case MESHDISPLAY_ELCOND:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_e, meshRect, h_e.mindim(), elC);
		break;

	case MESHDISPLAY_SACCUM:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n_e, meshRect, h_e.mindim(), S);
		break;

	case MESHDISPLAY_JSX:

		//pdisplay_vec_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n_e, meshRect, h_e.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(0));
		}
		break;

	case MESHDISPLAY_JSY:

		//pdisplay_vec_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n_e, meshRect, h_e.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(1));
		}
		break;

	case MESHDISPLAY_JSZ:

		//pdisplay_vec_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n_e, meshRect, h_e.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(2));
		}
		break;

	case MESHDISPLAY_TS:

		//pdisplay_vec_vec at maximumresolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n, meshRect, h.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinTorqueCUDA());
		}
		break;

	case MESHDISPLAY_TSI:

		//pdisplay_vec_vec at maximumresolution
		if (pMesh->pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n, meshRect, h.mindim(), dynamic_cast<STransport*>(pMesh->pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorqueCUDA(dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))));
		}
		break;

	case MESHDISPLAY_TEMPERATURE:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_t, meshRect, h_t.mindim(), Temp);
		break;

	case MESHDISPLAY_UDISP:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n_m, meshRect, h_m.mindim(), u_disp);
		break;

	case MESHDISPLAY_STRAINDIAG:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n_m, meshRect, h_m.mindim(), strain_diag);
		break;

	case MESHDISPLAY_STRAINODIAG:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n_m, meshRect, h_m.mindim(), strain_odiag);
		break;
	}
}

//return value of currently displayed mesh quantity at the given absolute position; the value is read directly from the storage VEC, not from the displayed PhysQ.
//Return an Any as the displayed quantity could be either a scalar or a vector.
Any MeshCUDA::GetDisplayedMeshValue(DBL3 abs_pos)
{
	DBL3 rel_pos = abs_pos - meshRect.s;

	switch (pMesh->displayedPhysicalQuantity) {

	case MESHDISPLAY_NONE:
		return (double)0.0;
		break;

	case MESHDISPLAY_MAGNETIZATION:
		return (*pdisplay_vec_vc_vec)[rel_pos];
		break;

	case MESHDISPLAY_MAGNETIZATION2:
		return (*pdisplay_vec_vc_vec)[rel_pos];
		break;

	case MESHDISPLAY_MAGNETIZATION12:
		return (*pdisplay_vec_vc_vec)[rel_pos];
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:
		return (*pdisplay_vec_vec)[rel_pos];
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:
		return (*pdisplay_vec_vec)[rel_pos];
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:
		return (*pdisplay_vec_vec)[rel_pos];
		break;

	case MESHDISPLAY_CURRDENSITY:
		return (*pdisplay_vec_vc_vec)[rel_pos];
		break;

	case MESHDISPLAY_VOLTAGE:
		return (*pdisplay_vec_vc_sca)[rel_pos];
		break;

	case MESHDISPLAY_ELCOND:
		return (*pdisplay_vec_vc_sca)[rel_pos];
		break;

	case MESHDISPLAY_SACCUM:
		return (*pdisplay_vec_vc_vec)[rel_pos];
		break;

	case MESHDISPLAY_JSX:
	case MESHDISPLAY_JSY:
	case MESHDISPLAY_JSZ:
		return (*pdisplay_vec_vec)[rel_pos];
		break;

	case MESHDISPLAY_TS:
	case MESHDISPLAY_TSI:
		return (*pdisplay_vec_vec)[rel_pos];
		break;

	case MESHDISPLAY_TEMPERATURE:
		return (*pdisplay_vec_vc_sca)[rel_pos];
		break;

	case MESHDISPLAY_UDISP:
		return (*pdisplay_vec_vc_vec)[rel_pos];
		break;

	case MESHDISPLAY_STRAINDIAG:
		return (*pdisplay_vec_vc_vec)[rel_pos];
		break;

	case MESHDISPLAY_STRAINODIAG:
		return (*pdisplay_vec_vc_vec)[rel_pos];
		break;

	case MESHDISPLAY_PARAMVAR:
		return pMesh->get_meshparam_s_scaling_value((PARAM_)pMesh->displayedParamVar, rel_pos, pMesh->pSMesh->GetStageTime());
		break;

	case MESHDISPLAY_ROUGHNESS:
		return dynamic_cast<Roughness*>(pMesh->pMod(MOD_ROUGHNESS))->GetRoughness()[rel_pos];
		break;

	case MESHDISPLAY_CUSTOM_VEC:
		if (pMesh->displayVEC_VEC.linear_size()) return pMesh->displayVEC_VEC[rel_pos];
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		if (pMesh->displayVEC_SCA.linear_size()) return pMesh->displayVEC_SCA[rel_pos];
		break;
	}

	return (double)0.0;
}

//return average value for currently displayed mesh quantity in the given relative rectangle
Any MeshCUDA::GetAverageDisplayedMeshValue(Rect rel_rect)
{
	switch (pMesh->displayedPhysicalQuantity) {

	case MESHDISPLAY_NONE:
		break;

	case MESHDISPLAY_MAGNETIZATION:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), M);
		return pdisplay_vec_vc_vec->average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_MAGNETIZATION2:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), M2);
		return pdisplay_vec_vc_vec->average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_MAGNETIZATION12:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), M);
		return pdisplay_vec_vc_vec->average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD:

		//pdisplay_vec_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), Heff);
		return pdisplay_vec_vec->average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:

		//pdisplay_vec_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), Heff2);
		return pdisplay_vec_vec->average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:

		//pdisplay_vec_vec at maximum resolution
		prepare_display(n, meshRect, h.mindim(), Heff);
		return pdisplay_vec_vec->average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_CURRDENSITY:

		//pdisplay_vec_vc_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n_e, meshRect, h_e.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetChargeCurrentCUDA());
			return pdisplay_vec_vc_vec->average_nonempty_omp(rel_rect);
		}
		break;

	case MESHDISPLAY_VOLTAGE:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_e, meshRect, h_e.mindim(), V);
		return pdisplay_vec_vc_sca->average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_ELCOND:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_e, meshRect, h_e.mindim(), elC);
		return pdisplay_vec_vc_sca->average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_SACCUM:

		//pdisplay_vec_vc_vec at maximum resolution
		prepare_display(n_e, meshRect, h_e.mindim(), S);
		return pdisplay_vec_vc_vec->average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_JSX:

		//pdisplay_vec_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n_e, meshRect, h_e.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(0));
			return pdisplay_vec_vec->average_nonempty_omp(rel_rect);
		}
		break;

	case MESHDISPLAY_JSY:

		//pdisplay_vec_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n_e, meshRect, h_e.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(1));
			return pdisplay_vec_vec->average_nonempty_omp(rel_rect);
		}
		break;

	case MESHDISPLAY_JSZ:

		//pdisplay_vec_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n_e, meshRect, h_e.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(2));
			return pdisplay_vec_vec->average_nonempty_omp(rel_rect);
		}
		break;

	case MESHDISPLAY_TS:

		//pdisplay_vec_vec at maximumresolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n, meshRect, h.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinTorqueCUDA());
			return pdisplay_vec_vec->average_nonempty_omp(rel_rect);
		}
		break;

	case MESHDISPLAY_TSI:

		//pdisplay_vec_vec at maximumresolution
		if (pMesh->pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n, meshRect, h.mindim(), dynamic_cast<STransport*>(pMesh->pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorqueCUDA(dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))));
			return pdisplay_vec_vec->average_nonempty_omp(rel_rect);
		}
		break;

	case MESHDISPLAY_TEMPERATURE:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_t, meshRect, h_t.mindim(), Temp);
		return pdisplay_vec_vc_sca->average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_UDISP:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_m, meshRect, h_m.mindim(), u_disp);
		return pdisplay_vec_vc_vec->average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_STRAINDIAG:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_m, meshRect, h_m.mindim(), strain_diag);
		return pdisplay_vec_vc_vec->average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_STRAINODIAG:

		//pdisplay_vec_vc_sca at maximum resolution
		prepare_display(n_m, meshRect, h_m.mindim(), strain_odiag);
		return pdisplay_vec_vc_vec->average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_ROUGHNESS:
		if (pMesh->IsModuleSet(MOD_ROUGHNESS)) {

			return dynamic_cast<Roughness*>(pMesh->pMod(MOD_ROUGHNESS))->GetRoughness().average_nonempty_omp(rel_rect);
		}
		break;

	case MESHDISPLAY_CUSTOM_VEC:
		return pMesh->displayVEC_VEC.average_nonempty_omp(rel_rect);
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		return pMesh->displayVEC_SCA.average_nonempty_omp(rel_rect);
		break;
	}

	return (double)0.0;
}

//copy aux_vec_sca in GPU memory to displayVEC in CPU memory
void MeshCUDA::copy_aux_vec_sca(VEC<double>& displayVEC)
{
	aux_vec_sca()->copy_to_cpuvec(displayVEC);
}

//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

//magnetization dynamics computation enabled
bool MeshCUDA::MComputation_Enabled(void)
{
	return pMesh->Heff.linear_size();
}

bool MeshCUDA::Magnetism_Enabled(void)
{ 
	return pMesh->M.linear_size();
}

//electrical conduction computation enabled
bool MeshCUDA::EComputation_Enabled(void)
{
	return pMesh->V.linear_size();
}

//thermal conduction computation enabled
bool MeshCUDA::TComputation_Enabled(void)
{
	return pMesh->Temp.linear_size();
}

//mechanical computation enabled
bool MeshCUDA::MechComputation_Enabled(void)
{
	return pMesh->u_disp.linear_size();
}

bool MeshCUDA::GInterface_Enabled(void)
{
	return (DBL2(pMesh->Gmix.get0()).norm() > 0);
}

//----------------------------------- OTHER MESH SHAPE CONTROL

//copy all meshes controlled using change_mesh_shape from cpu to gpu versions
BError MeshCUDA::copy_shapes_from_cpu(void)
{
	//Primary quantities are : M, elC, Temp, u_disp

	BError error(__FUNCTION__);
	
	bool success = true;

	//1a. shape magnetization
	if (M()->size_cpu().dim()) success &= M()->set_from_cpuvec(pMesh->M);

	//1b. shape magnetization for AFM meshes
	if (M2()->size_cpu().dim()) success &= M2()->set_from_cpuvec(pMesh->M2);

	//2. shape electrical conductivity
	if (elC()->size_cpu().dim()) success &= elC()->set_from_cpuvec(pMesh->elC);

	//3. shape temperature
	if (Temp()->size_cpu().dim()) success &= Temp()->set_from_cpuvec(pMesh->Temp);

	//4. shape mechanical properties
	if (u_disp()->size_cpu().dim()) success &= u_disp()->set_from_cpuvec(pMesh->u_disp);

	//if adding any more here also remember to edit change_mesh_shape

	if (!success) error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

//copy all meshes controlled using change_mesh_shape from gpu to cpu versions
BError MeshCUDA::copy_shapes_to_cpu(void)
{
	//Primary quantities are : M, elC, Temp

	BError error(__FUNCTION__);

	bool success = true;

	//1a. shape magnetization
	if (M()->size_cpu().dim()) success &= M()->set_cpuvec(pMesh->M);
	
	//1b. shape magnetization for AFM meshes
	if (M2()->size_cpu().dim()) success &= M2()->set_cpuvec(pMesh->M2);

	//2. shape electrical conductivity
	if (elC()->size_cpu().dim()) success &= elC()->set_cpuvec(pMesh->elC);

	//3. shape temperature
	if (Temp()->size_cpu().dim()) success &= Temp()->set_cpuvec(pMesh->Temp);

	//4. shape mechanical properties
	if (u_disp()->size_cpu().dim()) success &= u_disp()->set_cpuvec(pMesh->u_disp);

	//if adding any more here also remember to edit change_mesh_shape

	if (!success) error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

cu_obj<ManagedDiffEq_CommonCUDA>& MeshCUDA::Get_ManagedDiffEq_CommonCUDA(void)
{ 
	return pMesh->pSMesh->Get_ManagedDiffEq_CommonCUDA(); 
}

#endif