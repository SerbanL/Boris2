#include "stdafx.h"
#include "MeshCUDA.h"
#include "Mesh.h"
#include "SuperMesh.h"
#include "PhysQRep.h"
#include "BorisLib.h"
#include "OVF2_Handlers.h"

#if COMPILECUDA == 1

#include "ManagedDiffEqFMCUDA.h"

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

		elC()->copy_to_cpuvec(pMesh->elC);
		V()->copy_to_cpuvec(pMesh->V);
		E()->copy_to_cpuvec(pMesh->E);
		S()->copy_to_cpuvec(pMesh->S);

		//-----Thermal conduction properties

		Temp()->copy_to_cpuvec(pMesh->Temp);
		Temp_l()->copy_to_cpuvec(pMesh->Temp_l);

		//-----Mechanical properties

		u_disp()->copy_to_cpuvec(pMesh->u_disp);
		strain_diag()->copy_to_cpuvec(pMesh->strain_diag);
		strain_odiag()->copy_to_cpuvec(pMesh->strain_odiag);
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
		if ((MOD_)pMesh->Get_Module_Heff_Display() == MOD_ALL || (MOD_)pMesh->Get_Module_Heff_Display() == MOD_ERROR) {

			if (prepare_display(n, meshRect, detail_level, Heff)) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		else {

			MOD_ Module_Heff = (MOD_)pMesh->Get_ActualModule_Heff_Display();
			if (pMesh->IsModuleSet(Module_Heff)) {

				if (prepare_display(n, meshRect, detail_level, pMesh->pMod(Module_Heff)->Get_Module_HeffCUDA())) {

					//return PhysQ made from the cpu version of coarse mesh display.
					return PhysQ(pdisplay_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
				}
			}
		}
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:
		if ((MOD_)pMesh->Get_Module_Heff_Display() == MOD_ALL || (MOD_)pMesh->Get_Module_Heff_Display() == MOD_ERROR) {

			if (prepare_display(n, meshRect, detail_level, Heff2)) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		else {

			MOD_ Module_Heff = (MOD_)pMesh->Get_ActualModule_Heff_Display();
			if (pMesh->IsModuleSet(Module_Heff)) {

				if (prepare_display(n, meshRect, detail_level, pMesh->pMod(Module_Heff)->Get_Module_Heff2CUDA())) {

					//return PhysQ made from the cpu version of coarse mesh display.
					return PhysQ(pdisplay_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
				}
			}
		}
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:
		if ((MOD_)pMesh->Get_Module_Heff_Display() == MOD_ALL || (MOD_)pMesh->Get_Module_Heff_Display() == MOD_ERROR) {

			if (prepare_display(n, meshRect, detail_level, Heff, Heff2)) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, pdisplay2_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		else {

			MOD_ Module_Heff = (MOD_)pMesh->Get_ActualModule_Heff_Display();
			if (pMesh->IsModuleSet(Module_Heff)) {

				if (prepare_display(n, meshRect, detail_level, pMesh->pMod(Module_Heff)->Get_Module_HeffCUDA(), pMesh->pMod(Module_Heff)->Get_Module_Heff2CUDA())) {

					//return PhysQ made from the cpu version of coarse mesh display.
					return PhysQ(pdisplay_vec_vec, pdisplay2_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
				}
			}
		}
		break;

	case MESHDISPLAY_ENERGY:
	{
		MOD_ Module_Energy = (MOD_)pMesh->Get_ActualModule_Energy_Display();
		if (pMesh->IsModuleSet(Module_Energy)) {

			if (prepare_display(n, meshRect, detail_level, pMesh->pMod(Module_Energy)->Get_Module_EnergyCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_sca, physicalQuantity);
			}
		}
	}
		break;

	case MESHDISPLAY_ENERGY2:
	{
		MOD_ Module_Energy = (MOD_)pMesh->Get_ActualModule_Energy_Display();
		if (pMesh->IsModuleSet(Module_Energy)) {

			if (prepare_display(n, meshRect, detail_level, pMesh->pMod(Module_Energy)->Get_Module_Energy2CUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_sca, physicalQuantity);
			}
		}
	}
		break;

	case MESHDISPLAY_CURRDENSITY:

		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n_e, meshRect, detail_level, dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetChargeCurrentCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vc_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		else if (pMesh->IsModuleSet(MOD_TMR)) {

			if (prepare_display(n_e, meshRect, detail_level, dynamic_cast<TMR*>(pMesh->pMod(MOD_TMR))->GetChargeCurrentCUDA())) {

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
#ifdef MODULE_COMPILATION_TRANSPORT
		if (pMesh->pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n, meshRect, detail_level, 
				dynamic_cast<STransport*>(pMesh->pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorqueCUDA(dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))))) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, physicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
#endif
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
BError MeshCUDA::SaveOnScreenPhysicalQuantity(std::string fileName, std::string ovf2_dataType, MESHDISPLAY_ quantity)
{
	BError error(__FUNCTION__);

	OVF2 ovf2;

	switch ((quantity == MESHDISPLAY_NONE ? pMesh->displayedPhysicalQuantity : quantity)) {

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
		if ((MOD_)pMesh->Get_Module_Heff_Display() == MOD_ALL || (MOD_)pMesh->Get_Module_Heff_Display() == MOD_ERROR) {

			prepare_display(n, meshRect, h.mindim(), Heff);
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		else {

			MOD_ Module_Heff = (MOD_)pMesh->Get_ActualModule_Heff_Display();
			if (pMesh->IsModuleSet(Module_Heff)) {

				prepare_display(n, meshRect, h.mindim(), pMesh->pMod(Module_Heff)->Get_Module_HeffCUDA());
				error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
			}
		}
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:
		if ((MOD_)pMesh->Get_Module_Heff_Display() == MOD_ALL || (MOD_)pMesh->Get_Module_Heff_Display() == MOD_ERROR) {

			prepare_display(n, meshRect, h.mindim(), Heff2);
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		else {

			MOD_ Module_Heff = (MOD_)pMesh->Get_ActualModule_Heff_Display();
			if (pMesh->IsModuleSet(Module_Heff)) {

				prepare_display(n, meshRect, h.mindim(), pMesh->pMod(Module_Heff)->Get_Module_Heff2CUDA());
				error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
			}
		}
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:
		if ((MOD_)pMesh->Get_Module_Heff_Display() == MOD_ALL || (MOD_)pMesh->Get_Module_Heff_Display() == MOD_ERROR) {

			prepare_display(n, meshRect, h.mindim(), Heff);
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
		else {

			MOD_ Module_Heff = (MOD_)pMesh->Get_ActualModule_Heff_Display();
			if (pMesh->IsModuleSet(Module_Heff)) {

				prepare_display(n, meshRect, h.mindim(), pMesh->pMod(Module_Heff)->Get_Module_HeffCUDA());
				error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
			}
		}
		break;

	case MESHDISPLAY_ENERGY:
	{
		MOD_ Module_Energy = (MOD_)pMesh->Get_ActualModule_Energy_Display();
		if (pMesh->IsModuleSet(Module_Energy)) {

			prepare_display(n, meshRect, h.mindim(), pMesh->pMod(Module_Energy)->Get_Module_EnergyCUDA());
			error = ovf2.Write_OVF2_SCA(fileName, *pdisplay_vec_sca, ovf2_dataType);
		}
	}
		break;

	case MESHDISPLAY_ENERGY2:
	{
		MOD_ Module_Energy = (MOD_)pMesh->Get_ActualModule_Energy_Display();
		if (pMesh->IsModuleSet(Module_Energy)) {

			prepare_display(n, meshRect, h.mindim(), pMesh->pMod(Module_Energy)->Get_Module_Energy2CUDA());
			error = ovf2.Write_OVF2_SCA(fileName, *pdisplay_vec_sca, ovf2_dataType);
		}
	}
		break;

	case MESHDISPLAY_CURRDENSITY:

		//pdisplay_vec_vc_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n_e, meshRect, h_e.mindim(), dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetChargeCurrentCUDA());
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vc_vec, ovf2_dataType);
		}
		else if (pMesh->IsModuleSet(MOD_TMR)) {

			prepare_display(n_e, meshRect, h_e.mindim(), dynamic_cast<TMR*>(pMesh->pMod(MOD_TMR))->GetChargeCurrentCUDA());
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
#ifdef MODULE_COMPILATION_TRANSPORT
		//pdisplay_vec_vec at maximumresolution
		if (pMesh->pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && pMesh->IsModuleSet(MOD_TRANSPORT)) {

			prepare_display(n, meshRect, h.mindim(), dynamic_cast<STransport*>(pMesh->pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorqueCUDA(dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))));
			error = ovf2.Write_OVF2_VEC(fileName, *pdisplay_vec_vec, ovf2_dataType);
		}
#endif
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

//extract profile from focused mesh, from currently display mesh quantity, but reading directly from the quantity
//Displayed	mesh quantity can be scalar or a vector; pass in std::vector pointers, then check for nullptr to determine what type is displayed
//if do_average = true then build average and don't return anything, else return just a single-shot profile. If read_average = true then simply read out the internally stored averaged profile by assigning to pointer.
void MeshCUDA::GetPhysicalQuantityProfile(
	DBL3 start, DBL3 end, double step, DBL3 stencil, 
	std::vector<DBL3>*& pprofile_dbl3, std::vector<double>*& pprofile_dbl, 
	bool do_average, bool read_average, MESHDISPLAY_ quantity)
{
	size_t size = round((end - start).norm() / step) + 1;

	auto read_profile_vec = [&](void) -> void
	{
		if (pMesh->profile_storage_dbl3.size() != profile_storage_vec.size()) { pMesh->profile_storage_dbl3.resize(profile_storage_vec.size()); }
		profile_storage_vec.copy_to_vector(pMesh->profile_storage_dbl3);
		pprofile_dbl3 = &pMesh->profile_storage_dbl3;
	};

	auto read_profile_sca = [&](void) -> void
	{
		if (pMesh->profile_storage_dbl.size() != profile_storage_sca.size()) { pMesh->profile_storage_dbl.resize(profile_storage_sca.size()); }
		profile_storage_sca.copy_to_vector(pMesh->profile_storage_dbl);
		pprofile_dbl = &pMesh->profile_storage_dbl;
	};

	auto setup_profile_cuvecvc_vec = [&](cu_obj<cuVEC_VC<cuReal3>>& vec) -> void
	{
		if (profile_storage_vec.size() != size) { if (!profile_storage_vec.resize(size)) return; }
		if (do_average) {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil);
			average_mesh_profile(vec, pMesh->num_profile_averages);
		}
		else {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil, profile_storage_vec);
			if (pMesh->profile_storage_dbl3.size() != size) { pMesh->profile_storage_dbl3.resize(size); }
			profile_storage_vec.copy_to_vector(pMesh->profile_storage_dbl3);
			pprofile_dbl3 = &pMesh->profile_storage_dbl3;
		}
	};

	auto setup_profile_cuvec_vec = [&](cu_obj<cuVEC<cuReal3>>& vec) -> void
	{
		if (profile_storage_vec.size() != size) { if (!profile_storage_vec.resize(size)) return; }
		if (do_average) {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil);
			average_mesh_profile(vec, pMesh->num_profile_averages);
		}
		else {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil, profile_storage_vec);
			if (pMesh->profile_storage_dbl3.size() != size) { pMesh->profile_storage_dbl3.resize(size); }
			profile_storage_vec.copy_to_vector(pMesh->profile_storage_dbl3);
			pprofile_dbl3 = &pMesh->profile_storage_dbl3;
		}
	};

	auto setup_profile_cuvecvc_sca = [&](cu_obj<cuVEC_VC<cuBReal>>& vec) -> void
	{
		if (profile_storage_sca.size() != size) { if (!profile_storage_sca.resize(size)) return; }
		if (do_average) {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil);
			average_mesh_profile(vec, pMesh->num_profile_averages);
		}
		else {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil, profile_storage_sca);
			if (pMesh->profile_storage_dbl.size() != size) { pMesh->profile_storage_dbl.resize(size); }
			profile_storage_sca.copy_to_vector(pMesh->profile_storage_dbl);
			pprofile_dbl = &pMesh->profile_storage_dbl;
		}
	};

	auto setup_profile_cuvec_sca = [&](cu_obj<cuVEC<cuBReal>>& vec) -> void
	{
		if (profile_storage_sca.size() != size) { if (!profile_storage_sca.resize(size)) return; }
		if (do_average) {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil);
			average_mesh_profile(vec, pMesh->num_profile_averages);
		}
		else {

			vec()->extract_profile((cuReal3)start, (cuReal3)end, (cuBReal)step, (cuReal3)stencil, profile_storage_sca);
			if (pMesh->profile_storage_dbl.size() != size) { pMesh->profile_storage_dbl.resize(size); }
			profile_storage_sca.copy_to_vector(pMesh->profile_storage_dbl);
			pprofile_dbl = &pMesh->profile_storage_dbl;
		}
	};

	if (read_average) pMesh->num_profile_averages = 0;

	switch ((quantity == MESHDISPLAY_NONE ? pMesh->displayedPhysicalQuantity : quantity)) {

	default:
	case MESHDISPLAY_NONE:
		break;

	case MESHDISPLAY_MAGNETIZATION12:
	case MESHDISPLAY_MAGNETIZATION:
		if (read_average) { read_profile_vec(); return; }
		setup_profile_cuvecvc_vec(M);
		break;
		
	case MESHDISPLAY_MAGNETIZATION2:
		if (read_average) { read_profile_vec(); return; }
		setup_profile_cuvecvc_vec(M2);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:
	case MESHDISPLAY_EFFECTIVEFIELD:
		if (read_average) { read_profile_vec(); return; }
		if ((MOD_)pMesh->Get_Module_Heff_Display() == MOD_ALL || (MOD_)pMesh->Get_Module_Heff_Display() == MOD_ERROR) {

			setup_profile_cuvec_vec(Heff);
		}
		else {

			MOD_ Module_Heff = (MOD_)pMesh->Get_ActualModule_Heff_Display();
			if (pMesh->IsModuleSet(Module_Heff)) {

				setup_profile_cuvec_vec(pMesh->pMod(Module_Heff)->Get_Module_HeffCUDA());
			}
		}
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:
		if (read_average) { read_profile_vec(); return; }
		if ((MOD_)pMesh->Get_Module_Heff_Display() == MOD_ALL || (MOD_)pMesh->Get_Module_Heff_Display() == MOD_ERROR) {

			setup_profile_cuvec_vec(Heff2);
		}
		else {

			MOD_ Module_Heff = (MOD_)pMesh->Get_ActualModule_Heff_Display();
			if (pMesh->IsModuleSet(Module_Heff)) {

				setup_profile_cuvec_vec(pMesh->pMod(Module_Heff)->Get_Module_Heff2CUDA());
			}
		}
		break;

	case MESHDISPLAY_ENERGY:
	{
		if (read_average) { read_profile_sca(); return; }
		MOD_ Module_Energy = (MOD_)pMesh->Get_ActualModule_Energy_Display();
		if (pMesh->IsModuleSet(Module_Energy)) {

			setup_profile_cuvec_sca(pMesh->pMod(Module_Energy)->Get_Module_EnergyCUDA());
		}
	}
	break;

	case MESHDISPLAY_ENERGY2:
	{
		if (read_average) { read_profile_sca(); return; }
		MOD_ Module_Energy = (MOD_)pMesh->Get_ActualModule_Energy_Display();
		if (pMesh->IsModuleSet(Module_Energy)) {

			setup_profile_cuvec_sca(pMesh->pMod(Module_Energy)->Get_Module_Energy2CUDA());
		}
	}
	break;

	case MESHDISPLAY_CURRDENSITY:
		if (read_average) { read_profile_vec(); return; }
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			setup_profile_cuvecvc_vec(dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetChargeCurrentCUDA());
		}
		else if (pMesh->IsModuleSet(MOD_TMR)) {

			setup_profile_cuvecvc_vec(dynamic_cast<TMR*>(pMesh->pMod(MOD_TMR))->GetChargeCurrentCUDA());
		}
		break;

	case MESHDISPLAY_VOLTAGE:
		if (read_average) { read_profile_sca(); return; }
		setup_profile_cuvecvc_sca(V);
		break;

	case MESHDISPLAY_ELCOND:
		if (read_average) { read_profile_sca(); return; }
		setup_profile_cuvecvc_sca(elC);
		break;

	case MESHDISPLAY_SACCUM:
		if (read_average) { read_profile_vec(); return; }
		setup_profile_cuvecvc_vec(S);
		break;

	case MESHDISPLAY_JSX:
		if (read_average) { read_profile_vec(); return; }
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			setup_profile_cuvec_vec(dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(0));
		}
		break;

	case MESHDISPLAY_JSY:
		if (read_average) { read_profile_vec(); return; }
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			setup_profile_cuvec_vec(dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(1));
		}
		break;

	case MESHDISPLAY_JSZ:
		if (read_average) { read_profile_vec(); return; }
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			setup_profile_cuvec_vec(dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(2));
		}
		break;

	case MESHDISPLAY_TS:

		if (read_average) { read_profile_vec(); return; }
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			setup_profile_cuvec_vec(dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinTorqueCUDA());
		}
		break;

	case MESHDISPLAY_TSI:
#ifdef MODULE_COMPILATION_TRANSPORT
		if (read_average) { read_profile_vec(); return; }
		if (pMesh->pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && pMesh->IsModuleSet(MOD_TRANSPORT)) {

			setup_profile_cuvec_vec(dynamic_cast<STransport*>(pMesh->pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorqueCUDA(dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))));
		}
#endif
		break;

	case MESHDISPLAY_TEMPERATURE:
		if (read_average) { read_profile_sca(); return; }
		setup_profile_cuvecvc_sca(Temp);
		break;

	case MESHDISPLAY_UDISP:
		if (read_average) { read_profile_vec(); return; }
		setup_profile_cuvecvc_vec(u_disp);
		break;

	case MESHDISPLAY_STRAINDIAG:
		if (read_average) { read_profile_vec(); return; }
		setup_profile_cuvecvc_vec(strain_diag);
		break;

	case MESHDISPLAY_STRAINODIAG:
		if (read_average) { read_profile_vec(); return; }
		setup_profile_cuvecvc_vec(strain_odiag);
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

			pMesh->profile_storage_dbl = reinterpret_cast<VEC<double>*>(s_scaling)->extract_profile(start, end, step, stencil);
			pprofile_dbl = &pMesh->profile_storage_dbl;
		}
		else {

			pMesh->profile_storage_dbl3 = reinterpret_cast<VEC<DBL3>*>(s_scaling)->extract_profile(start, end, step, stencil);
			pprofile_dbl3 = &pMesh->profile_storage_dbl3;
		}
	}
	break;

	case MESHDISPLAY_ROUGHNESS:
		if (pMesh->IsModuleSet(MOD_ROUGHNESS)) {

			pMesh->profile_storage_dbl = dynamic_cast<Roughness*>(pMesh->pMod(MOD_ROUGHNESS))->GetRoughness().extract_profile(start, end, step, stencil);
			pprofile_dbl = &pMesh->profile_storage_dbl;
		}
		break;

	case MESHDISPLAY_CUSTOM_VEC:
		pMesh->profile_storage_dbl3 = pMesh->displayVEC_VEC.extract_profile(start, end, step, stencil);
		pprofile_dbl3 = &pMesh->profile_storage_dbl3;
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		pMesh->profile_storage_dbl = pMesh->displayVEC_SCA.extract_profile(start, end, step, stencil);
		pprofile_dbl = &pMesh->profile_storage_dbl;
		break;
	}
}

//return average value for currently displayed mesh quantity in the given relative rectangle
Any MeshCUDA::GetAverageDisplayedMeshValue(Rect rel_rect, MESHDISPLAY_ quantity)
{
	switch ((quantity == MESHDISPLAY_NONE ? pMesh->displayedPhysicalQuantity : quantity)) {

	default:
	case MESHDISPLAY_NONE:
		break;

	case MESHDISPLAY_MAGNETIZATION12:
	case MESHDISPLAY_MAGNETIZATION:
		return (DBL3)M()->average_nonempty(pMesh->M.linear_size(), (cuRect)rel_rect);
		break;

	case MESHDISPLAY_MAGNETIZATION2:
		return (DBL3)M2()->average_nonempty(pMesh->M2.linear_size(), (cuRect)rel_rect);
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:
	case MESHDISPLAY_EFFECTIVEFIELD:
		if ((MOD_)pMesh->Get_Module_Heff_Display() == MOD_ALL || (MOD_)pMesh->Get_Module_Heff_Display() == MOD_ERROR) {

			return (DBL3)Heff()->average_nonempty(pMesh->Heff.linear_size(), (cuRect)rel_rect);
		}
		else {

			MOD_ Module_Heff = (MOD_)pMesh->Get_ActualModule_Heff_Display();
			if (pMesh->IsModuleSet(Module_Heff)) {

				return (DBL3)pMesh->pMod(Module_Heff)->Get_Module_HeffCUDA()()->average_nonempty(pMesh->pMod(Module_Heff)->Get_Module_Heff().linear_size(), (cuRect)rel_rect);
			}
		}
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:
		if ((MOD_)pMesh->Get_Module_Heff_Display() == MOD_ALL || (MOD_)pMesh->Get_Module_Heff_Display() == MOD_ERROR) {

			return (DBL3)Heff2()->average_nonempty(pMesh->Heff2.linear_size(), (cuRect)rel_rect);
		}
		else {

			MOD_ Module_Heff = (MOD_)pMesh->Get_ActualModule_Heff_Display();
			if (pMesh->IsModuleSet(Module_Heff)) {

				return (DBL3)pMesh->pMod(Module_Heff)->Get_Module_Heff2CUDA()()->average_nonempty(pMesh->pMod(Module_Heff)->Get_Module_Heff2().linear_size(), (cuRect)rel_rect);
			}
		}
		break;

	case MESHDISPLAY_ENERGY:
	{
		MOD_ Module_Energy = (MOD_)pMesh->Get_ActualModule_Energy_Display();
		if (pMesh->IsModuleSet(Module_Energy)) {

			return (double)pMesh->pMod(Module_Energy)->Get_Module_EnergyCUDA()()->average_nonempty(pMesh->pMod(Module_Energy)->Get_Module_Energy().linear_size(), rel_rect);
		}
	}
	break;

	case MESHDISPLAY_ENERGY2:
	{
		MOD_ Module_Energy = (MOD_)pMesh->Get_ActualModule_Energy_Display();
		if (pMesh->IsModuleSet(Module_Energy)) {

			return (double)pMesh->pMod(Module_Energy)->Get_Module_Energy2CUDA()()->average_nonempty(pMesh->pMod(Module_Energy)->Get_Module_Energy2().linear_size(), rel_rect);
		}
	}
	break;

	case MESHDISPLAY_CURRDENSITY:

		//pdisplay_vec_vc_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {
			
			return dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetAverageChargeCurrent(rel_rect);
		}
		else if (pMesh->IsModuleSet(MOD_TMR)) {

			return dynamic_cast<TMR*>(pMesh->pMod(MOD_TMR))->GetAverageChargeCurrent(rel_rect);
		}
		break;

	case MESHDISPLAY_VOLTAGE:
		return (double)V()->average_nonempty(pMesh->V.linear_size(), (cuRect)rel_rect);
		break;

	case MESHDISPLAY_ELCOND:
		return (double)elC()->average_nonempty(pMesh->elC.linear_size(), (cuRect)rel_rect);
		break;

	case MESHDISPLAY_SACCUM:
		return (DBL3)S()->average_nonempty(pMesh->S.linear_size(), (cuRect)rel_rect);
		break;

	case MESHDISPLAY_JSX:

		//pdisplay_vec_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			return dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetAverageSpinCurrent(0, rel_rect);
		}
		break;

	case MESHDISPLAY_JSY:

		//pdisplay_vec_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			return dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetAverageSpinCurrent(1, rel_rect);
		}
		break;

	case MESHDISPLAY_JSZ:

		//pdisplay_vec_vec at maximum resolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			return dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetAverageSpinCurrent(2, rel_rect);
		}
		break;

	case MESHDISPLAY_TS:

		//pdisplay_vec_vec at maximumresolution
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			return dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetAverageSpinTorque(rel_rect);
		}
		break;

	case MESHDISPLAY_TSI:
#ifdef MODULE_COMPILATION_TRANSPORT
		//pdisplay_vec_vec at maximumresolution
		if (pMesh->pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && pMesh->IsModuleSet(MOD_TRANSPORT)) {

			//spin torque calculated internally in the Transport module, ready to be read out when needed
			dynamic_cast<STransport*>(pMesh->pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorque(dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT)));
			return dynamic_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetAverageInterfacialSpinTorque(rel_rect);
		}
#endif
		break;

	case MESHDISPLAY_TEMPERATURE:
		return (double)Temp()->average_nonempty(pMesh->Temp.linear_size(), (cuRect)rel_rect);
		break;

	case MESHDISPLAY_UDISP:
		return (DBL3)u_disp()->average_nonempty(pMesh->u_disp.linear_size(), (cuRect)rel_rect);
		break;

	case MESHDISPLAY_STRAINDIAG:
		return (DBL3)strain_diag()->average_nonempty(pMesh->strain_diag.linear_size(), (cuRect)rel_rect);
		break;

	case MESHDISPLAY_STRAINODIAG:
		return (DBL3)strain_odiag()->average_nonempty(pMesh->strain_odiag.linear_size(), (cuRect)rel_rect);
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

			return reinterpret_cast<VEC<double>*>(s_scaling)->average_nonempty_omp(rel_rect);
		}
		else {

			return reinterpret_cast<VEC<DBL3>*>(s_scaling)->average_nonempty_omp(rel_rect);
		}
	}
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

cu_obj<ManagedDiffEqFMCUDA>& MeshCUDA::Get_ManagedDiffEqFMCUDA(void)
{
	return dynamic_cast<DifferentialEquationFMCUDA*>(dynamic_cast<FMesh*>(pMesh)->Get_DifferentialEquation().Get_DifferentialEquationCUDA_ptr())->Get_ManagedDiffEqCUDA();
}

std::vector<DBL4>& MeshCUDA::get_tensorial_anisotropy(void)
{
	return pMesh->Kt;
}

std::vector<DBL4>& MeshCUDA::get_tensorial_anisotropy2(void)
{
	return pMesh->Kt2;
}

#endif