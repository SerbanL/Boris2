#include "stdafx.h"
#include "MeshCUDA.h"
#include "Mesh.h"
#include "SuperMesh.h"
#include "PhysQRep.h"

#if COMPILECUDA == 1

MeshCUDA::MeshCUDA(Mesh* pMesh) :
	MeshParamsCUDA(dynamic_cast<MeshParams*>(pMesh)),
	MeshDisplayCUDA(),
	meshRect(pMesh->meshRect),
	n(pMesh->n), h(pMesh->h),
	n_e(pMesh->n_e), h_e(pMesh->h_e),
	n_t(pMesh->n_t), h_t(pMesh->h_t)
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

	//-----Elastic properties
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

		//electrical potential - on n_e, h_e mesh
		V()->copy_to_cpuvec(pMesh->V);

		//electrical conductivity - on n_e, h_e mesh
		elC()->copy_to_cpuvec(pMesh->elC);

		//electrical field - on n_e, h_e mesh
		E()->copy_to_cpuvec(pMesh->E);

		//electrical current density - on n_e, h_e mesh
		S()->copy_to_cpuvec(pMesh->S);

		//-----Thermal conduction properties

		//temperature calculated by Heat module
		Temp()->copy_to_cpuvec(pMesh->Temp);

		//-----Elastic properties
	}
}

//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS

PhysQ MeshCUDA::FetchOnScreenPhysicalQuantity(double detail_level)
{
	switch (pMesh->displayedPhysicalQuantity) {

	case MESHDISPLAY_NONE:
		return PhysQ(meshRect, h, pMesh->displayedPhysicalQuantity);

	case MESHDISPLAY_MAGNETIZATION:
		
		if (prepare_display(n, meshRect, detail_level, M)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_vec, pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;

	case MESHDISPLAY_MAGNETIZATION2:

		if (prepare_display(n, meshRect, detail_level, M2)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_vec, pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;

	case MESHDISPLAY_MAGNETIZATION12:

		if (prepare_display(n, meshRect, detail_level, M, M2)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_vec, pdisplay2_vec_vc_vec, pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;
		
	case MESHDISPLAY_EFFECTIVEFIELD:
		
		if (prepare_display(n, meshRect, detail_level, Heff)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vec, pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;

	case MESHDISPLAY_EFFECTIVEFIELD2:

		if (prepare_display(n, meshRect, detail_level, Heff2)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vec, pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;

	case MESHDISPLAY_EFFECTIVEFIELD12:

		if (prepare_display(n, meshRect, detail_level, Heff, Heff2)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vec, pdisplay2_vec_vec, pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;
		
	case MESHDISPLAY_CURRDENSITY:

		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n_e, meshRect, detail_level, reinterpret_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetChargeCurrentCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vc_vec, pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		break;

	case MESHDISPLAY_VOLTAGE:

		if (prepare_display(n_e, meshRect, detail_level, V)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_sca, pMesh->displayedPhysicalQuantity);
		}
		break;

	case MESHDISPLAY_ELCOND:

		if (prepare_display(n_e, meshRect, detail_level, elC)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_sca, pMesh->displayedPhysicalQuantity);
		}
		break;

	case MESHDISPLAY_SACCUM:

		if (prepare_display(n_e, meshRect, detail_level, S)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_vec, pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
		}
		break;

	case MESHDISPLAY_JSX:
		
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n_e, meshRect, detail_level, reinterpret_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(0))) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		break;

	case MESHDISPLAY_JSY:

		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n_e, meshRect, detail_level, reinterpret_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(1))) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		break;

	case MESHDISPLAY_JSZ:

		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n_e, meshRect, detail_level, reinterpret_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(2))) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		break;

	case MESHDISPLAY_TS:

		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n, meshRect, detail_level, reinterpret_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinTorqueCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		break;

	case MESHDISPLAY_TSI:
		
		if (pMesh->pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n, meshRect, detail_level, 
				reinterpret_cast<STransport*>(pMesh->pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorqueCUDA(reinterpret_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))))) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
			}
		}
		break;

	case MESHDISPLAY_TEMPERATURE:
		
		if (prepare_display(n_t, meshRect, detail_level, Temp)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_sca, pMesh->displayedPhysicalQuantity);
		}
		break;

	case MESHDISPLAY_PARAMVAR:
	{
		void* s_scaling = pMesh->get_meshparam_s_scaling((PARAM_)pMesh->displayedParamVar);

		if (pMesh->is_s_scaling_scalar((PARAM_)pMesh->displayedParamVar))
			return PhysQ(reinterpret_cast<VEC<double>*>(s_scaling), pMesh->displayedPhysicalQuantity);
		else
			return PhysQ(reinterpret_cast<VEC<DBL3>*>(s_scaling), pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
	}
		break;

	case MESHDISPLAY_ROUGHNESS:
		if (pMesh->IsModuleSet(MOD_ROUGHNESS)) {

			return PhysQ(&reinterpret_cast<Roughness*>(pMesh->pMod(MOD_ROUGHNESS))->GetRoughness(), pMesh->displayedPhysicalQuantity);
		}
		break;

	case MESHDISPLAY_CUSTOM_VEC:
		return PhysQ(&pMesh->displayVEC_VEC, pMesh->displayedPhysicalQuantity, (VEC3REP_)pMesh->vec3rep);
		break;

	case MESHDISPLAY_CUSTOM_SCA:
		return PhysQ(&pMesh->displayVEC_SCA, pMesh->displayedPhysicalQuantity);
		break;
	}
	
	return PhysQ(meshRect, h, pMesh->displayedPhysicalQuantity);
}

//----------------------------------- MESH INFO GET/SET METHODS

int MeshCUDA::GetMeshType(void)
{
	return (int)pMesh->GetMeshType();
}

//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

//magnetization dynamics computation enabled
bool MeshCUDA::MComputation_Enabled(void)
{
	return pMesh->Heff.linear_size();
}

bool MeshCUDA::Magnetisation_Enabled(void)
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

bool MeshCUDA::GInterface_Enabled(void)
{
	return (DBL2(pMesh->Gmix.get0()).norm() > 0);
}

//check if the ODECommon::available flag is true (ode step solved)
bool MeshCUDA::CurrentTimeStepSolved(void)
{
	return pMesh->pSMesh->CurrentTimeStepSolved();
}

//check evaluation speedup flag in ODECommon
int MeshCUDA::EvaluationSpeedup(void)
{
	return pMesh->pSMesh->EvaluationSpeedup();
}

//check in ODECommon the type of field update we need to do depending on the ODE evaluation step
int MeshCUDA::Check_Step_Update(void)
{
	return pMesh->pSMesh->Check_Step_Update();
}

//----------------------------------- VALUE GETTERS

//get average magnetisation in given box (entire mesh if none specified)
cuReal3 MeshCUDA::GetAverageMagnetisation(cuRect rectangle)
{
	if (pMesh->M.linear_size()) return M()->average_nonempty(n.dim(), rectangle);
	else return cuReal3(0.0);
}

//get average magnetisation in given box (entire mesh if none specified); sub-lattice B
cuReal3 MeshCUDA::GetAverageMagnetisation2(cuRect rectangle)
{
	if (pMesh->M2.linear_size()) return M2()->average_nonempty(n.dim(), rectangle);
	else return cuReal3(0.0);
}

cuBReal MeshCUDA::GetAverageElectricalPotential(cuRect rectangle)
{
	if (pMesh->V.linear_size()) return V()->average_nonempty(n_e.dim(), rectangle);
	else return 0.0;
}

cuReal3 MeshCUDA::GetAverageSpinAccumulation(cuRect rectangle)
{
	if (pMesh->S.linear_size()) return S()->average_nonempty(n_e.dim(), rectangle);
	else return cuReal3(0.0);
}

cuBReal MeshCUDA::GetAverageElectricalConductivity(cuRect rectangle)
{
	if (pMesh->elC.linear_size()) return elC()->average_nonempty(n_e.dim(), rectangle);
	else return 0.0;
}

cuBReal MeshCUDA::GetAverageTemperature(cuRect rectangle)
{
	if (pMesh->Temp.linear_size()) return Temp()->average_nonempty(n_t.dim(), rectangle);
	else return pMesh->base_temperature;
}

//----------------------------------- OTHER MESH SHAPE CONTROL

//copy all meshes controlled using change_mesh_shape from cpu to gpu versions
BError MeshCUDA::copy_shapes_from_cpu(void)
{
	//Primary quantities are : M, elC, Temp

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

	//if adding any more here also remember to edit change_mesh_shape

	if (!success) error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

#endif