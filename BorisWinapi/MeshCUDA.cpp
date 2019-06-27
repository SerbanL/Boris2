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
	
	//effective field - sum total field of all the added modules
	if(!Heff()->set_from_cpuvec(pMesh->Heff)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//-----Electric conduction properties (Electron charge and spin Transport)

	//electrical potential - on n_e, h_e mesh
	if(!V()->set_from_cpuvec(pMesh->V)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//electrical conductivity - on n_e, h_e mesh
	if(!elC()->set_from_cpuvec(pMesh->elC)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

	//electrical current density - on n_e, h_e mesh
	if(!Jc()->set_from_cpuvec(pMesh->Jc)) error_on_create(BERROR_OUTOFGPUMEMORY_CRIT);

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

		//effective field - sum total field of all the added modules
		Heff()->copy_to_cpuvec(pMesh->Heff);

		//-----Electric conduction properties (Electron charge and spin Transport)

		//electrical potential - on n_e, h_e mesh
		V()->copy_to_cpuvec(pMesh->V);

		//electrical conductivity - on n_e, h_e mesh
		elC()->copy_to_cpuvec(pMesh->elC);

		//electrical current density - on n_e, h_e mesh
		Jc()->copy_to_cpuvec(pMesh->Jc);

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
			return PhysQ(pdisplay_vec_vc_vec, pMesh->displayedPhysicalQuantity);
		}
		break;
		
	case MESHDISPLAY_EFFECTIVEFIELD:
		
		if (prepare_display(n, meshRect, detail_level, Heff)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vec, pMesh->displayedPhysicalQuantity);
		}
		break;
		
	case MESHDISPLAY_CURRDENSITY:

		if (prepare_display(n_e, meshRect, detail_level, Jc)) {

			//return PhysQ made from the cpu version of coarse mesh display.
			return PhysQ(pdisplay_vec_vc_vec, pMesh->displayedPhysicalQuantity);
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
			return PhysQ(pdisplay_vec_vc_vec, pMesh->displayedPhysicalQuantity);
		}
		break;

	case MESHDISPLAY_JSX:
		
		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n_e, meshRect, detail_level, reinterpret_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(0))) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, pMesh->displayedPhysicalQuantity);
			}
		}
		break;

	case MESHDISPLAY_JSY:

		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n_e, meshRect, detail_level, reinterpret_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(1))) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, pMesh->displayedPhysicalQuantity);
			}
		}
		break;

	case MESHDISPLAY_JSZ:

		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n_e, meshRect, detail_level, reinterpret_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinCurrentCUDA(2))) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, pMesh->displayedPhysicalQuantity);
			}
		}
		break;

	case MESHDISPLAY_TS:

		if (pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n, meshRect, detail_level, reinterpret_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))->GetSpinTorqueCUDA())) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, pMesh->displayedPhysicalQuantity);
			}
		}
		break;

	case MESHDISPLAY_TSI:
		
		if (pMesh->pSMesh->IsSuperMeshModuleSet(MODS_STRANSPORT) && pMesh->IsModuleSet(MOD_TRANSPORT)) {

			if (prepare_display(n, meshRect, detail_level, 
				reinterpret_cast<STransport*>(pMesh->pSMesh->pSMod(MODS_STRANSPORT))->GetInterfacialSpinTorqueCUDA(reinterpret_cast<Transport*>(pMesh->pMod(MOD_TRANSPORT))))) {

				//return PhysQ made from the cpu version of coarse mesh display.
				return PhysQ(pdisplay_vec_vec, pMesh->displayedPhysicalQuantity);
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
			return PhysQ(reinterpret_cast<VEC<DBL3>*>(s_scaling), pMesh->displayedPhysicalQuantity);
	}
		break;

	case MESHDISPLAY_ROUGHNESS:
		if (pMesh->IsModuleSet(MOD_ROUGHNESS)) {

			return PhysQ(&reinterpret_cast<Roughness*>(pMesh->pMod(MOD_ROUGHNESS))->GetRoughness(), pMesh->displayedPhysicalQuantity);
		}
		break;
	}
	
	return PhysQ(meshRect, h, pMesh->displayedPhysicalQuantity);
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
	return (DBL2(pMesh->Gi).norm() > 0);
}

//check if the ODECommon::available flag is true (ode step solved)
bool MeshCUDA::CurrentTimeStepSolved(void)
{
	return pMesh->pSMesh->CurrentTimeStepSolved();
}

bool MeshCUDA::are_all_params_const(int numParams, ...)
{
	va_list params;

	va_start(params, numParams);

	for (int idx = 0; idx < numParams; idx++) {

		PARAM_ paramID = va_arg(params, PARAM_);

		if (pMesh->is_param_nonconst(paramID)) return false;
	}

	va_end(params);

	return true;
}

//----------------------------------- VALUE GETTERS

//get average magnetisation in given box (entire mesh if none specified)
cuReal3 MeshCUDA::GetAverageMagnetisation(cuRect rectangle)
{
	if (pMesh->M.linear_size()) return M()->average_nonempty(n.dim(), rectangle);
	else return cuReal3(0.0);
}

cuReal3 MeshCUDA::GetAverageChargeCurrentDensity(cuRect rectangle)
{
	if (pMesh->Jc.linear_size()) return Jc()->average_nonempty(n_e.dim(), rectangle);
	else return cuReal3(0.0);
}

cuReal MeshCUDA::GetAverageElectricalPotential(cuRect rectangle)
{
	if (pMesh->V.linear_size()) return V()->average_nonempty(n_e.dim(), rectangle);
	else return 0.0;
}

cuReal3 MeshCUDA::GetAverageSpinAccumulation(cuRect rectangle)
{
	if (pMesh->S.linear_size()) return S()->average_nonempty(n_e.dim(), rectangle);
	else return cuReal3(0.0);
}

cuReal MeshCUDA::GetAverageElectricalConductivity(cuRect rectangle)
{
	if (pMesh->elC.linear_size()) return elC()->average_nonempty(n_e.dim(), rectangle);
	else return 0.0;
}

cuReal MeshCUDA::GetAverageTemperature(cuRect rectangle)
{
	if (pMesh->Temp.linear_size()) return Temp()->average_nonempty(n_t.dim(), rectangle);
	else return pMesh->base_temperature;
}

//----------------------------------- OTHER MESH SHAPE CONTROL

//mask cells using bitmap image : white -> empty cells. black -> keep values. Apply mask up to given z depth number of cells depending on grayscale value (zDepth, all if 0).
BError MeshCUDA::applymask(cuReal zDepth_m, string fileName, function<vector<BYTE>(string, INT2)>& bitmap_loader)
{
	BError error;

	auto run_this = [](auto& VEC_quantity, auto default_value, cuReal zDepth_m, string fileName, function<vector<BYTE>(string, INT2)>& bitmap_loader) -> BError {

		BError error;

		cuINT3 cells = VEC_quantity()->size_cpu();
		cuReal3 cellsize = VEC_quantity()->cellsize_cpu();
		
		if (!VEC_quantity()->apply_bitmap_mask(bitmap_loader(fileName, INT2(cells.x, cells.y)), (int)round(zDepth_m / cellsize.z)))
			return error(BERROR_COULDNOTLOADFILE);
			
		return error;
	};

	error = change_mesh_shape(run_this, zDepth_m, fileName, bitmap_loader);

	if (!error) {

		//update mesh (and modules more importantly which will control any dependent VEC and VEC_VC quantites!)
		error = pMesh->pSMesh->UpdateConfiguration(UPDATECONFIG_MESHSHAPECHANGE);
	}

	return error;
}

//set cells to empty in given box (delete by setting entries to zero)
BError MeshCUDA::delrect(cuRect rectangle)
{
	BError error;

	auto run_this = [](auto& VEC_quantity, auto default_value, cuRect& rectangle) -> BError {

		VEC_quantity()->delrect(rectangle);
		return BError();
	};

	error = change_mesh_shape(run_this, rectangle);

	if (!error) {

		//update mesh (and modules more importantly which will control any dependent VEC and VEC_VC quantites!)
		error = pMesh->pSMesh->UpdateConfiguration(UPDATECONFIG_MESHSHAPECHANGE);
	}

	return error;
}

//set cells to non-empty in given box
BError MeshCUDA::setrect(cuRect rectangle)
{
	BError error;

	auto run_this = [](auto& VEC_quantity, auto default_value, cuRect& rectangle) -> BError {

		VEC_quantity()->setrect(rectangle, default_value);
		return BError();
	};

	error = change_mesh_shape(run_this, rectangle);

	if (!error) {

		//update mesh (and modules more importantly which will control any dependent VEC and VEC_VC quantites!)
		error = pMesh->pSMesh->UpdateConfiguration(UPDATECONFIG_MESHSHAPECHANGE);
	}

	return error;
}

//copy all meshes controlled using change_mesh_shape from cpu to gpu versions
BError MeshCUDA::copy_shapes_from_cpu(void)
{
	//Primary quantities are : M, elC, Temp

	BError error(__FUNCTION__);
	
	bool success = true;

	//1. shape magnetization
	if (M()->size_cpu().dim()) success &= M()->set_from_cpuvec(pMesh->M);

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

	//1. shape magnetization
	if (M()->size_cpu().dim()) success &= M()->set_cpuvec(pMesh->M);

	//2. shape electrical conductivity
	if (elC()->size_cpu().dim()) success &= elC()->set_cpuvec(pMesh->elC);

	//3. shape temperature
	if (Temp()->size_cpu().dim()) success &= Temp()->set_cpuvec(pMesh->Temp);

	//if adding any more here also remember to edit change_mesh_shape

	if (!success) error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

#endif