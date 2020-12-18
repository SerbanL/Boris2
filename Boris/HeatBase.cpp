#include "stdafx.h"
#include "HeatBase.h"

#ifdef MODULE_COMPILATION_HEAT

#include "MeshBase.h"

#if COMPILECUDA == 1
#include "HeatBaseCUDA.h"
#endif

HeatBase::HeatBase(MeshBase *pMeshBase_) :
	Q_equation({ "x", "y", "z", "t" })
{
	pMeshBase = pMeshBase_;
}

//-------------------Setters

void HeatBase::SetAmbientTemperature(double T_ambient_)
{
	if (IsZoP(T_ambient)) T_ambient = T_ambient_;

	SetRobinBoundaryConditions();
}

void HeatBase::SetAlphaBoundary(double alpha_boundary_)
{
	if (IsZoP(alpha_boundary_)) alpha_boundary = alpha_boundary_;

	SetRobinBoundaryConditions();
}

//Set Q_equation text equation object
BError HeatBase::SetQEquation(std::string equation_string, int step)
{
	BError error(CLASS_STR(HeatBase));

	DBL3 meshDim = pMeshBase->GetMeshDimensions();

	//set equation if not already set, or this is the first step in a stage
	if (!Q_equation.is_set() || step == 0) {

		//important to set user constants first : if these are required then the make_from_string call will fail. This may mean the user constants are not set when we expect them to.
		UpdateTEquationUserConstants(false);

		if (!Q_equation.make_from_string(equation_string, { {"Lx", meshDim.x}, {"Ly", meshDim.y}, {"Lz", meshDim.z}, {"Ss", 0} })) return error(BERROR_INCORRECTSTRING);
	}
	else {

		//equation set and not the first step : adjust Ss constant
		Q_equation.set_constant("Ss", step);
		UpdateTEquationUserConstants(false);
	}

	//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pModuleCUDA) error = dynamic_cast<HeatBaseCUDA*>(pModuleCUDA)->SetQEquation(Q_equation.get_scalar_fspec());
#endif

	return error;
}

//Update TEquation object with user constants values
void HeatBase::UpdateTEquationUserConstants(bool makeCuda)
{
	if (pMeshBase->userConstants.size()) {

		std::vector<std::pair<std::string, double>> constants(pMeshBase->userConstants.size());
		for (int idx = 0; idx < pMeshBase->userConstants.size(); idx++) {

			constants[idx] = { pMeshBase->userConstants.get_key_from_index(idx), pMeshBase->userConstants[idx] };
		}

		Q_equation.set_constants(constants);

		//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
		if (pModuleCUDA && makeCuda) dynamic_cast<HeatBaseCUDA*>(pModuleCUDA)->SetQEquation(Q_equation.get_scalar_fspec());
#endif
	}
}

//set insulating mesh sides flags. str can be "x", "-x", "y", "-y", "z", "-z"
void HeatBase::SetInsulatingSides(std::string literal, bool status)
{
	if (literal == "x") insulate_px = status;
	else if (literal == "-x") insulate_nx = status;
	else if (literal == "y") insulate_py = status;
	else if (literal == "-y") insulate_ny = status;
	else if (literal == "z") insulate_pz = status;
	else if (literal == "-z") insulate_nz = status;

	SetRobinBoundaryConditions();
}

void HeatBase::SetRobinBoundaryConditions(void)
{
	//set Robin boundary conditions
	DBL2 robin_values_px = DBL2(alpha_boundary, T_ambient) * (1 - insulate_px);
	DBL2 robin_values_nx = DBL2(alpha_boundary, T_ambient) * (1 - insulate_nx);
	DBL2 robin_values_py = DBL2(alpha_boundary, T_ambient) * (1 - insulate_py);
	DBL2 robin_values_ny = DBL2(alpha_boundary, T_ambient) * (1 - insulate_ny);
	DBL2 robin_values_pz = DBL2(alpha_boundary, T_ambient) * (1 - insulate_pz);
	DBL2 robin_values_nz = DBL2(alpha_boundary, T_ambient) * (1 - insulate_nz);

	pMeshBase->Temp.set_robin_conditions(DBL2(alpha_boundary, T_ambient), robin_values_px, robin_values_nx, robin_values_py, robin_values_ny, robin_values_pz, robin_values_nz);

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pMeshBase->pMeshBaseCUDA->Temp()->set_robin_conditions(cuReal2(alpha_boundary, T_ambient), robin_values_px, robin_values_nx, robin_values_py, robin_values_ny, robin_values_pz, robin_values_nz);
	}
#endif
}

//set temperature solver type
BError HeatBase::Set_TMType(TMTYPE_ tmtype_)
{
	BError error(CLASS_STR(HeatBase));

	if (tmtype_ == TMTYPE_DEFAULT) {

		tmtype = TMTYPE_1TM;
	}
	else {

		tmtype = tmtype_;
	}

	//Some temperature models are only appropriate for certain types of meshes
	switch (pMeshBase->GetMeshType()) {

	case MESH_INSULATOR:
	case MESH_DIPOLE:
		tmtype = TMTYPE_1TM;
		break;

		//Currently all other types of meshes can support the 2TM and 3TM models
	}

	error = UpdateConfiguration(UPDATECONFIG_HEAT_MODELTYPE);

	return error;
}

#endif