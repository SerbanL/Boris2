#include "stdafx.h"
#include "Heat.h"

#ifdef MODULE_HEAT

#include "Mesh.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

//-------------------Setters

void Heat::SetAmbientTemperature(double T_ambient_)
{
	if (IsZoP(T_ambient))
		T_ambient = T_ambient_;

	SetRobinBoundaryConditions();
}

void Heat::SetAlphaBoundary(double alpha_boundary_)
{
	if (IsZoP(alpha_boundary_))
		alpha_boundary = alpha_boundary_;

	SetRobinBoundaryConditions();
}

//set Temp uniformly to base temperature, unless a spatial variation is also specified through the cT mesh parameter
void Heat::SetBaseTemperature(double Temperature)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!pMesh->cT.is_sdep()) {

			reinterpret_cast<HeatCUDA*>(pModuleCUDA)->SetBaseTemperature(Temperature);
		}
		else {

			reinterpret_cast<HeatCUDA*>(pModuleCUDA)->SetBaseTemperature_Nonuniform(Temperature);
		}

		return;
	}
#endif

	if (!pMesh->cT.is_sdep()) {

		//uniform temperature
		pMesh->Temp.setnonempty(Temperature);

		if (pMesh->Temp_l.linear_size()) pMesh->Temp_l.setnonempty(Temperature);
	}
	else {

		//non-uniform temperature setting
#pragma omp parallel for
		for (int idx = 0; idx < pMesh->Temp.linear_size(); idx++) {

			if (pMesh->Temp.is_not_empty(idx)) {

				double cT = pMesh->cT;
				pMesh->update_parameters_tcoarse(idx, pMesh->cT, cT);

				pMesh->Temp[idx] = cT * Temperature;

				if (pMesh->Temp_l.linear_size()) {
					
					pMesh->Temp_l[idx] = cT * Temperature;
				}
			}
		}
	}
}

//Set Q_equation text equation object
BError Heat::SetQEquation(string equation_string, int step)
{
	BError error(CLASS_STR(Heat));

	DBL3 meshDim = pMesh->GetMeshDimensions();

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
	if (pModuleCUDA) error = reinterpret_cast<HeatCUDA*>(pModuleCUDA)->SetQEquation(Q_equation.get_scalar_fspec());
#endif

	return error;
}

//Update TEquation object with user constants values
void Heat::UpdateTEquationUserConstants(bool makeCuda)
{
	if (pMesh->userConstants.size()) {

		vector<pair<string, double>> constants(pMesh->userConstants.size());
		for (int idx = 0; idx < pMesh->userConstants.size(); idx++) {

			constants[idx] = { pMesh->userConstants.get_key_from_index(idx), pMesh->userConstants[idx] };
		}

		Q_equation.set_constants(constants);

		//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
		if (pModuleCUDA && makeCuda) reinterpret_cast<HeatCUDA*>(pModuleCUDA)->SetQEquation(Q_equation.get_scalar_fspec());
#endif
	}
}

//set insulating mesh sides flags. str can be "x", "-x", "y", "-y", "z", "-z"
void Heat::SetInsulatingSides(string literal, bool status)
{
	if (literal == "x") insulate_px = status;
	else if (literal == "-x") insulate_nx = status;
	else if (literal == "y") insulate_py = status;
	else if (literal == "-y") insulate_ny = status;
	else if (literal == "z") insulate_pz = status;
	else if (literal == "-z") insulate_nz = status;

	SetRobinBoundaryConditions();
}

bool Heat::GetInsulatingSide(string literal)
{
	if (literal == "x") return insulate_px;
	else if (literal == "-x") return insulate_nx;
	else if (literal == "y") return insulate_py;
	else if (literal == "-y") return insulate_ny;
	else if (literal == "z") return insulate_pz;
	else if (literal == "-z") return insulate_nz;

	return false;
}

//-------------------Others

//called by MoveMesh method in this mesh - move relevant transport quantities
void Heat::MoveMesh_Heat(double x_shift)
{
	double mesh_end_size = pMesh->meshRect.size().x * MOVEMESH_ENDRATIO;

	Rect shift_rect = Rect(pMesh->meshRect.s + DBL3(mesh_end_size, 0, 0), pMesh->meshRect.e - DBL3(mesh_end_size, 0, 0));

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pMesh->pMeshCUDA->Temp()->shift_x(pMesh->Temp.linear_size(), x_shift, shift_rect);

		if (pMesh->Temp_l.linear_size()) pMesh->pMeshCUDA->Temp_l()->shift_x(pMesh->Temp_l.linear_size(), x_shift, shift_rect);

		return;
	}
#endif

	pMesh->Temp.shift_x(x_shift, shift_rect);

	if (pMesh->Temp_l.linear_size()) pMesh->Temp_l.shift_x(x_shift, shift_rect);
}

void Heat::SetRobinBoundaryConditions(void)
{
	//set Robin boundary conditions
	DBL2 robin_values_px = DBL2(alpha_boundary, T_ambient) * (1 - insulate_px);
	DBL2 robin_values_nx = DBL2(alpha_boundary, T_ambient) * (1 - insulate_nx);
	DBL2 robin_values_py = DBL2(alpha_boundary, T_ambient) * (1 - insulate_py);
	DBL2 robin_values_ny = DBL2(alpha_boundary, T_ambient) * (1 - insulate_ny);
	DBL2 robin_values_pz = DBL2(alpha_boundary, T_ambient) * (1 - insulate_pz);
	DBL2 robin_values_nz = DBL2(alpha_boundary, T_ambient) * (1 - insulate_nz);

	pMesh->Temp.set_robin_conditions(DBL2(alpha_boundary, T_ambient), robin_values_px, robin_values_nx, robin_values_py, robin_values_ny, robin_values_pz, robin_values_nz);

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pMesh->pMeshCUDA->Temp()->set_robin_conditions(cuReal2(alpha_boundary, T_ambient), robin_values_px, robin_values_nx, robin_values_py, robin_values_ny, robin_values_pz, robin_values_nz);
	}
#endif
}

#endif