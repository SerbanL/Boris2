#include "stdafx.h"
#include "SuperMesh.h"

//--------------------------------------------------------- ELASTODYNAMICS SOLVER CONTROL : SuperMeshElastodynamics.cpp

//reset elastodynamics solver state
BError SuperMesh::Reset_ElSolver(void)
{
	BError error(__FUNCTION__);

	//all meshes
	for (int idx = 0; idx < pMesh.size(); idx++) {

		pMesh[idx]->CallModuleMethod(&MElastic::Reset_ElSolver);
	}

	return error;
}

//set text equation in given mesh for diagonal strain
BError SuperMesh::Set_Sd_Equation(std::string meshName, std::string text_equation)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	error = pMesh[meshName]->CallModuleMethod(&MElastic::Set_Sd_Equation, text_equation);
	//disable elastodynamics solver iteration
	CallModuleMethod(&SMElastic::set_el_dT, 0.0);

	return error;
}

//set text equation in given mesh for shear strain
BError SuperMesh::Set_Sod_Equation(std::string meshName, std::string text_equation)
{
	BError error(__FUNCTION__);

	if (!contains(meshName)) return error(BERROR_INCORRECTNAME);

	error = pMesh[meshName]->CallModuleMethod(&MElastic::Set_Sd_Equation, text_equation);
	//disable elastodynamics solver iteration
	CallModuleMethod(&SMElastic::set_el_dT, 0.0);

	return error;
}

//clear text equations in given mesh (or all meshes if meshName is empty)
BError SuperMesh::Clear_Sd_Sod_Equations(std::string meshName)
{
	BError error(__FUNCTION__);

	if (meshName.length() && !contains(meshName)) return error(BERROR_INCORRECTNAME);

	if (meshName.length()) {

		pMesh[meshName]->CallModuleMethod(&MElastic::Clear_Sd_Sod_Equations);
	}
	else {

		for (int idx = 0; idx < pMesh.size(); idx++) {

			pMesh[meshName]->CallModuleMethod(&MElastic::Clear_Sd_Sod_Equations);
		}
	}

	return error;
}

//add a fixed surface for elastodynamics solver (rectangle in absolute coordinates)
//can also specify a given face in a given mesh (face -x, x, -y, y, -z, or z)
BError SuperMesh::Add_Fixed_Surface(std::string meshName, std::string face, Rect surface_rect)
{
	BError error(__FUNCTION__);

	if (face.length() && !contains(meshName)) return error(BERROR_INCORRECTNAME);
	if (!IsSuperMeshModuleSet(MODS_SMELASTIC)) return error(BERROR_INCORRECTACTION);

	//get face rectangle from specified face in meshName
	if (face.length()) {

		if (face == "-x") surface_rect = pMesh[meshName]->meshRect.get_face_mx();
		else if (face == "x") surface_rect = pMesh[meshName]->meshRect.get_face_px();
		else if (face == "-y") surface_rect = pMesh[meshName]->meshRect.get_face_my();
		else if (face == "y") surface_rect = pMesh[meshName]->meshRect.get_face_py();
		else if (face == "-z") surface_rect = pMesh[meshName]->meshRect.get_face_mz();
		else if (face == "z") surface_rect = pMesh[meshName]->meshRect.get_face_pz();
		else return error(BERROR_INCORRECTNAME);
	}

	error = CallModuleMethod(&SMElastic::Add_Fixed_Surface, surface_rect);

	error = UpdateConfiguration(UPDATECONFIG_MELASTIC);

	return error;
}

//add a stress surface for elastodynamics solver (rectangle in absolute coordinates) with set vector equation
//can also specify a given face in a given mesh (face -x, x, -y, y, -z, or z)
BError SuperMesh::Add_Stress_Surface(std::string meshName, std::string face, Rect surface_rect, std::string equation)
{
	BError error(__FUNCTION__);

	if (face.length() && !contains(meshName)) return error(BERROR_INCORRECTNAME);
	if (!IsSuperMeshModuleSet(MODS_SMELASTIC)) return error(BERROR_INCORRECTACTION);

	//get face rectangle from specified face in meshName
	if (face.length()) {

		if (face == "-x") surface_rect = pMesh[meshName]->meshRect.get_face_mx();
		else if (face == "x") surface_rect = pMesh[meshName]->meshRect.get_face_px();
		else if (face == "-y") surface_rect = pMesh[meshName]->meshRect.get_face_my();
		else if (face == "y") surface_rect = pMesh[meshName]->meshRect.get_face_py();
		else if (face == "-z") surface_rect = pMesh[meshName]->meshRect.get_face_mz();
		else if (face == "z") surface_rect = pMesh[meshName]->meshRect.get_face_pz();
		else return error(BERROR_INCORRECTNAME);
	}

	error = CallModuleMethod(&SMElastic::Add_Stress_Surface, surface_rect, equation);

	error = UpdateConfiguration(UPDATECONFIG_MELASTIC);

	return error;
}

//Delete fixed or stressed surfaces with given index (-1 to delete all)
BError SuperMesh::Del_Fixed_Surface(int index)
{
	BError error(__FUNCTION__);

	if (!IsSuperMeshModuleSet(MODS_SMELASTIC)) return error(BERROR_INCORRECTACTION);

	if (index == -1) {

		while (CallModuleMethod(&SMElastic::Get_Num_Fixed_Surfaces)) {

			CallModuleMethod(&SMElastic::Del_Fixed_Surface, 0);
		}
	}
	else {

		if (index < CallModuleMethod(&SMElastic::Get_Num_Fixed_Surfaces)) CallModuleMethod(&SMElastic::Del_Fixed_Surface, index);
	}

	error = UpdateConfiguration(UPDATECONFIG_MELASTIC);

	return error;
}

BError SuperMesh::Del_Stress_Surface(int index)
{
	BError error(__FUNCTION__);

	if (!IsSuperMeshModuleSet(MODS_SMELASTIC)) return error(BERROR_INCORRECTACTION);

	if (index == -1) {

		while (CallModuleMethod(&SMElastic::Get_Num_Stress_Surfaces)) {

			CallModuleMethod(&SMElastic::Del_Stress_Surface, 0);
		}
	}
	else {

		if (index < CallModuleMethod(&SMElastic::Get_Num_Stress_Surfaces)) CallModuleMethod(&SMElastic::Del_Stress_Surface, index);
	}

	error = UpdateConfiguration(UPDATECONFIG_MELASTIC);

	return error;
}

//edit values of existing fixed or stress surface
BError SuperMesh::Edit_Fixed_Surface(int index, Rect rect)
{
	BError error(__FUNCTION__);

	if (!IsSuperMeshModuleSet(MODS_SMELASTIC)) return error(BERROR_INCORRECTACTION);

	if (index < CallModuleMethod(&SMElastic::Get_Num_Fixed_Surfaces)) CallModuleMethod(&SMElastic::Edit_Fixed_Surface, index, rect);

	error = UpdateConfiguration(UPDATECONFIG_MELASTIC);

	return error;
}

BError SuperMesh::Edit_Stress_Surface_Rectangle(int index, Rect rect)
{
	BError error(__FUNCTION__);

	if (!IsSuperMeshModuleSet(MODS_SMELASTIC)) return error(BERROR_INCORRECTACTION);

	if (index < CallModuleMethod(&SMElastic::Get_Num_Stress_Surfaces)) CallModuleMethod(&SMElastic::Edit_Stress_Surface_Rectangle, index, rect);

	error = UpdateConfiguration(UPDATECONFIG_MELASTIC);

	return error;
}

BError SuperMesh::Edit_Stress_Surface_Equation(int index, std::string equation)
{
	BError error(__FUNCTION__);

	if (!IsSuperMeshModuleSet(MODS_SMELASTIC)) return error(BERROR_INCORRECTACTION);

	if (index < CallModuleMethod(&SMElastic::Get_Num_Stress_Surfaces)) CallModuleMethod(&SMElastic::Edit_Stress_Surface_Equation, index, equation);

	error = UpdateConfiguration(UPDATECONFIG_MELASTIC);

	return error;
}