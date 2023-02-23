#include "stdafx.h"
#include "Atom_Mesh_Cubic.h"
#include "SuperMesh.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC

//compute topological charge density spatial dependence and have it available to display in Cust_S
//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
void Atom_Mesh_Cubic::Compute_TopoChargeDensity(void)
{
	if (M1.linear_size()) {

		displayVEC_SCA.resize(h, meshRect);

#if COMPILECUDA == 1
		if (paMeshCUDA) {

			//compute topological charge density spatial variation in MeshCUDA::aux_vec_sca
			paMeshCUDA->Compute_TopoChargeDensity();

			//now transfer it to displayVEC_SCA ready to display
			paMeshCUDA->copy_aux_vec_sca(displayVEC_SCA);
			return;
		}
#endif

#pragma omp parallel for
		for (int idx = 0; idx < M1.linear_size(); idx++) {

			if (M1.is_not_empty(idx)) {

				double Mnorm = M1[idx].norm();

				DBL33 M_grad = M1.grad_neu(idx);

				DBL3 dm_dx = M_grad.x / Mnorm;
				DBL3 dm_dy = M_grad.y / Mnorm;

				//divide by number of z cells (this is intended for 2D layers)
				displayVEC_SCA[idx] = (M1[idx] / Mnorm) * (dm_dx ^ dm_dy) * M1.h.x * M1.h.y / (4 * PI * M1.n.z);
			}
			else displayVEC_SCA[idx] = 0.0;
		}
	}
}

//return phase transition temperature (K) based on formula Tc = J*e*z/3kB
double Atom_Mesh_Cubic::Show_Transition_Temperature(void)
{
	return J * spinwave_factor * coordination_number / (3 * BOLTZMANN);
}

//return saturation magnetization (A/m) based on formula Ms = mu_s*n/a^3
double Atom_Mesh_Cubic::Show_Ms(void)
{
	return mu_s * (MUB / h.dim()) * atoms_per_cell;
}

//return exchange stiffness (J/m) based on formula A = J*n/2a
double Atom_Mesh_Cubic::Show_A(void)
{
	//formula applicable for a cubic cell
	double h_av = (h.x + h.y + h.z) / 3;

	return J * atoms_per_cell / (2 * h_av);
}

//return uniaxial anisotropy constant (J/m^3) based on formula K = k*n/a^3
double Atom_Mesh_Cubic::Show_Ku(void)
{
	return K1 * atoms_per_cell / h.dim();
}

//Fit domain wall along the x direction through centre of rectangle : fit the component which matches a tanh profile. Return centre position and width.
DBL2 Atom_Mesh_Cubic::FitDomainWall_X(Rect rectangle)
{
	if (rectangle.IsNull()) rectangle = meshRect;

#if COMPILECUDA == 1
	if (paMeshCUDA) {

		size_t size = round((rectangle.e.x - rectangle.s.x - h.x) / h.x) + 1;
		cu_arr<cuReal2>& cu_xy_data = dwPos.getcuda_xy_data_ref(size);
		
		switch (pSMesh->Get_DWPos_Component())
		{
		case 0:
			if (!paMeshCUDA->M1()->extract_profile_component_x(
				cuReal3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				cuReal3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.x, cuReal3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z),
				cu_xy_data)) return DBL2();
			break;

		case 1:
			if (!paMeshCUDA->M1()->extract_profile_component_y(
				cuReal3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				cuReal3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.x, cuReal3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z),
				cu_xy_data)) return DBL2();
			break;

		case 2:
			if (!paMeshCUDA->M1()->extract_profile_component_z(
				cuReal3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				cuReal3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.x, cuReal3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z),
				cu_xy_data)) return DBL2();
			break;

		default:
		case -1:
			if (!paMeshCUDA->M1()->extract_profile_component_max(
				cuReal3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				cuReal3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.x, cuReal3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z),
				cu_xy_data)) return DBL2();
			break;
		}

		return dwPos.FitDomainWallCUDA(rectangle.e.x - rectangle.s.x - h.x);
	}
	else {

		switch (pSMesh->Get_DWPos_Component())
		{
		case 0:
		{
			std::vector<DBL2>& xy_data = M1.extract_profile_component_x(
				DBL3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				DBL3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.x, DBL3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z));

			return dwPos.FitDomainWall(xy_data);
		}
			break;

		case 1:
		{
			std::vector<DBL2>& xy_data = M1.extract_profile_component_y(
				DBL3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				DBL3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.x, DBL3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z));

			return dwPos.FitDomainWall(xy_data);
		}
			break;

		case 2:
		{
			std::vector<DBL2>& xy_data = M1.extract_profile_component_z(
				DBL3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				DBL3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.x, DBL3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z));

			return dwPos.FitDomainWall(xy_data);
		}
			break;

		default:
		case -1:
		{
			std::vector<DBL2>& xy_data = M1.extract_profile_component_max(
				DBL3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				DBL3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.x, DBL3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z));

			return dwPos.FitDomainWall(xy_data);
		}
			break;
		}
	}
#else
	switch (pSMesh->Get_DWPos_Component())
	{
	case 0:
	{
		std::vector<DBL2>& xy_data = M1.extract_profile_component_x(
			DBL3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
			DBL3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
			h.x, DBL3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	case 1:
	{
		std::vector<DBL2>& xy_data = M1.extract_profile_component_y(
			DBL3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
			DBL3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
			h.x, DBL3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z));

		return dwPos.FitDomainWall(xy_data);
		}
	break;

	case 2:
	{
		std::vector<DBL2>& xy_data = M1.extract_profile_component_z(
			DBL3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
			DBL3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
			h.x, DBL3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	default:
	case -1:
	{
		std::vector<DBL2>& xy_data = M1.extract_profile_component_max(
			DBL3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
			DBL3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
			h.x, DBL3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;
	}
#endif
}

//Fit domain wall along the y direction through centre of rectangle : fit the component which matches a tanh profile. Return centre position and width.
DBL2 Atom_Mesh_Cubic::FitDomainWall_Y(Rect rectangle)
{
	if (rectangle.IsNull()) rectangle = meshRect;

#if COMPILECUDA == 1
	if (paMeshCUDA) {

		size_t size = round((rectangle.e.y - rectangle.s.y - h.y) / h.y) + 1;
		cu_arr<cuReal2>& cu_xy_data = dwPos.getcuda_xy_data_ref(size);

		switch (pSMesh->Get_DWPos_Component())
		{
		case 0:
			if (!paMeshCUDA->M1()->extract_profile_component_x(
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.y, cuReal3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z),
				cu_xy_data)) return DBL2();
			break;

		case 1:
			if (!paMeshCUDA->M1()->extract_profile_component_y(
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.y, cuReal3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z),
				cu_xy_data)) return DBL2();
			break;

		case 2:
			if (!paMeshCUDA->M1()->extract_profile_component_z(
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.y, cuReal3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z),
				cu_xy_data)) return DBL2();
			break;

		default:
		case -1:
			if (!paMeshCUDA->M1()->extract_profile_component_max(
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.y, cuReal3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z),
				cu_xy_data)) return DBL2();
			break;
		}

		return dwPos.FitDomainWallCUDA(rectangle.e.y - rectangle.s.y - h.y);
	}
	else {

		switch (pSMesh->Get_DWPos_Component())
		{
		case 0:
		{
			std::vector<DBL2>& xy_data = M1.extract_profile_component_x(
				DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.y, DBL3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		case 1:
		{
			std::vector<DBL2>& xy_data = M1.extract_profile_component_y(
				DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.y, DBL3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		case 2:
		{
			std::vector<DBL2>& xy_data = M1.extract_profile_component_z(
				DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.y, DBL3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		default:
		case -1:
		{
			std::vector<DBL2>& xy_data = M1.extract_profile_component_max(
				DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.y, DBL3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;
		}
	}
#else
	switch (pSMesh->Get_DWPos_Component())
	{
	case 0:
	{
		std::vector<DBL2>& xy_data = M1.extract_profile_component_x(
			DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
			DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
			h.y, DBL3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	case 1:
	{
		std::vector<DBL2>& xy_data = M1.extract_profile_component_y(
			DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
			DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
			h.y, DBL3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	case 2:
	{
		std::vector<DBL2>& xy_data = M1.extract_profile_component_z(
			DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
			DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
			h.y, DBL3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	default:
	case -1:
	{
		std::vector<DBL2>& xy_data = M1.extract_profile_component_max(
			DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
			DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
			h.y, DBL3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;
	}
#endif
}

//Fit domain wall along the z direction through centre of rectangle : fit the component which matches a tanh profile. Return centre position and width.
DBL2 Atom_Mesh_Cubic::FitDomainWall_Z(Rect rectangle)
{
	if (rectangle.IsNull()) rectangle = meshRect;

#if COMPILECUDA == 1
	if (paMeshCUDA) {

		size_t size = round((rectangle.e.z - rectangle.s.z - h.z) / h.z) + 1;
		cu_arr<cuReal2>& cu_xy_data = dwPos.getcuda_xy_data_ref(size);

		switch (pSMesh->Get_DWPos_Component())
		{
		case 0:
			if (!paMeshCUDA->M1()->extract_profile_component_x(
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
				h.z, cuReal3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z),
				cu_xy_data)) return DBL2();
			break;

		case 1:
			if (!paMeshCUDA->M1()->extract_profile_component_y(
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
				h.z, cuReal3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z),
				cu_xy_data)) return DBL2();
			break;

		case 2:
			if (!paMeshCUDA->M1()->extract_profile_component_z(
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
				h.z, cuReal3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z),
				cu_xy_data)) return DBL2();
			break;

		default:
		case -1:
			if (!paMeshCUDA->M1()->extract_profile_component_max(
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
				h.z, cuReal3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z),
				cu_xy_data)) return DBL2();
			break;
		}

		return dwPos.FitDomainWallCUDA(rectangle.e.z - rectangle.s.z - h.z);
	}
	else {

		switch (pSMesh->Get_DWPos_Component())
		{
		case 0:
		{
			std::vector<DBL2>& xy_data = M1.extract_profile_component_x(
				DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
				DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
				h.z, DBL3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		case 1:
		{
			std::vector<DBL2>& xy_data = M1.extract_profile_component_y(
				DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
				DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
				h.z, DBL3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		case 2:
		{
			std::vector<DBL2>& xy_data = M1.extract_profile_component_z(
				DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
				DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
				h.z, DBL3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		default:
		case -1:
		{
			std::vector<DBL2>& xy_data = M1.extract_profile_component_max(
				DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
				DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
				h.z, DBL3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;
		}
}
#else
	switch (pSMesh->Get_DWPos_Component())
	{
	case 0:
	{
		std::vector<DBL2>& xy_data = M1.extract_profile_component_x(
			DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
			DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
			h.z, DBL3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	case 1:
	{
		std::vector<DBL2>& xy_data = M1.extract_profile_component_y(
			DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
			DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
			h.z, DBL3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	case 2:
	{
		std::vector<DBL2>& xy_data = M1.extract_profile_component_z(
			DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
			DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
			h.z, DBL3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	default:
	case -1:
	{
		std::vector<DBL2>& xy_data = M1.extract_profile_component_max(
			DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
			DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
			h.z, DBL3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;
	}
#endif
}

#endif