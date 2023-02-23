#include "stdafx.h"
#include "Mesh_AntiFerromagnetic.h"
#include "SuperMesh.h"

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

//Fit domain wall along the x direction through centre of rectangle : fit the component which matches a tanh profile. Return centre position and width.
DBL2 AFMesh::FitDomainWall_X(Rect rectangle)
{
	if (rectangle.IsNull()) rectangle = meshRect;

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		size_t size = round((rectangle.e.x - rectangle.s.x - h.x) / h.x) + 1;
		cu_arr<cuReal2>& cu_xy_data = dwPos.getcuda_xy_data_ref(size);

		switch (pSMesh->Get_DWPos_Component())
		{
		case 0:
			if (!pMeshCUDA->M()->extract_profile_component_x(
				cuReal3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				cuReal3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.x, cuReal3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z),
				cu_xy_data)) return DBL2();
			break;

		case 1:
			if (!pMeshCUDA->M()->extract_profile_component_y(
				cuReal3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				cuReal3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.x, cuReal3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z),
				cu_xy_data)) return DBL2();
			break;

		case 2:
			if (!pMeshCUDA->M()->extract_profile_component_z(
				cuReal3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				cuReal3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.x, cuReal3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z),
				cu_xy_data)) return DBL2();
			break;

		default:
		case -1:
			if (!pMeshCUDA->M()->extract_profile_component_max(
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
			std::vector<DBL2>& xy_data = M.extract_profile_component_x(
				DBL3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				DBL3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.x, DBL3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		case 1:
		{
			std::vector<DBL2>& xy_data = M.extract_profile_component_y(
				DBL3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				DBL3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.x, DBL3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		case 2:
		{
			std::vector<DBL2>& xy_data = M.extract_profile_component_z(
				DBL3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				DBL3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.x, DBL3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		default:
		case -1:
		{
			std::vector<DBL2>& xy_data = M.extract_profile_component_max(
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
		std::vector<DBL2>& xy_data = M.extract_profile_component_x(
			DBL3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
			DBL3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
			h.x, DBL3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	case 1:
	{
		std::vector<DBL2>& xy_data = M.extract_profile_component_y(
			DBL3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
			DBL3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
			h.x, DBL3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	case 2:
	{
		std::vector<DBL2>& xy_data = M.extract_profile_component_z(
			DBL3(rectangle.s.x + h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
			DBL3(rectangle.e.x - h.x / 2, (rectangle.s.y + rectangle.e.y) / 2, (rectangle.s.z + rectangle.e.z) / 2),
			h.x, DBL3(h.x, rectangle.e.y - rectangle.s.y, rectangle.e.z - rectangle.s.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	default:
	case -1:
	{
		std::vector<DBL2>& xy_data = M.extract_profile_component_max(
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
DBL2 AFMesh::FitDomainWall_Y(Rect rectangle)
{
	if (rectangle.IsNull()) rectangle = meshRect;

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		size_t size = round((rectangle.e.y - rectangle.s.y - h.y) / h.y) + 1;
		cu_arr<cuReal2>& cu_xy_data = dwPos.getcuda_xy_data_ref(size);

		switch (pSMesh->Get_DWPos_Component())
		{
		case 0:
			if (!pMeshCUDA->M()->extract_profile_component_x(
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.y, cuReal3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z),
				cu_xy_data)) return DBL2();
			break;

		case 1:
			if (!pMeshCUDA->M()->extract_profile_component_y(
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.y, cuReal3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z),
				cu_xy_data)) return DBL2();
			break;

		case 2:
			if (!pMeshCUDA->M()->extract_profile_component_z(
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.y, cuReal3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z),
				cu_xy_data)) return DBL2();
			break;

		default:
		case -1:
			if (!pMeshCUDA->M()->extract_profile_component_max(
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
			std::vector<DBL2>& xy_data = M.extract_profile_component_x(
				DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.y, DBL3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		case 1:
		{
			std::vector<DBL2>& xy_data = M.extract_profile_component_y(
				DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.y, DBL3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		case 2:
		{
			std::vector<DBL2>& xy_data = M.extract_profile_component_z(
				DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
				h.y, DBL3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		default:
		case -1:
		{
			std::vector<DBL2>& xy_data = M.extract_profile_component_max(
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
		std::vector<DBL2>& xy_data = M.extract_profile_component_x(
			DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
			DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
			h.y, DBL3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	case 1:
	{
		std::vector<DBL2>& xy_data = M.extract_profile_component_y(
			DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
			DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
			h.y, DBL3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	case 2:
	{
		std::vector<DBL2>& xy_data = M.extract_profile_component_z(
			DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.s.y + h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
			DBL3((rectangle.s.x + rectangle.e.x) / 2, rectangle.e.y - h.y / 2, (rectangle.s.z + rectangle.e.z) / 2),
			h.y, DBL3(rectangle.e.x - rectangle.s.x, h.y, rectangle.e.z - rectangle.s.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	default:
	case -1:
	{
		std::vector<DBL2>& xy_data = M.extract_profile_component_max(
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
DBL2 AFMesh::FitDomainWall_Z(Rect rectangle)
{
	if (rectangle.IsNull()) rectangle = meshRect;

#if COMPILECUDA == 1
	if (pMeshCUDA) {

		size_t size = round((rectangle.e.z - rectangle.s.z - h.z) / h.z) + 1;
		cu_arr<cuReal2>& cu_xy_data = dwPos.getcuda_xy_data_ref(size);

		switch (pSMesh->Get_DWPos_Component())
		{
		case 0:
			if (!pMeshCUDA->M()->extract_profile_component_x(
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
				h.z, cuReal3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z),
				cu_xy_data)) return DBL2();
			break;

		case 1:
			if (!pMeshCUDA->M()->extract_profile_component_y(
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
				h.z, cuReal3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z),
				cu_xy_data)) return DBL2();
			break;

		case 2:
			if (!pMeshCUDA->M()->extract_profile_component_z(
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
				cuReal3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
				h.z, cuReal3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z),
				cu_xy_data)) return DBL2();
			break;

		default:
		case -1:
			if (!pMeshCUDA->M()->extract_profile_component_max(
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
			std::vector<DBL2>& xy_data = M.extract_profile_component_x(
				DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
				DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
				h.z, DBL3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		case 1:
		{
			std::vector<DBL2>& xy_data = M.extract_profile_component_y(
				DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
				DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
				h.z, DBL3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		case 2:
		{
			std::vector<DBL2>& xy_data = M.extract_profile_component_z(
				DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
				DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
				h.z, DBL3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z));

			return dwPos.FitDomainWall(xy_data);
		}
		break;

		default:
		case -1:
		{
			std::vector<DBL2>& xy_data = M.extract_profile_component_max(
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
		std::vector<DBL2>& xy_data = M.extract_profile_component_x(
			DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
			DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
			h.z, DBL3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	case 1:
	{
		std::vector<DBL2>& xy_data = M.extract_profile_component_y(
			DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
			DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
			h.z, DBL3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	case 2:
	{
		std::vector<DBL2>& xy_data = M.extract_profile_component_z(
			DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.s.z + h.z / 2),
			DBL3((rectangle.s.x + rectangle.e.x) / 2, (rectangle.s.y + rectangle.e.y) / 2, rectangle.e.z - h.z / 2),
			h.z, DBL3(rectangle.e.x - rectangle.s.x, rectangle.e.y - rectangle.s.y, h.z));

		return dwPos.FitDomainWall(xy_data);
	}
	break;

	default:
	case -1:
	{
		std::vector<DBL2>& xy_data = M.extract_profile_component_max(
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