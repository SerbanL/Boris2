#pragma once

#include "cuVEC_VC.h"

//CMBNDInfo describes a composite media boundary contact between 2 meshes of same type, used to calculate values at CMBND cells using boundary conditions
struct CMBNDInfoCUDA {

	//index of contacting meshes as : (secondary mesh, primary mesh)
	cuINT2 mesh_idx;

	//primary cellsize shift perpendicular to contact (from first primary cell-center position add half this DBL3 to reach the contact interface)
	cuReal3 hshift_primary;

	//secondary cellsize shift perpendicular to contact (from contact interface position add half this DBL3 to reach the first position to read value from on the contacting mesh, then add the full DBL3 to reach the second position)
	cuReal3 hshift_secondary;

	//weights depending on primary and secondary cellsizes perpendicular to the interface; in particular : i: max{h_primary, h_secondary} / h_secondary and j: max{h_primary, h_secondary} / h_primary
	cuReal2 weights;

	//integer cell shift on primary side (add full INT3 to reach second cell on primary side to read value from)
	cuINT3 cell_shift;

	//Box containing all cells on primary side of contact (though not all may be CMBND cells so check flag before setting boundary conditions)
	cuBox cells_box;

	__host__ void construct_cu_obj(void) {}
	
	__host__ void construct_cu_obj(const CMBNDInfoCUDA& copyThis)
	{
		assign_cu_obj(copyThis);
	}

	//assignment operator
	__host__ void assign_cu_obj(const CMBNDInfoCUDA& copyThis)
	{
		gpu_to_gpu(mesh_idx, copyThis.mesh_idx);
		gpu_to_gpu(hshift_primary, copyThis.hshift_primary);
		gpu_to_gpu(hshift_secondary, copyThis.hshift_secondary);
		gpu_to_gpu(weights, copyThis.weights);
		gpu_to_gpu(cell_shift, copyThis.cell_shift);
		gpu_to_gpu(cells_box, copyThis.cells_box);
	}

	__host__ void destruct_cu_obj(void) {}

	template <typename CMBNDInfo>
	__host__ void copy_from_CMBNDInfo(CMBNDInfo& cmbndInfo)
	{
		set_gpu_value(mesh_idx, (cuINT2)cmbndInfo.mesh_idx);

		set_gpu_value(hshift_primary, (cuReal3)cmbndInfo.hshift_primary);
		set_gpu_value(hshift_secondary, (cuReal3)cmbndInfo.hshift_secondary);

		set_gpu_value(weights, (cuReal2)cmbndInfo.weights);

		set_gpu_value(cell_shift, (cuINT3)cmbndInfo.cell_shift);
		set_gpu_value(cells_box, (cuBox)cmbndInfo.cells_box);
	}

	template <typename CMBNDInfo>
	CMBNDInfoCUDA(const CMBNDInfo& copyThis)
	{
		*this = copyThis;
	}

	template <typename CMBNDInfo>
	CMBNDInfoCUDA& operator=(const CMBNDInfo& copyThis)
	{
		mesh_idx = copyThis.mesh_idx;
		hshift_primary = copyThis.hshift_primary;
		hshift_secondary = copyThis.hshift_secondary;
		weights = copyThis.weights;
		cell_shift = copyThis.cell_shift;
		cells_box = copyThis.cells_box;

		return *this;
	}

	//check if primary is on +ve or -ve side of contact (e.g. top or bottom)
	__host__ __device__ bool IsPrimaryTop(void) { return (cell_shift >= cuINT3(0)); }
};