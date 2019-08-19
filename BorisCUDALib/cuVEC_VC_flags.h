#pragma once

#include "cuVEC_VC.h"

//--------------------------------------------IMPORTANT FLAG MANIPULATION METHODS

//from NF_DIRICHLET type flag and cell_idx return boundary value from one of the dirichlet vectors
template <typename VType>
__device__ VType cuVEC_VC<VType>::get_dirichlet_value(int dirichlet_flag, int cell_idx) const
{
	switch (dirichlet_flag) {

	case NF_DIRICHLETPX:
		return dirichlet_nx[((cell_idx / n.x) % n.y) + (cell_idx / (n.x*n.y))*n.y];

	case NF_DIRICHLETNX:
		return dirichlet_px[((cell_idx / n.x) % n.y) + (cell_idx / (n.x*n.y))*n.y];

	case NF_DIRICHLETPY:
		return dirichlet_ny[(cell_idx % n.x) + (cell_idx / (n.x*n.y))*n.x];

	case NF_DIRICHLETNY:
		return dirichlet_py[(cell_idx % n.x) + (cell_idx / (n.x*n.y))*n.x];

	case NF_DIRICHLETPZ:
		return dirichlet_nz[(cell_idx % n.x) + ((cell_idx / n.x) % n.y)*n.x];

	case NF_DIRICHLETNZ:
		return dirichlet_pz[(cell_idx % n.x) + ((cell_idx / n.x) % n.y)*n.x];
	}

	return VType();
}

template <typename VType>
__device__ VType cuVEC_VC<VType>::get_dirichlet_value(int dirichlet_flag, const cuINT3& ijk) const
{
	switch (dirichlet_flag) {

	case NF_DIRICHLETPX:
		return dirichlet_nx[ijk.j + ijk.k * n.y];

	case NF_DIRICHLETNX:
		return dirichlet_px[ijk.j + ijk.k * n.y];

	case NF_DIRICHLETPY:
		return dirichlet_ny[ijk.i + ijk.k * n.x];

	case NF_DIRICHLETNY:
		return dirichlet_py[ijk.i + ijk.k * n.x];

	case NF_DIRICHLETPZ:
		return dirichlet_nz[ijk.i + ijk.j * n.x];

	case NF_DIRICHLETNZ:
		return dirichlet_pz[ijk.i + ijk.j * n.x];
	}

	return VType();
}

//--------------------------------------------FLAG CHECKING

template <typename VType>
__host__ int cuVEC_VC<VType>::get_nonempty_cells_cpu(void)
{
	return get_gpu_value(nonempty_cells);
}

//check if all cells intersecting the rectangle (absolute coordinates) are empty
template <typename VType>
__device__ bool cuVEC_VC<VType>::is_empty(const cuRect& rectangle) const
{
	cuBox cells = box_from_rect_max(rectangle);

	for (int i = (cells.s.x >= 0 ? cells.s.x : 0); i < (cells.e.x <= n.x ? cells.e.x : n.x); i++) {

		for (int j = (cells.s.y >= 0 ? cells.s.y : 0); j < (cells.e.y <= n.y ? cells.e.y : n.y); j++) {
			
			for (int k = (cells.s.z >= 0 ? cells.s.z : 0); k < (cells.e.z <= n.z ? cells.e.z : n.z); k++) {

				if (is_not_empty(cuINT3(i, j, k))) return false;
			}
		}
	}

	return true;
}

//--------------------------------------------SET CELL FLAGS - EXTERNAL USE

//clear all dirichlet flags and vectors
template <typename VType>
__host__ void cuVEC_VC<VType>::clear_dirichlet_flags(void)
{
	gpu_and_managed(ngbrFlags, ~NF_DIRICHLET, get_ngbrFlags_size());

	set_dirichlet_size(0, NF_DIRICHLETPX);
	set_dirichlet_size(0, NF_DIRICHLETNX);
	set_dirichlet_size(0, NF_DIRICHLETPY);
	set_dirichlet_size(0, NF_DIRICHLETNY);
	set_dirichlet_size(0, NF_DIRICHLETPZ);
	set_dirichlet_size(0, NF_DIRICHLETNZ);
}

//clear all pbc flags : can also be achieved setting all flags to false in set_pbc but this one is more readable
template <typename VType>
__host__ void cuVEC_VC<VType>::clear_pbc(void)
{
	set_gpu_value(pbc_x, (int)0);
	set_gpu_value(pbc_y, (int)0);
	set_gpu_value(pbc_z, (int)0);

	gpu_and_managed(ngbrFlags, ~NF_PBC, get_ngbrFlags_size());
}

//set pbc conditions : setting any to false clears flags
template <typename VType>
__host__ void cuVEC_VC<VType>::set_pbc(bool pbc_x_, bool pbc_y_, bool pbc_z_)
{
	set_gpu_value(pbc_x, (int)pbc_x_);
	set_gpu_value(pbc_y, (int)pbc_y_);
	set_gpu_value(pbc_z, (int)pbc_z_);

	set_pbc_flags();
}

template <typename VType>
__host__ void cuVEC_VC<VType>::clear_cmbnd_flags(void)
{
	gpu_and_managed(ngbrFlags, ~NF_CMBND, get_ngbrFlags_size());
}

//clear all skip cell flags
template <typename VType>
__host__ void cuVEC_VC<VType>::clear_skipcells(void)
{
	gpu_and_managed(ngbrFlags, ~NF_SKIPCELL, get_ngbrFlags_size());
}

template <typename VType>
__host__ void cuVEC_VC<VType>::set_robin_conditions(cuReal2 robin_v_, cuReal2 robin_px_, cuReal2 robin_nx_, cuReal2 robin_py_, cuReal2 robin_ny_, cuReal2 robin_pz_, cuReal2 robin_nz_)
{
	set_robin(robin_px_, NF_ROBINPX);
	set_robin(robin_nx_, NF_ROBINNX);
	set_robin(robin_py_, NF_ROBINPY);
	set_robin(robin_ny_, NF_ROBINNY);
	set_robin(robin_pz_, NF_ROBINPZ);
	set_robin(robin_nz_, NF_ROBINNZ);
	set_robin(robin_v_, NF_ROBINV);

	//set flags
	set_robin_flags();
}

//clear all Robin boundary conditions and values
template <typename VType>
__host__ void cuVEC_VC<VType>::clear_robin_conditions(void)
{
	set_robin(cuReal2(), NF_ROBINPX);
	set_robin(cuReal2(), NF_ROBINNX);
	set_robin(cuReal2(), NF_ROBINPY);
	set_robin(cuReal2(), NF_ROBINNY);
	set_robin(cuReal2(), NF_ROBINPZ);
	set_robin(cuReal2(), NF_ROBINNZ);
	set_robin(cuReal2(), NF_ROBINV);

	//clear Robin flags
	gpu_and_managed(ngbrFlags, ~NF_ROBIN, get_ngbrFlags_size());
}