#pragma once

#include "cuVEC_VC.h"

//--------------------------------------------IMPORTANT FLAG MANIPULATION METHODS

//check if we need to use ngbrFlags2 (allocate memory etc.)
template <typename VType>
__host__ bool cuVEC_VC<VType>::use_extended_flags(void)
{
	//currently ngbrFlags2 used for:

	//ROBIN
	//DIRICHLET
	//HALO

	bool using_dirichlet = (
		get_dirichlet_size(NF2_DIRICHLETNX) 
		|| get_dirichlet_size(NF2_DIRICHLETPX) 
		|| get_dirichlet_size(NF2_DIRICHLETNY) 
		|| get_dirichlet_size(NF2_DIRICHLETPY) 
		|| get_dirichlet_size(NF2_DIRICHLETNZ) 
		|| get_dirichlet_size(NF2_DIRICHLETPZ));

	bool using_halo = (
		get_halo_size(NF2_HALOX)
		|| get_halo_size(NF2_HALOY)
		|| get_halo_size(NF2_HALOZ));

	bool using_robin = (
		get_gpu_value(robin_px) != cuReal2() 
		|| get_gpu_value(robin_nx) != cuReal2()
		|| get_gpu_value(robin_py) != cuReal2()
		|| get_gpu_value(robin_ny) != cuReal2()
		|| get_gpu_value(robin_pz) != cuReal2()
		|| get_gpu_value(robin_nz) != cuReal2()
		|| get_gpu_value(robin_v) != cuReal2());

	bool using_extended_flags_cpu = (using_dirichlet || using_robin || using_halo);

	//make sure ngbrFlags2 has the correct memory allocated only if currently empty
	if (using_extended_flags_cpu && isnullgpuptr(ngbrFlags2)) {

		cudaError_t error = gpu_alloc_managed(ngbrFlags2, get_ngbrFlags_size());
		if (error != cudaSuccess) return false;

		gpu_set_managed(ngbrFlags2, (int)0, get_ngbrFlags_size());
	}

	set_gpu_value(using_extended_flags, using_extended_flags_cpu);

	return using_extended_flags_cpu;
}

//from NF2_DIRICHLET type flag and cell_idx return boundary value from one of the dirichlet vectors
template <typename VType>
__device__ VType cuVEC_VC<VType>::get_dirichlet_value(int dirichlet_flag, int cell_idx) const
{
	switch (dirichlet_flag) {

	case NF2_DIRICHLETPX:
		return dirichlet_nx[((cell_idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (cell_idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y];

	case NF2_DIRICHLETNX:
		return dirichlet_px[((cell_idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y) + (cell_idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.y];

	case NF2_DIRICHLETPY:
		return dirichlet_ny[(cell_idx % cuVEC<VType>::n.x) + (cell_idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x];

	case NF2_DIRICHLETNY:
		return dirichlet_py[(cell_idx % cuVEC<VType>::n.x) + (cell_idx / (cuVEC<VType>::n.x*cuVEC<VType>::n.y))*cuVEC<VType>::n.x];

	case NF2_DIRICHLETPZ:
		return dirichlet_nz[(cell_idx % cuVEC<VType>::n.x) + ((cell_idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x];

	case NF2_DIRICHLETNZ:
		return dirichlet_pz[(cell_idx % cuVEC<VType>::n.x) + ((cell_idx / cuVEC<VType>::n.x) % cuVEC<VType>::n.y)*cuVEC<VType>::n.x];
	}

	return VType();
}

template <typename VType>
__device__ VType cuVEC_VC<VType>::get_dirichlet_value(int dirichlet_flag, const cuINT3& ijk) const
{
	switch (dirichlet_flag) {

	case NF2_DIRICHLETPX:
		return dirichlet_nx[ijk.j + ijk.k * cuVEC<VType>::n.y];

	case NF2_DIRICHLETNX:
		return dirichlet_px[ijk.j + ijk.k * cuVEC<VType>::n.y];

	case NF2_DIRICHLETPY:
		return dirichlet_ny[ijk.i + ijk.k * cuVEC<VType>::n.x];

	case NF2_DIRICHLETNY:
		return dirichlet_py[ijk.i + ijk.k * cuVEC<VType>::n.x];

	case NF2_DIRICHLETPZ:
		return dirichlet_nz[ijk.i + ijk.j * cuVEC<VType>::n.x];

	case NF2_DIRICHLETNZ:
		return dirichlet_pz[ijk.i + ijk.j * cuVEC<VType>::n.x];
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
	cuBox cells = cuVEC<VType>::box_from_rect_max(rectangle);

	for (int i = (cells.s.x >= 0 ? cells.s.x : 0); i < (cells.e.x <= cuVEC<VType>::n.x ? cells.e.x : cuVEC<VType>::n.x); i++) {

		for (int j = (cells.s.y >= 0 ? cells.s.y : 0); j < (cells.e.y <= cuVEC<VType>::n.y ? cells.e.y : cuVEC<VType>::n.y); j++) {
			
			for (int k = (cells.s.z >= 0 ? cells.s.z : 0); k < (cells.e.z <= cuVEC<VType>::n.z ? cells.e.z : cuVEC<VType>::n.z); k++) {

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
	//DIRICHLET flags used with extended ngbrFlags only.
	if (use_extended_flags()) {

		gpu_and_managed(ngbrFlags2, ~NF2_DIRICHLET, get_ngbrFlags2_size());

		set_dirichlet_size(0, NF2_DIRICHLETPX);
		set_dirichlet_size(0, NF2_DIRICHLETNX);
		set_dirichlet_size(0, NF2_DIRICHLETPY);
		set_dirichlet_size(0, NF2_DIRICHLETNY);
		set_dirichlet_size(0, NF2_DIRICHLETPZ);
		set_dirichlet_size(0, NF2_DIRICHLETNZ);

		//clear memory for extended ngbrFlags if now not used for anything else
		if (!use_extended_flags()) {

			gpu_free_managed(ngbrFlags2);
			set_gpu_value(using_extended_flags, false);
		}
	}
	else {

		gpu_free_managed(ngbrFlags2);
		set_gpu_value(using_extended_flags, false);
	}
}

//clear all halo flags and vectors
template <typename VType>
__host__ void cuVEC_VC<VType>::clear_halo_flags(void)
{
	//HALO flags used with extended ngbrFlags only.
	if (use_extended_flags()) {

		gpu_and_managed(ngbrFlags2, ~NF2_HALO, get_ngbrFlags2_size());

		set_halo_size(0, NF2_HALOX);
		set_halo_size(0, NF2_HALOY);
		set_halo_size(0, NF2_HALOZ);

		//clear memory for extended ngbrFlags if now not used for anything else
		if (!use_extended_flags()) {

			gpu_free_managed(ngbrFlags2);
			set_gpu_value(using_extended_flags, false);
		}
	}
	else {

		gpu_free_managed(ngbrFlags2);
		set_gpu_value(using_extended_flags, false);
	}
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
__host__ void cuVEC_VC<VType>::set_pbc(int pbc_x_, int pbc_y_, int pbc_z_)
{
	set_gpu_value(pbc_x, pbc_x_);
	set_gpu_value(pbc_y, pbc_y_);
	set_gpu_value(pbc_z, pbc_z_);

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
	set_robin(robin_px_, NF2_ROBINPX);
	set_robin(robin_nx_, NF2_ROBINNX);
	set_robin(robin_py_, NF2_ROBINPY);
	set_robin(robin_ny_, NF2_ROBINNY);
	set_robin(robin_pz_, NF2_ROBINPZ);
	set_robin(robin_nz_, NF2_ROBINNZ);
	set_robin(robin_v_, NF2_ROBINV);

	//set flags
	set_robin_flags();
}

//clear all Robin boundary conditions and values
template <typename VType>
__host__ void cuVEC_VC<VType>::clear_robin_conditions(void)
{
	set_robin(cuReal2(), NF2_ROBINPX);
	set_robin(cuReal2(), NF2_ROBINNX);
	set_robin(cuReal2(), NF2_ROBINPY);
	set_robin(cuReal2(), NF2_ROBINNY);
	set_robin(cuReal2(), NF2_ROBINPZ);
	set_robin(cuReal2(), NF2_ROBINNZ);
	set_robin(cuReal2(), NF2_ROBINV);

	if (use_extended_flags()) {

		//clear Robin flags
		gpu_and_managed(ngbrFlags, ~NF2_ROBIN, get_ngbrFlags2_size());
	}
	else {

		gpu_free_managed(ngbrFlags2);
		set_gpu_value(using_extended_flags, false);
	}
}