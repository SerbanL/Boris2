#pragma once

#include <vector>

#include "Funcs_Vectors.h"

#include "cuVEC_VC.h"

//--------------------------------------------MEMORY MANAGEMENT HELPER METHODS

template <typename VType>
//memory allocation for objects and initialize to default - only call at start of managed constructors
__host__ void cuVEC_VC<VType>::alloc_initialize_data(void)
{
	nullgpuptr(ngbrFlags);
	nullgpuptr(ngbrFlags2);
	nullgpuptr(dirichlet_px);
	nullgpuptr(dirichlet_nx);
	nullgpuptr(dirichlet_py);
	nullgpuptr(dirichlet_ny);
	nullgpuptr(dirichlet_pz);
	nullgpuptr(dirichlet_nz);

	nullgpuptr(halo_px);
	nullgpuptr(halo_nx);
	nullgpuptr(halo_py);
	nullgpuptr(halo_ny);
	nullgpuptr(halo_pz);
	nullgpuptr(halo_nz);

	set_ngbrFlags_size(0);
	set_gpu_value(using_extended_flags, false);
	set_gpu_value(nonempty_cells, (int)0);

	set_dirichlet_size(0, NF2_DIRICHLETPX);
	set_dirichlet_size(0, NF2_DIRICHLETNX);
	set_dirichlet_size(0, NF2_DIRICHLETPY);
	set_dirichlet_size(0, NF2_DIRICHLETNY);
	set_dirichlet_size(0, NF2_DIRICHLETPZ);
	set_dirichlet_size(0, NF2_DIRICHLETNZ);

	set_halo_size(0, NF2_HALOX);
	set_halo_size(0, NF2_HALOY);
	set_halo_size(0, NF2_HALOZ);

	set_robin(cuReal2(), NF2_ROBINPX);
	set_robin(cuReal2(), NF2_ROBINNX);
	set_robin(cuReal2(), NF2_ROBINPY);
	set_robin(cuReal2(), NF2_ROBINNY);
	set_robin(cuReal2(), NF2_ROBINPZ);
	set_robin(cuReal2(), NF2_ROBINNZ);
	set_robin(cuReal2(), NF2_ROBINV);

	set_gpu_value(shift_debt, cuReal3());

	set_gpu_value(pbc_x, (int)0);
	set_gpu_value(pbc_y, (int)0);
	set_gpu_value(pbc_z, (int)0);
}

//set size of ngbrFlags, allocating memory
template <typename VType>
__host__ cudaError_t cuVEC_VC<VType>::set_ngbrFlags_size(size_t size)
{
	cudaError_t error = cudaSuccess;

	//allocate memory for ngbrFlags (if size is 0 then allocate 1 element as 0 will fail; Still keep the real intended size in ngbrFlags_size though).
	error = gpu_alloc_managed(ngbrFlags, size);
	if (error != cudaSuccess) {

		//failed : nullptr
		nullgpuptr(ngbrFlags);
		size = 0;
	}
		
	//store ngbrFlags size
	set_gpu_value(ngbrFlags_size, size);

	return error;
}

//get ngbrFlags_size value in cpu memory
template <typename VType>
__host__ size_t cuVEC_VC<VType>::get_ngbrFlags_size(void) const
{
	return get_gpu_value(ngbrFlags_size);
}

//get ngbrFlags_size value in cpu memory
template <typename VType>
__host__ size_t cuVEC_VC<VType>::get_ngbrFlags2_size(void) const
{
	if (isnullgpuptr((int*&)ngbrFlags2)) return 0;
	else return get_gpu_value(ngbrFlags_size);
}

//set robin values using a robin id
template <typename VType>
__host__ void cuVEC_VC<VType>::set_robin(cuReal2 robin_value, int robin_id)
{
	switch (robin_id) {

	case NF2_ROBINPX:
		set_gpu_value(robin_px, robin_value);
		break;
		
	case NF2_ROBINNX:
		set_gpu_value(robin_nx, robin_value);
		break;
		
	case NF2_ROBINPY:
		set_gpu_value(robin_py, robin_value);
		break;
		
	case NF2_ROBINNY:
		set_gpu_value(robin_ny, robin_value);
		break;
		
	case NF2_ROBINPZ:
		set_gpu_value(robin_pz, robin_value);
		break;
		
	case NF2_ROBINNZ:
		set_gpu_value(robin_nz, robin_value);
		break;
		
	case NF2_ROBINV:
		set_gpu_value(robin_v, robin_value);
		break;
	}
}

//set dirichlet array size, allocating memory as required
template <typename VType>
__host__ cudaError_t cuVEC_VC<VType>::set_dirichlet_size(size_t size, int dirichlet_id)
{
	//same pattern as in set_ngbrFlags_size method so read comments there
	cudaError_t error = cudaSuccess;
	//do not reallocate if array already has requested size
	if (get_dirichlet_size(dirichlet_id) == size) return error;

	switch (dirichlet_id) {

	case NF2_DIRICHLETPX:
		if (size > 0) error = gpu_alloc_managed(dirichlet_px, size);
		else error = gpu_free_managed(dirichlet_px);
		if (error != cudaSuccess) {

			nullgpuptr(dirichlet_px);
			size = 0;
		}
		set_gpu_value(dirichlet_px_size, size);
		break;
		
	case NF2_DIRICHLETNX:
		if (size > 0) error = gpu_alloc_managed(dirichlet_nx, size);
		else error = gpu_free_managed(dirichlet_nx);
		if (error != cudaSuccess) {

			nullgpuptr(dirichlet_nx);
			size = 0;
		}
		set_gpu_value(dirichlet_nx_size, size);
		break;
		
	case NF2_DIRICHLETPY:
		if (size > 0) error = gpu_alloc_managed(dirichlet_py, size);
		else error = gpu_free_managed(dirichlet_py);
		if (error != cudaSuccess) {

			nullgpuptr(dirichlet_py);
			size = 0;
		}
		set_gpu_value(dirichlet_py_size, size);
		break;
		
	case NF2_DIRICHLETNY:
		if (size > 0) error = gpu_alloc_managed(dirichlet_ny, size);
		else error = gpu_free_managed(dirichlet_ny);
		if (error != cudaSuccess) {

			nullgpuptr(dirichlet_ny);
			size = 0;
		}
		set_gpu_value(dirichlet_ny_size, size);
		break;
		
	case NF2_DIRICHLETPZ:
		if (size > 0) error = gpu_alloc_managed(dirichlet_pz, size);
		else error = gpu_free_managed(dirichlet_pz);
		if (error != cudaSuccess) {

			nullgpuptr(dirichlet_pz);
			size = 0;
		}
		set_gpu_value(dirichlet_pz_size, size);
		break;
		
	case NF2_DIRICHLETNZ:
		if (size > 0) error = gpu_alloc_managed(dirichlet_nz, size);
		else error = gpu_free_managed(dirichlet_nz);
		if (error != cudaSuccess) {

			nullgpuptr(dirichlet_nz);
			size = 0;
		}
		set_gpu_value(dirichlet_nz_size, size);
		break;
	}

	return error;
}

//get dirichlet array size
template <typename VType>
__host__ size_t cuVEC_VC<VType>::get_dirichlet_size(int dirichlet_id) const
{
	switch (dirichlet_id) {

	case NF2_DIRICHLETPX:
		return get_gpu_value(dirichlet_px_size);
	
	case NF2_DIRICHLETNX:
		return get_gpu_value(dirichlet_nx_size);
	
	case NF2_DIRICHLETPY:
		return get_gpu_value(dirichlet_py_size);
	
	case NF2_DIRICHLETNY:
		return get_gpu_value(dirichlet_ny_size);
	
	case NF2_DIRICHLETPZ:
		return get_gpu_value(dirichlet_pz_size);
	
	case NF2_DIRICHLETNZ:
		return get_gpu_value(dirichlet_nz_size);
	}

	return 0;
}

//set halo array size, allocating memory as required
//halo_id is one of NF2_HALOX, NF2_HALOY, NF2_HALOZ, since both n and p halos are allocated for a given axis.
//also accepts individual halo flags for halo_id
template <typename VType>
__host__ cudaError_t cuVEC_VC<VType>::set_halo_size(size_t size, int halo_id)
{
	//same pattern as in set_ngbrFlags_size method so read comments there
	cudaError_t error = cudaSuccess;
	
	//do not reallocate if array already has requested size
	if (get_halo_size(halo_id) == size) return error;
	
	switch (halo_id) {

	case NF2_HALOPX:
	case NF2_HALONX:
	case NF2_HALOX:
		if (size > 0) {

			error = gpu_alloc_managed(halo_px, size);
			error = gpu_alloc_managed(halo_nx, size);
		}
		else {

			error = gpu_free_managed(halo_px);
			error = gpu_free_managed(halo_nx);
		}
		if (error != cudaSuccess) {

			nullgpuptr(halo_px);
			nullgpuptr(halo_nx);
			size = 0;
		}
		set_gpu_value(halo_x_size, size);
		break;

	case NF2_HALOPY:
	case NF2_HALONY:
	case NF2_HALOY:
		if (size > 0) {

			error = gpu_alloc_managed(halo_py, size);
			error = gpu_alloc_managed(halo_ny, size);
		}
		else {

			error = gpu_free_managed(halo_py);
			error = gpu_free_managed(halo_ny);
		}
		if (error != cudaSuccess) {

			nullgpuptr(halo_py);
			nullgpuptr(halo_ny);
			size = 0;
		}
		set_gpu_value(halo_y_size, size);
		break;

	case NF2_HALOPZ:
	case NF2_HALONZ:
	case NF2_HALOZ:
		if (size > 0) {

			error = gpu_alloc_managed(halo_pz, size);
			error = gpu_alloc_managed(halo_nz, size);
		}
		else {

			error = gpu_free_managed(halo_pz);
			error = gpu_free_managed(halo_nz);
		}
		if (error != cudaSuccess) {

			nullgpuptr(halo_pz);
			nullgpuptr(halo_nz);
			size = 0;
		}
		set_gpu_value(halo_z_size, size);
		break;
	}

	return error;
}

//get halo array size
//halo_id is one of NF2_HALOX, NF2_HALOY, NF2_HALOZ, since both n and p halos are allocated for a given axis.
//also accepts individual halo flags for halo_id
template <typename VType>
__host__ size_t cuVEC_VC<VType>::get_halo_size(int halo_id) const
{
	switch (halo_id) {

	case NF2_HALOPX:
	case NF2_HALONX:
	case NF2_HALOX:
		return get_gpu_value(halo_x_size);

	case NF2_HALOPY:
	case NF2_HALONY:
	case NF2_HALOY:
		return get_gpu_value(halo_y_size);

	case NF2_HALOPZ:
	case NF2_HALONZ:
	case NF2_HALOZ:
		return get_gpu_value(halo_z_size);
	}

	return 0;
}

//--------------------------------------------CONSTRUCTORS : cu_obj "managed constructors" only. Real constructors are never called since you should never make a real instance of a cuVEC_VC.

//void constructor
template <typename VType>
__host__ void cuVEC_VC<VType>::construct_cu_obj(void)
{
	cuVEC<VType>::construct_cu_obj();

	alloc_initialize_data();
}

//construct to given number of cells : n_ is in cpu memory
template <typename VType>
__host__ void cuVEC_VC<VType>::construct_cu_obj(const cuSZ3& n_)
{
	cuVEC<VType>::construct_cu_obj(n_);

	alloc_initialize_data();

	if (get_gpu_value(cuVEC<VType>::n).dim()) {

		cudaError_t error = set_ngbrFlags_size(get_gpu_value(cuVEC<VType>::n).dim());
		if (error != cudaSuccess) cuVEC<VType>::set_n(cuSZ3());
		else gpu_set_managed(ngbrFlags, NF_NOTEMPTY, get_gpu_value(cuVEC<VType>::n).dim());
	}
}

//construct to given dimensions : h_ and rect_ are in cpu memory
template <typename VType>
__host__ void cuVEC_VC<VType>::construct_cu_obj(const cuReal3& h_, const cuRect& rect_)
{
	cuVEC<VType>::construct_cu_obj(h_, rect_);

	alloc_initialize_data();

	if (get_gpu_value(cuVEC<VType>::n).dim()) {

		cudaError_t error = set_ngbrFlags_size(get_gpu_value(cuVEC<VType>::n).dim());
		if (error != cudaSuccess) cuVEC<VType>::set_n(cuSZ3());
		else {

			gpu_set_managed(ngbrFlags, NF_NOTEMPTY, get_gpu_value(cuVEC<VType>::n).dim());
			set_ngbrFlags();
		}
	}
}

//construct to given dimensions and initialize to given value : h_ and rect_ are in cpu memory
template <typename VType>
__host__ void cuVEC_VC<VType>::construct_cu_obj(const cuReal3& h_, const cuRect& rect_, VType value)
{
	cuVEC<VType>::construct_cu_obj(h_, rect_, value);

	alloc_initialize_data();

	if (get_gpu_value(cuVEC<VType>::n).dim()) {

		cudaError_t error = set_ngbrFlags_size(get_gpu_value(cuVEC<VType>::n).dim());
		if (error != cudaSuccess) cuVEC<VType>::set_n(cuSZ3());
		else {

			gpu_set_managed(ngbrFlags, NF_NOTEMPTY, get_gpu_value(cuVEC<VType>::n).dim());
			set_ngbrFlags();
		}
	}
}

//copy constructor
template <typename VType>
__host__ void cuVEC_VC<VType>::construct_cu_obj(const cuVEC_VC<VType>& copyThis)
{
	cuVEC<VType>::construct_cu_obj();
	
	alloc_initialize_data();

	assign_cu_obj(copyThis);
}

//assignment operator
template <typename VType>
__host__ void cuVEC_VC<VType>::assign_cu_obj(const cuVEC_VC<VType>& copyThis)
{
	//copy n
	gpu_to_gpu(cuVEC<VType>::n, copyThis.n);
	//copy h
	gpu_to_gpu(cuVEC<VType>::h, copyThis.h);
	//copy rect
	gpu_to_gpu(cuVEC<VType>::rect, copyThis.rect);

	//allocate memory and copy array
	cudaError_t error = allocate_quantity(get_gpu_value(cuVEC<VType>::n));
	if (error == cudaSuccess) {

		gpu_to_gpu_managed(cuVEC<VType>::quantity, copyThis.quantity, get_gpu_value(cuVEC<VType>::n).dim());
	}
	else {

		cuVEC<VType>::set_h(cuReal3());
		cuVEC<VType>::set_rect(cuRect());
	}

	//transfer info might not be valid any more
	cuVEC<VType>::transfer.clear_transfer_data();

	//copy ngbrFlags_size
	gpu_to_gpu(ngbrFlags_size, copyThis.ngbrFlags_size);

	//allocate memory and copy ngbrFlags
	error = set_ngbrFlags_size(get_gpu_value(cuVEC<VType>::n).dim());
	if (error == cudaSuccess) {

		gpu_to_gpu_managed(ngbrFlags, copyThis.ngbrFlags, get_ngbrFlags_size());
	}

	//copy ngbrFlags2 if needed (else free and null it)
	if (isnullgpuptr((int*&)copyThis.ngbrFlags2)) {

		gpu_free_managed(ngbrFlags2);
		set_gpu_value(using_extended_flags, false);
	}
	else {

		error = gpu_alloc_managed(ngbrFlags2, get_ngbrFlags_size());
		if (error != cudaSuccess) {

			//failed : nullptr
			nullgpuptr(ngbrFlags2);
			set_gpu_value(using_extended_flags, false);
		}
		else {

			gpu_to_gpu_managed(ngbrFlags2, copyThis.ngbrFlags2, get_ngbrFlags_size());
			set_gpu_value(using_extended_flags, true);
		}
	}

	//copy nonempty_cells
	gpu_to_gpu(nonempty_cells, copyThis.nonempty_cells);

	//copy dirichlet sizes
	gpu_to_gpu(dirichlet_px_size, copyThis.dirichlet_px_size);
	gpu_to_gpu(dirichlet_nx_size, copyThis.dirichlet_nx_size);
	gpu_to_gpu(dirichlet_py_size, copyThis.dirichlet_py_size);
	gpu_to_gpu(dirichlet_ny_size, copyThis.dirichlet_ny_size);
	gpu_to_gpu(dirichlet_pz_size, copyThis.dirichlet_pz_size);
	gpu_to_gpu(dirichlet_nz_size, copyThis.dirichlet_nz_size);

	//copy halo sizes
	gpu_to_gpu(halo_x_size, copyThis.halo_x_size);
	gpu_to_gpu(halo_y_size, copyThis.halo_y_size);
	gpu_to_gpu(halo_z_size, copyThis.halo_z_size);

	//copy robin values
	gpu_to_gpu(robin_px, copyThis.robin_px);
	gpu_to_gpu(robin_nx, copyThis.robin_nx);
	gpu_to_gpu(robin_py, copyThis.robin_py);
	gpu_to_gpu(robin_ny, copyThis.robin_ny);
	gpu_to_gpu(robin_pz, copyThis.robin_pz);
	gpu_to_gpu(robin_nz, copyThis.robin_nz);
	gpu_to_gpu(robin_v, copyThis.robin_v);

	//copy shift debt
	gpu_to_gpu(shift_debt, copyThis.shift_debt);

	//copy dirichlet values arrays
	size_t size;
	
	size = get_dirichlet_size(NF2_DIRICHLETPX);
	if (size) {

		error = set_dirichlet_size(size, NF2_DIRICHLETPX);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(dirichlet_px, copyThis.dirichlet_px, size);
		}
	}

	size = get_dirichlet_size(NF2_DIRICHLETNX);
	if (size) {

		error = set_dirichlet_size(size, NF2_DIRICHLETNX);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(dirichlet_nx, copyThis.dirichlet_nx, size);
		}
	}

	size = get_dirichlet_size(NF2_DIRICHLETPY);
	if (size) {

		error = set_dirichlet_size(size, NF2_DIRICHLETPY);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(dirichlet_py, copyThis.dirichlet_py, size);
		}
	}

	size = get_dirichlet_size(NF2_DIRICHLETNY);
	if (size) {

		error = set_dirichlet_size(size, NF2_DIRICHLETNY);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(dirichlet_ny, copyThis.dirichlet_ny, size);
		}
	}

	size = get_dirichlet_size(NF2_DIRICHLETPZ);
	if (size) {

		error = set_dirichlet_size(size, NF2_DIRICHLETPZ);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(dirichlet_pz, copyThis.dirichlet_pz, size);
		}
	}

	size = get_dirichlet_size(NF2_DIRICHLETNZ);
	if (size) {

		error = set_dirichlet_size(size, NF2_DIRICHLETNZ);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(dirichlet_nz, copyThis.dirichlet_nz, size);
		}
	}

	//copy halo values arrays

	size = get_halo_size(NF2_HALOX);
	if (size) {

		error = set_halo_size(size, NF2_HALOX);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(halo_px, copyThis.halo_px, size);
			gpu_to_gpu_managed(halo_nx, copyThis.halo_nx, size);
		}
	}

	size = get_halo_size(NF2_HALOY);
	if (size) {

		error = set_halo_size(size, NF2_HALOY);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(halo_py, copyThis.halo_py, size);
			gpu_to_gpu_managed(halo_ny, copyThis.halo_ny, size);
		}
	}

	size = get_halo_size(NF2_HALOZ);
	if (size) {

		error = set_halo_size(size, NF2_HALOZ);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(halo_pz, copyThis.halo_pz, size);
			gpu_to_gpu_managed(halo_nz, copyThis.halo_nz, size);
		}
	}

	//copy pbc settings
	gpu_to_gpu(pbc_x, copyThis.pbc_x);
	gpu_to_gpu(pbc_y, copyThis.pbc_y);
	gpu_to_gpu(pbc_z, copyThis.pbc_z);
}

//destructor
template <typename VType>
__host__ void cuVEC_VC<VType>::destruct_cu_obj(void)
{
	gpu_free_managed(ngbrFlags);
	
	gpu_free_managed(ngbrFlags2);
	set_gpu_value(using_extended_flags, false);
	
	gpu_free_managed(dirichlet_px);
	gpu_free_managed(dirichlet_nx);
	gpu_free_managed(dirichlet_py);
	gpu_free_managed(dirichlet_ny);
	gpu_free_managed(dirichlet_pz);
	gpu_free_managed(dirichlet_nz);

	gpu_free_managed(halo_px);
	gpu_free_managed(halo_nx);
	gpu_free_managed(halo_py);
	gpu_free_managed(halo_ny);
	gpu_free_managed(halo_pz);
	gpu_free_managed(halo_nz);

	cuVEC<VType>::destruct_cu_obj();
}

//--------------------------------------------COPY FROM VEC_VC

//copy everything from a VEC_VC - type must be convertible. Return false if failed (memory could not be allocated)
template <typename VType>
template <typename cpuVEC_VC>
__host__ bool cuVEC_VC<VType>::set_from_cpuvec(cpuVEC_VC& vec_vc)
{
	//-----------

	//things to copy:
	//n, h, rect
	//quantity
	//ngbrFlags (and set ngbrFlags_size)
	//nonempty_cells
	//dirichlet vectors (and set sizes)
	//robin values
	//shift debt

	//-----------

	cuSZ3 cpu_n = vec_vc.n;
	cuReal3 cpu_h = vec_vc.h;
	cuRect cpu_rect = vec_vc.rect;

	//-----------

	cuVEC<VType>::set_from_cpuvec(vec_vc);

	//----------- VEC part now done, onto the VEC_VC part

	cudaError_t error;
	
	//copy ngbrFlags
	if(get_ngbrFlags_size() != cpu_n.dim()) {

		//first allocate memory if sizes don't match
		error = gpu_alloc_managed(ngbrFlags, cpu_n.dim());
		if (error != cudaSuccess) return false; //not enough memory so abort

		//set ngbrFlags_size
		set_gpu_value(ngbrFlags_size, cpu_n.dim());
	}

	//now copy
	cpu_to_gpu_managed(ngbrFlags, vec_vc.ngbrFlags_ref().data(), cpu_n.dim());

	//copy ngbrFlags2
	if (vec_vc.ngbrFlags2_ref().size()) {

		//first allocate memory
		error = gpu_alloc_managed(ngbrFlags2, cpu_n.dim());
		if (error != cudaSuccess) return false; //not enough memory so abort

		set_gpu_value(using_extended_flags, true);

		//now copy
		cpu_to_gpu_managed(ngbrFlags2, vec_vc.ngbrFlags2_ref().data(), cpu_n.dim());

	}
	else {

		//ngbrFlags2 not used
		gpu_free_managed(ngbrFlags2);
		set_gpu_value(using_extended_flags, false);
	}

	//-----------

	//copy nonempty_cells
	set_gpu_value(nonempty_cells, vec_vc.get_nonempty_cells());

	//-----------

	//copy dirichlet vectors
	
	size_t dirichlet_size;

	//---NF2_DIRICHLETPX
	
	//allocate gpu memory if sizes don't match
	dirichlet_size = vec_vc.dirichlet_px_ref().size();
	if (get_dirichlet_size(NF2_DIRICHLETPX) != dirichlet_size) {

		error = gpu_alloc_managed(dirichlet_px, dirichlet_size);
		if (error != cudaSuccess) return false; //not enough memory so abort

		//set size in gpu memory
		set_gpu_value(dirichlet_px_size, dirichlet_size);
	}

	//copy over
	cpu_to_gpu_managed(dirichlet_px, vec_vc.dirichlet_px_ref().data(), dirichlet_size);

	//---NF2_DIRICHLETNX

	//allocate gpu memory if sizes don't match
	dirichlet_size = vec_vc.dirichlet_nx_ref().size();
	if (get_dirichlet_size(NF2_DIRICHLETNX) != dirichlet_size) {

		error = gpu_alloc_managed(dirichlet_nx, dirichlet_size);
		if (error != cudaSuccess) return false; //not enough memory so abort

		//set size in gpu memory
		set_gpu_value(dirichlet_nx_size, dirichlet_size);
	}

	//copy over
	cpu_to_gpu_managed(dirichlet_nx, vec_vc.dirichlet_nx_ref().data(), dirichlet_size);
	
	//---NF2_DIRICHLETPY

	//allocate gpu memory if sizes don't match
	dirichlet_size = vec_vc.dirichlet_py_ref().size();
	if (get_dirichlet_size(NF2_DIRICHLETPY) != dirichlet_size) {

		error = gpu_alloc_managed(dirichlet_py, dirichlet_size);
		if (error != cudaSuccess) return false; //not enough memory so abort

		//set size in gpu memory
		set_gpu_value(dirichlet_py_size, dirichlet_size);
	}

	//copy over
	cpu_to_gpu_managed(dirichlet_py, vec_vc.dirichlet_py_ref().data(), dirichlet_size);

	//---NF2_DIRICHLETNY

	//allocate gpu memory if sizes don't match
	dirichlet_size = vec_vc.dirichlet_ny_ref().size();
	if (get_dirichlet_size(NF2_DIRICHLETNY) != dirichlet_size) {

		error = gpu_alloc_managed(dirichlet_ny, dirichlet_size);
		if (error != cudaSuccess) return false; //not enough memory so abort

		//set size in gpu memory
		set_gpu_value(dirichlet_ny_size, dirichlet_size);
	}

	//copy over
	cpu_to_gpu_managed(dirichlet_ny, vec_vc.dirichlet_ny_ref().data(), dirichlet_size);

	//---NF2_DIRICHLETPZ

	//allocate gpu memory if sizes don't match
	dirichlet_size = vec_vc.dirichlet_pz_ref().size();
	if (get_dirichlet_size(NF2_DIRICHLETPZ) != dirichlet_size) {

		error = gpu_alloc_managed(dirichlet_pz, dirichlet_size);
		if (error != cudaSuccess) return false; //not enough memory so abort

		//set size in gpu memory
		set_gpu_value(dirichlet_pz_size, dirichlet_size);
	}

	//copy over
	cpu_to_gpu_managed(dirichlet_pz, vec_vc.dirichlet_pz_ref().data(), dirichlet_size);

	//---NF2_DIRICHLETNZ

	//allocate gpu memory if sizes don't match
	dirichlet_size = vec_vc.dirichlet_nz_ref().size();
	if (get_dirichlet_size(NF2_DIRICHLETNZ) != dirichlet_size) {

		error = gpu_alloc_managed(dirichlet_nz, dirichlet_size);
		if (error != cudaSuccess) return false; //not enough memory so abort

		//set size in gpu memory
		set_gpu_value(dirichlet_nz_size, dirichlet_size);
	}

	//copy over
	cpu_to_gpu_managed(dirichlet_nz, vec_vc.dirichlet_nz_ref().data(), dirichlet_size);
	
	//-----------

	//halo not defined for VEC_VC
	//don't do anything about halos here, these would be managed at policy class level for multi-device management

	//-----------

	//copy robin values
	set_gpu_value(robin_px, (cuReal2)vec_vc.robin_px_ref());
	set_gpu_value(robin_nx, (cuReal2)vec_vc.robin_nx_ref());
	set_gpu_value(robin_py, (cuReal2)vec_vc.robin_py_ref());
	set_gpu_value(robin_ny, (cuReal2)vec_vc.robin_ny_ref());
	set_gpu_value(robin_pz, (cuReal2)vec_vc.robin_pz_ref());
	set_gpu_value(robin_nz, (cuReal2)vec_vc.robin_nz_ref());
	set_gpu_value(robin_v, (cuReal2)vec_vc.robin_v_ref());

	//copy shift debt
	set_gpu_value(shift_debt, (cuReal3)vec_vc.shift_debt_ref());

	//copy pbc parameters
	set_gpu_value(pbc_x, (int)vec_vc.pbc_x_ref());
	set_gpu_value(pbc_y, (int)vec_vc.pbc_y_ref());
	set_gpu_value(pbc_z, (int)vec_vc.pbc_z_ref());

	return true;
}

//--------------------------------------------COPY TO VEC_VC

//copy everything to a VEC_VC - type must be convertible. Return false if failed (memory could not be allocated)
template <typename VType>
template <typename cpuVEC_VC>
__host__ bool cuVEC_VC<VType>::set_cpuvec(cpuVEC_VC& vec_vc)
{
	//things to copy:
	//n, h, rect
	//quantity
	//ngbrFlags
	//nonempty_cells
	//dirichlet vectors
	//robin values
	//shift debt

	//----------- copy to VEC part (n, h, rect, quantity)

	cuVEC<VType>::set_cpuvec(vec_vc);

	//-----------

	SZ3 cpu_n = get_gpu_value(cuVEC<VType>::n);
	DBL3 cpu_h = get_gpu_value(cuVEC<VType>::h);
	Rect cpu_rect = get_gpu_value(cuVEC<VType>::rect);

	//-----------

	//copy ngbrFlags
	if (vec_vc.ngbrFlags_ref().size() != cpu_n.dim()) {

		//first allocate memory if sizes don't match
		if (!malloc_vector(vec_vc.ngbrFlags_ref(), cpu_n.dim())) return false;		//couldn't allocate memory
	}

	//now copy
	gpu_to_cpu_managed(vec_vc.ngbrFlags_ref().data(), ngbrFlags, cpu_n.dim());

	//copy ngbrFlags2 if needed
	if (!isnullgpuptr(ngbrFlags2)) {

		//first allocate memory
		if (!malloc_vector(vec_vc.ngbrFlags2_ref(), cpu_n.dim())) return false;		//couldn't allocate memory

		//now copy
		gpu_to_cpu_managed(vec_vc.ngbrFlags2_ref().data(), ngbrFlags2, cpu_n.dim());
	}
	else {
		
		//ngbrFlags2 not used
		vec_vc.ngbrFlags2_ref().clear();
		vec_vc.ngbrFlags2_ref().shrink_to_fit();
	}

	//-----------

	//copy nonempty_cells
	vec_vc.nonempty_cells_ref() = get_gpu_value(nonempty_cells);

	//-----------
	
	//Set dirichlet vectors
	
	size_t dirichlet_size;

	//---NF2_DIRICHLETPX

	dirichlet_size = get_dirichlet_size(NF2_DIRICHLETPX);

	//copy dirichlet array to cpu memory
	if (!malloc_vector(vec_vc.dirichlet_px_ref(), dirichlet_size)) return false;		//couldn't allocate memory	
	gpu_to_cpu_managed(vec_vc.dirichlet_px_ref().data(), dirichlet_px, dirichlet_size);

	//---NF2_DIRICHLETNX

	dirichlet_size = get_dirichlet_size(NF2_DIRICHLETNX);

	//copy dirichlet array to cpu memory
	if (!malloc_vector(vec_vc.dirichlet_nx_ref(), dirichlet_size)) return false;		//couldn't allocate memory	
	gpu_to_cpu_managed(vec_vc.dirichlet_nx_ref().data(), dirichlet_nx, dirichlet_size);

	//---NF2_DIRICHLETPY

	dirichlet_size = get_dirichlet_size(NF2_DIRICHLETPY);

	//copy dirichlet array to cpu memory
	if (!malloc_vector(vec_vc.dirichlet_py_ref(), dirichlet_size)) return false;		//couldn't allocate memory	
	gpu_to_cpu_managed(vec_vc.dirichlet_py_ref().data(), dirichlet_py, dirichlet_size);

	//---NF2_DIRICHLETNY

	dirichlet_size = get_dirichlet_size(NF2_DIRICHLETNY);

	//copy dirichlet array to cpu memory
	if (!malloc_vector(vec_vc.dirichlet_ny_ref(), dirichlet_size)) return false;		//couldn't allocate memory	
	gpu_to_cpu_managed(vec_vc.dirichlet_ny_ref().data(), dirichlet_ny, dirichlet_size);

	//---NF2_DIRICHLETPZ

	dirichlet_size = get_dirichlet_size(NF2_DIRICHLETPZ);

	//copy dirichlet array to cpu memory
	if (!malloc_vector(vec_vc.dirichlet_pz_ref(), dirichlet_size)) return false;		//couldn't allocate memory	
	gpu_to_cpu_managed(vec_vc.dirichlet_pz_ref().data(), dirichlet_pz, dirichlet_size);

	//---NF2_DIRICHLETNZ

	dirichlet_size = get_dirichlet_size(NF2_DIRICHLETNZ);

	//copy dirichlet array to cpu memory
	if (!malloc_vector(vec_vc.dirichlet_nz_ref(), dirichlet_size)) return false;		//couldn't allocate memory	
	gpu_to_cpu_managed(vec_vc.dirichlet_nz_ref().data(), dirichlet_nz, dirichlet_size);

	//-----------

	//halo not defined for VEC_VC

	//-----------

	//copy robin values
	vec_vc.robin_px_ref() = get_gpu_value(robin_px);
	vec_vc.robin_nx_ref() = get_gpu_value(robin_nx);
	vec_vc.robin_py_ref() = get_gpu_value(robin_py);
	vec_vc.robin_ny_ref() = get_gpu_value(robin_ny);
	vec_vc.robin_pz_ref() = get_gpu_value(robin_pz);
	vec_vc.robin_nz_ref() = get_gpu_value(robin_nz);
	vec_vc.robin_v_ref() = get_gpu_value(robin_v);
	
	//copy shift debt
	vec_vc.shift_debt_ref() = get_gpu_value(shift_debt);

	//copy pbc parameters
	vec_vc.pbc_x_ref() = get_gpu_value(pbc_x);
	vec_vc.pbc_y_ref() = get_gpu_value(pbc_y);
	vec_vc.pbc_z_ref() = get_gpu_value(pbc_z);

	return true;
}

//faster version of set_from_cpuvec, where it is assumed the cpu vec already has the same sizes as this cuVEC_VC
template <typename VType>
template <typename cpuVEC_VC>
__host__ bool cuVEC_VC<VType>::copy_from_cpuvec(cpuVEC_VC& vec_vc)
{
	//copy quantity
	if (!cuVEC<VType>::copy_from_cpuvec(vec_vc)) return false;

	//copy shift debt
	set_gpu_value(shift_debt, (cuReal3)vec_vc.shift_debt_ref());

	//copy pbc parameters
	set_gpu_value(pbc_x, (int)vec_vc.pbc_x_ref());
	set_gpu_value(pbc_y, (int)vec_vc.pbc_y_ref());
	set_gpu_value(pbc_z, (int)vec_vc.pbc_z_ref());

	//copy ngbrFlags
	cudaError_t error = cpu_to_gpu_managed(ngbrFlags, vec_vc.ngbrFlags_ref().data(), vec_vc.linear_size());
	if (error != cudaSuccess) return false;
	
	//copy ngbrFlags2 if needed
	if (vec_vc.ngbrFlags2_ref().size()) {

		//first allocate memory
		error = gpu_alloc_managed(ngbrFlags2, vec_vc.linear_size());
		if (error != cudaSuccess) return false; //not enough memory so abort

		set_gpu_value(using_extended_flags, true);

		//now copy
		cpu_to_gpu_managed(ngbrFlags2, vec_vc.ngbrFlags2_ref().data(), vec_vc.linear_size());
	}
	else {

		//ngbrFlags2 not used
		gpu_free_managed(ngbrFlags2);
		set_gpu_value(using_extended_flags, false);
	}
	
	return true;
}

//faster version of set_cpuvec, where it is assumed the cpu vec already has the same sizes as this cuVEC_VC
template <typename VType>
template <typename cpuVEC_VC>
__host__ bool cuVEC_VC<VType>::copy_to_cpuvec(cpuVEC_VC& vec_vc)
{
	//copy quantity
	if (!cuVEC<VType>::copy_to_cpuvec(vec_vc)) return false;

	//copy shift debt
	vec_vc.shift_debt_ref() = get_gpu_value(shift_debt);

	//copy pbc parameters
	vec_vc.pbc_x_ref() = get_gpu_value(pbc_x);
	vec_vc.pbc_y_ref() = get_gpu_value(pbc_y);
	vec_vc.pbc_z_ref() = get_gpu_value(pbc_z);

	//now copy ngbrFlags
	cudaError_t error = gpu_to_cpu_managed(vec_vc.ngbrFlags_ref().data(), ngbrFlags, vec_vc.linear_size());
	if (error != cudaSuccess) return false;
	
	//copy ngbrFlags2 if needed
	if (!isnullgpuptr(ngbrFlags2)) {

		//first allocate memory
		if (!malloc_vector(vec_vc.ngbrFlags2_ref(), vec_vc.linear_size())) return false;		//couldn't allocate memory

		//now copy
		gpu_to_cpu_managed(vec_vc.ngbrFlags2_ref().data(), ngbrFlags2, vec_vc.linear_size());
	}
	else {

		//ngbrFlags2 not used
		vec_vc.ngbrFlags2_ref().clear();
		vec_vc.ngbrFlags2_ref().shrink_to_fit();
	}
	
	return true;
}

//copy flags only from vec_vc, where sizes must match
template <typename VType>
template <typename cpuVEC_VC>
__host__ bool cuVEC_VC<VType>::copyflags_from_cpuvec(cpuVEC_VC& vec_vc)
{
	//copy pbc parameters : need to copy these since they are associated with setting flags
	set_gpu_value(pbc_x, (int)vec_vc.pbc_x_ref());
	set_gpu_value(pbc_y, (int)vec_vc.pbc_y_ref());
	set_gpu_value(pbc_z, (int)vec_vc.pbc_z_ref());

	std::vector<int>& ngbrFlags_cpu = vec_vc.ngbrFlags_ref();

	cudaError_t error = cpu_to_gpu_managed(ngbrFlags, ngbrFlags_cpu.data(), ngbrFlags_cpu.size());
	if (error != cudaSuccess) return false;

	//copy ngbrFlags2 if needed
	if (vec_vc.ngbrFlags2_ref().size()) {

		//first allocate memory
		error = gpu_alloc_managed(ngbrFlags2, ngbrFlags_cpu.size());
		if (error != cudaSuccess) return false; //not enough memory so abort

		set_gpu_value(using_extended_flags, true);

		//now copy
		cpu_to_gpu_managed(ngbrFlags2, vec_vc.ngbrFlags2_ref().data(), ngbrFlags_cpu.size());
	}
	else {

		//ngbrFlags2 not used
		gpu_free_managed(ngbrFlags2);
		set_gpu_value(using_extended_flags, false);
	}

	return true;
}
