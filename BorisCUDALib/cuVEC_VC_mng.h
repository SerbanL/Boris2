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
	nullgpuptr(dirichlet_px);
	nullgpuptr(dirichlet_nx);
	nullgpuptr(dirichlet_py);
	nullgpuptr(dirichlet_ny);
	nullgpuptr(dirichlet_pz);
	nullgpuptr(dirichlet_nz);

	set_ngbrFlags_size(0);
	set_gpu_value(nonempty_cells, (int)0);

	set_dirichlet_size(0, NF_DIRICHLETPX);
	set_dirichlet_size(0, NF_DIRICHLETNX);
	set_dirichlet_size(0, NF_DIRICHLETPY);
	set_dirichlet_size(0, NF_DIRICHLETNY);
	set_dirichlet_size(0, NF_DIRICHLETPZ);
	set_dirichlet_size(0, NF_DIRICHLETNZ);

	set_robin(cuReal2(), NF_ROBINPX);
	set_robin(cuReal2(), NF_ROBINNX);
	set_robin(cuReal2(), NF_ROBINPY);
	set_robin(cuReal2(), NF_ROBINNY);
	set_robin(cuReal2(), NF_ROBINPZ);
	set_robin(cuReal2(), NF_ROBINNZ);
	set_robin(cuReal2(), NF_ROBINV);

	set_gpu_value(shift_debt, cuReal3());

	set_gpu_value(aSOR_damping, (cuBReal)1.0);
	set_gpu_value(aSOR_lasterror, (cuBReal)0.0);
	set_gpu_value(aSOR_lastgrad, (cuBReal)0.0);

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

//set robin values using a robin id
template <typename VType>
__host__ void cuVEC_VC<VType>::set_robin(cuReal2 robin_value, int robin_id)
{
	switch (robin_id) {

	case NF_ROBINPX:
		set_gpu_value(robin_px, robin_value);
		break;
		
	case NF_ROBINNX:
		set_gpu_value(robin_nx, robin_value);
		break;
		
	case NF_ROBINPY:
		set_gpu_value(robin_py, robin_value);
		break;
		
	case NF_ROBINNY:
		set_gpu_value(robin_ny, robin_value);
		break;
		
	case NF_ROBINPZ:
		set_gpu_value(robin_pz, robin_value);
		break;
		
	case NF_ROBINNZ:
		set_gpu_value(robin_nz, robin_value);
		break;
		
	case NF_ROBINV:
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

	switch (dirichlet_id) {

	case NF_DIRICHLETPX:
		error = gpu_alloc_managed(dirichlet_px, size);
		if (error != cudaSuccess) {

			nullgpuptr(dirichlet_px);
			size = 0;
		}
		set_gpu_value(dirichlet_px_size, size);
		break;
		
	case NF_DIRICHLETNX:
		error = gpu_alloc_managed(dirichlet_nx, size);
		if (error != cudaSuccess) {

			nullgpuptr(dirichlet_nx);
			size = 0;
		}
		set_gpu_value(dirichlet_nx_size, size);
		break;
		
	case NF_DIRICHLETPY:
		error = gpu_alloc_managed(dirichlet_py, size);
		if (error != cudaSuccess) {

			nullgpuptr(dirichlet_py);
			size = 0;
		}
		set_gpu_value(dirichlet_py_size, size);
		break;
		
	case NF_DIRICHLETNY:
		error = gpu_alloc_managed(dirichlet_ny, size);
		if (error != cudaSuccess) {

			nullgpuptr(dirichlet_ny);
			size = 0;
		}
		set_gpu_value(dirichlet_ny_size, size);
		break;
		
	case NF_DIRICHLETPZ:
		error = gpu_alloc_managed(dirichlet_pz, size);
		if (error != cudaSuccess) {

			nullgpuptr(dirichlet_pz);
			size = 0;
		}
		set_gpu_value(dirichlet_pz_size, size);
		break;
		
	case NF_DIRICHLETNZ:
		error = gpu_alloc_managed(dirichlet_nz, size);
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

	case NF_DIRICHLETPX:
		return get_gpu_value(dirichlet_px_size);
	
	case NF_DIRICHLETNX:
		return get_gpu_value(dirichlet_nx_size);
	
	case NF_DIRICHLETPY:
		return get_gpu_value(dirichlet_py_size);
	
	case NF_DIRICHLETNY:
		return get_gpu_value(dirichlet_ny_size);
	
	case NF_DIRICHLETPZ:
		return get_gpu_value(dirichlet_pz_size);
	
	case NF_DIRICHLETNZ:
		return get_gpu_value(dirichlet_nz_size);
	}

	return 0;
}

//--------------------------------------------CONSTRUCTORS : cu_obj "managed constructors" only. Real constructors are never called since you should never make a real instance of a cuVEC_VC.

//void constructor
template <typename VType>
__host__ void cuVEC_VC<VType>::construct_cu_obj(void)
{
	cuVEC::construct_cu_obj();

	alloc_initialize_data();
}

//construct to given number of cells : n_ is in cpu memory
template <typename VType>
__host__ void cuVEC_VC<VType>::construct_cu_obj(const cuSZ3& n_)
{
	cuVEC::construct_cu_obj(n_);

	alloc_initialize_data();

	if (get_gpu_value(n).dim()) {

		cudaError_t error = set_ngbrFlags_size(get_gpu_value(n).dim());
		if (error != cudaSuccess) set_n(cuSZ3());
		else gpu_set_managed(ngbrFlags, NF_NOTEMPTY, get_gpu_value(n).dim());
	}
}

//construct to given dimensions : h_ and rect_ are in cpu memory
template <typename VType>
__host__ void cuVEC_VC<VType>::construct_cu_obj(const cuReal3& h_, const cuRect& rect_)
{
	cuVEC::construct_cu_obj(h_, rect_);

	alloc_initialize_data();

	if (get_gpu_value(n).dim()) {

		cudaError_t error = set_ngbrFlags_size(get_gpu_value(n).dim());
		if (error != cudaSuccess) set_n(cuSZ3());
		else {

			gpu_set_managed(ngbrFlags, NF_NOTEMPTY, get_gpu_value(n).dim());
			set_ngbrFlags();
		}
	}
}

//construct to given dimensions and initialize to given value : h_ and rect_ are in cpu memory
template <typename VType>
__host__ void cuVEC_VC<VType>::construct_cu_obj(const cuReal3& h_, const cuRect& rect_, VType value)
{
	cuVEC::construct_cu_obj(h_, rect_, value);

	alloc_initialize_data();

	if (get_gpu_value(n).dim()) {

		cudaError_t error = set_ngbrFlags_size(get_gpu_value(n).dim());
		if (error != cudaSuccess) set_n(cuSZ3());
		else {

			gpu_set_managed(ngbrFlags, NF_NOTEMPTY, get_gpu_value(n).dim());
			set_ngbrFlags();
		}
	}
}

//copy constructor
template <typename VType>
__host__ void cuVEC_VC<VType>::construct_cu_obj(const cuVEC_VC<VType>& copyThis)
{
	cuVEC::construct_cu_obj();
	
	alloc_initialize_data();

	assign_cu_obj(copyThis);
}

//assignment operator
template <typename VType>
__host__ void cuVEC_VC<VType>::assign_cu_obj(const cuVEC_VC<VType>& copyThis)
{
	//copy n
	gpu_to_gpu(n, copyThis.n);
	//copy h
	gpu_to_gpu(h, copyThis.h);
	//copy rect
	gpu_to_gpu(rect, copyThis.rect);

	//allocate memory and copy array
	cudaError_t error = allocate_quantity(get_gpu_value(n));
	if (error == cudaSuccess) {

		gpu_to_gpu_managed(quantity, copyThis.quantity, get_gpu_value(n).dim());
	}
	else {

		set_h(cuReal3());
		set_rect(cuRect());
	}

	//transfer info might not be valid any more
	transfer.clear_transfer_data();

	//copy ngbrFlags_size
	gpu_to_gpu(ngbrFlags_size, copyThis.ngbrFlags_size);

	//allocate memory and copy ngbrFlags
	error = set_ngbrFlags_size(get_gpu_value(n).dim());
	if (error == cudaSuccess) {

		gpu_to_gpu_managed(ngbrFlags, copyThis.ngbrFlags, get_ngbrFlags_size());
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
	
	size = get_dirichlet_size(NF_DIRICHLETPX);
	if (size) {

		error = set_dirichlet_size(size, NF_DIRICHLETPX);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(dirichlet_px, copyThis.dirichlet_px, size);
		}
	}

	size = get_dirichlet_size(NF_DIRICHLETNX);
	if (size) {

		error = set_dirichlet_size(size, NF_DIRICHLETNX);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(dirichlet_nx, copyThis.dirichlet_nx, size);
		}
	}

	size = get_dirichlet_size(NF_DIRICHLETPY);
	if (size) {

		error = set_dirichlet_size(size, NF_DIRICHLETPY);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(dirichlet_py, copyThis.dirichlet_py, size);
		}
	}

	size = get_dirichlet_size(NF_DIRICHLETNY);
	if (size) {

		error = set_dirichlet_size(size, NF_DIRICHLETNY);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(dirichlet_ny, copyThis.dirichlet_ny, size);
		}
	}

	size = get_dirichlet_size(NF_DIRICHLETPZ);
	if (size) {

		error = set_dirichlet_size(size, NF_DIRICHLETPZ);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(dirichlet_pz, copyThis.dirichlet_pz, size);
		}
	}

	size = get_dirichlet_size(NF_DIRICHLETNZ);
	if (size) {

		error = set_dirichlet_size(size, NF_DIRICHLETNZ);
		if (error == cudaSuccess) {

			gpu_to_gpu_managed(dirichlet_nz, copyThis.dirichlet_nz, size);
		}
	}

	//copy aSOR parameters
	gpu_go_gpu(aSOR_lasterror, copyThis.aSOR_lasterror);
	gpu_go_gpu(aSOR_damping, copyThis.aSOR_damping);
	gpu_go_gpu(aSOR_lastgrad, copyThis.aSOR_lastgrad);

	//copy pbc settings
	gpu_go_gpu(pbc_x, copyThis.pbc_x);
	gpu_go_gpu(pbc_y, copyThis.pbc_y);
	gpu_go_gpu(pbc_z, copyThis.pbc_z);
}

//destructor
template <typename VType>
__host__ void cuVEC_VC<VType>::destruct_cu_obj(void)
{
	gpu_free_managed(ngbrFlags);
	gpu_free_managed(dirichlet_px);
	gpu_free_managed(dirichlet_nx);
	gpu_free_managed(dirichlet_py);
	gpu_free_managed(dirichlet_ny);
	gpu_free_managed(dirichlet_pz);
	gpu_free_managed(dirichlet_nz);

	cuVEC::destruct_cu_obj();
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
	//aSOR parameters

	//-----------

	cuSZ3 cpu_n = vec_vc.n;
	cuReal3 cpu_h = vec_vc.h;
	cuRect cpu_rect = vec_vc.rect;

	//-----------

	cuVEC::set_from_cpuvec(vec_vc);

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

	//-----------

	//copy nonempty_cells
	set_gpu_value(nonempty_cells, vec_vc.get_nonempty_cells());

	//-----------

	//copy dirichlet vectors
	
	size_t dirichlet_size;

	//---NF_DIRICHLETPX
	
	//allocate gpu memory if sizes don't match
	dirichlet_size = vec_vc.dirichlet_px_ref().size();
	if (get_dirichlet_size(NF_DIRICHLETPX) != dirichlet_size) {

		error = gpu_alloc_managed(dirichlet_px, dirichlet_size);
		if (error != cudaSuccess) return false; //not enough memory so abort

		//set size in gpu memory
		set_gpu_value(dirichlet_px_size, dirichlet_size);
	}

	//copy over
	cpu_to_gpu_managed(dirichlet_px, vec_vc.dirichlet_px_ref().data(), dirichlet_size);

	//---NF_DIRICHLETNX

	//allocate gpu memory if sizes don't match
	dirichlet_size = vec_vc.dirichlet_nx_ref().size();
	if (get_dirichlet_size(NF_DIRICHLETNX) != dirichlet_size) {

		error = gpu_alloc_managed(dirichlet_nx, dirichlet_size);
		if (error != cudaSuccess) return false; //not enough memory so abort

		//set size in gpu memory
		set_gpu_value(dirichlet_nx_size, dirichlet_size);
	}

	//copy over
	cpu_to_gpu_managed(dirichlet_nx, vec_vc.dirichlet_nx_ref().data(), dirichlet_size);
	
	//---NF_DIRICHLETPY

	//allocate gpu memory if sizes don't match
	dirichlet_size = vec_vc.dirichlet_py_ref().size();
	if (get_dirichlet_size(NF_DIRICHLETPY) != dirichlet_size) {

		error = gpu_alloc_managed(dirichlet_py, dirichlet_size);
		if (error != cudaSuccess) return false; //not enough memory so abort

		//set size in gpu memory
		set_gpu_value(dirichlet_py_size, dirichlet_size);
	}

	//copy over
	cpu_to_gpu_managed(dirichlet_py, vec_vc.dirichlet_py_ref().data(), dirichlet_size);

	//---NF_DIRICHLETNY

	//allocate gpu memory if sizes don't match
	dirichlet_size = vec_vc.dirichlet_ny_ref().size();
	if (get_dirichlet_size(NF_DIRICHLETNY) != dirichlet_size) {

		error = gpu_alloc_managed(dirichlet_ny, dirichlet_size);
		if (error != cudaSuccess) return false; //not enough memory so abort

		//set size in gpu memory
		set_gpu_value(dirichlet_ny_size, dirichlet_size);
	}

	//copy over
	cpu_to_gpu_managed(dirichlet_ny, vec_vc.dirichlet_ny_ref().data(), dirichlet_size);

	//---NF_DIRICHLETPZ

	//allocate gpu memory if sizes don't match
	dirichlet_size = vec_vc.dirichlet_pz_ref().size();
	if (get_dirichlet_size(NF_DIRICHLETPZ) != dirichlet_size) {

		error = gpu_alloc_managed(dirichlet_pz, dirichlet_size);
		if (error != cudaSuccess) return false; //not enough memory so abort

		//set size in gpu memory
		set_gpu_value(dirichlet_pz_size, dirichlet_size);
	}

	//copy over
	cpu_to_gpu_managed(dirichlet_pz, vec_vc.dirichlet_pz_ref().data(), dirichlet_size);

	//---NF_DIRICHLETNZ

	//allocate gpu memory if sizes don't match
	dirichlet_size = vec_vc.dirichlet_nz_ref().size();
	if (get_dirichlet_size(NF_DIRICHLETNZ) != dirichlet_size) {

		error = gpu_alloc_managed(dirichlet_nz, dirichlet_size);
		if (error != cudaSuccess) return false; //not enough memory so abort

		//set size in gpu memory
		set_gpu_value(dirichlet_nz_size, dirichlet_size);
	}

	//copy over
	cpu_to_gpu_managed(dirichlet_nz, vec_vc.dirichlet_nz_ref().data(), dirichlet_size);
	
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

	//copy aSOR parameters
	set_gpu_value(aSOR_damping, (cuBReal)vec_vc.aSOR_damping_ref());

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
	//aSOR parameters

	//----------- copy to VEC part (n, h, rect, quantity)

	cuVEC::set_cpuvec(vec_vc);

	//-----------

	SZ3 cpu_n = get_gpu_value(n);
	DBL3 cpu_h = get_gpu_value(h);
	Rect cpu_rect = get_gpu_value(rect);

	//-----------

	//copy ngbrFlags
	if (vec_vc.ngbrFlags_ref().size() != cpu_n.dim()) {

		//first allocate memory if sizes don't match
		if (!malloc_vector(vec_vc.ngbrFlags_ref(), cpu_n.dim())) return false;		//couldn't allocate memory
	}

	//now copy
	gpu_to_cpu_managed(vec_vc.ngbrFlags_ref().data(), ngbrFlags, cpu_n.dim());

	//-----------

	//copy nonempty_cells
	vec_vc.nonempty_cells_ref() = get_gpu_value(nonempty_cells);

	//-----------
	
	//Set dirichlet vectors
	
	size_t dirichlet_size;

	//---NF_DIRICHLETPX

	dirichlet_size = get_dirichlet_size(NF_DIRICHLETPX);

	//copy dirichlet array to cpu memory
	if (!malloc_vector(vec_vc.dirichlet_px_ref(), dirichlet_size)) return false;		//couldn't allocate memory	
	gpu_to_cpu_managed(vec_vc.dirichlet_px_ref().data(), dirichlet_px, dirichlet_size);

	//---NF_DIRICHLETNX

	dirichlet_size = get_dirichlet_size(NF_DIRICHLETNX);

	//copy dirichlet array to cpu memory
	if (!malloc_vector(vec_vc.dirichlet_nx_ref(), dirichlet_size)) return false;		//couldn't allocate memory	
	gpu_to_cpu_managed(vec_vc.dirichlet_nx_ref().data(), dirichlet_nx, dirichlet_size);

	//---NF_DIRICHLETPY

	dirichlet_size = get_dirichlet_size(NF_DIRICHLETPY);

	//copy dirichlet array to cpu memory
	if (!malloc_vector(vec_vc.dirichlet_py_ref(), dirichlet_size)) return false;		//couldn't allocate memory	
	gpu_to_cpu_managed(vec_vc.dirichlet_py_ref().data(), dirichlet_py, dirichlet_size);

	//---NF_DIRICHLETNY

	dirichlet_size = get_dirichlet_size(NF_DIRICHLETNY);

	//copy dirichlet array to cpu memory
	if (!malloc_vector(vec_vc.dirichlet_ny_ref(), dirichlet_size)) return false;		//couldn't allocate memory	
	gpu_to_cpu_managed(vec_vc.dirichlet_ny_ref().data(), dirichlet_ny, dirichlet_size);

	//---NF_DIRICHLETPZ

	dirichlet_size = get_dirichlet_size(NF_DIRICHLETPZ);

	//copy dirichlet array to cpu memory
	if (!malloc_vector(vec_vc.dirichlet_pz_ref(), dirichlet_size)) return false;		//couldn't allocate memory	
	gpu_to_cpu_managed(vec_vc.dirichlet_pz_ref().data(), dirichlet_pz, dirichlet_size);

	//---NF_DIRICHLETNZ

	dirichlet_size = get_dirichlet_size(NF_DIRICHLETNZ);

	//copy dirichlet array to cpu memory
	if (!malloc_vector(vec_vc.dirichlet_nz_ref(), dirichlet_size)) return false;		//couldn't allocate memory	
	gpu_to_cpu_managed(vec_vc.dirichlet_nz_ref().data(), dirichlet_nz, dirichlet_size);

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

	//copy aSOR parameters
	vec_vc.aSOR_damping_ref() = get_gpu_value(aSOR_damping);

	//copy pbc parameters
	vec_vc.pbc_x_ref() = get_gpu_value(pbc_x);
	vec_vc.pbc_y_ref() = get_gpu_value(pbc_y);
	vec_vc.pbc_z_ref() = get_gpu_value(pbc_z);

	return true;
}

//faster version of set_from_cpuvec, where it is assumed the cpu vec already has the same sizes as this cuVEC : only quantity is copied.
template <typename VType>
template <typename cpuVEC_VC>
__host__ bool cuVEC_VC<VType>::copy_from_cpuvec(cpuVEC_VC& vec_vc)
{
	//copy quantity
	if (!cuVEC::copy_from_cpuvec(vec_vc)) return false;

	//copy shift debt
	set_gpu_value(shift_debt, (cuReal3)vec_vc.shift_debt_ref());

	//copy aSOR parameters
	set_gpu_value(aSOR_damping, (cuBReal)vec_vc.aSOR_damping_ref());

	//copy pbc parameters
	set_gpu_value(pbc_x, (int)vec_vc.pbc_x_ref());
	set_gpu_value(pbc_y, (int)vec_vc.pbc_y_ref());
	set_gpu_value(pbc_z, (int)vec_vc.pbc_z_ref());

	//now copy ngbrFlags
	cudaError_t error = cpu_to_gpu_managed(ngbrFlags, vec_vc.ngbrFlags_ref().data(), vec_vc.linear_size());
	if (error != cudaSuccess) return false;
	else return true;
}

//faster version of set_cpuvec, where it is assumed the cpu vec already has the same sizes as this cuVEC_VC : only quantity and ngbrFlags are copied.
template <typename VType>
template <typename cpuVEC_VC>
__host__ bool cuVEC_VC<VType>::copy_to_cpuvec(cpuVEC_VC& vec_vc)
{
	//copy quantity
	if (!cuVEC::copy_to_cpuvec(vec_vc)) return false;

	//copy shift debt
	vec_vc.shift_debt_ref() = get_gpu_value(shift_debt);

	//copy aSOR parameters
	vec_vc.aSOR_damping_ref() = get_gpu_value(aSOR_damping);

	//copy pbc parameters
	vec_vc.pbc_x_ref() = get_gpu_value(pbc_x);
	vec_vc.pbc_y_ref() = get_gpu_value(pbc_y);
	vec_vc.pbc_z_ref() = get_gpu_value(pbc_z);

	//now copy ngbrFlags
	cudaError_t error = gpu_to_cpu_managed(vec_vc.ngbrFlags_ref().data(), ngbrFlags, vec_vc.linear_size());
	if (error != cudaSuccess) return false;
	else return true;
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
	else return true;
}
