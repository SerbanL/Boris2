#pragma once

#include "cuVEC.h"
#include "launchers.h"

//--------------------------------------------MEMORY MANAGEMENT HELPER METHODS

//memory allocation for objects and initialize to default - only call at start of managed constructors
template <typename VType>
__host__ void cuVEC<VType>::alloc_initialize_data(void)
{	
	nullgpuptr(quantity);
	nullgpuptr(aux_block_values);

	nullgpuptr(line_profile);
	nullgpuptr(line_profile_component_x);
	nullgpuptr(line_profile_component_y);
	nullgpuptr(line_profile_component_z);
	nullgpuptr(line_profile_avpoints);
	set_gpu_value(line_profile_component_size, (size_t)0);

	nullgpuptr(histogram);
	set_gpu_value(histogram_size, (size_t)0);

	nullgpuptr(histogram_preaverage);
	set_gpu_value(histogram_preaverage_size, (size_t)0);

	set_n(cuSZ3());
	set_h(cuReal3());
	set_rect(cuRect());

	transfer.construct_cu_obj();
	transfer2.construct_cu_obj();
}

//--------------------------------------------HELPER METHODS

template <typename VType>
__host__ void cuVEC<VType>::SetMeshRect(void)
{ 
	cuRect cpu_rect_value = cuRect(get_gpu_value(rect).s, get_gpu_value(rect).s + (get_gpu_value(h) & get_gpu_value(n)));
	set_gpu_value(rect, cpu_rect_value);
}

template <typename VType>
__host__ cuSZ3 cuVEC<VType>::get_n_from_h_and_rect(const cuReal3& h_, const cuRect& rect_)
{
	//calculate new n from rect and current h
	cuSZ3 new_n = cu_round(rect_ / h_);
	if (new_n.x < 1) new_n.x = 1;
	if (new_n.y < 1) new_n.y = 1;
	if (new_n.z < 1) new_n.z = 1;

	return new_n;
}

//from current rectangle and h value set n. h may also need to be adjusted since n must be an integer
template <typename VType>
__host__ bool cuVEC<VType>::set_n_adjust_h(void)
{
	//make sure the cellsize divides the mesh rectangle

	//calculate new n from rect and current h

	cuSZ3 new_n = get_n_from_h_and_rect(get_gpu_value(h), get_gpu_value(rect));

	//adjust h for new n (but save it first in case memory allocation fails)
	cuReal3 h_save = get_gpu_value(h);
	set_h(get_gpu_value(rect) / new_n);

	//now set n - either map to new sizes (if current object is not empty) or allocate memory for new object
	if (!resize(new_n)) {

		//failed : go back to previous h (n unchanged)
		set_h(h_save);
		return false;
	}
	else return true;
}

template <typename VType>
__host__ bool cuVEC<VType>::allocate_profile_component_memory(size_t size)
{
	if (size == get_gpu_value(line_profile_component_size)) return true;

	auto fail_memory_allocation = [&]() {

		gpu_free_managed(line_profile);
		gpu_free_managed(line_profile_component_x);
		gpu_free_managed(line_profile_component_y);
		gpu_free_managed(line_profile_component_z);
		gpu_free_managed(line_profile_avpoints);
		set_gpu_value(line_profile_component_size, (size_t)0);
	};

	//try to allocate required memory for profile component array
	cudaError_t error = gpu_alloc_managed(line_profile, size);
	if (error != cudaSuccess) { fail_memory_allocation(); return false; }
	error = gpu_alloc_managed(line_profile_component_x, size);
	if (error != cudaSuccess) { fail_memory_allocation(); return false; }
	error = gpu_alloc_managed(line_profile_component_y, size);
	if (error != cudaSuccess) { fail_memory_allocation(); return false; }
	error = gpu_alloc_managed(line_profile_component_z, size);
	if (error != cudaSuccess) { fail_memory_allocation(); return false; }
	error = gpu_alloc_managed(line_profile_avpoints, size);
	if (error != cudaSuccess) { fail_memory_allocation(); return false; }

	set_gpu_value(line_profile_component_size, size);
	return true;
}

//memory management for histogram array : attempt to resize to new size if not already exactly given size
template <typename VType>
__host__ bool cuVEC<VType>::allocate_histogram_memory(size_t histogram_size_cpu, size_t histogram_preaverage_size_cpu)
{
	if (histogram_size_cpu == get_gpu_value(histogram_size) && histogram_preaverage_size_cpu == get_gpu_value(histogram_preaverage_size)) return true;

	auto fail_memory_allocation = [&]() {

		gpu_free_managed(histogram);
		set_gpu_value(histogram_size, (size_t)0);

		gpu_free_managed(histogram_preaverage);
		set_gpu_value(histogram_preaverage_size, (size_t)0);
	};

	//try to allocate required memory for profile component array
	cudaError_t error = cudaSuccess;
		
	if (histogram_size_cpu != get_gpu_value(histogram_size)) {

		if (histogram_size_cpu > 0) {

			error = gpu_alloc_managed(histogram, histogram_size_cpu);
			if (error != cudaSuccess) { fail_memory_allocation(); return false; }
		}
		else gpu_free_managed(histogram);

		set_gpu_value(histogram_size, histogram_size_cpu);
	}
	
	if (histogram_preaverage_size_cpu != get_gpu_value(histogram_preaverage_size)) {

		if (histogram_preaverage_size_cpu > 0) {

			error = gpu_alloc_managed(histogram_preaverage, histogram_preaverage_size_cpu);
			if (error != cudaSuccess) { fail_memory_allocation(); return false; }
		}
		else gpu_free_managed(histogram_preaverage);

		set_gpu_value(histogram_preaverage_size, histogram_preaverage_size_cpu);
	}

	return true;
}

//--------------------------------------------GET/SET FROM/TO GPU MEMORY

template <typename VType>
__host__ cudaError cuVEC<VType>::allocate_quantity(cuSZ3 new_n)
{
	if (new_n == get_gpu_value(n)) return cudaSuccess;

	//try to allocate required memory for quantity
	cudaError_t error = gpu_alloc_managed(quantity, new_n.dim());

	//must also adjust space for aux_block_values
	if (error == cudaSuccess) error = gpu_alloc_managed(aux_block_values, (new_n.dim() + CUDATHREADS) / CUDATHREADS);

	//copy n value if memory allocated correctly
	if (error == cudaSuccess) {

		set_n(new_n);
	}
	//failed : set nullptr
	else {

		nullgpuptr(quantity);
		nullgpuptr(aux_block_values);
		set_n(cuSZ3());
	}

	return error;
}

//set n in gpu memory from cpu memory value - can take both l-values and r-values
template <typename VType>
__host__ void cuVEC<VType>::set_n(const cuSZ3& n_)
{
	set_gpu_value(n, n_);
}

//set h in gpu memory from cpu memory value - can take both l-values and r-values
template <typename VType>
__host__ void cuVEC<VType>::set_h(const cuReal3& h_)
{
	set_gpu_value(h, h_);
}

//set rect gpu memory from cpu memory value - can take both l-values and r-values
template <typename VType>
__host__ void cuVEC<VType>::set_rect(const cuRect& rect_)
{
	set_gpu_value(rect, rect_);
}

//--------------------------------------------CONSTRUCTORS : cu_obj "managed constructors" only. Real constructors are never called since you should never make a real instance of a cuVEC.

//void constructor
template <typename VType>
__host__ void cuVEC<VType>::construct_cu_obj(void)
{	
	alloc_initialize_data();
}

//construct to given number of cells : n_ is in cpu memory
template <typename VType>
__host__ void cuVEC<VType>::construct_cu_obj(const cuSZ3& n_)
{
	alloc_initialize_data();

	cudaError_t error = allocate_quantity(n_);
	if (error == cudaSuccess) gpu_set_managed(quantity, VType(), n_.dim());
}
	
//construct to given dimensions : h_ and rect_ are in cpu memory
template <typename VType>
__host__ void cuVEC<VType>::construct_cu_obj(const cuReal3& h_, const cuRect& rect_)
{
	alloc_initialize_data();

	set_h(h_);
	set_rect(rect_);

	if (!set_n_adjust_h()) {

		set_h(cuReal3());
		set_rect(cuRect());
	}
	else {

		//initialize to value
		gpu_set_managed(quantity, VType(), get_gpu_value(n).dim());
	}
}

//construct to given dimensions and initialize to given value : h_ and rect_ are in cpu memory
template <typename VType>
__host__ void cuVEC<VType>::construct_cu_obj(const cuReal3& h_, const cuRect& rect_, VType value)
{
	alloc_initialize_data();

	set_h(h_);
	set_rect(rect_);

	if (!set_n_adjust_h()) {

		set_h(cuReal3());
		set_rect(cuRect());
	}
	else {

		//initialize to value
		gpu_set_managed(quantity, value, get_gpu_value(n).dim());
	}
}

//copy constructor
template <typename VType>
__host__ void cuVEC<VType>::construct_cu_obj(const cuVEC<VType>& copyThis)
{
	alloc_initialize_data();

	assign_cu_obj(copyThis);
}

//assignment operator
template <typename VType>
__host__ void cuVEC<VType>::assign_cu_obj(const cuVEC<VType>& copyThis)
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
	transfer2.clear_transfer_data();
}

template <typename VType>
__host__ void cuVEC<VType>::destruct_cu_obj(void)
{
	gpu_free_managed(quantity);
	gpu_free_managed(aux_block_values);
	
	gpu_free_managed(line_profile);
	gpu_free_managed(line_profile_component_x);
	gpu_free_managed(line_profile_component_y);
	gpu_free_managed(line_profile_component_z);
	gpu_free_managed(line_profile_avpoints);
	set_gpu_value(line_profile_component_size, (size_t)0);

	gpu_free_managed(histogram);
	set_gpu_value(histogram_size, (size_t)0);

	gpu_free_managed(histogram_preaverage);
	set_gpu_value(histogram_preaverage_size, (size_t)0);

	transfer.destruct_cu_obj();
	transfer2.destruct_cu_obj();
}

//--------------------------------------------SIZING

template <typename VType>
__host__ bool cuVEC<VType>::resize(cuSZ3 new_n)
{
	if (new_n == get_gpu_value(n)) return true;

	//transfer info might not be valid any more
	transfer.clear_transfer_data();
	transfer2.clear_transfer_data();

	//current zero size : set new size and rect
	if (get_gpu_value(n) == cuSZ3(0)) {

		//check if memory for new size can be allocated - it might still fail even after passing this check
		if (new_n.dim() <= get_gpu_value(n).dim() || new_n.dim() - get_gpu_value(n).dim() <= cudaMemGetFree() / sizeof(VType)) {

			cudaError_t error = allocate_quantity(new_n);
			if (error != cudaSuccess) return false;

			//initialize to default
			gpu_set_managed(quantity, VType(), new_n.dim());

			SetMeshRect();
			return true;
		}
	}
	else {

		//remap to new size and set rect
		if (mapmesh_newdims(new_n)) {

			SetMeshRect();
			return true;
		}
	}
	
	//if here then couldn't allocate memory : fail and previous size maintained
	return false;
}

template <typename VType>
__host__ bool cuVEC<VType>::resize(cuReal3 new_h, cuRect new_rect)
{
	if (new_h == get_gpu_value(h) && new_rect == get_gpu_value(rect)) return true;

	//transfer info might not be valid any more
	transfer.clear_transfer_data();
	transfer2.clear_transfer_data();

	//save h and rect in case we cannot resize
	cuReal3 save_h = get_gpu_value(h);
	cuRect save_rect = get_gpu_value(rect);

	//set new h and rect for now
	set_h(new_h);
	set_rect(new_rect);

	if (!set_n_adjust_h()) {

		//failed : go back to previous dimensions
		set_h(save_h);
		set_rect(save_rect);

		return false;
	}
	else return true;
}

template <typename VType>
__host__ bool cuVEC<VType>::assign(cuSZ3 new_n, VType value)
{ 
	//transfer info might not be valid any more
	transfer.clear_transfer_data();
	transfer2.clear_transfer_data();

	//is there enough memory?
	if (new_n.dim() <= get_gpu_value(n).dim() || new_n.dim() - get_gpu_value(n).dim() <= cudaMemGetFree() / sizeof(VType)) {

		//maybe...
		cudaError_t error = allocate_quantity(new_n);
		if (error != cudaSuccess) return false;

		gpu_set_managed(quantity, value, get_gpu_value(n).dim());

		SetMeshRect();
		return true;
	}
	else return false;
}

//works like resize(h_, rect_) but sets given value also
template <typename VType>
__host__ bool cuVEC<VType>::assign(cuReal3 new_h, cuRect new_rect, VType value)
{ 
	//transfer info might not be valid any more
	transfer.clear_transfer_data();
	transfer2.clear_transfer_data();

	//calculate new n from rect and current h
	cuSZ3 new_n = get_n_from_h_and_rect(new_h, new_rect);

	//is there enough memory?
	if (new_n.dim() <= get_gpu_value(n).dim() || new_n.dim() - get_gpu_value(n).dim() <= cudaMemGetFree() / sizeof(VType)) {

		//maybe...
		cudaError_t error = allocate_quantity(new_n);
		if (error != cudaSuccess) return false;

		gpu_set_managed(quantity, value, get_gpu_value(n).dim());

		//now set dimensions and set value
		set_h(new_rect / new_n);
		set_rect(new_rect);

		return true;
	}
	else return false;
}

template <typename VType>
__host__ void cuVEC<VType>::clear(void)
{
	transfer.clear_transfer_data();
	transfer2.clear_transfer_data();

	set_n(cuSZ3());
	gpu_free_managed(quantity);
	gpu_free_managed(aux_block_values);

	gpu_free_managed(line_profile);
	gpu_free_managed(line_profile_component_x);
	gpu_free_managed(line_profile_component_y);
	gpu_free_managed(line_profile_component_z);
	gpu_free_managed(line_profile_avpoints);
	set_gpu_value(line_profile_component_size, (size_t)0);

	gpu_free_managed(histogram);
	set_gpu_value(histogram_size, (size_t)0);

	gpu_free_managed(histogram_preaverage);
	set_gpu_value(histogram_preaverage_size, (size_t)0);

	SetMeshRect(); 
}

//set rect start (i.e. shift the entire rectangle to align with given absolute starting coordinates
template <typename VType>
__host__ void cuVEC<VType>::set_rect_start(const cuReal3& rect_start)
{
	cuRect rect_cpu = get_gpu_value(rect);
	rect_cpu += (rect_start - rect_cpu.s);
	set_rect(rect_cpu);
}

template <typename VType>
__host__ void cuVEC<VType>::shift_rect_start(const cuReal3& shift)
{
	cuRect rect_cpu = get_gpu_value(rect);
	rect_cpu += shift;
	set_rect(rect_cpu);
}

//--------------------------------------------COPY FROM VEC

//copy everything from a VEC - type must be convertible. Return false if failed (memory could not be allocated)
template <typename VType>
template <typename cpuVEC>
__host__ bool cuVEC<VType>::set_from_cpuvec(cpuVEC& vec)
{
	//transfer info might not be valid any more
	transfer.clear_transfer_data();
	transfer2.clear_transfer_data();

	//-----------

	//things to copy:
	//n, h, rect
	//quantity

	cuSZ3 cpu_n = vec.n;
	cuReal3 cpu_h = vec.h;
	cuRect cpu_rect = vec.rect;

	//-----------
	
	cudaError_t error = cudaSuccess;

	//first allocate memory if sizes don't match
	error = allocate_quantity(cpu_n);
	if (error != cudaSuccess) return false;
	
	//now copy
	error = cpu_to_gpu_managed(quantity, vec.data(), vec.linear_size());
	if (error != cudaSuccess) return false;

	//-----------

	//copy sizes : n, h, rect
	set_gpu_value(n, cpu_n);
	set_gpu_value(h, cpu_h);
	set_gpu_value(rect, cpu_rect);
	
	return true;
}

//--------------------------------------------COPY TO VEC

//copy everything to a VEC - type must be convertible. Return false if failed (memory could not be allocated)
template <typename VType>
template <typename cpuVEC>
__host__ bool cuVEC<VType>::set_cpuvec(cpuVEC& vec)
{
	//-----------

	//things to copy:
	//n, h, rect
	//quantity

	//-----------

	SZ3 cpu_n = get_gpu_value(n);
	DBL3 cpu_h = get_gpu_value(h);
	Rect cpu_rect = get_gpu_value(rect);

	//-----------

	cudaError_t error = cudaSuccess;

	//First allocate memory
	if (!malloc_vector(vec.quantity_ref(), cpu_n.dim())) return false;		//couldn't allocate memory

	//set sizes
	vec.n = cpu_n;
	vec.h = cpu_h;
	vec.rect = cpu_rect;

	//now copy quantity
	error = gpu_to_cpu_managed(vec.quantity_ref().data(), quantity, cpu_n.dim());
	if (error != cudaSuccess) return false;

	return true;
}

//faster version of set_from_cpuvec, where it is assumed the cpu vec already has the same sizes as this cuVEC : only quantity is copied.
template <typename VType>
template <typename cpuVEC>
__host__ bool cuVEC<VType>::copy_from_cpuvec(cpuVEC& vec)
{
	cudaError_t error = cpu_to_gpu_managed(quantity, vec.data(), vec.linear_size());
	if (error != cudaSuccess) return false;

	return true;
}

//faster version of set_cpuvec, where it is assumed the cpu vec already has the same sizes as this cuVEC : only quantity is copied.
template <typename VType>
template <typename cpuVEC>
__host__ bool cuVEC<VType>::copy_to_cpuvec(cpuVEC& vec)
{
	cudaError_t error = gpu_to_cpu_managed(vec.data(), quantity, vec.linear_size());
	if (error != cudaSuccess) return false;

	return true;
}

template <typename VType>
template <typename SType>
__host__ bool cuVEC<VType>::copy_from_vector(std::vector<SType>& vec)
{
	if (vec.size() != get_gpu_value(n).dim()) return false;

	cudaError_t error = cpu_to_gpu_managed(quantity, vec.data(), vec.size());
	if (error != cudaSuccess) return false;

	else return true;
}

template <typename VType>
template <typename SType>
__host__ bool cuVEC<VType>::copy_to_vector(std::vector<SType>& vec)
{
	if (vec.size() != get_gpu_value(n).dim()) return false;

	cudaError_t error = gpu_to_cpu_managed(vec.data(), quantity, vec.size());
	if (error != cudaSuccess) return false;

	return true;
}
