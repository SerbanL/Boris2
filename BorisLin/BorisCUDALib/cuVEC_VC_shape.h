#pragma once

#include "cuVEC_VC.h"

//--------------------------------------------SIZING

//resize and set shape using linked vec
template <typename VType>
template <typename LVType>
__host__ bool cuVEC_VC<VType>::resize(cuSZ3 new_n, cuVEC_VC<LVType>& linked_vec)
{
	if (new_n != get_gpu_value(cuVEC<VType>::n)) {

		//memory reserved so map flags to new size before changing the n value
		cudaError_t error = resize_ngbrFlags(new_n);
		if (error != cudaSuccess) return false;

		if (!cuVEC<VType>::resize(new_n)) return false;
	}

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec.n, linked_vec.h, linked_vec.rect, linked_vec.ngbrFlags);

	return true;
}

//resize but keep shape
template <typename VType>
__host__ bool cuVEC_VC<VType>::resize(cuSZ3 new_n)
{
	if (new_n != get_gpu_value(cuVEC<VType>::n)) {

		//memory reserved so map flags to new size before changing the n value
		cudaError_t error = resize_ngbrFlags(new_n);
		if (error != cudaSuccess) return false;

		if (!cuVEC<VType>::resize(new_n)) return false;

		//all good, finish off by setting flags.
		set_ngbrFlags();
	}

	return true;
}

//resize and set shape using linked vec
template <typename VType>
template <typename LVType>
__host__ bool cuVEC_VC<VType>::resize(cuReal3 new_h, cuRect new_rect, cuVEC_VC<LVType>& linked_vec)
{
	if (new_h != get_gpu_value(cuVEC<VType>::h) || new_rect != get_gpu_value(cuVEC<VType>::rect)) {

		cuSZ3 new_n = cuVEC<VType>::get_n_from_h_and_rect(new_h, new_rect);

		//memory reserved so map flags to new size before changing the n value
		cudaError_t error = resize_ngbrFlags(new_n);
		if (error != cudaSuccess) return false;

		if (!cuVEC<VType>::resize(new_h, new_rect)) return false;
	}

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec.n, linked_vec.h, linked_vec.rect, linked_vec.ngbrFlags_ref());

	return true;
}

//resize but keep shape
template <typename VType>
__host__ bool cuVEC_VC<VType>::resize(cuReal3 new_h, cuRect new_rect)
{
	if (new_h != get_gpu_value(cuVEC<VType>::h) || new_rect != get_gpu_value(cuVEC<VType>::rect)) {

		cuSZ3 new_n = cuVEC<VType>::get_n_from_h_and_rect(new_h, new_rect);

		//memory reserved so map flags to new size before changing the n value
		cudaError_t error = resize_ngbrFlags(new_n);
		if (error != cudaSuccess) return false;

		if (!cuVEC<VType>::resize(new_h, new_rect)) return false;

		//all good, finish off by setting flags.
		set_ngbrFlags();
	}

	return true;
}

//set value and shape from linked vec
template <typename VType>
template <typename LVType>
__host__ bool cuVEC_VC<VType>::assign(cuSZ3 new_n, VType value, cuVEC_VC<LVType>& linked_vec)
{
	if (new_n != get_gpu_value(cuVEC<VType>::n)) {

		//memory reserved so map flags to new size before changing the n value
		cudaError_t error = resize_ngbrFlags(new_n);
		if (error != cudaSuccess) return false;
	}
	
	if (!cuVEC<VType>::assign(new_n, value)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec.n, linked_vec.h, linked_vec.rect, linked_vec.ngbrFlags);

	return true;
}

//set value but keep shape - empty cells will retain zero value : i.e. set value everywhere but in empty cells
template <typename VType>
__host__ bool cuVEC_VC<VType>::assign(cuSZ3 new_n, VType value)
{
	if (new_n != get_gpu_value(cuVEC<VType>::n)) {

		//memory reserved so map flags to new size before changing the n value
		cudaError_t error = resize_ngbrFlags(new_n);
		if (error != cudaSuccess) return false;

	}

	if (!cuVEC<VType>::assign(new_n, value)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags();

	return true;
}

//set value and shape from linked vec
template <typename VType>
template <typename LVType>
__host__ bool cuVEC_VC<VType>::assign(cuReal3 new_h, cuRect new_rect, VType value, cuVEC_VC<LVType>& linked_vec)
{
	if (new_h != get_gpu_value(cuVEC<VType>::h) || new_rect != get_gpu_value(cuVEC<VType>::rect)) {

		cuSZ3 new_n = cuVEC<VType>::get_n_from_h_and_rect(new_h, new_rect);

		//memory reserved so map flags to new size before changing the n value
		cudaError_t error = resize_ngbrFlags(new_n);
		if (error != cudaSuccess) return false;
	}

	if (!cuVEC<VType>::assign(new_h, new_rect, value)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec.n, linked_vec.h, linked_vec.rect, linked_vec.ngbrFlags_ref());

	return true;
}

//set value but keep shape - empty cells will retain zero value : i.e. set value everywhere but in empty cells
template <typename VType>
__host__ bool cuVEC_VC<VType>::assign(cuReal3 new_h, cuRect new_rect, VType value)
{
	if (new_h != get_gpu_value(cuVEC<VType>::h) || new_rect != get_gpu_value(cuVEC<VType>::rect)) {

		cuSZ3 new_n = cuVEC<VType>::get_n_from_h_and_rect(new_h, new_rect);

		//memory reserved so map flags to new size before changing the n value
		cudaError_t error = resize_ngbrFlags(new_n);
		if (error != cudaSuccess) return false;
	}

	if (!cuVEC<VType>::assign(new_h, new_rect, value)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags();

	return true;
}

template <typename VType>
__host__ void cuVEC_VC<VType>::clear(void)
{
	cuVEC<VType>::clear();
	set_ngbrFlags_size(0);

	gpu_free_managed(ngbrFlags2);
	set_gpu_value(using_extended_flags, false);

	set_gpu_value(shift_debt, cuReal3());

	clear_dirichlet_flags();
	clear_halo_flags();

	set_gpu_value(pbc_x, (int)0);
	set_gpu_value(pbc_y, (int)0);
	set_gpu_value(pbc_z, (int)0);
}

//--------------------------------------------MULTIPLE ENTRIES SETTERS - SHAPE CHANGERS

//set value in rectangle (i.e. in cells intersecting the rectangle), where the rectangle is relative to this VEC's rectangle - all cells become non-empty cells irrespective of value set
template <typename VType>
__host__ void cuVEC_VC<VType>::setrect(cuRect rectangle, VType value)
{
	if (!get_gpu_value(cuVEC<VType>::rect).intersects(rectangle + get_gpu_value(cuVEC<VType>::rect).s)) return;

	cuBox cells_box = cuVEC<VType>::box_from_rect_max_cpu(rectangle + get_gpu_value(cuVEC<VType>::rect).s);

	setbox(cells_box, value);
}