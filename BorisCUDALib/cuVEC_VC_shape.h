#pragma once

#include "cuVEC_VC.h"

//--------------------------------------------SIZING

//resize and set shape using linked vec
template <typename VType>
template <typename LVType>
__host__ bool cuVEC_VC<VType>::resize(cuSZ3 new_n, cuVEC_VC<LVType>& linked_vec)
{
	//memory reserved so map flags to new size before changing the n value
	cudaError_t error = resize_ngbrFlags(new_n);
	if (error != cudaSuccess) return false;

	if (!cuVEC::resize(new_n)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec.n, linked_vec.h, linked_vec.rect, linked_vec.ngbrFlags);

	return true;
}

//resize but keep shape
template <typename VType>
__host__ bool cuVEC_VC<VType>::resize(cuSZ3 new_n)
{
	//memory reserved so map flags to new size before changing the n value
	cudaError_t error = resize_ngbrFlags(new_n);
	if (error != cudaSuccess) return false;

	if (!cuVEC::resize(new_n)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags();

	return true;
}

//resize and set shape using linked vec
template <typename VType>
template <typename LVType>
__host__ bool cuVEC_VC<VType>::resize(cuReal3 new_h, cuRect new_rect, cuVEC_VC<LVType>& linked_vec)
{
	cuSZ3 new_n = get_n_from_h_and_rect(new_h, new_rect);

	//memory reserved so map flags to new size before changing the n value
	cudaError_t error = resize_ngbrFlags(new_n);
	if (error != cudaSuccess) return false;

	if (!cuVEC::resize(new_h, new_rect)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec.n, linked_vec.h, linked_vec.rect, linked_vec.ngbrFlags_ref());

	return true;
}

//resize but keep shape
template <typename VType>
__host__ bool cuVEC_VC<VType>::resize(cuReal3 new_h, cuRect new_rect)
{
	cuSZ3 new_n = get_n_from_h_and_rect(new_h, new_rect);

	//memory reserved so map flags to new size before changing the n value
	cudaError_t error = resize_ngbrFlags(new_n);
	if (error != cudaSuccess) return false;

	if (!cuVEC::resize(new_h, new_rect)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags();

	return true;
}

//set value and shape from linked vec
template <typename VType>
template <typename LVType>
__host__ bool cuVEC_VC<VType>::assign(cuSZ3 new_n, VType value, cuVEC_VC<LVType>& linked_vec)
{
	//memory reserved so map flags to new size before changing the n value
	cudaError_t error = resize_ngbrFlags(new_n);
	if (error != cudaSuccess) return false;

	if (!cuVEC::assign(new_n, value)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec.n, linked_vec.h, linked_vec.rect, linked_vec.ngbrFlags);

	return true;
}

//set value but keep shape - empty cells will retain zero value : i.e. set value everywhere but in empty cells
template <typename VType>
__host__ bool cuVEC_VC<VType>::assign(cuSZ3 new_n, VType value)
{
	//memory reserved so map flags to new size before changing the n value
	cudaError_t error = resize_ngbrFlags(new_n);
	if (error != cudaSuccess) return false;

	if (!cuVEC::assign(new_n, value)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags();

	return true;
}

//set value and shape from linked vec
template <typename VType>
template <typename LVType>
__host__ bool cuVEC_VC<VType>::assign(cuReal3 new_h, cuRect new_rect, VType value, cuVEC_VC<LVType>& linked_vec)
{
	cuSZ3 new_n = get_n_from_h_and_rect(new_h, new_rect);

	//memory reserved so map flags to new size before changing the n value
	cudaError_t error = resize_ngbrFlags(new_n);
	if (error != cudaSuccess) return false;

	if (!cuVEC::assign(new_h, new_rect, value)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec.n, linked_vec.h, linked_vec.rect, linked_vec.ngbrFlags_ref());

	return true;
}

//set value but keep shape - empty cells will retain zero value : i.e. set value everywhere but in empty cells
template <typename VType>
__host__ bool cuVEC_VC<VType>::assign(cuReal3 new_h, cuRect new_rect, VType value)
{
	cuSZ3 new_n = get_n_from_h_and_rect(new_h, new_rect);

	//memory reserved so map flags to new size before changing the n value
	cudaError_t error = resize_ngbrFlags(new_n);
	if (error != cudaSuccess) return false;

	if (!cuVEC::assign(new_h, new_rect, value)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags();

	return true;
}

template <typename VType>
__host__ void cuVEC_VC<VType>::clear(void)
{
	cuVEC::clear();
	set_ngbrFlags_size(0);

	set_gpu_value(shift_debt, cuReal3());

	clear_dirichlet_flags();

	set_gpu_value(aSOR_damping, (cuReal)1.0);
}

//--------------------------------------------MULTIPLE ENTRIES SETTERS - SHAPE CHANGERS

//set value in rectangle (i.e. in cells intersecting the rectangle), where the rectangle is relative to this VEC's rectangle - all cells become non-empty cells irrespective of value set
template <typename VType>
__host__ void cuVEC_VC<VType>::setrect(cuRect rectangle, VType value)
{
	if (!get_gpu_value(rect).intersects(rectangle + get_gpu_value(rect).s)) return;

	cuBox cells_box = box_from_rect_max_cpu(rectangle + get_gpu_value(rect).s);

	setbox(cells_box, value);
}