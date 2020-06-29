#pragma once

#include "VEC_VC.h"

//--------------------------------------------GET MIN-MAX by MAGNITUDE

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax(const Box& box) const
{
	VEC<VType>::magnitude_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * VEC<VType>::n.x + (k + box.s.k) * VEC<VType>::n.x * VEC<VType>::n.y;

		if (idx < 0 || idx >= VEC<VType>::n.dim() || !(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

		VEC<VType>::magnitude_reduction.reduce_minmax(GetMagnitude(VEC<VType>::quantity[idx]));
	}

	return VEC<VType>::magnitude_reduction.minmax();
}

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax(const Rect& rectangle) const
{
	//if empty rectangle then use entire mesh
	if (rectangle.IsNull()) return get_minmax(Box(VEC<VType>::n));

	//... otherwise rectangle must intersect with this mesh
	if (!VEC<VType>::rect.intersects(rectangle + VEC<VType>::rect.s)) return VAL2<PType>();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return get_minmax(VEC<VType>::box_from_rect_max(rectangle + VEC<VType>::rect.s));
}

//--------------------------------------------GET MIN-MAX COMPONENT X

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax_component_x(const Box& box) const
{
	VEC<VType>::magnitude_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * VEC<VType>::n.x + (k + box.s.k) * VEC<VType>::n.x * VEC<VType>::n.y;

		if (idx < 0 || idx >= VEC<VType>::n.dim() || !(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

		VEC<VType>::magnitude_reduction.reduce_minmax(VEC<VType>::quantity[idx].x);
	}

	return VEC<VType>::magnitude_reduction.minmax();
}

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax_component_x(const Rect& rectangle) const
{
	//if empty rectangle then use entire mesh
	if (rectangle.IsNull()) return get_minmax_component_x(Box(VEC<VType>::n));

	//... otherwise rectangle must intersect with this mesh
	if (!VEC<VType>::rect.intersects(rectangle + VEC<VType>::rect.s)) return VAL2<PType>();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return get_minmax_component_x(VEC<VType>::box_from_rect_max(rectangle + VEC<VType>::rect.s));
}

//--------------------------------------------GET MIN-MAX COMPONENT Y

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax_component_y(const Box& box) const
{
	VEC<VType>::magnitude_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * VEC<VType>::n.x + (k + box.s.k) * VEC<VType>::n.x * VEC<VType>::n.y;

		if (idx < 0 || idx >= VEC<VType>::n.dim() || !(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

		VEC<VType>::magnitude_reduction.reduce_minmax(VEC<VType>::quantity[idx].y);
	}

	return VEC<VType>::magnitude_reduction.minmax();
}

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax_component_y(const Rect& rectangle) const
{
	//if empty rectangle then use entire mesh
	if (rectangle.IsNull()) return get_minmax_component_y(Box(VEC<VType>::n));

	//... otherwise rectangle must intersect with this mesh
	if (!VEC<VType>::rect.intersects(rectangle + VEC<VType>::rect.s)) return VAL2<PType>();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return get_minmax_component_y(VEC<VType>::box_from_rect_max(rectangle + VEC<VType>::rect.s));
}

//--------------------------------------------GET MIN-MAX COMPONENT Z

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax_component_z(const Box& box) const
{
	VEC<VType>::magnitude_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * VEC<VType>::n.x + (k + box.s.k) * VEC<VType>::n.x * VEC<VType>::n.y;

		if (idx < 0 || idx >= VEC<VType>::n.dim() || !(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

		VEC<VType>::magnitude_reduction.reduce_minmax(VEC<VType>::quantity[idx].z);
	}

	return VEC<VType>::magnitude_reduction.minmax();
}

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax_component_z(const Rect& rectangle) const
{
	//if empty rectangle then use entire mesh
	if (rectangle.IsNull()) return get_minmax_component_z(Box(VEC<VType>::n));

	//... otherwise rectangle must intersect with this mesh
	if (!VEC<VType>::rect.intersects(rectangle + VEC<VType>::rect.s)) return VAL2<PType>();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return get_minmax_component_z(VEC<VType>::box_from_rect_max(rectangle + VEC<VType>::rect.s));
}