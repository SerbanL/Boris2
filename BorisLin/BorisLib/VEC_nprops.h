#pragma once

#include "VEC.h"

//--------------------------------------------GET MIN-MAX by MAGNITUDE

template <typename VType>
template <typename PType>
VAL2<PType> VEC<VType>::get_minmax(const Box& box) const
{
	magnitude_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * n.x + (k + box.s.k) * n.x * n.y;

		if (idx < 0 || idx >= n.dim()) continue;

		magnitude_reduction.reduce_minmax(GetMagnitude(quantity[idx]));
	}

	return magnitude_reduction.minmax();
}

template <typename VType>
template <typename PType>
VAL2<PType> VEC<VType>::get_minmax(const Rect& rectangle) const
{
	//if empty rectangle then use entire mesh
	if (rectangle.IsNull()) return get_minmax(Box(n));

	//... otherwise rectangle must intersect with this mesh
	if (!rect.intersects(rectangle + rect.s)) return VAL2<PType>();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return get_minmax(box_from_rect_max(rectangle + rect.s));
}

//--------------------------------------------GET MIN-MAX COMPONENT X

template <typename VType>
template <typename PType>
VAL2<PType> VEC<VType>::get_minmax_component_x(const Box& box) const
{
	magnitude_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * n.x + (k + box.s.k) * n.x * n.y;

		if (idx < 0 || idx >= n.dim()) continue;

		magnitude_reduction.reduce_minmax(quantity[idx].x);
	}

	return magnitude_reduction.minmax();
}

template <typename VType>
template <typename PType>
VAL2<PType> VEC<VType>::get_minmax_component_x(const Rect& rectangle) const
{
	//if empty rectangle then use entire mesh
	if (rectangle.IsNull()) return get_minmax_component_x(Box(n));

	//... otherwise rectangle must intersect with this mesh
	if (!rect.intersects(rectangle + rect.s)) return VAL2<PType>();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return get_minmax_component_x(box_from_rect_max(rectangle + rect.s));
}

//--------------------------------------------GET MIN-MAX COMPONENT Y

template <typename VType>
template <typename PType>
VAL2<PType> VEC<VType>::get_minmax_component_y(const Box& box) const
{
	magnitude_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * n.x + (k + box.s.k) * n.x * n.y;

		if (idx < 0 || idx >= n.dim()) continue;

		magnitude_reduction.reduce_minmax(quantity[idx].y);
	}

	return magnitude_reduction.minmax();
}

template <typename VType>
template <typename PType>
VAL2<PType> VEC<VType>::get_minmax_component_y(const Rect& rectangle) const
{
	//if empty rectangle then use entire mesh
	if (rectangle.IsNull()) return get_minmax_component_y(Box(n));

	//... otherwise rectangle must intersect with this mesh
	if (!rect.intersects(rectangle + rect.s)) return VAL2<PType>();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return get_minmax_component_y(box_from_rect_max(rectangle + rect.s));
}

//--------------------------------------------GET MIN-MAX COMPONENT Z

template <typename VType>
template <typename PType>
VAL2<PType> VEC<VType>::get_minmax_component_z(const Box& box) const
{
	magnitude_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * n.x + (k + box.s.k) * n.x * n.y;

		if (idx < 0 || idx >= n.dim()) continue;

		magnitude_reduction.reduce_minmax(quantity[idx].z);
	}

	return magnitude_reduction.minmax();
}

template <typename VType>
template <typename PType>
VAL2<PType> VEC<VType>::get_minmax_component_z(const Rect& rectangle) const
{
	//if empty rectangle then use entire mesh
	if (rectangle.IsNull()) return get_minmax_component_z(Box(n));

	//... otherwise rectangle must intersect with this mesh
	if (!rect.intersects(rectangle + rect.s)) return VAL2<PType>();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return get_minmax_component_z(box_from_rect_max(rectangle + rect.s));
}