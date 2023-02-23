#pragma once

#include "VEC_VC.h"

//--------------------------------------------AVERAGING OPERATIONS

template <typename VType>
VType VEC_VC<VType>::average_nonempty(const Box& box) const
{
	VType av = VType();
	int count = 0;

	for (int j = (box.s.y >= 0 ? box.s.y : 0); j < (box.e.y <= VEC<VType>::n.y ? box.e.y : VEC<VType>::n.y); j++) {
		for (int k = (box.s.z >= 0 ? box.s.z : 0); k < (box.e.z <= VEC<VType>::n.z ? box.e.z : VEC<VType>::n.z); k++) {
			for (int i = (box.s.x >= 0 ? box.s.x : 0); i < (box.e.x <= VEC<VType>::n.x ? box.e.x : VEC<VType>::n.x); i++) {

				int idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

				//get running average
				if (ngbrFlags[idx] & NF_NOTEMPTY) {

					count++;
					av = (av * (count - 1) + VEC<VType>::quantity[idx]) / count;
				}
			}
		}
	}

	return av;
}

//average over non-empty cells intersecting given rectangle (relative to this VEC's VEC<VType>::rect)
template <typename VType>
VType VEC_VC<VType>::average_nonempty(const Rect& rectangle) const
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average_nonempty(Box(VEC<VType>::n));

	//if rect start and end point are the same, then just read single value
	if (rectangle.s == rectangle.e && VEC<VType>::rect.contains(rectangle.s)) {

		return (*this)[rectangle.s];
	}

	//... otherwise rectangle must intersect with this mesh
	if (!VEC<VType>::rect.intersects(rectangle + VEC<VType>::rect.s)) return VType();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return average_nonempty(VEC<VType>::box_from_rect_max(rectangle + VEC<VType>::rect.s));
}

template <typename VType>
VType VEC_VC<VType>::average_nonempty_omp(const Box& box) const
{
	VEC<VType>::reduction.new_average_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * VEC<VType>::n.x + (k + box.s.k) * VEC<VType>::n.x * VEC<VType>::n.y;

		//you shouldn't call this with a box that is not contained in Box(VEC<VType>::n), but this stops any bad memory accesses
		if (idx < 0 || idx >= VEC<VType>::n.dim()) continue;

		//only include non-empty cells
		if (ngbrFlags[idx] & NF_NOTEMPTY) {

			VEC<VType>::reduction.reduce_average(VEC<VType>::quantity[idx]);
		}
	}

	return VEC<VType>::reduction.average();
}

//average over non-empty cells intersecting given rectangle (relative to this VEC's VEC<VType>::rect)
template <typename VType>
VType VEC_VC<VType>::average_nonempty_omp(const Rect& rectangle) const
{
	//if empty rectangle then average over the entire mesh
	if (rectangle.IsNull()) return average_nonempty_omp(Box(VEC<VType>::n));

	//... otherwise rectangle must intersect with this mesh
	if (!VEC<VType>::rect.intersects(rectangle + VEC<VType>::rect.s)) return VType();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return average_nonempty_omp(VEC<VType>::box_from_rect_max(rectangle + VEC<VType>::rect.s));
}

//SUMMING FUNCTIONS

template <typename VType>
VType VEC_VC<VType>::sum_nonempty_omp(const Box& box) const
{
	VEC<VType>::reduction.new_sum_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * VEC<VType>::n.x + (k + box.s.k) * VEC<VType>::n.x * VEC<VType>::n.y;

		//you shouldn't call this with a box that is not contained in Box(VEC<VType>::n), but this stops any bad memory accesses
		if (idx < 0 || idx >= VEC<VType>::n.dim()) continue;

		//only include non-empty cells
		if (ngbrFlags[idx] & NF_NOTEMPTY) {

			VEC<VType>::reduction.reduce_sum(VEC<VType>::quantity[idx]);
		}
	}

	return VEC<VType>::reduction.sum();
}

//summ all non-empty cells intersecting given rectangle (relative to this VEC's VEC<VType>::rect)
template <typename VType>
VType VEC_VC<VType>::sum_nonempty_omp(const Rect& rectangle) const
{
	//if empty rectangle then average over the entire mesh
	if (rectangle.IsNull()) return sum_nonempty_omp(Box(VEC<VType>::n));

	//... otherwise rectangle must intersect with this mesh
	if (!VEC<VType>::rect.intersects(rectangle + VEC<VType>::rect.s)) return VType();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return sum_nonempty_omp(VEC<VType>::box_from_rect_max(rectangle + VEC<VType>::rect.s));
}
