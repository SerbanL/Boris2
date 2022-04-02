#pragma once

#include "VEC.h"

//--------------------------------------------AVERAGING OPERATIONS

template <typename VType>
VType VEC<VType>::average(const Box& box) const
{
	VType av = VType();
	int count = 0;

	for (int i = (box.s.x >= 0 ? box.s.x : 0); i < (box.e.x <= n.x ? box.e.x : n.x); i++) {
		for (int j = (box.s.y >= 0 ? box.s.y : 0); j < (box.e.y <= n.y ? box.e.y : n.y); j++) {
			for (int k = (box.s.z >= 0 ? box.s.z : 0); k < (box.e.z <= n.z ? box.e.z : n.z); k++) {

				//get running average
				count++;
				av = (av * (count - 1) + quantity[i + j * n.x + k * n.x*n.y]) / count;
			}
		}
	}

	return av;
}

//average over given rectangle (relative to this VEC's rect)
template <typename VType>
VType VEC<VType>::average(const Rect& rectangle) const
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average(Box(n));

	//if rect start and end point are the same, then just read single value
	if (rectangle.s == rectangle.e && rect.contains(rectangle.s)) {

		return (*this)[rectangle.s];
	}

	//... otherwise rectangle must intersect with this mesh
	if (!rect.intersects(rectangle + rect.s)) return VType();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return average(box_from_rect_max(rectangle + rect.s));
}

template <typename VType>
VType VEC<VType>::average_omp(const Box& box) const
{
	reduction.new_average_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * n.x + (k + box.s.k) * n.x * n.y;

		if (idx < 0 || idx >= n.dim()) continue;

		reduction.reduce_average(quantity[idx]);
	}

	return reduction.average();
}

template <typename VType>
VType VEC<VType>::average_omp(const Rect& rectangle) const
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average_omp(Box(n));

	//... otherwise rectangle must intersect with this mesh
	if (!rect.intersects(rectangle + rect.s)) return VType();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return average_omp(box_from_rect_max(rectangle + rect.s));
}

template <typename VType>
VType VEC<VType>::average_nonempty(const Box& box) const
{
	VType av = VType();
	int count = 0;

	for (int i = (box.s.x >= 0 ? box.s.x : 0); i < (box.e.x <= n.x ? box.e.x : n.x); i++) {
		for (int j = (box.s.y >= 0 ? box.s.y : 0); j < (box.e.y <= n.y ? box.e.y : n.y); j++) {
			for (int k = (box.s.z >= 0 ? box.s.z : 0); k < (box.e.z <= n.z ? box.e.z : n.z); k++) {

				//get running average
				VType value = quantity[i + j * n.x + k * n.x*n.y];

				if (IsNZ(GetMagnitude(value))) {

					count++;
					av = (av * (count - 1) + value) / count;
				}
			}
		}
	}

	return av;
}

template <typename VType>
VType VEC<VType>::average_nonempty(const Rect& rectangle) const
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average_nonempty(Box(n));

	//... otherwise rectangle must intersect with this mesh
	if (!rect.intersects(rectangle + rect.s)) return VType();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return average_nonempty(box_from_rect_max(rectangle + rect.s));
}

template <typename VType>
VType VEC<VType>::average_nonempty_omp(const Box& box) const
{
	reduction.new_average_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * n.x + (k + box.s.k) * n.x * n.y;

		if (idx < 0 || idx >= n.dim()) continue;

		VType value = quantity[idx];

		//only include non-empty cells
		if (IsNZ(GetMagnitude(value))) {

			reduction.reduce_average(value);
		}
	}

	return reduction.average();
}

template <typename VType>
VType VEC<VType>::average_nonempty_omp(const Rect& rectangle) const
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average_nonempty_omp(Box(n));

	//... otherwise rectangle must intersect with this mesh
	if (!rect.intersects(rectangle + rect.s)) return VType();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return average_nonempty_omp(box_from_rect_max(rectangle + rect.s));
}

template <typename VType>
VType VEC<VType>::weighted_average(const DBL3& coord, const DBL3& stencil) const
{
	//smoothedValue must be initialized to zero
	VType smoothedValue = VType();

	//1. find indexes of lower-left and uper-right cells intersecting with the stencil (i.e. find cells which contain its lower-left corner and upper-right corner)
	//Note, it could be a single cell containing both corners : either the stencil is strictly included in it or coincides with it

	//positions of lower-left and uper-right corners of cell relative to containing mesh
	DBL3 pos_ll = coord - stencil / 2;
	DBL3 pos_ur = pos_ll + stencil;

	INT3 idx_ll = floor(pos_ll / h);
	INT3 idx_ur = ceil(pos_ur / h);

	//2. obtain weighted average value from cells intersecting with the stencil (the larger the distance of the cell centre from coord, the lower the weight)
	double d_max = GetMagnitude(stencil / 2 + h / 2);						//this is the maximum possible distance (a cell at this distance will get a weight of zero)	
	double d_total = 0;														//total reciprocal distance

	for (int ii = (idx_ll.i >= 0 ? idx_ll.i : 0); ii < (idx_ur.i < n.x ? idx_ur.i : n.x); ii++) {
		for (int jj = (idx_ll.j >= 0 ? idx_ll.j : 0); jj < (idx_ur.j < n.y ? idx_ur.j : n.y); jj++) {
			for (int kk = (idx_ll.k >= 0 ? idx_ll.k : 0); kk < (idx_ur.k < n.z ? idx_ur.k : n.z); kk++) {

				//find reciprocal distance for each cell : this is its weight * total reciprocal distance (to divide by d_total at the end)
				double d_recip = d_max - get_distance(coord, DBL3((ii + 0.5)*h.x, (jj + 0.5)*h.y, (kk + 0.5)*h.z));
				d_total += d_recip;

				smoothedValue += d_recip * quantity[ii + jj * n.x + kk * n.x*n.y];
			}
		}
	}

	//finally divide by total reciprocal distance to finish off the averaging
	if (d_total) smoothedValue /= d_total;

	return smoothedValue;
}

template <typename VType>
VType VEC<VType>::weighted_average(const INT3& ijk, const DBL3& cs) const
{
	if (cs == h) return quantity[ijk.i + ijk.j * n.x + ijk.k * n.x * n.y];

	//smoothedValue must be initialized to zero
	VType smoothedValue = VType();

	//1. find indexes of lower-left and uper-right cells intersecting with the stencil (i.e. find cells which contain its lower-left corner and upper-right corner)
	//Note, it could be a single cell containing both corners : either the stencil is strictly included in it or coincides with it

	//positions of lower-left and uper-right corners of cell relative to containing mesh
	DBL3 pos_ll = cs & ijk;
	DBL3 pos_ur = pos_ll + cs;

	INT3 idx_ll = floor(pos_ll / h);
	INT3 idx_ur = ceil(pos_ur / h);

	DBL3 coord = pos_ll + (cs / 2);

	//2. obtain weighted average value from cells intersecting with the stencil (the larger the distance of the cell centre from coord, the lower the weight)
	double d_max = GetMagnitude(cs / 2 + h / 2);						//this is the maximum possible distance (a cell at this distance will get a weight of zero)	
	double d_total = 0;														//total reciprocal distance

	for (int ii = (idx_ll.i >= 0 ? idx_ll.i : 0); ii < (idx_ur.i < n.x ? idx_ur.i : n.x); ii++) {
		for (int jj = (idx_ll.j >= 0 ? idx_ll.j : 0); jj < (idx_ur.j < n.y ? idx_ur.j : n.y); jj++) {
			for (int kk = (idx_ll.k >= 0 ? idx_ll.k : 0); kk < (idx_ur.k < n.z ? idx_ur.k : n.z); kk++) {

				//find reciprocal distance for each cell : this is its weight * total reciprocal distance (to divide by d_total at the end)
				double d_recip = d_max - get_distance(coord, DBL3((ii + 0.5)*h.x, (jj + 0.5)*h.y, (kk + 0.5)*h.z));
				d_total += d_recip;

				smoothedValue += d_recip * quantity[ii + jj * n.x + kk * n.x*n.y];
			}
		}
	}

	//finally divide by total reciprocal distance to finish off the averaging
	if (d_total) smoothedValue /= d_total;

	return smoothedValue;
}