#pragma once

#include "VEC_VC.h"

//--------------------------------------------MULTIPLE ENTRIES SETTERS - OTHERS

template <typename VType>
void VEC_VC<VType>::setnonempty(VType value)
{
#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		if (ngbrFlags[idx] & NF_NOTEMPTY) {

			quantity[idx] = value;
		}
	}
}

//set value in non-empty cells only in given rectangle (relative coordinates)
template <typename VType>
void VEC_VC<VType>::setrectnonempty(const Rect& rectangle, VType value)
{
	Box box = box_from_rect_max(rectangle + rect.s);

#pragma omp parallel for
	for (int j = (box.s.y >= 0 ? box.s.y : 0); j < (box.e.y <= n.y ? box.e.y : n.y); j++) {
		for (int k = (box.s.z >= 0 ? box.s.z : 0); k < (box.e.z <= n.z ? box.e.z : n.z); k++) {
			for (int i = (box.s.x >= 0 ? box.s.x : 0); i < (box.e.x <= n.x ? box.e.x : n.x); i++) {

				int idx = i + j * n.x + k * n.x*n.y;

				if (ngbrFlags[idx] & NF_NOTEMPTY) {

					quantity[idx] = value;
				}
			}
		}
	}
}

//re-normalize all non-zero values to have the new magnitude (multiply by new_norm and divide by current magnitude)
template <typename VType>
template <typename PType>
void VEC_VC<VType>::renormalize(PType new_norm)
{
#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		PType curr_norm = GetMagnitude(quantity[idx]);

		if ((ngbrFlags[idx] & NF_NOTEMPTY) && IsNZ(curr_norm)) {

			quantity[idx] *= new_norm / curr_norm;
		}
	}
}

//copy values from copy_this but keep current dimensions - if necessary map values from copy_this to local dimensions; from flags only copy the shape but not the boundary condition values or anything else - these are reset
template <typename VType>
void VEC_VC<VType>::copy_values(const VEC_VC<VType>& copy_this)
{
	//copy values
	VEC<VType>::copy_values(copy_this);

	//copy shape

	//first clear current flags
	ngbrFlags.assign(n.dim(), 0);

	//now map shape from copy_this.ngbrFlags to ngbrFlags
	SZ3 source_n = copy_this.n;
	DBL3 sourceIdx = (DBL3)source_n / n;

#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		int _x = (int)floor((idx % n.x) * sourceIdx.x);
		int _y = (int)floor(((idx / n.x) % n.y) * sourceIdx.y);
		int _z = (int)floor((idx / (n.x*n.y)) * sourceIdx.z);

		if (copy_this.ngbrFlags[_x + _y * source_n.x + _z * (source_n.x*source_n.y)] & NF_NOTEMPTY)
			ngbrFlags[idx] = NF_NOTEMPTY;
	}

	//recalculate neighbor flags
	set_ngbrFlags();
}

//scale all stored values by the given constant
template <typename VType>
void VEC_VC<VType>::scale_values(double constant)
{
#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		quantity[idx] *= constant;
	}
}

//shift all the values in this VEC by the given delta (units same as h)
template <typename VType>
void VEC_VC<VType>::shift(const DBL3& delta, const Rect& shift_rect)
{
	Box shift_box = box_from_rect_min(shift_rect);

	int i_start, j_start, k_start, i_end, j_end, k_end, i_delta, j_delta, k_delta;

	if (delta.x < 0) { i_start = shift_box.s.x; i_end = shift_box.e.x; i_delta = 1; }
	else { i_start = shift_box.e.x - 1; i_end = shift_box.s.x - 1; i_delta = -1; }

	if (delta.y < 0) { j_start = shift_box.s.y; j_end = shift_box.e.y; j_delta = 1; }
	else { j_start = shift_box.e.y - 1; j_end = shift_box.s.y - 1; j_delta = -1; }

	if (delta.z < 0) { k_start = shift_box.s.z; k_end = shift_box.e.z; k_delta = 1; }
	else { k_start = shift_box.e.z - 1; k_end = shift_box.s.z - 1; k_delta = -1; }

	for (int k = k_start; k != k_end; k += k_delta) {
		for (int j = j_start; j != j_end; j += j_delta) {
			for (int i = i_start; i != i_end; i += i_delta) {

				int cell_idx = i + j * n.x + k * n.x*n.y;

				DBL3 position = DBL3(h.x * ((double)i + 0.5), h.y * ((double)j + 0.5), h.z * ((double)k + 0.5)) - delta;

				if (rect.contains(position) && (ngbrFlags[cell_idx] & NF_NOTEMPTY)) {

					quantity[cell_idx] = weighted_average(position, h);
				}
			}
		}
	}
}

//shift all the values in this VEC by the given delta (units same as h)
template <typename VType>
void VEC_VC<VType>::shift_keepmag(const DBL3& delta, const Rect& shift_rect)
{
	Box shift_box = box_from_rect_min(shift_rect);

	int i_start, j_start, k_start, i_end, j_end, k_end, i_delta, j_delta, k_delta;

	if (delta.x < 0) { i_start = shift_box.s.x; i_end = shift_box.e.x; i_delta = 1; }
	else { i_start = shift_box.e.x - 1; i_end = shift_box.s.x - 1; i_delta = -1; }

	if (delta.y < 0) { j_start = shift_box.s.y; j_end = shift_box.e.y; j_delta = 1; }
	else { j_start = shift_box.e.y - 1; j_end = shift_box.s.y - 1; j_delta = -1; }

	if (delta.z < 0) { k_start = shift_box.s.z; k_end = shift_box.e.z; k_delta = 1; }
	else { k_start = shift_box.e.z - 1; k_end = shift_box.s.z - 1; k_delta = -1; }

	for (int k = k_start; k != k_end; k += k_delta) {
		for (int j = j_start; j != j_end; j += j_delta) {
			for (int i = i_start; i != i_end; i += i_delta) {

				int cell_idx = i + j * n.x + k * n.x*n.y;

				DBL3 position = DBL3(h.x * ((double)i + 0.5), h.y * ((double)j + 0.5), h.z * ((double)k + 0.5)) - delta;

				if (rect.contains(position) && (ngbrFlags[cell_idx] & NF_NOTEMPTY)) {

					//for vectorial quantities (e.g. DBL3) the operator * is a scalar product.
					double old_magnitude_squared = (double)(quantity[cell_idx] * quantity[cell_idx]);

					VType new_value = weighted_average(position, h);

					double new_magnitude_squared = (double)(new_value * new_value);

					if (new_magnitude_squared) {

						quantity[cell_idx] = new_value * sqrt(old_magnitude_squared / new_magnitude_squared);
					}
				}
			}
		}
	}
}

//shift all the values in this VEC by the given delta (units same as h)
template <typename VType>
void VEC_VC<VType>::shift_x(double delta, const Rect& shift_rect)
{
	if (fabs(shift_debt.x + delta) < h.x) {

		//total shift not enough : bank it and return
		shift_debt.x += delta;
		return;
	}

	//only shift an integer number of cells : there might be a sub-cellsize remainder so just bank it to be used next time
	int cells_shift = (int)((shift_debt.x + delta) / h.x);
	shift_debt.x -= h.x * cells_shift - delta;

	Box shift_box = box_from_rect_min(shift_rect);

	if (cells_shift < 0) {

		for (int i = shift_box.s.x; i < shift_box.e.x + cells_shift; i++) {
#pragma omp parallel for
			for (int j = shift_box.s.y; j < shift_box.e.y; j++) {
				for (int k = shift_box.s.z; k < shift_box.e.z; k++) {

					int cell_idx = i + j * n.x + k * n.x*n.y;
					int shift_cell_idx = cell_idx - cells_shift;

					if ((ngbrFlags[cell_idx] & NF_NOTEMPTY) && (ngbrFlags[shift_cell_idx] & NF_NOTEMPTY)) {
						
						quantity[cell_idx] = quantity[shift_cell_idx];
					}
				}
			}
		}
	}
	else {

		for (int i = shift_box.e.x - 1; i >= shift_box.s.x + cells_shift; i--) {
#pragma omp parallel for
			for (int j = shift_box.s.y; j < shift_box.e.y; j++) {
				for (int k = shift_box.s.z; k < shift_box.e.z; k++) {

					int cell_idx = i + j * n.x + k * n.x*n.y;
					int shift_cell_idx = cell_idx - cells_shift;

					if ((ngbrFlags[cell_idx] & NF_NOTEMPTY) && (ngbrFlags[shift_cell_idx] & NF_NOTEMPTY)) {

						quantity[cell_idx] = quantity[shift_cell_idx];
					}
				}
			}
		}
	}
}

//shift all the values in this VEC by the given delta (units same as h)
template <typename VType>
void VEC_VC<VType>::shift_y(double delta, const Rect& shift_rect)
{
	if (fabs(shift_debt.y + delta) < h.y) {

		//total shift not enough : bank it and return
		shift_debt.y += delta;
		return;
	}

	//only shift an integer number of cells : there might be a sub-cellsize remainder so just bank it to be used next time
	int cells_shift = (int)((shift_debt.y + delta) / h.y);
	shift_debt.y -= h.y * cells_shift - delta;

	Box shift_box = box_from_rect_min(shift_rect);

	if (cells_shift < 0) {

		for (int j = shift_box.s.y; j < shift_box.e.y + cells_shift; j++) {
#pragma omp parallel for
			for (int i = shift_box.s.x; i < shift_box.e.x; i++) {
				for (int k = shift_box.s.z; k < shift_box.e.z; k++) {

					int cell_idx = i + j * n.x + k * n.x*n.y;
					int shift_cell_idx = cell_idx - cells_shift * n.x;

					if ((ngbrFlags[cell_idx] & NF_NOTEMPTY) && (ngbrFlags[shift_cell_idx] & NF_NOTEMPTY)) {

						quantity[cell_idx] = quantity[shift_cell_idx];
					}
				}
			}
		}
	}
	else {

		for (int j = shift_box.e.y - 1; j >= shift_box.s.y + cells_shift; j--) {
#pragma omp parallel for
			for (int i = shift_box.s.x; i < shift_box.e.x; i++) {
				for (int k = shift_box.s.z; k < shift_box.e.z; k++) {

					int cell_idx = i + j * n.x + k * n.x*n.y;
					int shift_cell_idx = cell_idx - cells_shift * n.x;

					if ((ngbrFlags[cell_idx] & NF_NOTEMPTY) && (ngbrFlags[shift_cell_idx] & NF_NOTEMPTY)) {

						quantity[cell_idx] = quantity[shift_cell_idx];
					}
				}
			}
		}
	}
}

//--------------------------------------------OPERATIONS

template <typename VType>
VType VEC_VC<VType>::average_nonempty(const Box& box) const
{
	VType av = VType();
	int count = 0;

	for (int j = (box.s.y >= 0 ? box.s.y : 0); j < (box.e.y <= n.y ? box.e.y : n.y); j++) {
		for (int k = (box.s.z >= 0 ? box.s.z : 0); k < (box.e.z <= n.z ? box.e.z : n.z); k++) {
			for (int i = (box.s.x >= 0 ? box.s.x : 0); i < (box.e.x <= n.x ? box.e.x : n.x); i++) {

				int idx = i + j * n.x + k * n.x*n.y;

				//get running average
				if (ngbrFlags[idx] & NF_NOTEMPTY) {

					count++;
					av = (av * (count - 1) + quantity[idx]) / count;
				}
			}
		}
	}

	return av;
}

//average over non-empty cells intersecting given rectangle (relative to this VEC's rect)
template <typename VType>
VType VEC_VC<VType>::average_nonempty(const Rect& rectangle) const
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average_nonempty(Box(n));

	//... otherwise rectangle must intersect with this mesh
	if (!rect.intersects(rectangle + rect.s)) return VType();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return average_nonempty(box_from_rect_max(rectangle + rect.s));
}

template <typename VType>
VType VEC_VC<VType>::average_nonempty_omp(const Box& box) const
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

		//you shouldn't call this with a box that is not contained in Box(n), but this stops any bad memory accesses
		if (idx < 0 || idx >= n.dim()) continue;

		int tn = omp_get_thread_num();

		//only include non-empty cells
		if (ngbrFlags[idx] & NF_NOTEMPTY) {

			reduction.reduce_average(quantity[idx]);
		}
	}

	return reduction.average();;
}

//average over non-empty cells intersecting given rectangle (relative to this VEC's rect)
template <typename VType>
VType VEC_VC<VType>::average_nonempty_omp(const Rect& rectangle) const
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average_nonempty_omp(Box(n));

	//... otherwise rectangle must intersect with this mesh
	if (!rect.intersects(rectangle + rect.s)) return VType();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return average_nonempty_omp(box_from_rect_max(rectangle + rect.s));
}

