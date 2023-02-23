#pragma once

#include "VEC_VC.h"

//--------------------------------------------MULTIPLE ENTRIES SETTERS - OTHERS

template <typename VType>
void VEC_VC<VType>::setnonempty(VType value)
{
#pragma omp parallel for
	for (int idx = 0; idx < VEC<VType>::n.dim(); idx++) {

		if (ngbrFlags[idx] & NF_NOTEMPTY) {

			VEC<VType>::quantity[idx] = value;
		}
	}
}

//set value in non-empty cells only in given rectangle (relative coordinates)
template <typename VType>
void VEC_VC<VType>::setrectnonempty(const Rect& rectangle, VType value)
{
	Box box = VEC<VType>::box_from_rect_max(rectangle + VEC<VType>::rect.s);

#pragma omp parallel for
	for (int j = (box.s.y >= 0 ? box.s.y : 0); j < (box.e.y <= VEC<VType>::n.y ? box.e.y : VEC<VType>::n.y); j++) {
		for (int k = (box.s.z >= 0 ? box.s.z : 0); k < (box.e.z <= VEC<VType>::n.z ? box.e.z : VEC<VType>::n.z); k++) {
			for (int i = (box.s.x >= 0 ? box.s.x : 0); i < (box.e.x <= VEC<VType>::n.x ? box.e.x : VEC<VType>::n.x); i++) {

				int idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

				if (ngbrFlags[idx] & NF_NOTEMPTY) {

					VEC<VType>::quantity[idx] = value;
				}
			}
		}
	}
}

//set value in solid object only containing relpos
template <typename VType>
void VEC_VC<VType>::setobject(VType value, DBL3 relpos)
{
	//use a simple serial algorithm for filling a solid object
	//it's fast enough even for very large meshes (almost instantaneous in terms of response to user input) and since it's not meant to be used at runtime a parallel algorithm is not really necessary

	int start_idx = VEC<VType>::position_to_cellidx(relpos);
	if (is_empty(start_idx)) return;

	//allocate memory in chunks when needed, to avoid having to allocate memory too often
	const int memory_chunk = 10000;

	//keep track of marked cells here so we can restore ngbrFlags state at the end
	std::vector<int> marked_cells(memory_chunk);
	int num_marked_cells = 0;

	std::vector<int> array1, array2(memory_chunk);
	int num_previous_cells = 0;
	int num_current_cells = 0;

	VEC<VType>::quantity[start_idx] = value;
	ngbrFlags[start_idx] &= ~NF_NOTEMPTY;
	marked_cells[num_marked_cells++] = start_idx;
	array2[num_current_cells++] = start_idx;

	while (num_current_cells > 0) {

		num_previous_cells = num_current_cells;

		//reset current cells to zero before going through previous cells so we can recount number of current cells
		num_current_cells = 0;

		//array1 must always have the list of previously marked cells, and array2 will store the marked cells for next iteration, so swap their storage every time
		array1.swap(array2);

		//go through cells marked last time and find new cells to mark
		//NOTE: this could be parallelized in theory by working with sets of non-interacting sub-grids (and move memory allocation outside of for loop)
		//e.g. rather than getting indexes from array1 only, you could use several of them, each storing non-interacting cell indexes
		for (int idx = 0; idx < num_previous_cells; idx++) {

			//make sure marked_cells and array for next iteration has enough memory allocated
			if (num_marked_cells + 6 > marked_cells.size()) marked_cells.resize(marked_cells.size() + memory_chunk);
			if (num_current_cells + 6 > array2.size()) array2.resize(array2.size() + memory_chunk);

			//cell index for which we consider its neighbors
			int pidx = array1[idx];

			//must be a neighbor, and not already marked (NF_NOTEMPTY)
			if ((ngbrFlags[pidx] & NF_NPX) && (ngbrFlags[pidx + 1] & NF_NOTEMPTY)) {

				//set value and mark it as already marked
				VEC<VType>::quantity[pidx + 1] = value;
				ngbrFlags[pidx + 1] &= ~NF_NOTEMPTY;
				//save it so we can restore the NF_NOTEMPTY flags (all cells in which we set value are not actually empty)
				marked_cells[num_marked_cells++] = pidx + 1;
				//save it in number of cells to use on next iteration, which will be used to generate new neighbors
				array2[num_current_cells++] = pidx + 1;
			}

			if ((ngbrFlags[pidx] & NF_NNX) && (ngbrFlags[pidx - 1] & NF_NOTEMPTY)) {

				VEC<VType>::quantity[pidx - 1] = value;
				ngbrFlags[pidx - 1] &= ~NF_NOTEMPTY;
				marked_cells[num_marked_cells++] = pidx - 1;
				array2[num_current_cells++] = pidx - 1;
			}

			if ((ngbrFlags[pidx] & NF_NPY) && (ngbrFlags[pidx + VEC<VType>::n.x] & NF_NOTEMPTY)) {

				VEC<VType>::quantity[pidx + VEC<VType>::n.x] = value;
				ngbrFlags[pidx + VEC<VType>::n.x] &= ~NF_NOTEMPTY;
				marked_cells[num_marked_cells++] = pidx + VEC<VType>::n.x;
				array2[num_current_cells++] = pidx + VEC<VType>::n.x;
			}

			if ((ngbrFlags[pidx] & NF_NNY) && (ngbrFlags[pidx - VEC<VType>::n.x] & NF_NOTEMPTY)) {

				VEC<VType>::quantity[pidx - VEC<VType>::n.x] = value;
				ngbrFlags[pidx - VEC<VType>::n.x] &= ~NF_NOTEMPTY;
				marked_cells[num_marked_cells++] = pidx - VEC<VType>::n.x;
				array2[num_current_cells++] = pidx - VEC<VType>::n.x;
			}

			if ((ngbrFlags[pidx] & NF_NPZ) && (ngbrFlags[pidx + VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY)) {

				VEC<VType>::quantity[pidx + VEC<VType>::n.x*VEC<VType>::n.y] = value;
				ngbrFlags[pidx + VEC<VType>::n.x*VEC<VType>::n.y] &= ~NF_NOTEMPTY;
				marked_cells[num_marked_cells++] = pidx + VEC<VType>::n.x*VEC<VType>::n.y;
				array2[num_current_cells++] = pidx + VEC<VType>::n.x*VEC<VType>::n.y;
			}

			if ((ngbrFlags[pidx] & NF_NNZ) && (ngbrFlags[pidx - VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY)) {

				VEC<VType>::quantity[pidx - VEC<VType>::n.x*VEC<VType>::n.y] = value;
				ngbrFlags[pidx - VEC<VType>::n.x*VEC<VType>::n.y] &= ~NF_NOTEMPTY;
				marked_cells[num_marked_cells++] = pidx - VEC<VType>::n.x*VEC<VType>::n.y;
				array2[num_current_cells++] = pidx - VEC<VType>::n.x*VEC<VType>::n.y;
			}
		}
	} 
	
	//restore cells marked with NF_NOTEMPTY
#pragma omp parallel for
	for (int idx = 0; idx < num_marked_cells; idx++) {

		ngbrFlags[marked_cells[idx]] |= NF_NOTEMPTY;
	}
}

//re-normalize all non-zero values to have the new magnitude (multiply by new_norm and divide by current magnitude)
template <typename VType>
template <typename PType>
void VEC_VC<VType>::renormalize(PType new_norm)
{
#pragma omp parallel for
	for (int idx = 0; idx < VEC<VType>::n.dim(); idx++) {

		PType curr_norm = GetMagnitude(VEC<VType>::quantity[idx]);

		if ((ngbrFlags[idx] & NF_NOTEMPTY) && IsNZ(curr_norm)) {

			VEC<VType>::quantity[idx] *= new_norm / curr_norm;
		}
	}
}

//copy values from copy_this but keep current dimensions - if necessary map values from copy_this to local dimensions; from flags only copy the shape but not the boundary condition values or anything else - these are reset
template <typename VType>
void VEC_VC<VType>::copy_values(const VEC_VC<VType>& copy_this, Rect dstRect, Rect srcRect)
{
	//copy values
	VEC<VType>::copy_values(copy_this, dstRect, srcRect);

	//copy shape

	if (dstRect.IsNull()) dstRect = VEC<VType>::rect - VEC<VType>::rect.s;
	if (srcRect.IsNull()) srcRect = copy_this.rect - copy_this.rect.s;

	Box cells_box_dst = VEC<VType>::box_from_rect_max(dstRect + VEC<VType>::rect.s);
	Box cells_box_src = copy_this.box_from_rect_max(srcRect + copy_this.rect.s);

	SZ3 dst_n = cells_box_dst.size();
	SZ3 src_n = cells_box_src.size();

	DBL3 sourceIdx = (DBL3)src_n / dst_n;

	//now map shape from copy_this.ngbrFlags to ngbrFlags

#pragma omp parallel for
	for (int j = 0; j < dst_n.j; j++) {
		for (int k = 0; k < dst_n.k; k++) {
			for (int i = 0; i < dst_n.i; i++) {

				int idx_box_dst = i + j * dst_n.x + k * dst_n.x*dst_n.y;

				int _x = (int)floor(i * sourceIdx.x);
				int _y = (int)floor(j * sourceIdx.y);
				int _z = (int)floor(k * sourceIdx.z);

				int idx_out = (i + cells_box_dst.s.i) + (j + cells_box_dst.s.j) * VEC<VType>::n.x + (k + cells_box_dst.s.k) * VEC<VType>::n.x*VEC<VType>::n.y;
				int idx_in = (_x + cells_box_src.s.i) + (_y + cells_box_src.s.j) * copy_this.n.x + (_z + cells_box_src.s.k) * (copy_this.n.x*copy_this.n.y);

				if (idx_out < VEC<VType>::n.dim() && idx_in < copy_this.n.dim()) {

					if (copy_this.is_empty(idx_in)) mark_empty(idx_out);
					else mark_not_empty(idx_out);
				}
			}
		}
	}

	//recalculate neighbor flags
	set_ngbrFlags();
}

//copy values from copy_this but keep current dimensions - if necessary map values from copy_this to local dimensions. Points with zero values are set as empty.
template <typename VType>
void VEC_VC<VType>::copy_values(const VEC<VType>& copy_this, Rect dstRect, Rect srcRect)
{
	VEC<VType>::copy_values(copy_this, dstRect, srcRect);

	//mark zero points as empty before setting flags; otherwise they are not empty -> sets shape
#pragma omp parallel for
	for (int idx = 0; idx < VEC<VType>::n.dim(); idx++) {

		if (VEC<VType>::is_empty(idx)) mark_empty(idx);
		else mark_not_empty(idx);
	}

	//recalculate neighbor flags : current shape is maintained.
	set_ngbrFlags();
}

//shift all the values in this VEC by the given delta (units same as VEC<VType>::h)
template <typename VType>
void VEC_VC<VType>::shift(const DBL3& delta, const Rect& shift_rect)
{
	Box shift_box = VEC<VType>::box_from_rect_min(shift_rect);

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

				int cell_idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

				DBL3 position = DBL3(VEC<VType>::h.x * ((double)i + 0.5), VEC<VType>::h.y * ((double)j + 0.5), VEC<VType>::h.z * ((double)k + 0.5)) - delta;

				if (VEC<VType>::rect.contains(position) && (ngbrFlags[cell_idx] & NF_NOTEMPTY)) {

					VEC<VType>::quantity[cell_idx] = weighted_average(position, VEC<VType>::h);
				}
			}
		}
	}
}

//shift all the values in this VEC by the given delta (units same as VEC<VType>::h)
template <typename VType>
void VEC_VC<VType>::shift_keepmag(const DBL3& delta, const Rect& shift_rect)
{
	Box shift_box = VEC<VType>::box_from_rect_min(shift_rect);

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

				int cell_idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

				DBL3 position = DBL3(VEC<VType>::h.x * ((double)i + 0.5), VEC<VType>::h.y * ((double)j + 0.5), VEC<VType>::h.z * ((double)k + 0.5)) - delta;

				if (VEC<VType>::rect.contains(position) && (ngbrFlags[cell_idx] & NF_NOTEMPTY)) {

					//for vectorial quantities (e.g. DBL3) the operator * is a scalar product.
					double old_magnitude_squared = (double)(VEC<VType>::quantity[cell_idx] * VEC<VType>::quantity[cell_idx]);

					VType new_value = weighted_average(position, VEC<VType>::h);

					double new_magnitude_squared = (double)(new_value * new_value);

					if (new_magnitude_squared) {

						VEC<VType>::quantity[cell_idx] = new_value * sqrt(old_magnitude_squared / new_magnitude_squared);
					}
				}
			}
		}
	}
}

//shift all the values in this VEC by the given delta (units same as VEC<VType>::h)
template <typename VType>
void VEC_VC<VType>::shift_x(double delta, const Rect& shift_rect)
{
	if (fabs(shift_debt.x + delta) < VEC<VType>::h.x) {

		//total shift not enough : bank it and return
		shift_debt.x += delta;
		return;
	}

	//only shift an integer number of cells : there might be a sub-cellsize remainder so just bank it to be used next time
	int cells_shift = (int)((shift_debt.x + delta) / VEC<VType>::h.x);
	shift_debt.x -= VEC<VType>::h.x * cells_shift - delta;

	Box shift_box = VEC<VType>::box_from_rect_min(shift_rect);

	if (cells_shift < 0) {

		for (int i = shift_box.s.x; i < shift_box.e.x + cells_shift; i++) {
#pragma omp parallel for
			for (int j = shift_box.s.y; j < shift_box.e.y; j++) {
				for (int k = shift_box.s.z; k < shift_box.e.z; k++) {

					int cell_idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;
					int shift_cell_idx = cell_idx - cells_shift;

					if ((ngbrFlags[cell_idx] & NF_NOTEMPTY) && (ngbrFlags[shift_cell_idx] & NF_NOTEMPTY)) {
						
						VEC<VType>::quantity[cell_idx] = VEC<VType>::quantity[shift_cell_idx];
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

					int cell_idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;
					int shift_cell_idx = cell_idx - cells_shift;

					if ((ngbrFlags[cell_idx] & NF_NOTEMPTY) && (ngbrFlags[shift_cell_idx] & NF_NOTEMPTY)) {

						VEC<VType>::quantity[cell_idx] = VEC<VType>::quantity[shift_cell_idx];
					}
				}
			}
		}
	}
}

//shift all the values in this VEC by the given delta (units same as VEC<VType>::h)
template <typename VType>
void VEC_VC<VType>::shift_y(double delta, const Rect& shift_rect)
{
	if (fabs(shift_debt.y + delta) < VEC<VType>::h.y) {

		//total shift not enough : bank it and return
		shift_debt.y += delta;
		return;
	}

	//only shift an integer number of cells : there might be a sub-cellsize remainder so just bank it to be used next time
	int cells_shift = (int)((shift_debt.y + delta) / VEC<VType>::h.y);
	shift_debt.y -= VEC<VType>::h.y * cells_shift - delta;

	Box shift_box = VEC<VType>::box_from_rect_min(shift_rect);

	if (cells_shift < 0) {

		for (int j = shift_box.s.y; j < shift_box.e.y + cells_shift; j++) {
#pragma omp parallel for
			for (int i = shift_box.s.x; i < shift_box.e.x; i++) {
				for (int k = shift_box.s.z; k < shift_box.e.z; k++) {

					int cell_idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;
					int shift_cell_idx = cell_idx - cells_shift * VEC<VType>::n.x;

					if ((ngbrFlags[cell_idx] & NF_NOTEMPTY) && (ngbrFlags[shift_cell_idx] & NF_NOTEMPTY)) {

						VEC<VType>::quantity[cell_idx] = VEC<VType>::quantity[shift_cell_idx];
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

					int cell_idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;
					int shift_cell_idx = cell_idx - cells_shift * VEC<VType>::n.x;

					if ((ngbrFlags[cell_idx] & NF_NOTEMPTY) && (ngbrFlags[shift_cell_idx] & NF_NOTEMPTY)) {

						VEC<VType>::quantity[cell_idx] = VEC<VType>::quantity[shift_cell_idx];
					}
				}
			}
		}
	}
}