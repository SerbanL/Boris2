#pragma once

#include "VEC_VC.h"

//--------------------------------------------IMPORTANT FLAG MANIPULATION METHODS

//check if we need to use ngbrFlags2 (allocate memory etc.)
template <typename VType>
bool VEC_VC<VType>::use_extended_flags(void)
{
	//currently ngbrFlags2 used for:

	//ROBIN
	//DIRICHLET

	bool using_dirichlet = (dirichlet_nx.size() || dirichlet_px.size() || dirichlet_ny.size() || dirichlet_py.size() || dirichlet_nz.size() || dirichlet_pz.size());

	bool using_robin = (robin_px != DBL2() || robin_nx != DBL2() || robin_py != DBL2() || robin_ny != DBL2() || robin_pz != DBL2() || robin_nz != DBL2() || robin_v != DBL2());

	bool using_extended_flags = (using_dirichlet || using_robin);

	//make sure ngbrFlags2 has the correct memory allocated only if currently empty
	if (using_extended_flags && !ngbrFlags2.size()) {

		malloc_vector(ngbrFlags2, VEC<VType>::n.dim());
		ngbrFlags2.assign(VEC<VType>::n.dim(), 0);
	}

	return using_extended_flags;
}

template <typename VType>
void VEC_VC<VType>::resize_ngbrFlags(SZ3 new_n)
{
	if (!ngbrFlags.size()) {

		//if the VEC is being resized from a zero size then set solid shape : memory has been reserved for ngbrFlags so this cannot fail
		ngbrFlags.assign(new_n.dim(), NF_NOTEMPTY);
	}
	else {

		if (new_n != VEC<VType>::n) {

			//VEC is being resized from non-zero size : map current shape to new size

			//save old flags
			std::vector<int> old_ngbrFlags = ngbrFlags;

			//resize flags to new size and zero : memory has been reserved for ngbrFlags so this cannot fail 
			ngbrFlags.assign(new_n.dim(), 0);

			//now map shape from ngbrFlags_swapspace to ngbrFlags
			DBL3 sourceIdx = (DBL3)VEC<VType>::n / new_n;

#pragma omp parallel for
			for (int idx = 0; idx < new_n.dim(); idx++) {

				int _x = (int)floor((idx % new_n.x) * sourceIdx.x);
				int _y = (int)floor(((idx / new_n.x) % new_n.y) * sourceIdx.y);
				int _z = (int)floor((idx / (new_n.x*new_n.y)) * sourceIdx.z);

				if (old_ngbrFlags[_x + _y * VEC<VType>::n.x + _z * (VEC<VType>::n.x*VEC<VType>::n.y)] & NF_NOTEMPTY)
					ngbrFlags[idx] = NF_NOTEMPTY;
			}
		}
		else {

			//VEC is not being resized : clear everything but the shape

#pragma omp parallel for
			for (int idx = 0; idx < VEC<VType>::n.dim(); idx++) {

				if (ngbrFlags[idx] & NF_NOTEMPTY) ngbrFlags[idx] = NF_NOTEMPTY;
				else ngbrFlags[idx] = 0;
			}
		}
	}

	//clear ngbrFlags2 if in use
	if (use_extended_flags()) {

		ngbrFlags2.assign(new_n.dim(), 0);
	}
}

//initialization method for neighbor flags : set flags at size VEC<VType>::n, counting neighbors etc. Use current shape in ngbrFlags
template <typename VType>
void VEC_VC<VType>::set_ngbrFlags(void)
{
	if (VEC<VType>::rect.IsNull() || VEC<VType>::h == DBL3()) return;

	//clear any shift debt as mesh has been resized so not valid anymore
	shift_debt = DBL3();

	//dirichlet flags will be cleared from ngbrFlags, so also clear the dirichlet vectors
	clear_dirichlet_flags();

	//also count non-empty points
	int cellsCount = 0;

	//1. Count all the neighbors

#pragma omp parallel for reduction(+:cellsCount)
	for (int j = 0; j < VEC<VType>::n.y; j++) {
		for (int k = 0; k < VEC<VType>::n.z; k++) {
			for (int i = 0; i < VEC<VType>::n.x; i++) {

				int idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

				if (ngbrFlags[idx] & NF_NOTEMPTY) {

					cellsCount++;

					if (i < VEC<VType>::n.x - 1) {

						//on-axis
						if (ngbrFlags[idx + 1] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NPX; }

						//off-axis (z slice : xy)
						if (j < VEC<VType>::n.y - 1) {

							if (ngbrFlags[idx + 1 + VEC<VType>::n.x] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XY_PXPY; }
						}

						if (j > 0) {

							if (ngbrFlags[idx + 1 - VEC<VType>::n.x] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XY_PXNY; }
						}
					}

					if (i > 0) {

						//on-axis
						if (ngbrFlags[idx - 1] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NNX; }

						//off-axis (z slice : xy)
						if (j < VEC<VType>::n.y - 1) {

							if (ngbrFlags[idx - 1 + VEC<VType>::n.x] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XY_NXPY; }
						}

						if (j > 0) {

							if (ngbrFlags[idx - 1 - VEC<VType>::n.x] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XY_NXNY; }
						}
					}

					if (j < VEC<VType>::n.y - 1) {

						//on-axis
						if (ngbrFlags[idx + VEC<VType>::n.x] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NPY; }

						//off-axis (x slice : yz)
						if (k < VEC<VType>::n.z - 1) {

							if (ngbrFlags[idx + VEC<VType>::n.x + VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_YZ_PYPZ; }
						}

						if (k > 0) {

							if (ngbrFlags[idx + VEC<VType>::n.x - VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_YZ_PYNZ; }
						}
					}

					if (j > 0) {

						//on-axis
						if (ngbrFlags[idx - VEC<VType>::n.x] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NNY; }

						//off-axis (x slice : yz)
						if (k < VEC<VType>::n.z - 1) {

							if (ngbrFlags[idx - VEC<VType>::n.x + VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_YZ_NYPZ; }
						}

						if (k > 0) {

							if (ngbrFlags[idx - VEC<VType>::n.x - VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_YZ_NYNZ; }
						}
					}

					if (k < VEC<VType>::n.z - 1) {

						//on-axis
						if (ngbrFlags[idx + VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NPZ; }

						//off-axis (y slice : xz)
						if (i < VEC<VType>::n.x - 1) {

							if (ngbrFlags[idx + 1 + VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XZ_PXPZ; }
						}

						if (i > 0) {

							if (ngbrFlags[idx - 1 + VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XZ_NXPZ; }
						}
					}

					if (k > 0) {

						//on-axis
						if (ngbrFlags[idx - VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NNZ; }

						//off-axis (y slice : xz)
						if (i < VEC<VType>::n.x - 1) {

							if (ngbrFlags[idx + 1 - VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XZ_PXNZ; }
						}

						if (i > 0) {

							if (ngbrFlags[idx - 1 - VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XZ_NXNZ; }
						}
					}

					//calculate mixed second order differentials off-axis stencils
					bool right_side, left_side, center_column;

					//XY plane
					right_side = ((int(ngbrFlags[idx] & NF_XY_PXPY) + int(ngbrFlags[idx] & NF_XY_PXNY) + int(ngbrFlags[idx] & NF_NPX)) >= 2);
					left_side = ((int(ngbrFlags[idx] & NF_XY_NXPY) + int(ngbrFlags[idx] & NF_XY_NXNY) + int(ngbrFlags[idx] & NF_NNX)) >= 2);
					center_column = ((ngbrFlags[idx] & NF_NPY) || (ngbrFlags[idx] & NF_NNY));

					if ((right_side && left_side) || (right_side && center_column) || (center_column && left_side)) ngbrFlags[idx] |= NF_XY_OASTENCIL;

					//XZ plane
					right_side = ((int(ngbrFlags[idx] & NF_XZ_PXPZ) + int(ngbrFlags[idx] & NF_XZ_PXNZ) + int(ngbrFlags[idx] & NF_NPX)) >= 2);
					left_side = ((int(ngbrFlags[idx] & NF_XZ_NXPZ) + int(ngbrFlags[idx] & NF_XZ_NXNZ) + int(ngbrFlags[idx] & NF_NNX)) >= 2);
					center_column = ((ngbrFlags[idx] & NF_NPZ) || (ngbrFlags[idx] & NF_NNZ));

					if ((right_side && left_side) || (right_side && center_column) || (center_column && left_side)) ngbrFlags[idx] |= NF_XZ_OASTENCIL;

					//YZ plane
					right_side = ((int(ngbrFlags[idx] & NF_YZ_PYPZ) + int(ngbrFlags[idx] & NF_YZ_PYNZ) + int(ngbrFlags[idx] & NF_NPY)) >= 2);
					left_side = ((int(ngbrFlags[idx] & NF_YZ_NYPZ) + int(ngbrFlags[idx] & NF_YZ_NYNZ) + int(ngbrFlags[idx] & NF_NNY)) >= 2);
					//center_column same as above

					if ((right_side && left_side) || (right_side && center_column) || (center_column && left_side)) ngbrFlags[idx] |= NF_YZ_OASTENCIL;
				}
				else {

					VEC<VType>::quantity[idx] = VType();
				}
			}
		}
	}

	//2. set Robin flags depending on set conditions
	set_robin_flags();

	//3. set pbc flags depending on set conditions and currently calculated flags
	set_pbc_flags();

	nonempty_cells = cellsCount;
}


//initialization method for neighbor flags : set flags at size VEC<VType>::n, counting neighbors etc.
//Set empty cell values using information in linked_vec (keep same shape) - this must have same rectangle
template <typename VType>
template <typename LVType>
void VEC_VC<VType>::set_ngbrFlags(const VEC_VC<LVType>& linked_vec)
{
	//linked_vec must have same mesh rect
	if (linked_vec.rect != VEC<VType>::rect || VEC<VType>::rect.IsNull() || VEC<VType>::h == DBL3()) return;

	//copy shape from linked_vec
#pragma omp parallel for
	for (int idx = 0; idx < VEC<VType>::n.dim(); idx++) {

		if (!linked_vec.is_empty(VEC<VType>::get_cellrect(idx))) mark_not_empty(idx);
		else mark_empty(idx);
	}

	//now continue with set_ngbrFlags
	set_ngbrFlags();
}

//from NF2_DIRICHLET type flag and cell_idx return boundary value from one of the dirichlet vectors
template <typename VType>
VType VEC_VC<VType>::get_dirichlet_value(int dirichlet_flag, int cell_idx) const
{
	switch (dirichlet_flag) {

	case NF2_DIRICHLETPX:
		return dirichlet_nx[((cell_idx / VEC<VType>::n.x) % VEC<VType>::n.y) + (cell_idx / (VEC<VType>::n.x*VEC<VType>::n.y))*VEC<VType>::n.y];

	case NF2_DIRICHLETNX:
		return dirichlet_px[((cell_idx / VEC<VType>::n.x) % VEC<VType>::n.y) + (cell_idx / (VEC<VType>::n.x*VEC<VType>::n.y))*VEC<VType>::n.y];

	case NF2_DIRICHLETPY:
		return dirichlet_ny[(cell_idx % VEC<VType>::n.x) + (cell_idx / (VEC<VType>::n.x*VEC<VType>::n.y))*VEC<VType>::n.x];

	case NF2_DIRICHLETNY:
		return dirichlet_py[(cell_idx % VEC<VType>::n.x) + (cell_idx / (VEC<VType>::n.x*VEC<VType>::n.y))*VEC<VType>::n.x];

	case NF2_DIRICHLETPZ:
		return dirichlet_nz[(cell_idx % VEC<VType>::n.x) + ((cell_idx / VEC<VType>::n.x) % VEC<VType>::n.y)*VEC<VType>::n.x];

	case NF2_DIRICHLETNZ:
		return dirichlet_pz[(cell_idx % VEC<VType>::n.x) + ((cell_idx / VEC<VType>::n.x) % VEC<VType>::n.y)*VEC<VType>::n.x];
	}

	return VType();
}

//set robin flags from robin values and shape. Doesn't affect any other flags. Call from set_ngbrFlags after counting neighbors, and after setting robin values
template <typename VType>
void VEC_VC<VType>::set_robin_flags(void)
{
	//ROBIN flags used with extended ngbrFlags only
	if (use_extended_flags()) {

#pragma omp parallel for
		for (int j = 0; j < VEC<VType>::n.y; j++) {
			for (int k = 0; k < VEC<VType>::n.z; k++) {
				for (int i = 0; i < VEC<VType>::n.x; i++) {

					int idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

					//first clear any robin flags already set
					ngbrFlags2[idx] &= ~NF2_ROBIN;

					if (ngbrFlags[idx] & NF_NOTEMPTY) {

						//neighbors
						if (i < VEC<VType>::n.x - 1) {

							if (!(ngbrFlags[idx + 1] & NF_NOTEMPTY))
								//inner cell next to a void cell on the -x side
								if (IsNZ(robin_v.i) && i > 0 && ngbrFlags[idx - 1] & NF_NOTEMPTY) { ngbrFlags2[idx] |= (NF2_ROBINNX + NF2_ROBINV); }
						}
						//surface cell on the -x side of the surface
						else if (IsNZ(robin_px.i) && i > 0 && (ngbrFlags[idx - 1] & NF_NOTEMPTY)) ngbrFlags2[idx] |= NF2_ROBINNX;

						if (i > 0) {

							if (!(ngbrFlags[idx - 1] & NF_NOTEMPTY))
								if (IsNZ(robin_v.i) && i < VEC<VType>::n.x - 1 && ngbrFlags[idx + 1] & NF_NOTEMPTY) { ngbrFlags2[idx] |= (NF2_ROBINPX + NF2_ROBINV); }
						}
						else if (IsNZ(robin_nx.i) && i < VEC<VType>::n.x - 1 && (ngbrFlags[idx + 1] & NF_NOTEMPTY)) ngbrFlags2[idx] |= NF2_ROBINPX;

						if (j < VEC<VType>::n.y - 1) {

							if (!(ngbrFlags[idx + VEC<VType>::n.x] & NF_NOTEMPTY))
								if (IsNZ(robin_v.i) && j > 0 && (ngbrFlags[idx - VEC<VType>::n.x] & NF_NOTEMPTY)) { ngbrFlags2[idx] |= (NF2_ROBINNY + NF2_ROBINV); }
						}
						else if (IsNZ(robin_py.i) && j > 0 && (ngbrFlags[idx - VEC<VType>::n.x] & NF_NOTEMPTY)) ngbrFlags2[idx] |= NF2_ROBINNY;

						if (j > 0) {

							if (!(ngbrFlags[idx - VEC<VType>::n.x] & NF_NOTEMPTY))
								if (IsNZ(robin_v.i) && j < VEC<VType>::n.y - 1 && (ngbrFlags[idx + VEC<VType>::n.x] & NF_NOTEMPTY)) { ngbrFlags2[idx] |= (NF2_ROBINPY + NF2_ROBINV); }
						}
						else if (IsNZ(robin_ny.i) && j < VEC<VType>::n.y - 1 && (ngbrFlags[idx + VEC<VType>::n.x] & NF_NOTEMPTY)) ngbrFlags2[idx] |= NF2_ROBINPY;

						if (k < VEC<VType>::n.z - 1) {

							if (!(ngbrFlags[idx + VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY))
								if (IsNZ(robin_v.i) && k > 0 && (ngbrFlags[idx - VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY)) { ngbrFlags2[idx] |= (NF2_ROBINNZ + NF2_ROBINV); }
						}
						else if (IsNZ(robin_pz.i) && k > 0 && (ngbrFlags[idx - VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY)) ngbrFlags2[idx] |= NF2_ROBINNZ;

						if (k > 0) {

							if (!(ngbrFlags[idx - VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY))
								if (IsNZ(robin_v.i) && k < VEC<VType>::n.z - 1 && (ngbrFlags[idx + VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY)) { ngbrFlags2[idx] |= (NF2_ROBINPZ + NF2_ROBINV); }
						}
						else if (IsNZ(robin_nz.i) && k < VEC<VType>::n.z - 1 && (ngbrFlags[idx + VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY)) ngbrFlags2[idx] |= NF2_ROBINPZ;
					}
				}
			}
		}
	}
}

//--------------------------------------------FLAG CHECKING

//check if all cells intersecting the rectangle (absolute coordinates) are empty
template <typename VType>
bool VEC_VC<VType>::is_empty(const Rect& rectangle) const
{
	Box cells = VEC<VType>::box_from_rect_max(rectangle);

	for (int i = cells.s.x; i < cells.e.x; i++) {
		for (int j = cells.s.y; j < cells.e.y; j++) {
			for (int k = cells.s.z; k < cells.e.z; k++) {

				if (is_not_empty(INT3(i, j, k))) return false;
			}
		}
	}

	return true;
}

//check if all cells intersecting the rectangle (absolute coordinates) are not empty
template <typename VType>
bool VEC_VC<VType>::is_not_empty(const Rect& rectangle) const
{
	Box cells = VEC<VType>::box_from_rect_max(rectangle);

	for (int i = cells.s.x; i < cells.e.x; i++) {
		for (int j = cells.s.y; j < cells.e.y; j++) {
			for (int k = cells.s.z; k < cells.e.z; k++) {

				if (is_empty(INT3(i, j, k))) return false;
			}
		}
	}

	return true;
}

//--------------------------------------------SET CELL FLAGS - EXTERNAL USE

//set dirichlet boundary conditions from surface_rect (must be a rectangle intersecting with one of the surfaces of this mesh) and value
////return false on memory allocation failure only, otherwise return true even if surface_rect was not valid
template <typename VType>
bool VEC_VC<VType>::set_dirichlet_conditions(const Rect& surface_rect, VType value)
{
	if (!VEC<VType>::rect.intersects(surface_rect)) return true;

	Rect intersection = VEC<VType>::rect.get_intersection(surface_rect);
	//if (!intersection.IsPlane()) return true;

	auto set_dirichlet_value = [&](Rect intersection, VType value, int flag_value) -> void {

		INT3 start = VEC<VType>::box_from_rect_max(intersection).s;
		INT3 end = VEC<VType>::box_from_rect_max(intersection).e;

		//set cells
#pragma omp parallel for
		for (int j = (start.y >= 0 ? start.y : 0); j < (end.y < VEC<VType>::n.y ? end.y : VEC<VType>::n.y); j++) {
			for (int k = (start.z >= 0 ? start.z : 0); k < (end.z < VEC<VType>::n.z ? end.z : VEC<VType>::n.z); k++) {
				for (int i = (start.x >= 0 ? start.x : 0); i < (end.x < VEC<VType>::n.x ? end.x : VEC<VType>::n.x); i++) {

					int cell_idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

					switch (flag_value) {

					case NF2_DIRICHLETPX:
						if ((cell_idx % VEC<VType>::n.x) != 0 || !(ngbrFlags[cell_idx] & NF_NPX)) continue;
						dirichlet_nx[((cell_idx / VEC<VType>::n.x) % VEC<VType>::n.y) + (cell_idx / (VEC<VType>::n.x*VEC<VType>::n.y)) * VEC<VType>::n.y] = value;
						break;

					case NF2_DIRICHLETNX:
						if ((cell_idx % VEC<VType>::n.x) != (VEC<VType>::n.x - 1) || !(ngbrFlags[cell_idx] & NF_NNX)) continue;
						dirichlet_px[((cell_idx / VEC<VType>::n.x) % VEC<VType>::n.y) + (cell_idx / (VEC<VType>::n.x*VEC<VType>::n.y)) * VEC<VType>::n.y] = value;
						break;

					case NF2_DIRICHLETPY:
						if (((cell_idx / VEC<VType>::n.x) % VEC<VType>::n.y) != 0 || !(ngbrFlags[cell_idx] & NF_NPY)) continue;
						dirichlet_ny[(cell_idx % VEC<VType>::n.x) + (cell_idx / (VEC<VType>::n.x*VEC<VType>::n.y)) * VEC<VType>::n.x] = value;
						break;

					case NF2_DIRICHLETNY:
						if (((cell_idx / VEC<VType>::n.x) % VEC<VType>::n.y) != (VEC<VType>::n.y - 1) || !(ngbrFlags[cell_idx] & NF_NNY)) continue;
						dirichlet_py[(cell_idx % VEC<VType>::n.x) + (cell_idx / (VEC<VType>::n.x*VEC<VType>::n.y)) * VEC<VType>::n.x] = value;
						break;

					case NF2_DIRICHLETPZ:
						if ((cell_idx / (VEC<VType>::n.x*VEC<VType>::n.y)) != 0 || !(ngbrFlags[cell_idx] & NF_NPZ)) continue;
						dirichlet_nz[(cell_idx % VEC<VType>::n.x) + ((cell_idx / VEC<VType>::n.x) % VEC<VType>::n.y) * VEC<VType>::n.x] = value;
						break;

					case NF2_DIRICHLETNZ:
						if ((cell_idx / (VEC<VType>::n.x*VEC<VType>::n.y)) != (VEC<VType>::n.z - 1) || !(ngbrFlags[cell_idx] & NF_NNZ)) continue;
						dirichlet_pz[(cell_idx % VEC<VType>::n.x) + ((cell_idx / VEC<VType>::n.x) % VEC<VType>::n.y) * VEC<VType>::n.x] = value;
						break;
					}

					ngbrFlags2[cell_idx] |= flag_value;
				}
			}
		}
	};

	//y-z plane
	if (IsZ(intersection.s.x - intersection.e.x)) {

		//on lower x side
		if (IsZ(VEC<VType>::rect.s.x - intersection.s.x)) {

			if (!dirichlet_nx.size()) {

				if (!malloc_vector(dirichlet_nx, VEC<VType>::n.y*VEC<VType>::n.z)) return false;

				//DIRICHLET flags used with extended ngbrFlags so make sure it has memory allocated now
				use_extended_flags();
			}

			intersection.e.x += VEC<VType>::h.x;
			set_dirichlet_value(intersection, value, NF2_DIRICHLETPX);
		}
		//on upper x side
		else if (IsZ(VEC<VType>::rect.e.x - intersection.s.x)) {

			if (!dirichlet_px.size()) {

				if (!malloc_vector(dirichlet_px, VEC<VType>::n.y*VEC<VType>::n.z)) return false;

				//DIRICHLET flags used with extended ngbrFlags so make sure it has memory allocated now
				use_extended_flags();
			}

			intersection.s.x -= VEC<VType>::h.x;
			set_dirichlet_value(intersection, value, NF2_DIRICHLETNX);
		}
	}
	//x-z plane
	else if (IsZ(intersection.s.y - intersection.e.y)) {

		//on lower y side
		if (IsZ(VEC<VType>::rect.s.y - intersection.s.y)) {

			if (!dirichlet_ny.size()) {

				if (!malloc_vector(dirichlet_ny, VEC<VType>::n.x*VEC<VType>::n.z)) return false;

				//DIRICHLET flags used with extended ngbrFlags so make sure it has memory allocated now
				use_extended_flags();
			}

			intersection.e.y += VEC<VType>::h.y;
			set_dirichlet_value(intersection, value, NF2_DIRICHLETPY);
		}
		//on upper y side
		else if (IsZ(VEC<VType>::rect.e.y - intersection.s.y)) {

			if (!dirichlet_py.size()) {

				if (!malloc_vector(dirichlet_py, VEC<VType>::n.x*VEC<VType>::n.z)) return false;

				//DIRICHLET flags used with extended ngbrFlags so make sure it has memory allocated now
				use_extended_flags();
			}

			intersection.s.y -= VEC<VType>::h.y;
			set_dirichlet_value(intersection, value, NF2_DIRICHLETNY);
		}
	}
	//x-y plane
	else if (IsZ(intersection.s.z - intersection.e.z)) {

		//on lower z side
		if (IsZ(VEC<VType>::rect.s.z - intersection.s.z)) {

			if (!dirichlet_nz.size()) {

				if (!malloc_vector(dirichlet_nz, VEC<VType>::n.x*VEC<VType>::n.y)) return false;

				//DIRICHLET flags used with extended ngbrFlags so make sure it has memory allocated now
				use_extended_flags();
			}

			intersection.e.z += VEC<VType>::h.z;
			set_dirichlet_value(intersection, value, NF2_DIRICHLETPZ);
		}
		//on upper z side
		else if (IsZ(VEC<VType>::rect.e.z - intersection.s.z)) {

			if (!dirichlet_pz.size()) {

				if (!malloc_vector(dirichlet_pz, VEC<VType>::n.x*VEC<VType>::n.y)) return false;

				//DIRICHLET flags used with extended ngbrFlags so make sure it has memory allocated now
				use_extended_flags();
			}

			intersection.s.z -= VEC<VType>::h.z;
			set_dirichlet_value(intersection, value, NF2_DIRICHLETNZ);
		}
	}

	return true;
}

//clear all dirichlet flags and vectors
template <typename VType>
void VEC_VC<VType>::clear_dirichlet_flags(void)
{
	//DIRICHLET flags used with extended ngbrFlags only.
	if (use_extended_flags()) {

#pragma omp parallel for
		for (int idx = 0; idx < (int)ngbrFlags2.size(); idx++) {

			ngbrFlags2[idx] &= ~NF2_DIRICHLET;
		}

		dirichlet_px.clear();
		dirichlet_nx.clear();
		dirichlet_py.clear();
		dirichlet_ny.clear();
		dirichlet_pz.clear();
		dirichlet_nz.clear();

		dirichlet_px.shrink_to_fit();
		dirichlet_nx.shrink_to_fit();
		dirichlet_py.shrink_to_fit();
		dirichlet_ny.shrink_to_fit();
		dirichlet_pz.shrink_to_fit();
		dirichlet_nz.shrink_to_fit();

		//clear memory for extended ngbrFlags if now not used for anything else
		if (!use_extended_flags()) {

			ngbrFlags2.clear();
			ngbrFlags2.shrink_to_fit();
		}
	}
	else {

		ngbrFlags2.clear();
		ngbrFlags2.shrink_to_fit();
	}
}

//set pbc flags depending on set conditions and currently calculated flags - ngbrFlags must already be calculated before using this
template <typename VType>
void VEC_VC<VType>::set_pbc_flags(void)
{
#pragma omp parallel for
	for (int i = 0; i < VEC<VType>::n.x; i++) {
		for (int j = 0; j < VEC<VType>::n.y; j++) {
			for (int k = 0; k < VEC<VType>::n.z; k++) {

				int idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

				//skip empty cells
				if (!(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

				//first clear pbc in this cell before recalculating pbc flags
				ngbrFlags[idx] &= ~NF_PBC;

				//-x side, and on the +x side there is a non-empty cell : set pbc
				if (pbc_x) {

					//-x side, and on the +x side there is a non-empty cell : set pbc
					if (i == 0 && (ngbrFlags[VEC<VType>::n.x - 1 + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_PBCX;

					//+x side, and on the -x side there is a non-empty cell : set pbc
					if (i == VEC<VType>::n.x - 1 && (ngbrFlags[j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_PBCX;
				}

				if (pbc_y) {

					//-y side, and on the +y side there is a non-empty cell : set pbc
					if (j == 0 && (ngbrFlags[i + (VEC<VType>::n.y - 1) * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_PBCY;

					//+y side, and on the -y side there is a non-empty cell : set pbc
					if (j == VEC<VType>::n.y - 1 && (ngbrFlags[i + k * VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_PBCY;
				}

				if (pbc_z) {

					//-z side, and on the +z side there is a non-empty cell : set pbc
					if (k == 0 && (ngbrFlags[i + j * VEC<VType>::n.x + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_PBCZ;

					//+z side, and on the -z side there is a non-empty cell : set pbc
					if (k == VEC<VType>::n.z - 1 && (ngbrFlags[i + j * VEC<VType>::n.x] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_PBCZ;
				}
			}
		}
	}
}

//set pbc conditions : setting any to false clears flags
template <typename VType>
void VEC_VC<VType>::set_pbc(bool pbc_x_, bool pbc_y_, bool pbc_z_)
{
	pbc_x = pbc_x_;
	pbc_y = pbc_y_;
	pbc_z = pbc_z_;

	set_pbc_flags();
}

//clear all pbc flags
template <typename VType>
void VEC_VC<VType>::clear_pbc(void)
{
	pbc_x = 0;
	pbc_y = 0;
	pbc_z = 0;

#pragma omp parallel for
	for (int idx = 0; idx < (int)ngbrFlags.size(); idx++) {

		ngbrFlags[idx] &= ~NF_PBC;
	}
}

template <typename VType>
void VEC_VC<VType>::clear_cmbnd_flags(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < (int)ngbrFlags.size(); idx++)
		ngbrFlags[idx] &= ~NF_CMBND;
}

//mark cells included in this rectangle (absolute coordinates) to be skipped during some computations
template <typename VType>
void VEC_VC<VType>::set_skipcells(const Rect& rectangle, bool status)
{
	Box cells = VEC<VType>::box_from_rect_max(rectangle);

	for (int i = (cells.s.x >= 0 ? cells.s.x : 0); i < (cells.e.x <= VEC<VType>::n.x ? cells.e.x : VEC<VType>::n.x); i++) {
		for (int j = (cells.s.y >= 0 ? cells.s.y : 0); j < (cells.e.y <= VEC<VType>::n.y ? cells.e.y : VEC<VType>::n.y); j++) {
			for (int k = (cells.s.z >= 0 ? cells.s.z : 0); k < (cells.e.z <= VEC<VType>::n.z ? cells.e.z : VEC<VType>::n.z); k++) {

				if(status) ngbrFlags[i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y] |= NF_SKIPCELL;
				else ngbrFlags[i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y] &= ~NF_SKIPCELL;
			}
		}
	}
}

//clear all skip cell flags
template <typename VType>
void VEC_VC<VType>::clear_skipcells(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < (int)ngbrFlags.size(); idx++)
		ngbrFlags[idx] &= ~NF_SKIPCELL;
}

template <typename VType>
void VEC_VC<VType>::set_robin_conditions(DBL2 robin_v_, DBL2 robin_px_, DBL2 robin_nx_, DBL2 robin_py_, DBL2 robin_ny_, DBL2 robin_pz_, DBL2 robin_nz_)
{
	robin_px = robin_px_;
	robin_nx = robin_nx_;
	robin_py = robin_py_;
	robin_ny = robin_ny_;
	robin_pz = robin_pz_;
	robin_nz = robin_nz_;

	robin_v = robin_v_;

	//set flags
	set_robin_flags();
}

//clear all Robin boundary conditions and values
template <typename VType>
void VEC_VC<VType>::clear_robin_conditions(void)
{
	robin_px = DBL2();
	robin_nx = DBL2();
	robin_py = DBL2();
	robin_ny = DBL2();
	robin_pz = DBL2();
	robin_nz = DBL2();

	robin_v = DBL2();

	if (use_extended_flags()) {

		//clear Robin flags
#pragma omp parallel for
		for (int idx = 0; idx < (int)ngbrFlags2.size(); idx++) {

			ngbrFlags2[idx] &= ~NF2_ROBIN;
		}
	}
	else {

		ngbrFlags2.clear();
		ngbrFlags2.shrink_to_fit();
	}
}
