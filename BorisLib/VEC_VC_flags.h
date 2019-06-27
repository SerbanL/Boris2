#pragma once

#include "VEC_VC.h"

//--------------------------------------------IMPORTANT FLAG MANIPULATION METHODS

template <typename VType>
void VEC_VC<VType>::resize_ngbrFlags(SZ3 new_n)
{
	if (!ngbrFlags.size()) {

		//if the VEC is being resized from a zero size then set solid shape : memory has been reserved for ngbrFlags so this cannot fail
		ngbrFlags.assign(new_n.dim(), NF_NOTEMPTY);
	}
	else {

		if (new_n != n) {

			//VEC is being resized from non-zero size : map current shape to new size

			//save old flags
			std::vector<int> old_ngbrFlags = ngbrFlags;

			//resize flags to new size and zero : memory has been reserved for ngbrFlags so this cannot fail 
			ngbrFlags.assign(new_n.dim(), 0);

			//now map shape from ngbrFlags_swapspace to ngbrFlags
			DBL3 sourceIdx = (DBL3)n / new_n;

#pragma omp parallel for
			for (int idx = 0; idx < new_n.dim(); idx++) {

				int _x = (int)floor((idx % new_n.x) * sourceIdx.x);
				int _y = (int)floor(((idx / new_n.x) % new_n.y) * sourceIdx.y);
				int _z = (int)floor((idx / (new_n.x*new_n.y)) * sourceIdx.z);

				if (old_ngbrFlags[_x + _y * n.x + _z * (n.x*n.y)] & NF_NOTEMPTY)
					ngbrFlags[idx] = NF_NOTEMPTY;
			}
		}
		else {

			//VEC is not being resized : clear everything but the shape

#pragma omp parallel for
			for (int idx = 0; idx < n.dim(); idx++) {

				if (ngbrFlags[idx] & NF_NOTEMPTY) ngbrFlags[idx] = NF_NOTEMPTY;
				else ngbrFlags[idx] = 0;
			}
		}
	}
}


//initialization method for neighbor flags : set flags at size n, counting neighbors etc. Use current shape in ngbrFlags
template <typename VType>
void VEC_VC<VType>::set_ngbrFlags(void)
{
	if (rect.IsNull() || h == DBL3()) return;

	//clear any shift debt as mesh has been resized so not valid anymore
	shift_debt = DBL3();

	//dirichlet flags will be cleared from ngbrFlags, so also clear the dirichlet vectors
	clear_dirichlet_flags();

	//also count non-empty points
	int cellsCount = 0;

	//1. Count all the neighbors

#pragma omp parallel for reduction(+:cellsCount)
	for (int j = 0; j < n.y; j++) {
		for (int k = 0; k < n.z; k++) {
			for (int i = 0; i < n.x; i++) {

				int idx = i + j * n.x + k * n.x*n.y;

				if (ngbrFlags[idx] & NF_NOTEMPTY) {

					cellsCount++;

					//neighbors
					if (i < n.x - 1) {

						if (ngbrFlags[idx + 1] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NPX; }
					}

					if (i > 0) {

						if (ngbrFlags[idx - 1] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NNX; }
					}

					if (j < n.y - 1) {

						if (ngbrFlags[idx + n.x] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NPY; }
					}

					if (j > 0) {

						if (ngbrFlags[idx - n.x] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NNY; }
					}

					if (k < n.z - 1) {

						if (ngbrFlags[idx + n.x*n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NPZ; }
					}

					if (k > 0) {

						if (ngbrFlags[idx - n.x*n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NNZ; }
					}

					if ((ngbrFlags[idx] & NF_NPX) && (ngbrFlags[idx] & NF_NNX)) ngbrFlags[idx] |= NF_BOTHX;
					if ((ngbrFlags[idx] & NF_NPY) && (ngbrFlags[idx] & NF_NNY)) ngbrFlags[idx] |= NF_BOTHY;
					if ((ngbrFlags[idx] & NF_NPZ) && (ngbrFlags[idx] & NF_NNZ)) ngbrFlags[idx] |= NF_BOTHZ;
				}
				else {

					quantity[idx] = VType();
				}
			}
		}
	}

	//2. set Robin flags depending on set conditions
	set_robin_flags();

	nonempty_cells = cellsCount;
}


//initialization method for neighbor flags : set flags at size n, counting neighbors etc.
//Set empty cell values using information in linked_vec (keep same shape) - this must have same rectangle
template <typename VType>
template <typename LVType>
void VEC_VC<VType>::set_ngbrFlags(VEC_VC<LVType>& linked_vec)
{
	//linked_vec must have same mesh rect
	if (linked_vec.rect != rect || rect.IsNull() || h == DBL3()) return;

	//clear any shift debt as mesh has been resized so not valid anymore
	shift_debt = DBL3();

	//dirichlet flags will be cleared from ngbrFlags, so also clear the dirichlet vectors
	clear_dirichlet_flags();

	//also count non-empty points
	int cellsCount = 0;

	//1. Count all the neighbors

#pragma omp parallel for reduction(+:cellsCount)
	for (int j = 0; j < n.y; j++) {
		for (int k = 0; k < n.z; k++) {
			for (int i = 0; i < n.x; i++) {

				int idx = i + j * n.x + k * n.x*n.y;

				if (!linked_vec.is_empty(get_cellrect(INT3(i, j, k)))) {

					//mark cell as not empty
					mark_not_empty(idx);
					cellsCount++;

					if (i < n.x - 1) {

						if (!linked_vec.is_empty(get_cellrect(INT3(i + 1, j, k)))) { ngbrFlags[idx] |= NF_NPX; }
					}

					if (i > 0) {

						if (!linked_vec.is_empty(get_cellrect(INT3(i - 1, j, k)))) { ngbrFlags[idx] |= NF_NNX; }
					}

					if (j < n.y - 1) {

						if (!linked_vec.is_empty(get_cellrect(INT3(i, j + 1, k)))) { ngbrFlags[idx] |= NF_NPY; }
					}

					if (j > 0) {

						if (!linked_vec.is_empty(get_cellrect(INT3(i, j - 1, k)))) { ngbrFlags[idx] |= NF_NNY; }
					}

					if (k < n.z - 1) {

						if (!linked_vec.is_empty(get_cellrect(INT3(i, j, k + 1)))) { ngbrFlags[idx] |= NF_NPZ; }
					}

					if (k > 0) {

						if (!linked_vec.is_empty(get_cellrect(INT3(i, j, k - 1)))) { ngbrFlags[idx] |= NF_NNZ; }
					}

					if ((ngbrFlags[idx] & NF_NPX) && (ngbrFlags[idx] & NF_NNX)) ngbrFlags[idx] |= NF_BOTHX;
					if ((ngbrFlags[idx] & NF_NPY) && (ngbrFlags[idx] & NF_NNY)) ngbrFlags[idx] |= NF_BOTHY;
					if ((ngbrFlags[idx] & NF_NPZ) && (ngbrFlags[idx] & NF_NNZ)) ngbrFlags[idx] |= NF_BOTHZ;
				}
				else {

					//if linked quantity is empty then also empty the quantity in this VEC. Note, the linked vec could actually be the same vec - in this case the shape is not changed since the NF_NOTEMPTY flags have already been mapped to the new size
					mark_empty(idx);
				}
			}
		}
	}

	//2. set Robin flags depending on set conditions
	set_robin_flags();

	nonempty_cells = cellsCount;
}

//from NF_DIRICHLET type flag and cell_idx return boundary value from one of the dirichlet vectors
template <typename VType>
VType VEC_VC<VType>::get_dirichlet_value(int dirichlet_flag, int cell_idx) const
{
	switch (dirichlet_flag) {

	case NF_DIRICHLETPX:
		return dirichlet_nx[((cell_idx / n.x) % n.y) + (cell_idx / (n.x*n.y))*n.y];

	case NF_DIRICHLETNX:
		return dirichlet_px[((cell_idx / n.x) % n.y) + (cell_idx / (n.x*n.y))*n.y];

	case NF_DIRICHLETPY:
		return dirichlet_ny[(cell_idx % n.x) + (cell_idx / (n.x*n.y))*n.x];

	case NF_DIRICHLETNY:
		return dirichlet_py[(cell_idx % n.x) + (cell_idx / (n.x*n.y))*n.x];

	case NF_DIRICHLETPZ:
		return dirichlet_nz[(cell_idx % n.x) + ((cell_idx / n.x) % n.y)*n.x];

	case NF_DIRICHLETNZ:
		return dirichlet_pz[(cell_idx % n.x) + ((cell_idx / n.x) % n.y)*n.x];
	}

	return VType();
}

//set robin flags from robin values and shape. Doesn't affect any other flags. Call from set_ngbrFlags after counting neighbors, and after setting robin values
template <typename VType>
void VEC_VC<VType>::set_robin_flags(void)
{
#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {
		for (int k = 0; k < n.z; k++) {
			for (int i = 0; i < n.x; i++) {

				int idx = i + j * n.x + k * n.x*n.y;

				if (ngbrFlags[idx] & NF_NOTEMPTY) {

					//neighbors
					if (i < n.x - 1) {

						if (!(ngbrFlags[idx + 1] & NF_NOTEMPTY))
							//inner cell next to a void cell on the -x side
							if (IsNZ(robin_v.i) && i > 0 && ngbrFlags[idx - 1] & NF_NOTEMPTY) { ngbrFlags[idx] |= (NF_ROBINNX + NF_ROBINV); }
					}
					//surface cell on the -x side of the surface
					else if (IsNZ(robin_px.i) && i > 0 && (ngbrFlags[idx - 1] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_ROBINNX;

					if (i > 0) {

						if (!(ngbrFlags[idx - 1] & NF_NOTEMPTY))
							if (IsNZ(robin_v.i) && i < n.x - 1 && ngbrFlags[idx + 1] & NF_NOTEMPTY) { ngbrFlags[idx] |= (NF_ROBINPX + NF_ROBINV); }
					}
					else if (IsNZ(robin_nx.i) && i < n.x - 1 && (ngbrFlags[idx + 1] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_ROBINPX;

					if (j < n.y - 1) {

						if (!(ngbrFlags[idx + n.x] & NF_NOTEMPTY))
							if (IsNZ(robin_v.i) && j > 0 && (ngbrFlags[idx - n.x] & NF_NOTEMPTY)) { ngbrFlags[idx] |= (NF_ROBINNY + NF_ROBINV); }
					}
					else if (IsNZ(robin_py.i) && j > 0 && (ngbrFlags[idx - n.x] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_ROBINNY;

					if (j > 0) {

						if (!(ngbrFlags[idx - n.x] & NF_NOTEMPTY))
							if (IsNZ(robin_v.i) && j < n.y - 1 && (ngbrFlags[idx + n.x] & NF_NOTEMPTY)) { ngbrFlags[idx] |= (NF_ROBINPY + NF_ROBINV); }
					}
					else if (IsNZ(robin_ny.i) && j < n.y - 1 && (ngbrFlags[idx + n.x] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_ROBINPY;

					if (k < n.z - 1) {

						if (!(ngbrFlags[idx + n.x*n.y] & NF_NOTEMPTY))
							if (IsNZ(robin_v.i) && k > 0 && (ngbrFlags[idx - n.x*n.y] & NF_NOTEMPTY)) { ngbrFlags[idx] |= (NF_ROBINNZ + NF_ROBINV); }
					}
					else if (IsNZ(robin_pz.i) && k > 0 && (ngbrFlags[idx - n.x*n.y] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_ROBINNZ;

					if (k > 0) {

						if (!(ngbrFlags[idx - n.x*n.y] & NF_NOTEMPTY))
							if (IsNZ(robin_v.i) && k < n.z - 1 && (ngbrFlags[idx + n.x*n.y] & NF_NOTEMPTY)) { ngbrFlags[idx] |= (NF_ROBINPZ + NF_ROBINV); }
					}
					else if (IsNZ(robin_nz.i) && k < n.z - 1 && (ngbrFlags[idx + n.x*n.y] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_ROBINPZ;
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
	Box cells = box_from_rect_max(rectangle);

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
	Box cells = box_from_rect_max(rectangle);

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
	if (!rect.intersects(surface_rect)) return true;

	Rect intersection = rect.get_intersection(surface_rect);
	if (!intersection.IsPlane()) return true;

	auto set_dirichlet_value = [&](Rect intersection, VType value, int flag_value) -> void {

		INT3 start = box_from_rect_max(intersection).s;
		INT3 end = box_from_rect_max(intersection).e;

		//set cells
#pragma omp parallel for
		for (int j = (start.y >= 0 ? start.y : 0); j < (end.y < n.y ? end.y : n.y); j++) {
			for (int k = (start.z >= 0 ? start.z : 0); k < (end.z < n.z ? end.z : n.z); k++) {
				for (int i = (start.x >= 0 ? start.x : 0); i < (end.x < n.x ? end.x : n.x); i++) {

					int cell_idx = i + j * n.x + k * n.x*n.y;

					switch (flag_value) {

					case NF_DIRICHLETPX:
						if ((cell_idx % n.x) != 0 || !(ngbrFlags[cell_idx] & NF_NPX)) continue;
						dirichlet_nx[((cell_idx / n.x) % n.y) + (cell_idx / (n.x*n.y)) * n.y] = value;
						break;

					case NF_DIRICHLETNX:
						if ((cell_idx % n.x) != (n.x - 1) || !(ngbrFlags[cell_idx] & NF_NNX)) continue;
						dirichlet_px[((cell_idx / n.x) % n.y) + (cell_idx / (n.x*n.y)) * n.y] = value;
						break;

					case NF_DIRICHLETPY:
						if (((cell_idx / n.x) % n.y) != 0 || !(ngbrFlags[cell_idx] & NF_NPY)) continue;
						dirichlet_ny[(cell_idx % n.x) + (cell_idx / (n.x*n.y)) * n.x] = value;
						break;

					case NF_DIRICHLETNY:
						if (((cell_idx / n.x) % n.y) != (n.y - 1) || !(ngbrFlags[cell_idx] & NF_NNY)) continue;
						dirichlet_py[(cell_idx % n.x) + (cell_idx / (n.x*n.y)) * n.x] = value;
						break;

					case NF_DIRICHLETPZ:
						if ((cell_idx / (n.x*n.y)) != 0 || !(ngbrFlags[cell_idx] & NF_NPZ)) continue;
						dirichlet_nz[(cell_idx % n.x) + ((cell_idx / n.x) % n.y) * n.x] = value;
						break;

					case NF_DIRICHLETNZ:
						if ((cell_idx / (n.x*n.y)) != (n.z - 1) || !(ngbrFlags[cell_idx] & NF_NNZ)) continue;
						dirichlet_pz[(cell_idx % n.x) + ((cell_idx / n.x) % n.y) * n.x] = value;
						break;
					}

					ngbrFlags[cell_idx] |= flag_value;
				}
			}
		}
	};

	//y-z plane
	if (IsZ(intersection.s.x - intersection.e.x)) {

		//on lower x side
		if (IsZ(rect.s.x - intersection.s.x)) {

			if (!dirichlet_nx.size()) {

				if (!malloc_vector(dirichlet_nx, n.y*n.z)) return false;
			}

			intersection.e.x += h.x;
			set_dirichlet_value(intersection, value, NF_DIRICHLETPX);
		}
		//on upper x side
		else if (IsZ(rect.e.x - intersection.s.x)) {

			if (!dirichlet_px.size()) {

				if (!malloc_vector(dirichlet_px, n.y*n.z)) return false;
			}

			intersection.s.x -= h.x;
			set_dirichlet_value(intersection, value, NF_DIRICHLETNX);
		}
	}
	//x-z plane
	else if (IsZ(intersection.s.y - intersection.e.y)) {

		//on lower y side
		if (IsZ(rect.s.y - intersection.s.y)) {

			if (!dirichlet_ny.size()) {

				if (!malloc_vector(dirichlet_ny, n.x*n.z)) return false;
			}

			intersection.e.y += h.y;
			set_dirichlet_value(intersection, value, NF_DIRICHLETPY);
		}
		//on upper y side
		else if (IsZ(rect.e.y - intersection.s.y)) {

			if (!dirichlet_py.size()) {

				if (!malloc_vector(dirichlet_py, n.x*n.z)) return false;
			}

			intersection.s.y -= h.y;
			set_dirichlet_value(intersection, value, NF_DIRICHLETNY);
		}
	}
	//x-y plane
	else if (IsZ(intersection.s.z - intersection.e.z)) {

		//on lower z side
		if (IsZ(rect.s.z - intersection.s.z)) {

			if (!dirichlet_nz.size()) {

				if (!malloc_vector(dirichlet_nz, n.x*n.y)) return false;
			}

			intersection.e.z += h.z;
			set_dirichlet_value(intersection, value, NF_DIRICHLETPZ);
		}
		//on upper z side
		else if (IsZ(rect.e.z - intersection.s.z)) {

			if (!dirichlet_pz.size()) {

				if (!malloc_vector(dirichlet_pz, n.x*n.y)) return false;
			}

			intersection.s.z -= h.z;
			set_dirichlet_value(intersection, value, NF_DIRICHLETNZ);
		}
	}

	return true;
}

//clear all dirichlet flags and vectors
template <typename VType>
void VEC_VC<VType>::clear_dirichlet_flags(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < (int)ngbrFlags.size(); idx++) {

		ngbrFlags[idx] &= ~NF_DIRICHLET;
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
	Box cells = box_from_rect_max(rectangle);

	for (int i = (cells.s.x >= 0 ? cells.s.x : 0); i < (cells.e.x <= n.x ? cells.e.x : n.x); i++) {
		for (int j = (cells.s.y >= 0 ? cells.s.y : 0); j < (cells.e.y <= n.y ? cells.e.y : n.y); j++) {
			for (int k = (cells.s.z >= 0 ? cells.s.z : 0); k < (cells.e.z <= n.z ? cells.e.z : n.z); k++) {

				if(status) ngbrFlags[i + j * n.x + k * n.x*n.y] |= NF_SKIPCELL;
				else ngbrFlags[i + j * n.x + k * n.x*n.y] &= ~NF_SKIPCELL;
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

	//clear Robin flags
#pragma omp parallel for
	for (int idx = 0; idx < (int)ngbrFlags.size(); idx++)
		ngbrFlags[idx] &= ~NF_ROBIN;
}
