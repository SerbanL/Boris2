#pragma once

#include "cuVEC_VC.h"
#include "cuFuncs_Math.h"
#include "launchers.h"
#include "Reduction.cuh"

//------------------------------------------------------------------- RESIZE NGBRFLAGS

//ngbrFlags size kept : just erase every flag but the shape
__global__ static void resize_ngbrFlags_samesize_kernel(size_t size, int*& ngbrFlags)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		if (ngbrFlags[idx] & NF_NOTEMPTY) ngbrFlags[idx] = NF_NOTEMPTY;
		else ngbrFlags[idx] = 0;
	}
}

//map ngbrFlags to new size keeping only the shape in old_ngbrFlags
__global__ static void resize_ngbrFlags_newsize_kernel(cuSZ3 new_n, cuSZ3 old_n, int*& ngbrFlags, int* old_ngbrFlags)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//now also transfer mesh values to new dimensions
	cuReal3 sourceIdx = (cuReal3)old_n / new_n;

	if (idx < new_n.dim()) {

		int _x = (int)floor((idx % new_n.x) * sourceIdx.x);
		int _y = (int)floor(((idx / new_n.x) % new_n.y) * sourceIdx.y);
		int _z = (int)floor((idx / (new_n.x*new_n.y)) * sourceIdx.z);

		if (old_ngbrFlags[_x + _y * int(old_n.x) + _z * (int(old_n.x*old_n.y))] & NF_NOTEMPTY) ngbrFlags[idx] = NF_NOTEMPTY;
		else ngbrFlags[idx] = 0;
	}
}

template cudaError_t cuVEC_VC<float>::resize_ngbrFlags(cuSZ3 new_n);
template cudaError_t cuVEC_VC<double>::resize_ngbrFlags(cuSZ3 new_n);

template cudaError_t cuVEC_VC<cuFLT3>::resize_ngbrFlags(cuSZ3 new_n);
template cudaError_t cuVEC_VC<cuDBL3>::resize_ngbrFlags(cuSZ3 new_n);

template cudaError_t cuVEC_VC<cuFLT33>::resize_ngbrFlags(cuSZ3 new_n);
template cudaError_t cuVEC_VC<cuDBL33>::resize_ngbrFlags(cuSZ3 new_n);

//set size of ngbrFlags to new_n also mapping shape from current size to new size (if current zero size set solid shape). Memory must be reserved in ngbrFlags to guarantee success. Also n should still have the old value : call this before changing it.
template <typename VType>
__host__ cudaError_t cuVEC_VC<VType>::resize_ngbrFlags(cuSZ3 new_n)
{
	cudaError_t error = cudaSuccess;

	if (!get_ngbrFlags_size()) {

		//if the VEC is being resized from a zero size then set solid shape
		error = set_ngbrFlags_size(new_n.dim());
		if (error == cudaSuccess) gpu_set_managed(ngbrFlags, NF_NOTEMPTY, new_n.dim());
	}
	else {

		if (new_n != get_gpu_value(cuVEC<VType>::n)) {

			//VEC is being resized from non-zero size : map current shape to new size

			//save old flags
			int* old_ngbrFlags = nullptr;

			error = gpu_alloc(old_ngbrFlags, get_gpu_value(cuVEC<VType>::n).dim());
			if (error != cudaSuccess) return error;

			error = gpu_to_gpu_managed2nd(old_ngbrFlags, ngbrFlags, get_gpu_value(cuVEC<VType>::n).dim());
			if (error != cudaSuccess) {

				gpu_free(old_ngbrFlags);
				return error;
			}

			//resize flags to new size and zero
			error = set_ngbrFlags_size(new_n.dim());
			if (error == cudaSuccess) gpu_set_managed(ngbrFlags, 0, new_n.dim());
			else {

				gpu_free(old_ngbrFlags);
				return error;
			}

			//now map shape from ngbrFlags_swapspace to ngbrFlags
			resize_ngbrFlags_newsize_kernel <<< (new_n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (new_n, get_gpu_value(cuVEC<VType>::n), ngbrFlags, old_ngbrFlags);

			//free temporary
			gpu_free(old_ngbrFlags);
		}
		else {

			//VEC is not being resized : clear everything but the shape
			resize_ngbrFlags_samesize_kernel <<< (new_n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (new_n.dim(), ngbrFlags);
		}
	}

	//clear ngbrFlags2 if in use
	if (use_extended_flags()) {

		gpu_alloc_managed(ngbrFlags2, new_n.dim());
		gpu_set_managed(ngbrFlags2, (int)0, new_n.dim());
		set_gpu_value(using_extended_flags, true);
	}

	return error;
}

//------------------------------------------------------------------- SET NGBRFLAGS void

//ngbrFlags size kept : just erase every flag but the shape
template <typename VType>
__global__ void set_ngbrFlags_kernel(const cuSZ3& n, int*& ngbrFlags, VType*& quantity, int& nonempty_cells)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	if (idx < n.dim()) {

		if (ngbrFlags[idx] & NF_NOTEMPTY) {

			atomicAdd(&nonempty_cells, 1);

			if (ijk.i < n.x - 1) {

				//on-axis
				if (ngbrFlags[idx + 1] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NPX; }
				else { ngbrFlags[idx] &= ~NF_NPX; }

				//off-axis (z slice : xy)
				if (ijk.j < n.y - 1) {

					if (ngbrFlags[idx + 1 + n.x] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XY_PXPY; }
					else { ngbrFlags[idx] &= ~NF_XY_PXPY; }
				}

				if (ijk.j > 0) {

					if (ngbrFlags[idx + 1 - n.x] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XY_PXNY; }
					else { ngbrFlags[idx] &= ~NF_XY_PXNY; }
				}
			}

			if (ijk.i > 0) {

				//on-axis
				if (ngbrFlags[idx - 1] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NNX; }
				else { ngbrFlags[idx] &= ~NF_NNX; }

				//off-axis (z slice : xy)
				if (ijk.j < n.y - 1) {

					if (ngbrFlags[idx - 1 + n.x] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XY_NXPY; }
					else { ngbrFlags[idx] &= ~NF_XY_NXPY; }
				}

				if (ijk.j > 0) {

					if (ngbrFlags[idx - 1 - n.x] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XY_NXNY; }
					else { ngbrFlags[idx] &= ~NF_XY_NXNY; }
				}
			}

			if (ijk.j < n.y - 1) {

				//on-axis
				if (ngbrFlags[idx + n.x] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NPY; }
				else { ngbrFlags[idx] &= ~NF_NPY; }

				//off-axis (x slice : yz)
				if (ijk.k < n.z - 1) {

					if (ngbrFlags[idx + n.x + n.x*n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_YZ_PYPZ; }
					else { ngbrFlags[idx] &= ~NF_YZ_PYPZ; }
				}

				if (ijk.k > 0) {

					if (ngbrFlags[idx + n.x - n.x*n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_YZ_PYNZ; }
					else { ngbrFlags[idx] &= ~NF_YZ_PYNZ; }
				}
			}

			if (ijk.j > 0) {

				//on-axis
				if (ngbrFlags[idx - n.x] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NNY; }
				else { ngbrFlags[idx] &= ~NF_NNY; }

				//off-axis (x slice : yz)
				if (ijk.k < n.z - 1) {

					if (ngbrFlags[idx - n.x + n.x*n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_YZ_NYPZ; }
					else { ngbrFlags[idx] &= ~NF_YZ_NYPZ; }
				}

				if (ijk.k > 0) {

					if (ngbrFlags[idx - n.x - n.x*n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_YZ_NYNZ; }
					else { ngbrFlags[idx] &= ~NF_YZ_NYNZ; }
				}
			}

			if (ijk.k < n.z - 1) {

				//on-axis
				if (ngbrFlags[idx + n.x*n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NPZ; }
				else { ngbrFlags[idx] &= ~NF_NPZ; }

				//off-axis (y slice : xz)
				if (ijk.i < n.x - 1) {

					if (ngbrFlags[idx + 1 + n.x*n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XZ_PXPZ; }
					else { ngbrFlags[idx] &= ~NF_XZ_PXPZ; }
				}

				if (ijk.i > 0) {

					if (ngbrFlags[idx - 1 + n.x*n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XZ_NXPZ; }
					else { ngbrFlags[idx] &= ~NF_XZ_NXPZ; }
				}
			}

			if (ijk.k > 0) {

				//on-axis
				if (ngbrFlags[idx - n.x*n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_NNZ; }
				else { ngbrFlags[idx] &= ~NF_NNZ; }

				//off-axis (y slice : xz)
				if (ijk.i < n.x - 1) {

					if (ngbrFlags[idx + 1 - n.x*n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XZ_PXNZ; }
					else { ngbrFlags[idx] &= ~NF_XZ_PXNZ; }
				}

				if (ijk.i > 0) {

					if (ngbrFlags[idx - 1 - n.x*n.y] & NF_NOTEMPTY) { ngbrFlags[idx] |= NF_XZ_NXNZ; }
					else { ngbrFlags[idx] &= ~NF_XZ_NXNZ; }
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

			quantity[idx] = VType();
		}
	}
}

template void cuVEC_VC<float>::set_ngbrFlags(void);
template void cuVEC_VC<double>::set_ngbrFlags(void);

template void cuVEC_VC<cuFLT3>::set_ngbrFlags(void);
template void cuVEC_VC<cuDBL3>::set_ngbrFlags(void);

template void cuVEC_VC<cuFLT33>::set_ngbrFlags(void);
template void cuVEC_VC<cuDBL33>::set_ngbrFlags(void);

//initialization method for neighbor flags : set flags at size n, counting neighbors etc. Use current shape in ngbrFlags
template <typename VType>
__host__ void cuVEC_VC<VType>::set_ngbrFlags(void)
{
	if (get_gpu_value(cuVEC<VType>::rect).IsNull() || get_gpu_value(cuVEC<VType>::h) == cuReal3()) return;

	//clear any shift debt as mesh has been resized so not valid anymore
	set_gpu_value(shift_debt, cuReal3());

	//dirichlet flags will be cleared from ngbrFlags, so also clear the dirichlet vectors
	clear_dirichlet_flags();

	//1. Count all the neighbors
	set_gpu_value(nonempty_cells, (int)0);
	set_ngbrFlags_kernel <<< (get_ngbrFlags_size() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, ngbrFlags, cuVEC<VType>::quantity, nonempty_cells);

	//2. set Robin flags depending on set conditions
	set_robin_flags();

	//3. set pbc flags depending on set conditions and currently calculated flags
	set_pbc_flags();
}

//------------------------------------------------------------------- SET NGBRFLAGS linked cuVEC

//copy shape from linked vec
template <typename VType>
__global__ void set_ngbrFlags_copylinkedshape_kernel(const cuSZ3& n, const cuReal3& h, const cuRect& rect, int*& ngbrFlags, const cuSZ3& linked_n, const cuReal3& linked_h, const cuRect& linked_rect, int*& linked_ngbrFlags)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	auto get_cellrect = [](const cuReal3& h, const cuRect& rect, const cuINT3& ijk) -> cuRect {

		return cuRect(rect.s + (h & ijk), rect.s + (h & ijk) + h);
	};

	auto cellidx_from_position = [](const cuSZ3& n, const cuReal3& h, const cuRect& rect, const cuReal3& absolute_position) -> cuINT3 {

		cuINT3 ijk = cuINT3(
			(int)cu_floor_epsilon((absolute_position.x - rect.s.x) / h.x),
			(int)cu_floor_epsilon((absolute_position.y - rect.s.y) / h.y),
			(int)cu_floor_epsilon((absolute_position.z - rect.s.z) / h.z));

		if (ijk.i < 0) ijk.i = 0;
		if (ijk.j < 0) ijk.j = 0;
		if (ijk.k < 0) ijk.k = 0;

		if (ijk.i > n.x) ijk.i = n.x;
		if (ijk.j > n.y) ijk.j = n.y;
		if (ijk.k > n.z) ijk.k = n.z;

		return ijk;
	};

	auto box_from_rect_max = [&cellidx_from_position](const cuSZ3& n, const cuReal3& h, const cuRect& rect, const cuRect& rectangle) -> cuBox {

		if (!rectangle.intersects(rect)) return cuBox();

		//get start point. this will be limited to 0 to n (inclusive)
		cuINT3 start = cellidx_from_position(n, h, rect, rectangle.s);

		//the Rect could be a plane rectangle on one of the surfaces of this mesh, so adjust start point for this
		if (start.x >= n.x) start.x = n.x - 1;
		if (start.y >= n.y) start.y = n.y - 1;
		if (start.z >= n.z) start.z = n.z - 1;

		//get end point. this will be limited to 0 to n (inclusive)
		cuINT3 end = cellidx_from_position(n, h, rect, rectangle.e);

		cuReal3 snap = (h & end) + rect.s;

		//add 1 since end point must be included in the box, unless the rectangle end point is already at the end of a cell
		if (cuIsNZ(snap.x - rectangle.e.x) && end.x < n.x) end.x++;
		if (cuIsNZ(snap.y - rectangle.e.y) && end.y < n.y) end.y++;
		if (cuIsNZ(snap.z - rectangle.e.z) && end.z < n.z) end.z++;

		return cuBox(start, end);
	};

	auto is_empty = [&box_from_rect_max](const cuSZ3& n, const cuReal3& h, const cuRect& rect, int* ngbrFlags, const cuRect& rectangle) -> bool {

		cuBox cells = box_from_rect_max(n, h, rect, rectangle);

		for (int i = (cells.s.x >= 0 ? cells.s.x : 0); i < (cells.e.x <= n.x ? cells.e.x : n.x); i++) {

			for (int j = (cells.s.y >= 0 ? cells.s.y : 0); j < (cells.e.y <= n.y ? cells.e.y : n.y); j++) {

				for (int k = (cells.s.z >= 0 ? cells.s.z : 0); k < (cells.e.z <= n.z ? cells.e.z : n.z); k++) {

					cuINT3 ijk = cuINT3(i, j, k);
					if (ngbrFlags[ijk.i + ijk.j*n.x + ijk.k*n.x*n.y] & NF_NOTEMPTY) return false;
				}
			}
		}

		return true;
	};

	if (idx < n.dim()) {

		cuRect cellRect = get_cellrect(h, rect, ijk);

		if (!is_empty(linked_n, linked_h, linked_rect, linked_ngbrFlags, cellRect)) {

			//mark cell as not empty
			ngbrFlags[idx] |= NF_NOTEMPTY;
		}
		else {

			//if linked quantity is empty then also empty the quantity in this VEC. Note, the linked vec could actually be the same vec - in this case the shape is not changed since the NF_NOTEMPTY flags have already been mapped to the new size
			ngbrFlags[idx] &= ~NF_NOTEMPTY;
		}
	}
}

template void cuVEC_VC<float>::set_ngbrFlags(const cuSZ3& linked_n, const cuReal3& linked_h, const cuRect& linked_rect, int*& linked_ngbrFlags);
template void cuVEC_VC<double>::set_ngbrFlags(const cuSZ3& linked_n, const cuReal3& linked_h, const cuRect& linked_rect, int*& linked_ngbrFlags);

template void cuVEC_VC<cuFLT3>::set_ngbrFlags(const cuSZ3& linked_n, const cuReal3& linked_h, const cuRect& linked_rect, int*& linked_ngbrFlags);
template void cuVEC_VC<cuDBL3>::set_ngbrFlags(const cuSZ3& linked_n, const cuReal3& linked_h, const cuRect& linked_rect, int*& linked_ngbrFlags);

template void cuVEC_VC<cuFLT33>::set_ngbrFlags(const cuSZ3& linked_n, const cuReal3& linked_h, const cuRect& linked_rect, int*& linked_ngbrFlags);
template void cuVEC_VC<cuDBL33>::set_ngbrFlags(const cuSZ3& linked_n, const cuReal3& linked_h, const cuRect& linked_rect, int*& linked_ngbrFlags);

//initialization method for neighbor flags : set flags at size n, counting neighbors etc. Use current shape in ngbrFlags
template <typename VType>
__host__ void cuVEC_VC<VType>::set_ngbrFlags(const cuSZ3& linked_n, const cuReal3& linked_h, const cuRect& linked_rect, int*& linked_ngbrFlags)
{
	if (get_gpu_value(cuVEC<VType>::rect).IsNull() || get_gpu_value(cuVEC<VType>::h) == cuReal3()) return;

	//copy shape from linked vec (linked_ngbrFlags)
	set_ngbrFlags_copylinkedshape_kernel<VType> <<< (get_ngbrFlags_size() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, cuVEC<VType>::h, cuVEC<VType>::rect, ngbrFlags, linked_n, linked_h, linked_rect, linked_ngbrFlags);

	//now continue with set_ngbrFlags as normal
	set_ngbrFlags();
}

//------------------------------------------------------------------- SET DIRICHLET CONDITIONS

template <typename VType>
__global__ void set_dirichlet_conditions_kernel(
	const cuSZ3& n, int*& ngbrFlags, int*& ngbrFlags2,
	VType*& dirichlet_px, VType*& dirichlet_nx, VType*& dirichlet_py, VType*& dirichlet_ny, VType*& dirichlet_pz, VType*& dirichlet_nz,
	int dirichlet_flag, cuBox flags_box, VType value)
{
	int cell_idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(cell_idx % n.x, (cell_idx / n.x) % n.y, cell_idx / (n.x*n.y));

	if (cell_idx < n.dim() && flags_box.Contains(ijk)) {

		switch (dirichlet_flag) {

		case NF2_DIRICHLETPX:
			if ((cell_idx % n.x) != 0 || !(ngbrFlags[cell_idx] & NF_NPX)) return;
			dirichlet_nx[((cell_idx / n.x) % n.y) + (cell_idx / (n.x*n.y)) * n.y] = value;
			break;

		case NF2_DIRICHLETNX:
			if ((cell_idx % n.x) != (n.x - 1) || !(ngbrFlags[cell_idx] & NF_NNX)) return;
			dirichlet_px[((cell_idx / n.x) % n.y) + (cell_idx / (n.x*n.y)) * n.y] = value;
			break;

		case NF2_DIRICHLETPY:
			if (((cell_idx / n.x) % n.y) != 0 || !(ngbrFlags[cell_idx] & NF_NPY)) return;
			dirichlet_ny[(cell_idx % n.x) + (cell_idx / (n.x*n.y)) * n.x] = value;
			break;

		case NF2_DIRICHLETNY:
			if (((cell_idx / n.x) % n.y) != (n.y - 1) || !(ngbrFlags[cell_idx] & NF_NNY)) return;
			dirichlet_py[(cell_idx % n.x) + (cell_idx / (n.x*n.y)) * n.x] = value;
			break;

		case NF2_DIRICHLETPZ:
			if ((cell_idx / (n.x*n.y)) != 0 || !(ngbrFlags[cell_idx] & NF_NPZ)) return;
			dirichlet_nz[(cell_idx % n.x) + ((cell_idx / n.x) % n.y) * n.x] = value;
			break;

		case NF2_DIRICHLETNZ:
			if ((cell_idx / (n.x*n.y)) != (n.z - 1) || !(ngbrFlags[cell_idx] & NF_NNZ)) return;
			dirichlet_pz[(cell_idx % n.x) + ((cell_idx / n.x) % n.y) * n.x] = value;
			break;
		}

		ngbrFlags2[cell_idx] |= dirichlet_flag;
	}
}

template bool cuVEC_VC<float>::set_dirichlet_conditions(cuRect surface_rect, float value);
template bool cuVEC_VC<double>::set_dirichlet_conditions(cuRect surface_rect, double value);

template bool cuVEC_VC<cuFLT3>::set_dirichlet_conditions(cuRect surface_rect, cuFLT3 value);
template bool cuVEC_VC<cuDBL3>::set_dirichlet_conditions(cuRect surface_rect, cuDBL3 value);

template bool cuVEC_VC<cuFLT33>::set_dirichlet_conditions(cuRect surface_rect, cuFLT33 value);
template bool cuVEC_VC<cuDBL33>::set_dirichlet_conditions(cuRect surface_rect, cuDBL33 value);

//set dirichlet boundary conditions from surface_rect (must be a rectangle intersecting with one of the surfaces of this mesh) and value
////return false on memory allocation failure only, otherwise return true even if surface_rect was not valid
template <typename VType>
__host__ bool cuVEC_VC<VType>::set_dirichlet_conditions(cuRect surface_rect, VType value)
{
	cuReal3 n_ = get_gpu_value(cuVEC<VType>::n);
	cuReal3 h_ = get_gpu_value(cuVEC<VType>::h);
	cuRect rect_ = get_gpu_value(cuVEC<VType>::rect);

	if (!rect_.intersects(surface_rect)) return true;

	cuRect intersection = rect_.get_intersection(surface_rect);
	// if (!intersection.IsPlane()) return true;

	//y-z plane
	if (cuIsZ(intersection.s.x - intersection.e.x)) {

		//on lower x side
		if (cuIsZ(rect_.s.x - intersection.s.x)) {

			if (!get_dirichlet_size(NF2_DIRICHLETNX)) {

				if (set_dirichlet_size(n_.y*n_.z, NF2_DIRICHLETNX) != cudaSuccess) return false;

				//DIRICHLET flags used with extended ngbrFlags so make sure it has memory allocated now
				use_extended_flags();
			}

			intersection.e.x += h_.x;
			cuBox flags_box = cuBox(cuVEC<VType>::box_from_rect_max_cpu(intersection).s, cuVEC<VType>::box_from_rect_max_cpu(intersection).e);
			set_dirichlet_conditions_kernel <<< (n_.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> 
				(cuVEC<VType>::n, ngbrFlags, ngbrFlags2, dirichlet_px, dirichlet_nx, dirichlet_py, dirichlet_ny, dirichlet_pz, dirichlet_nz, NF2_DIRICHLETPX, flags_box, value);
		}
		//on upper x side
		else if (cuIsZ(rect_.e.x - intersection.s.x)) {

			if (!get_dirichlet_size(NF2_DIRICHLETPX)) {

				if (set_dirichlet_size(n_.y*n_.z, NF2_DIRICHLETPX) != cudaSuccess) return false;

				//DIRICHLET flags used with extended ngbrFlags so make sure it has memory allocated now
				use_extended_flags();
			}

			intersection.s.x -= h_.x;
			cuBox flags_box = cuBox(cuVEC<VType>::box_from_rect_max_cpu(intersection).s, cuVEC<VType>::box_from_rect_max_cpu(intersection).e);
			set_dirichlet_conditions_kernel <<< (n_.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
				(cuVEC<VType>::n, ngbrFlags, ngbrFlags2, dirichlet_px, dirichlet_nx, dirichlet_py, dirichlet_ny, dirichlet_pz, dirichlet_nz, NF2_DIRICHLETNX, flags_box, value);
		}
	}
	//x-z plane
	else if (cuIsZ(intersection.s.y - intersection.e.y)) {

		//on lower y side
		if (cuIsZ(rect_.s.y - intersection.s.y)) {

			if (!get_dirichlet_size(NF2_DIRICHLETNY)) {

				if (set_dirichlet_size(n_.x*n_.z, NF2_DIRICHLETNY) != cudaSuccess) return false;

				//DIRICHLET flags used with extended ngbrFlags so make sure it has memory allocated now
				use_extended_flags();
			}
			
			intersection.e.y += h_.y;
			cuBox flags_box = cuBox(cuVEC<VType>::box_from_rect_max_cpu(intersection).s, cuVEC<VType>::box_from_rect_max_cpu(intersection).e);
			set_dirichlet_conditions_kernel <<< (n_.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> 
				(cuVEC<VType>::n, ngbrFlags, ngbrFlags2, dirichlet_px, dirichlet_nx, dirichlet_py, dirichlet_ny, dirichlet_pz, dirichlet_nz, NF2_DIRICHLETPY, flags_box, value);
		}
		//on upper y side
		else if (cuIsZ(rect_.e.y - intersection.s.y)) {

			if (!get_dirichlet_size(NF2_DIRICHLETPY)) {

				if (set_dirichlet_size(n_.x*n_.z, NF2_DIRICHLETPY) != cudaSuccess) return false;

				//DIRICHLET flags used with extended ngbrFlags so make sure it has memory allocated now
				use_extended_flags();
			}
			
			intersection.s.y -= h_.y;
			cuBox flags_box = cuBox(cuVEC<VType>::box_from_rect_max_cpu(intersection).s, cuVEC<VType>::box_from_rect_max_cpu(intersection).e);
			set_dirichlet_conditions_kernel <<< (n_.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> 
				(cuVEC<VType>::n, ngbrFlags, ngbrFlags2, dirichlet_px, dirichlet_nx, dirichlet_py, dirichlet_ny, dirichlet_pz, dirichlet_nz, NF2_DIRICHLETNY, flags_box, value);
		}
	}
	//x-y plane
	else if (cuIsZ(intersection.s.z - intersection.e.z)) {

		//on lower z side
		if (cuIsZ(rect_.s.z - intersection.s.z)) {

			if (!get_dirichlet_size(NF2_DIRICHLETNZ)) {

				if (set_dirichlet_size(n_.x*n_.y, NF2_DIRICHLETNZ) != cudaSuccess) return false;

				//DIRICHLET flags used with extended ngbrFlags so make sure it has memory allocated now
				use_extended_flags();
			}
			
			intersection.e.z += h_.z;
			cuBox flags_box = cuBox(cuVEC<VType>::box_from_rect_max_cpu(intersection).s, cuVEC<VType>::box_from_rect_max_cpu(intersection).e);
			set_dirichlet_conditions_kernel <<< (n_.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> 
				(cuVEC<VType>::n, ngbrFlags, ngbrFlags2, dirichlet_px, dirichlet_nx, dirichlet_py, dirichlet_ny, dirichlet_pz, dirichlet_nz, NF2_DIRICHLETPZ, flags_box, value);
		}
		//on upper z side
		else if (cuIsZ(rect_.e.z - intersection.s.z)) {

			if (!get_dirichlet_size(NF2_DIRICHLETPZ)) {

				if (set_dirichlet_size(n_.x*n_.y, NF2_DIRICHLETPZ) != cudaSuccess) return false;

				//DIRICHLET flags used with extended ngbrFlags so make sure it has memory allocated now
				use_extended_flags();
			}
			
			intersection.s.z -= h_.z;
			cuBox flags_box = cuBox(cuVEC<VType>::box_from_rect_max_cpu(intersection).s, cuVEC<VType>::box_from_rect_max_cpu(intersection).e);
			set_dirichlet_conditions_kernel <<< (n_.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> 
				(cuVEC<VType>::n, ngbrFlags, ngbrFlags2, dirichlet_px, dirichlet_nx, dirichlet_py, dirichlet_ny, dirichlet_pz, dirichlet_nz, NF2_DIRICHLETNZ, flags_box, value);
		}
	}

	return true;
}

//------------------------------------------------------------------- SET ROBIN FLAGS

__global__ static void set_robin_flags_kernel(
	const cuSZ3& n, int*& ngbrFlags, int*& ngbrFlags2,
	const cuReal2& robin_px, const cuReal2& robin_nx,
	const cuReal2& robin_py, const cuReal2& robin_ny,
	const cuReal2& robin_pz, const cuReal2& robin_nz,
	const cuReal2& robin_v)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	if (idx < n.dim()) {

		//first clear any robin flags already set
		ngbrFlags2[idx] &= ~NF2_ROBIN;

		if (ngbrFlags[idx] & NF_NOTEMPTY) {

			//neighbors
			if (ijk.i < n.x - 1) {

				if (!(ngbrFlags[idx + 1] & NF_NOTEMPTY))
					//inner cell next to a void cell on the -x side
					if (cuIsNZ(robin_v.i) && ijk.i > 0 && ngbrFlags[idx - 1] & NF_NOTEMPTY) { ngbrFlags2[idx] |= (NF2_ROBINNX + NF2_ROBINV); }
			}
			//surface cell on the -x side of the surface
			else if (cuIsNZ(robin_px.i) && ijk.i > 0 && (ngbrFlags[idx - 1] & NF_NOTEMPTY)) ngbrFlags2[idx] |= NF2_ROBINNX;

			if (ijk.i > 0) {

				if (!(ngbrFlags[idx - 1] & NF_NOTEMPTY))
					if (cuIsNZ(robin_v.i) && ijk.i < n.x - 1 && ngbrFlags[idx + 1] & NF_NOTEMPTY) { ngbrFlags2[idx] |= (NF2_ROBINPX + NF2_ROBINV); }
			}
			else if (cuIsNZ(robin_nx.i) && ijk.i < n.x - 1 && (ngbrFlags[idx + 1] & NF_NOTEMPTY)) ngbrFlags2[idx] |= NF2_ROBINPX;

			if (ijk.j < n.y - 1) {

				if (!(ngbrFlags[idx + n.x] & NF_NOTEMPTY))
					if (cuIsNZ(robin_v.i) && ijk.j > 0 && (ngbrFlags[idx - n.x] & NF_NOTEMPTY)) { ngbrFlags2[idx] |= (NF2_ROBINNY + NF2_ROBINV); }
			}
			else if (cuIsNZ(robin_py.i) && ijk.j > 0 && (ngbrFlags[idx - n.x] & NF_NOTEMPTY)) ngbrFlags2[idx] |= NF2_ROBINNY;

			if (ijk.j > 0) {

				if (!(ngbrFlags[idx - n.x] & NF_NOTEMPTY))
					if (cuIsNZ(robin_v.i) && ijk.j < n.y - 1 && (ngbrFlags[idx + n.x] & NF_NOTEMPTY)) { ngbrFlags2[idx] |= (NF2_ROBINPY + NF2_ROBINV); }
			}
			else if (cuIsNZ(robin_ny.i) && ijk.j < n.y - 1 && (ngbrFlags[idx + n.x] & NF_NOTEMPTY)) ngbrFlags2[idx] |= NF2_ROBINPY;

			if (ijk.k < n.z - 1) {

				if (!(ngbrFlags[idx + n.x*n.y] & NF_NOTEMPTY))
					if (cuIsNZ(robin_v.i) && ijk.k > 0 && (ngbrFlags[idx - n.x*n.y] & NF_NOTEMPTY)) { ngbrFlags2[idx] |= (NF2_ROBINNZ + NF2_ROBINV); }
			}
			else if (cuIsNZ(robin_pz.i) && ijk.k > 0 && (ngbrFlags[idx - n.x*n.y] & NF_NOTEMPTY)) ngbrFlags2[idx] |= NF2_ROBINNZ;

			if (ijk.k > 0) {

				if (!(ngbrFlags[idx - n.x*n.y] & NF_NOTEMPTY))
					if (cuIsNZ(robin_v.i) && ijk.k < n.z - 1 && (ngbrFlags[idx + n.x*n.y] & NF_NOTEMPTY)) { ngbrFlags2[idx] |= (NF2_ROBINPZ + NF2_ROBINV); }
			}
			else if (cuIsNZ(robin_nz.i) && ijk.k < n.z - 1 && (ngbrFlags[idx + n.x*n.y] & NF_NOTEMPTY)) ngbrFlags2[idx] |= NF2_ROBINPZ;
		}
	}
}

template void cuVEC_VC<float>::set_robin_flags(void);
template void cuVEC_VC<double>::set_robin_flags(void);

template void cuVEC_VC<cuFLT3>::set_robin_flags(void);
template void cuVEC_VC<cuDBL3>::set_robin_flags(void);

template void cuVEC_VC<cuFLT33>::set_robin_flags(void);
template void cuVEC_VC<cuDBL33>::set_robin_flags(void);

//set robin flags from robin values and shape. Doesn't affect any other flags. Call from set_ngbrFlags after counting neighbors, and after setting robin values
template <typename VType>
__host__ void cuVEC_VC<VType>::set_robin_flags(void)
{
	if (use_extended_flags()) {

		set_robin_flags_kernel << < (get_gpu_value(cuVEC<VType>::n).dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuVEC<VType>::n, ngbrFlags, ngbrFlags2, robin_px, robin_nx, robin_py, robin_ny, robin_pz, robin_nz, robin_v);
	}
}

//------------------------------------------------------------------- SET PBC FLAGS

__global__ static void set_pbc_flags_kernel(
	const cuSZ3& n, int*& ngbrFlags,
	const int& pbc_x, const int& pbc_y, const int& pbc_z)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	if (idx < n.dim()) {

		if (ngbrFlags[idx] & NF_NOTEMPTY) {

			//first clear pbc in this cell before recalculating pbc flags
			ngbrFlags[idx] &= ~NF_PBC;

			//-x side, and on the +x side there is a non-empty cell : set pbc
			if (pbc_x) {

				if (ijk.i == 0 && (ngbrFlags[n.x - 1 + ijk.j * n.x + ijk.k * n.x*n.y] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_PBCX;

				//+x side, and on the -x side there is a non-empty cell : set pbc
				if (ijk.i == n.x - 1 && (ngbrFlags[ijk.j * n.x + ijk.k * n.x*n.y] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_PBCX;
			}

			if (pbc_y) {

				//-y side, and on the +y side there is a non-empty cell : set pbc
				if (ijk.j == 0 && (ngbrFlags[ijk.i + (n.y - 1) * n.x + ijk.k * n.x*n.y] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_PBCY;

				//+y side, and on the -y side there is a non-empty cell : set pbc
				if (ijk.j == n.y - 1 && (ngbrFlags[ijk.i + ijk.k * n.x*n.y] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_PBCY;
			}

			if (pbc_z) {
			
				//-z side, and on the +z side there is a non-empty cell : set pbc
				if (ijk.k == 0 && (ngbrFlags[ijk.i + ijk.j * n.x + (n.z - 1) * n.x*n.y] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_PBCZ;

				//+z side, and on the -z side there is a non-empty cell : set pbc
				if (ijk.k == n.z - 1 && (ngbrFlags[ijk.i + ijk.j * n.x] & NF_NOTEMPTY)) ngbrFlags[idx] |= NF_PBCZ;
			}
		}
	}
}

template void cuVEC_VC<float>::set_pbc_flags(void);
template void cuVEC_VC<double>::set_pbc_flags(void);

template void cuVEC_VC<cuFLT3>::set_pbc_flags(void);
template void cuVEC_VC<cuDBL3>::set_pbc_flags(void);

template void cuVEC_VC<cuFLT33>::set_pbc_flags(void);
template void cuVEC_VC<cuDBL33>::set_pbc_flags(void);

//set pbc flags depending on set conditions and currently calculated flags - ngbrFlags must already be calculated before using this
template <typename VType>
__host__ void cuVEC_VC<VType>::set_pbc_flags(void)
{
	set_pbc_flags_kernel <<< (get_gpu_value(cuVEC<VType>::n).dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, ngbrFlags, pbc_x, pbc_y, pbc_z);
}

//------------------------------------------------------------------- SET SKIP CELLS FLAGS

__global__ static void set_skipcells_kernel(const cuSZ3& n, int*& ngbrFlags, cuBox cells, bool status)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	if (idx < n.dim() && cells.Contains(ijk)) {

		if (status) ngbrFlags[idx] |= NF_SKIPCELL;
		else ngbrFlags[idx]  &= ~NF_SKIPCELL;
	}
}

template void cuVEC_VC<float>::set_skipcells(cuRect rectangle, bool status);
template void cuVEC_VC<double>::set_skipcells(cuRect rectangle, bool status);

template void cuVEC_VC<cuFLT3>::set_skipcells(cuRect rectangle, bool status);
template void cuVEC_VC<cuDBL3>::set_skipcells(cuRect rectangle, bool status);

template void cuVEC_VC<cuFLT33>::set_skipcells(cuRect rectangle, bool status);
template void cuVEC_VC<cuDBL33>::set_skipcells(cuRect rectangle, bool status);

//mark cells included in this rectangle (absolute coordinates) to be skipped during some computations
template <typename VType>
__host__ void cuVEC_VC<VType>::set_skipcells(cuRect rectangle, bool status)
{
	cuBox cells = cuVEC<VType>::box_from_rect_max_cpu(rectangle);

	set_skipcells_kernel <<< (get_gpu_value(cuVEC<VType>::n).dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, ngbrFlags, cells, status);
}