#pragma once

#include "cuVEC_VC.h"
#include "cuFuncs_Math.h"
#include "launchers.h"
#include "cuVEC_aux.cuh"

template <typename VType = cuBReal>
__global__ void normalize_error_kernel(cuBReal& max_error, cuBReal& max_val)
{
	if (threadIdx.x == 0) {

		if (max_val) max_error /= max_val;
	}
}

//---------------------------------------------------- LAPLACE

//-------------------- GLOBAL RED

template <typename VType>
__global__ void IterateLaplace_SOR_red_kernel(cuVEC_VC<VType>& vec, cuBReal& damping)
{
	vec.IterateLaplace_SOR_red(damping);
}

//-------------------- GLOBAL BLACK

template <typename VType>
__global__ void IterateLaplace_SOR_black_kernel(cuVEC_VC<VType>& vec, cuBReal& damping, cuBReal& max_error, cuBReal& max_val)
{
	vec.IterateLaplace_SOR_black(damping, max_error, max_val);
}

//-------------------- DEVICE RED

template <typename VType>
__device__ void cuVEC_VC<VType>::IterateLaplace_SOR_red(cuBReal damping)
{
	//this method must be called with half-size : arr_size / 2 = n.dim() / 2, i.e. <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
	//double idx : idx values will now take on even values
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) * 2;
	
	//ijk coordinates corresponding to idx
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));
	
	//for red-black passes must adjust the even idx values so we keep to a 3D checkerboard pattern
	//"red" : start at 0, "black" : start at 1
	//The following rules can be easily checked (remember i is the column number, j is the row number, k is the plane number)
	
	//For red squares:
	//If n.x is even : nudge idx by +1 on a) odd rows and even planes, and b) even rows and odd planes
	//If n.x is odd and n.y is even, nudge idx by +1 on odd planes only
	//If n.x is odd and n.y is odd, no nudge is needed : idx is automatically on the right 3D checkerboard pattern.

	//For black squares:
	//If n.x is even : nudge idx by +1 on a) even rows and even planes, and b) odd rows and odd planes
	//If n.x is odd and n.y is even, nudge idx by +1 on even planes only
	//If n.x is odd and n.y is odd, nudge by +1 everywhere : idx is automatically on the red 3D checkerboard pattern. (so nudge by 1 to find black checkerboard pattern).
	
	if (n.x % 2 == 1) { 
		if (n.y % 2 == 0) {

			//nx is odd and ny is even : for red squares nudge on odd planes only
			idx += (int)(n.z % 2);
		}
		//else : nx is odd and ny is odd, no nudge is needed
	}
	else {
		
		//nx is even : for red squares nudge on a) odd rows and even planes, and b) even rows and odd planes

		//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
		bool red_nudge = (((ijk.j % 2) == 1 && (ijk.k % 2) == 0) || (((ijk.j % 2) == 0 && (ijk.k % 2) == 1)));

		idx += (int)red_nudge;
	}

	//calculate new value only in non-empty cells; also skip if indicated as a composite media boundary condition cell
	bool calculate_idx = idx < n.dim() && (!(ngbrFlags[idx] & NF_CMBND) && (ngbrFlags[idx] & NF_NOTEMPTY));

	if (calculate_idx) {

		//get maximum cell side
		cuBReal h_max = cu_maximum(h.x, h.y, h.z);

		//get weights
		cuBReal w_x = (h_max / h.x) * (h_max / h.x);
		cuBReal w_y = (h_max / h.y) * (h_max / h.y);
		cuBReal w_z = (h_max / h.z) * (h_max / h.z);

		VType weighted_sum = VType();
		cuBReal total_weight = 0;

		//x direction
		if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

			total_weight += 2 * w_x;
			weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			total_weight += 6 * w_x;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

				weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
			}
			else {

				weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRX) {

			total_weight += w_x;

			if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * quantity[idx + 1];
			else						 weighted_sum += w_x * quantity[idx - 1];
		}

		//y direction
		if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

			total_weight += 2 * w_y;
			weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			total_weight += 6 * w_y;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

				weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
			}
			else {

				weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRY) {

			total_weight += w_y;

			if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * quantity[idx + n.x];
			else						 weighted_sum += w_y * quantity[idx - n.x];
		}

		//z direction
		if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

			total_weight += 2 * w_z;
			weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

			total_weight += 6 * w_z;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

				weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
			}
			else {

				weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRZ) {

			total_weight += w_z;

			if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * quantity[idx + n.x*n.y];
			else						 weighted_sum += w_z * quantity[idx - n.x*n.y];
		}

		//old_value = quantity[idx];
		quantity[idx] = quantity[idx] * (1.0 - damping) + damping * (weighted_sum / total_weight);
	}
}

//-------------------- DEVICE BLACK

template <typename VType>
__device__ void cuVEC_VC<VType>::IterateLaplace_SOR_black(cuBReal damping, cuBReal& max_error, cuBReal& max_val)
{
	//this method must be called with half-size : arr_size / 2 = n.dim() / 2, i.e. <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
	//double idx : idx values will now take on even values
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) * 2;

	//ijk coordinates corresponding to idx
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	//for red-black passes must adjust the even idx values so we keep to a 3D checkerboard pattern
	//"red" : start at 0, "black" : start at 1
	//The following rules can be easily checked (remember i is the column number, j is the row number, k is the plane number)

	//For red squares:
	//If n.x is even : nudge idx by +1 on a) odd rows and even planes, and b) even rows and odd planes
	//If n.x is odd and n.y is even, nudge idx by +1 on odd planes only
	//If n.x is odd and n.y is odd, no nudge is needed : idx is automatically on the right 3D checkerboard pattern.

	//For black squares:
	//If n.x is even : nudge idx by +1 on a) even rows and even planes, and b) odd rows and odd planes
	//If n.x is odd and n.y is even, nudge idx by +1 on even planes only
	//If n.x is odd and n.y is odd, nudge by +1 everywhere : idx is automatically on the red 3D checkerboard pattern. (so nudge by 1 to find black checkerboard pattern).

	if (n.x % 2 == 1) {
		if (n.y % 2 == 0) {

			//nx is odd and ny is even : for black squares nudge on even planes only
			idx += (int)(n.z % 2 == 0);
		}
		else {

			//nx is odd and ny is odd, nudge everything by 1 for black squares
			idx++;
		}
	}
	else {
	
		//nx is even : for black squares nudge on a) even rows and even planes, and b) odd rows and odd planes

		//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
		bool black_nudge = (((ijk.j % 2) == 0 && (ijk.k % 2) == 0) || (((ijk.j % 2) == 1 && (ijk.k % 2) == 1)));

		idx += (int)black_nudge;
	}

	//calculate new value only in non-empty cells; also skip if indicated as a composite media boundary condition cell
	bool calculate_idx = idx < n.dim() && (!(ngbrFlags[idx] & NF_CMBND) && (ngbrFlags[idx] & NF_NOTEMPTY));

	VType old_value = VType();

	if (calculate_idx) {

		//get maximum cell side
		cuBReal h_max = cu_maximum(h.x, h.y, h.z);

		//get weights
		cuBReal w_x = (h_max / h.x) * (h_max / h.x);
		cuBReal w_y = (h_max / h.y) * (h_max / h.y);
		cuBReal w_z = (h_max / h.z) * (h_max / h.z);

		VType weighted_sum = VType();
		cuBReal total_weight = 0;

		//x direction
		if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

			total_weight += 2 * w_x;
			weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			total_weight += 6 * w_x;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

				weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
			}
			else {

				weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRX) {

			total_weight += w_x;

			if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * quantity[idx + 1];
			else						 weighted_sum += w_x * quantity[idx - 1];
		}

		//y direction
		if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

			total_weight += 2 * w_y;
			weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			total_weight += 6 * w_y;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

				weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
			}
			else {

				weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRY) {

			total_weight += w_y;

			if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * quantity[idx + n.x];
			else						 weighted_sum += w_y * quantity[idx - n.x];
		}

		//z direction
		if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

			total_weight += 2 * w_z;
			weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

			total_weight += 6 * w_z;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

				weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
			}
			else {

				weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRZ) {

			total_weight += w_z;

			if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * quantity[idx + n.x*n.y];
			else						 weighted_sum += w_z * quantity[idx - n.x*n.y];
		}

		old_value = quantity[idx];
		quantity[idx] = quantity[idx] * (1.0 - damping) + damping * (weighted_sum / total_weight);
	}

	reduction_delta(idx, n.dim(), quantity, old_value, max_error, calculate_idx);
	reduction_delta(idx, n.dim(), quantity, VType(), max_val, calculate_idx);
}

//-------------------- LAUNCHER

template <typename VType>
__host__ void cuVEC_VC<VType>::IterateLaplace_SOR(size_t arr_size, cuBReal& damping, cuBReal& max_error, cuBReal& max_val)
{
	IterateLaplace_SOR_red_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, damping);

	IterateLaplace_SOR_black_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, damping, max_error, max_val);
}

//---------------------------------------------------- POISSON

//-------------------- GLOBAL RED

template <typename VType, typename Class_Poisson_RHS>
__global__ void IteratePoisson_SOR_red_kernel(cuVEC_VC<VType>& vec, Class_Poisson_RHS& obj, cuBReal& damping)
{
	vec.IteratePoisson_SOR_red(obj, damping);
}

//-------------------- GLOBAL BLACK

template <typename VType, typename Class_Poisson_RHS>
__global__ void IteratePoisson_SOR_black_kernel(cuVEC_VC<VType>& vec, Class_Poisson_RHS& obj, cuBReal& damping, cuBReal& max_error, cuBReal& max_val)
{
	vec.IteratePoisson_SOR_black(obj, damping, max_error, max_val);
}

//-------------------- DEVICE RED

template <typename VType>
template <typename Class_Poisson_RHS>
__device__ void cuVEC_VC<VType>::IteratePoisson_SOR_red(Class_Poisson_RHS& obj, cuBReal damping)
{
	//this method must be called with half-size : arr_size / 2 = n.dim() / 2, i.e. <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
	//double idx : idx values will now take on even values
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) * 2;

	//ijk coordinates corresponding to idx
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	//for red-black passes must adjust the even idx values so we keep to a 3D checkerboard pattern
	//"red" : start at 0, "black" : start at 1
	//The following rules can be easily checked (remember i is the column number, j is the row number, k is the plane number)

	//For red squares:
	//If n.x is even : nudge idx by +1 on a) odd rows and even planes, and b) even rows and odd planes
	//If n.x is odd and n.y is even, nudge idx by +1 on odd planes only
	//If n.x is odd and n.y is odd, no nudge is needed : idx is automatically on the right 3D checkerboard pattern.

	//For black squares:
	//If n.x is even : nudge idx by +1 on a) even rows and even planes, and b) odd rows and odd planes
	//If n.x is odd and n.y is even, nudge idx by +1 on even planes only
	//If n.x is odd and n.y is odd, nudge by +1 everywhere : idx is automatically on the red 3D checkerboard pattern. (so nudge by 1 to find black checkerboard pattern).

	if (n.x % 2 == 1) {
		if (n.y % 2 == 0) {

			//nx is odd and ny is even : for red squares nudge on odd planes only
			idx += (int)(n.z % 2);
		}
		//else : nx is odd and ny is odd, no nudge is needed
	}
	else {

		//nx is even : for red squares nudge on a) odd rows and even planes, and b) even rows and odd planes

		//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
		bool red_nudge = (((ijk.j % 2) == 1 && (ijk.k % 2) == 0) || (((ijk.j % 2) == 0 && (ijk.k % 2) == 1)));

		idx += (int)red_nudge;
	}

	//calculate new value only in non-empty cells; also skip if indicated as a composite media boundary condition cell
	bool calculate_idx = idx < n.dim() && (!(ngbrFlags[idx] & NF_CMBND) && (ngbrFlags[idx] & NF_NOTEMPTY));

	if (calculate_idx) {

		//get maximum cell side
		cuBReal h_max = cu_maximum(h.x, h.y, h.z);

		//get weights
		cuBReal w_x = (h_max / h.x) * (h_max / h.x);
		cuBReal w_y = (h_max / h.y) * (h_max / h.y);
		cuBReal w_z = (h_max / h.z) * (h_max / h.z);

		VType weighted_sum = VType();
		cuBReal total_weight = 0;

		//x direction
		if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

			total_weight += 2 * w_x;
			weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			total_weight += 6 * w_x;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

				weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
			}
			else {

				weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRX) {

			total_weight += w_x;

			if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * quantity[idx + 1];
			else						 weighted_sum += w_x * quantity[idx - 1];
		}

		//y direction
		if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

			total_weight += 2 * w_y;
			weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			total_weight += 6 * w_y;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

				weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
			}
			else {

				weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRY) {

			total_weight += w_y;

			if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * quantity[idx + n.x];
			else						 weighted_sum += w_y * quantity[idx - n.x];
		}

		//z direction
		if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

			total_weight += 2 * w_z;
			weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
		}
		else if (using_extended_flags && ngbrFlags2[idx] & NF2_DIRICHLETZ) {

			total_weight += 6 * w_z;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

				weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
			}
			else {

				weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRZ) {

			total_weight += w_z;

			if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * quantity[idx + n.x*n.y];
			else						 weighted_sum += w_z * quantity[idx - n.x*n.y];
		}

		quantity[idx] = quantity[idx] * (1.0 - damping) + damping * ((weighted_sum - h_max*h_max * obj.Poisson_RHS(idx)) / total_weight);
	}
}

//-------------------- DEVICE BLACK

template <typename VType>
template <typename Class_Poisson_RHS>
__device__ void cuVEC_VC<VType>::IteratePoisson_SOR_black(Class_Poisson_RHS& obj, cuBReal damping, cuBReal& max_error, cuBReal& max_val)
{
	//this method must be called with half-size : arr_size / 2 = n.dim() / 2, i.e. <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
	//double idx : idx values will now take on even values
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) * 2;

	//ijk coordinates corresponding to idx
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	//for red-black passes must adjust the even idx values so we keep to a 3D checkerboard pattern
	//"red" : start at 0, "black" : start at 1
	//The following rules can be easily checked (remember i is the column number, j is the row number, k is the plane number)

	//For red squares:
	//If n.x is even : nudge idx by +1 on a) odd rows and even planes, and b) even rows and odd planes
	//If n.x is odd and n.y is even, nudge idx by +1 on odd planes only
	//If n.x is odd and n.y is odd, no nudge is needed : idx is automatically on the right 3D checkerboard pattern.

	//For black squares:
	//If n.x is even : nudge idx by +1 on a) even rows and even planes, and b) odd rows and odd planes
	//If n.x is odd and n.y is even, nudge idx by +1 on even planes only
	//If n.x is odd and n.y is odd, nudge by +1 everywhere : idx is automatically on the red 3D checkerboard pattern. (so nudge by 1 to find black checkerboard pattern).

	if (n.x % 2 == 1) {
		if (n.y % 2 == 0) {

			//nx is odd and ny is even : for black squares nudge on even planes only
			idx += (int)(n.z % 2 == 0);
		}
		else {

			//nx is odd and ny is odd, nudge everything by 1 for black squares
			idx++;
		}
	}
	else {

		//nx is even : for black squares nudge on a) even rows and even planes, and b) odd rows and odd planes

		//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
		bool black_nudge = (((ijk.j % 2) == 0 && (ijk.k % 2) == 0) || (((ijk.j % 2) == 1 && (ijk.k % 2) == 1)));

		idx += (int)black_nudge;
	}

	//calculate new value only in non-empty cells; also skip if indicated as a composite media boundary condition cell
	bool calculate_idx = idx < n.dim() && (!(ngbrFlags[idx] & NF_CMBND) && (ngbrFlags[idx] & NF_NOTEMPTY));

	VType old_value = VType();

	if (calculate_idx) {

		//get maximum cell side
		cuBReal h_max = cu_maximum(h.x, h.y, h.z);

		//get weights
		cuBReal w_x = (h_max / h.x) * (h_max / h.x);
		cuBReal w_y = (h_max / h.y) * (h_max / h.y);
		cuBReal w_z = (h_max / h.z) * (h_max / h.z);

		VType weighted_sum = VType();
		cuBReal total_weight = 0;

		//x direction
		if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

			total_weight += 2 * w_x;
			weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			total_weight += 6 * w_x;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

				weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
			}
			else {

				weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRX) {

			total_weight += w_x;

			if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * quantity[idx + 1];
			else						 weighted_sum += w_x * quantity[idx - 1];
		}

		//y direction
		if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

			total_weight += 2 * w_y;
			weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			total_weight += 6 * w_y;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

				weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
			}
			else {

				weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRY) {

			total_weight += w_y;

			if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * quantity[idx + n.x];
			else						 weighted_sum += w_y * quantity[idx - n.x];
		}

		//z direction
		if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

			total_weight += 2 * w_z;
			weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

			total_weight += 6 * w_z;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

				weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
			}
			else {

				weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRZ) {

			total_weight += w_z;

			if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * quantity[idx + n.x*n.y];
			else						 weighted_sum += w_z * quantity[idx - n.x*n.y];
		}

		old_value = quantity[idx];
		quantity[idx] = quantity[idx] * (1.0 - damping) + damping * ((weighted_sum - h_max * h_max * obj.Poisson_RHS(idx)) / total_weight);
	}

	reduction_delta(idx, n.dim(), quantity, old_value, max_error, calculate_idx);
	reduction_delta(idx, n.dim(), quantity, VType(), max_val, calculate_idx);
}

//-------------------- LAUNCHER

template <typename VType>
template <typename Class_Poisson_RHS>
__host__ void cuVEC_VC<VType>::IteratePoisson_SOR(size_t arr_size, Class_Poisson_RHS& obj, cuBReal& damping, cuBReal& max_error, cuBReal& max_val)
{
	IteratePoisson_SOR_red_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, obj, damping);

	IteratePoisson_SOR_black_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, obj, damping, max_error, max_val);
}

//---------------------------------------------------- POISSON with non-homogeneous Neumann boundaries

//-------------------- GLOBAL RED

template <typename VType, typename Class_Poisson_NNeu>
__global__ void IteratePoisson_NNeu_SOR_red_kernel(cuVEC_VC<VType>& vec, Class_Poisson_NNeu& obj, cuBReal& damping)
{
	vec.IteratePoisson_NNeu_SOR_red(obj, damping);
}

//-------------------- GLOBAL BLACK

template <typename VType, typename Class_Poisson_NNeu>
__global__ void IteratePoisson_NNeu_SOR_black_kernel(cuVEC_VC<VType>& vec, Class_Poisson_NNeu& obj, cuBReal& damping, cuBReal& max_error, cuBReal& max_val)
{
	vec.IteratePoisson_NNeu_SOR_black(obj, damping, max_error, max_val);
}

//-------------------- DEVICE RED

template <typename VType>
template <typename Class_Poisson_NNeu>
__device__ void cuVEC_VC<VType>::IteratePoisson_NNeu_SOR_red(Class_Poisson_NNeu& obj, cuBReal damping)
{
	//this method must be called with half-size : arr_size / 2 = n.dim() / 2, i.e. <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
	//double idx : idx values will now take on even values
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) * 2;

	//ijk coordinates corresponding to idx
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	//for red-black passes must adjust the even idx values so we keep to a 3D checkerboard pattern
	//"red" : start at 0, "black" : start at 1
	//The following rules can be easily checked (remember i is the column number, j is the row number, k is the plane number)

	//For red squares:
	//If n.x is even : nudge idx by +1 on a) odd rows and even planes, and b) even rows and odd planes
	//If n.x is odd and n.y is even, nudge idx by +1 on odd planes only
	//If n.x is odd and n.y is odd, no nudge is needed : idx is automatically on the right 3D checkerboard pattern.

	//For black squares:
	//If n.x is even : nudge idx by +1 on a) even rows and even planes, and b) odd rows and odd planes
	//If n.x is odd and n.y is even, nudge idx by +1 on even planes only
	//If n.x is odd and n.y is odd, nudge by +1 everywhere : idx is automatically on the red 3D checkerboard pattern. (so nudge by 1 to find black checkerboard pattern).

	if (n.x % 2 == 1) {
		if (n.y % 2 == 0) {

			//nx is odd and ny is even : for red squares nudge on odd planes only
			idx += (int)(n.z % 2);
		}
		//else : nx is odd and ny is odd, no nudge is needed
	}
	else {

		//nx is even : for red squares nudge on a) odd rows and even planes, and b) even rows and odd planes

		//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
		bool red_nudge = (((ijk.j % 2) == 1 && (ijk.k % 2) == 0) || (((ijk.j % 2) == 0 && (ijk.k % 2) == 1)));

		idx += (int)red_nudge;
	}

	//calculate new value only in non-empty cells; also skip if indicated as a composite media boundary condition cell
	bool calculate_idx = idx < n.dim() && (!(ngbrFlags[idx] & NF_CMBND) && (ngbrFlags[idx] & NF_NOTEMPTY));

	if (calculate_idx) {

		//get maximum cell side
		cuBReal h_max = cu_maximum(h.x, h.y, h.z);

		//get weights
		cuBReal w_x = (h_max / h.x) * (h_max / h.x);
		cuBReal w_y = (h_max / h.y) * (h_max / h.y);
		cuBReal w_z = (h_max / h.z) * (h_max / h.z);

		VType weighted_sum = VType();
		cuBReal total_weight = 0;

		//x direction
		if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

			total_weight += 2 * w_x;
			weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			total_weight += 6 * w_x;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

				weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
			}
			else {

				weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRX) {

			total_weight += w_x;

			if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * (quantity[idx + 1] - obj.bdiff(idx).x * h.x);
			else						 weighted_sum += w_x * (quantity[idx - 1] + obj.bdiff(idx).x * h.x);
		}

		//y direction
		if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

			total_weight += 2 * w_y;
			weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			total_weight += 6 * w_y;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

				weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
			}
			else {

				weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRY) {

			total_weight += w_y;

			if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * (quantity[idx + n.x] - obj.bdiff(idx).y * h.y);
			else						 weighted_sum += w_y * (quantity[idx - n.x] + obj.bdiff(idx).y * h.y);
		}

		//z direction
		if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

			total_weight += 2 * w_z;
			weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

			total_weight += 6 * w_z;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

				weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
			}
			else {

				weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRZ) {

			total_weight += w_z;

			if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * (quantity[idx + n.x*n.y] - obj.bdiff(idx).z * h.z);
			else						 weighted_sum += w_z * (quantity[idx - n.x*n.y] + obj.bdiff(idx).z * h.z);
		}

		quantity[idx] = quantity[idx] * (1.0 - damping) + damping * ((weighted_sum - h_max * h_max * obj.Poisson_RHS(idx)) / total_weight);
	}
}

//-------------------- DEVICE BLACK

template <typename VType>
template <typename Class_Poisson_NNeu>
__device__ void cuVEC_VC<VType>::IteratePoisson_NNeu_SOR_black(Class_Poisson_NNeu& obj, cuBReal damping, cuBReal& max_error, cuBReal& max_val)
{
	//this method must be called with half-size : arr_size / 2 = n.dim() / 2, i.e. <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
	//double idx : idx values will now take on even values
	int idx = (blockIdx.x * blockDim.x + threadIdx.x) * 2;

	//ijk coordinates corresponding to idx
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	//for red-black passes must adjust the even idx values so we keep to a 3D checkerboard pattern
	//"red" : start at 0, "black" : start at 1
	//The following rules can be easily checked (remember i is the column number, j is the row number, k is the plane number)

	//For red squares:
	//If n.x is even : nudge idx by +1 on a) odd rows and even planes, and b) even rows and odd planes
	//If n.x is odd and n.y is even, nudge idx by +1 on odd planes only
	//If n.x is odd and n.y is odd, no nudge is needed : idx is automatically on the right 3D checkerboard pattern.

	//For black squares:
	//If n.x is even : nudge idx by +1 on a) even rows and even planes, and b) odd rows and odd planes
	//If n.x is odd and n.y is even, nudge idx by +1 on even planes only
	//If n.x is odd and n.y is odd, nudge by +1 everywhere : idx is automatically on the red 3D checkerboard pattern. (so nudge by 1 to find black checkerboard pattern).

	if (n.x % 2 == 1) {
		if (n.y % 2 == 0) {

			//nx is odd and ny is even : for black squares nudge on even planes only
			idx += (int)(n.z % 2 == 0);
		}
		else {

			//nx is odd and ny is odd, nudge everything by 1 for black squares
			idx++;
		}
	}
	else {

		//nx is even : for black squares nudge on a) even rows and even planes, and b) odd rows and odd planes

		//red_nudge = true for odd rows and even planes or for even rows and odd planes - have to keep index on the checkerboard pattern
		bool black_nudge = (((ijk.j % 2) == 0 && (ijk.k % 2) == 0) || (((ijk.j % 2) == 1 && (ijk.k % 2) == 1)));

		idx += (int)black_nudge;
	}

	//calculate new value only in non-empty cells; also skip if indicated as a composite media boundary condition cell
	bool calculate_idx = idx < n.dim() && (!(ngbrFlags[idx] & NF_CMBND) && (ngbrFlags[idx] & NF_NOTEMPTY));

	VType old_value = VType();

	if (calculate_idx) {

		//get maximum cell side
		cuBReal h_max = cu_maximum(h.x, h.y, h.z);

		//get weights
		cuBReal w_x = (h_max / h.x) * (h_max / h.x);
		cuBReal w_y = (h_max / h.y) * (h_max / h.y);
		cuBReal w_z = (h_max / h.z) * (h_max / h.z);

		VType weighted_sum = VType();
		cuBReal total_weight = 0;

		//x direction
		if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

			total_weight += 2 * w_x;
			weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

			total_weight += 6 * w_x;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

				weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
			}
			else {

				weighted_sum += w_x * (4 * get_dirichlet_value(NF2_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRX) {

			total_weight += w_x;

			if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * (quantity[idx + 1] - obj.bdiff(idx).x * h.x);
			else						 weighted_sum += w_x * (quantity[idx - 1] + obj.bdiff(idx).x * h.x);
		}

		//y direction
		if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

			total_weight += 2 * w_y;
			weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

			total_weight += 6 * w_y;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

				weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
			}
			else {

				weighted_sum += w_y * (4 * get_dirichlet_value(NF2_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRY) {

			total_weight += w_y;

			if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * (quantity[idx + n.x] - obj.bdiff(idx).y * h.y);
			else						 weighted_sum += w_y * (quantity[idx - n.x] + obj.bdiff(idx).y * h.y);
		}

		//z direction
		if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

			total_weight += 2 * w_z;
			weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
		}
		else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

			total_weight += 6 * w_z;

			if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

				weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
			}
			else {
				
				weighted_sum += w_z * (4 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
			}
		}
		else if (ngbrFlags[idx] & NF_NGBRZ) {

			total_weight += w_z;

			if (ngbrFlags[idx] & NF_NPZ) weighted_sum += w_z * (quantity[idx + n.x*n.y] - obj.bdiff(idx).z * h.z);
			else						 weighted_sum += w_z * (quantity[idx - n.x*n.y] + obj.bdiff(idx).z * h.z);
		}

		old_value = quantity[idx];
		quantity[idx] = quantity[idx] * (1.0 - damping) + damping * ((weighted_sum - h_max * h_max * obj.Poisson_RHS(idx)) / total_weight);
	}

	reduction_delta(idx, n.dim(), quantity, old_value, max_error, calculate_idx);
	reduction_delta(idx, n.dim(), quantity, VType(), max_val, calculate_idx);
}

//-------------------- LAUNCHER

template <typename VType>
template <typename Class_Poisson_NNeu>
__host__ void cuVEC_VC<VType>::IteratePoisson_NNeu_SOR(size_t arr_size, Class_Poisson_NNeu& obj, cuBReal& damping, cuBReal& max_error, cuBReal& max_val)
{
	IteratePoisson_NNeu_SOR_red_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, obj, damping);

	IteratePoisson_NNeu_SOR_black_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, obj, damping, max_error, max_val);
}