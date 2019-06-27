#pragma once

#include "cuVEC_VC.h"
#include "cuFuncs_Math.h"
#include "launchers.h"
#include "cuVEC_aux.cuh"

template <typename VType = cuReal>
__global__ void normalize_error_kernel(cuReal& max_error, cuReal& max_val)
{
	if (threadIdx.x == 0) {

		if (max_val) max_error /= max_val;
	}
}

//---------------------------------------------------- LAPLACE

//-------------------- GLOBAL RED

template <typename VType>
__global__ void IterateLaplace_SOR_red_kernel(cuVEC_VC<VType>& vec, cuReal& damping)
{
	vec.IterateLaplace_SOR_red(damping);
}

//-------------------- GLOBAL BLACK

template <typename VType>
__global__ void IterateLaplace_SOR_black_kernel(cuVEC_VC<VType>& vec, cuReal& damping, cuReal& max_error, cuReal& max_val)
{
	vec.IterateLaplace_SOR_black(damping, max_error, max_val);
}

//-------------------- DEVICE RED

template <typename VType>
__device__ void cuVEC_VC<VType>::IterateLaplace_SOR_red(cuReal damping)
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

	//VType old_value = VType();

	if (calculate_idx) {

		//get maximum cell side
		cuReal h_max = cu_maximum(h.x, h.y, h.z);

		//get weights
		cuReal w_x = (h_max / h.x) * (h_max / h.x);
		cuReal w_y = (h_max / h.y) * (h_max / h.y);
		cuReal w_z = (h_max / h.z) * (h_max / h.z);

		VType weighted_sum = VType();
		cuReal total_weight = 0;

		//x direction
		if (ngbrFlags[idx] & NF_BOTHX) {

			total_weight += 2 * w_x;
			weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETX) {

			total_weight += 6 * w_x;

			if (ngbrFlags[idx] & NF_DIRICHLETPX) weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
			else								 weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
		}
		else if (ngbrFlags[idx] & NF_NGBRX) {

			total_weight += w_x;

			if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * quantity[idx + 1];
			else						 weighted_sum += w_x * quantity[idx - 1];
		}

		//y direction
		if (ngbrFlags[idx] & NF_BOTHY) {

			total_weight += 2 * w_y;
			weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETY) {

			total_weight += 6 * w_y;

			if (ngbrFlags[idx] & NF_DIRICHLETPY) weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
			else								 weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
		}
		else if (ngbrFlags[idx] & NF_NGBRY) {

			total_weight += w_y;

			if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * quantity[idx + n.x];
			else						 weighted_sum += w_y * quantity[idx - n.x];
		}

		//z direction
		if (ngbrFlags[idx] & NF_BOTHZ) {

			total_weight += 2 * w_z;
			weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

			total_weight += 6 * w_z;

			if (ngbrFlags[idx] & NF_DIRICHLETPZ) weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
			else								 weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
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
__device__ void cuVEC_VC<VType>::IterateLaplace_SOR_black(cuReal damping, cuReal& max_error, cuReal& max_val)
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
		cuReal h_max = cu_maximum(h.x, h.y, h.z);

		//get weights
		cuReal w_x = (h_max / h.x) * (h_max / h.x);
		cuReal w_y = (h_max / h.y) * (h_max / h.y);
		cuReal w_z = (h_max / h.z) * (h_max / h.z);

		VType weighted_sum = VType();
		cuReal total_weight = 0;

		//x direction
		if (ngbrFlags[idx] & NF_BOTHX) {

			total_weight += 2 * w_x;
			weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETX) {

			total_weight += 6 * w_x;

			if (ngbrFlags[idx] & NF_DIRICHLETPX) weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
			else								 weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
		}
		else if (ngbrFlags[idx] & NF_NGBRX) {

			total_weight += w_x;

			if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * quantity[idx + 1];
			else						 weighted_sum += w_x * quantity[idx - 1];
		}

		//y direction
		if (ngbrFlags[idx] & NF_BOTHY) {

			total_weight += 2 * w_y;
			weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETY) {

			total_weight += 6 * w_y;

			if (ngbrFlags[idx] & NF_DIRICHLETPY) weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
			else								 weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
		}
		else if (ngbrFlags[idx] & NF_NGBRY) {

			total_weight += w_y;

			if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * quantity[idx + n.x];
			else						 weighted_sum += w_y * quantity[idx - n.x];
		}

		//z direction
		if (ngbrFlags[idx] & NF_BOTHZ) {

			total_weight += 2 * w_z;
			weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

			total_weight += 6 * w_z;

			if (ngbrFlags[idx] & NF_DIRICHLETPZ) weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
			else								 weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
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
__host__ void cuVEC_VC<VType>::IterateLaplace_SOR(size_t arr_size, cuReal& damping, cuReal& max_error, cuReal& max_val)
{
	IterateLaplace_SOR_red_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, damping);

	IterateLaplace_SOR_black_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, damping, max_error, max_val);
}

//---------------------------------------------------- POISSON

//-------------------- GLOBAL RED

template <typename VType, typename Class_Poisson_RHS>
__global__ void IteratePoisson_SOR_red_kernel(cuVEC_VC<VType>& vec, Class_Poisson_RHS& obj, cuReal& damping)
{
	vec.IteratePoisson_SOR_red(obj, damping);
}

//-------------------- GLOBAL BLACK

template <typename VType, typename Class_Poisson_RHS>
__global__ void IteratePoisson_SOR_black_kernel(cuVEC_VC<VType>& vec, Class_Poisson_RHS& obj, cuReal& damping, cuReal& max_error, cuReal& max_val)
{
	vec.IteratePoisson_SOR_black(obj, damping, max_error, max_val);
}

//-------------------- DEVICE RED

template <typename VType>
template <typename Class_Poisson_RHS>
__device__ void cuVEC_VC<VType>::IteratePoisson_SOR_red(Class_Poisson_RHS& obj, cuReal damping)
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

	//VType old_value = VType();

	if (calculate_idx) {

		//get maximum cell side
		cuReal h_max = cu_maximum(h.x, h.y, h.z);

		//get weights
		cuReal w_x = (h_max / h.x) * (h_max / h.x);
		cuReal w_y = (h_max / h.y) * (h_max / h.y);
		cuReal w_z = (h_max / h.z) * (h_max / h.z);

		VType weighted_sum = VType();
		cuReal total_weight = 0;

		//x direction
		if (ngbrFlags[idx] & NF_BOTHX) {

			total_weight += 2 * w_x;
			weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETX) {

			total_weight += 6 * w_x;

			if (ngbrFlags[idx] & NF_DIRICHLETPX) weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
			else								 weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
		}
		else if (ngbrFlags[idx] & NF_NGBRX) {

			total_weight += w_x;

			if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * quantity[idx + 1];
			else						 weighted_sum += w_x * quantity[idx - 1];
		}

		//y direction
		if (ngbrFlags[idx] & NF_BOTHY) {

			total_weight += 2 * w_y;
			weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETY) {

			total_weight += 6 * w_y;

			if (ngbrFlags[idx] & NF_DIRICHLETPY) weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
			else								 weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
		}
		else if (ngbrFlags[idx] & NF_NGBRY) {

			total_weight += w_y;

			if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * quantity[idx + n.x];
			else						 weighted_sum += w_y * quantity[idx - n.x];
		}

		//z direction
		if (ngbrFlags[idx] & NF_BOTHZ) {

			total_weight += 2 * w_z;
			weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

			total_weight += 6 * w_z;

			if (ngbrFlags[idx] & NF_DIRICHLETPZ) weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
			else								 weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
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
__device__ void cuVEC_VC<VType>::IteratePoisson_SOR_black(Class_Poisson_RHS& obj, cuReal damping, cuReal& max_error, cuReal& max_val)
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
		cuReal h_max = cu_maximum(h.x, h.y, h.z);

		//get weights
		cuReal w_x = (h_max / h.x) * (h_max / h.x);
		cuReal w_y = (h_max / h.y) * (h_max / h.y);
		cuReal w_z = (h_max / h.z) * (h_max / h.z);

		VType weighted_sum = VType();
		cuReal total_weight = 0;

		//x direction
		if (ngbrFlags[idx] & NF_BOTHX) {

			total_weight += 2 * w_x;
			weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETX) {

			total_weight += 6 * w_x;

			if (ngbrFlags[idx] & NF_DIRICHLETPX) weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
			else								 weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
		}
		else if (ngbrFlags[idx] & NF_NGBRX) {

			total_weight += w_x;

			if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * quantity[idx + 1];
			else						 weighted_sum += w_x * quantity[idx - 1];
		}

		//y direction
		if (ngbrFlags[idx] & NF_BOTHY) {

			total_weight += 2 * w_y;
			weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETY) {

			total_weight += 6 * w_y;

			if (ngbrFlags[idx] & NF_DIRICHLETPY) weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
			else								 weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
		}
		else if (ngbrFlags[idx] & NF_NGBRY) {

			total_weight += w_y;

			if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * quantity[idx + n.x];
			else						 weighted_sum += w_y * quantity[idx - n.x];
		}

		//z direction
		if (ngbrFlags[idx] & NF_BOTHZ) {

			total_weight += 2 * w_z;
			weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

			total_weight += 6 * w_z;

			if (ngbrFlags[idx] & NF_DIRICHLETPZ) weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
			else								 weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
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
__host__ void cuVEC_VC<VType>::IteratePoisson_SOR(size_t arr_size, Class_Poisson_RHS& obj, cuReal& damping, cuReal& max_error, cuReal& max_val)
{
	IteratePoisson_SOR_red_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, obj, damping);

	IteratePoisson_SOR_black_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, obj, damping, max_error, max_val);
}

//---------------------------------------------------- POISSON with adaptive SOR

//-------------------- GLOBAL

template <typename VType>
__global__ void adjust_aSOR_damping_kernel(cuVEC_VC<VType>& vec, bool start_iters, cuReal err_limit, cuReal& error, cuReal& max_val)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if(idx == 0) vec.adjust_aSOR_damping(start_iters, err_limit, error, max_val);
}

//-------------------- GLOBAL RED

template <typename VType, typename Class_Poisson_RHS>
__global__ void IteratePoisson_aSOR_red_kernel(cuVEC_VC<VType>& vec, Class_Poisson_RHS& obj)
{
	vec.IteratePoisson_SOR_red(obj, vec.aSOR_damping);
}

//-------------------- GLOBAL BLACK

template <typename VType, typename Class_Poisson_RHS>
__global__ void IteratePoisson_aSOR_black_kernel(cuVEC_VC<VType>& vec, Class_Poisson_RHS& obj, cuReal& max_error, cuReal& max_val)
{
	vec.IteratePoisson_SOR_black(obj, vec.aSOR_damping, max_error, max_val);
}

//-------------------- DEVICE

template <typename VType>
__device__ void cuVEC_VC<VType>::adjust_aSOR_damping(bool start_iters, cuReal err_limit, cuReal& error, cuReal& max_val)
{
#define ASOR_SPIKEVAL	0.4F
#define ASOR_EXPONENT	2.1F
#define ASOR_BIAS	0.02F
#define ASOR_NUDGE	1.018F
#define	ASOR_MINDAMPING	0.2F
#define	ASOR_MAXDAMPING	2.0F

	cuReal norm_error = (max_val > 0 ? error / max_val : error);

	//prepare start of a sequence of iterations but don't adjust damping at the start
	if (start_iters) {

		aSOR_lasterror = norm_error;
		aSOR_lastgrad = 0.0;
		return;
	}

	//adjust damping
	cuReal grad_lnerror = (log(norm_error) - log(aSOR_lasterror)) / aSOR_damping;

	//apply full adjustment mechanism only if error is above threshold : below this cannot apply the normal mechanism due to "numerical noise"
	if (norm_error > err_limit) {

		//positive gradient - should decrease damping
		if (grad_lnerror >= 0) {

			//avoid spikes : do nothing if simple spike detected
			if (aSOR_lastgrad <= 0 && grad_lnerror > ASOR_SPIKEVAL) {

				//save parameters from this iteration
				aSOR_lasterror = norm_error;
				aSOR_lastgrad = grad_lnerror;
				return;
			}

			//decrease damping using formula : larger g results in bigger decrease
			aSOR_damping *= exp(-grad_lnerror * ASOR_EXPONENT - ASOR_BIAS);
		}
		//negative g - might be able to do better by increasing damping
		else {

			aSOR_damping *= ASOR_NUDGE;
		}
	}
	else {

		//error is below threshold, but don't want to be stuck with a low damping value : give it a gentle increase
		aSOR_damping *= ASOR_NUDGE;
	}

	//make sure damping is within bounds
	if (aSOR_damping < ASOR_MINDAMPING) aSOR_damping = ASOR_MINDAMPING;
	if (aSOR_damping > ASOR_MAXDAMPING) aSOR_damping = ASOR_MAXDAMPING;

	//save parameters from this iteration
	aSOR_lasterror = norm_error;
	aSOR_lastgrad = grad_lnerror;
}

//-------------------- LAUNCHER

template <typename VType>
template <typename Class_Poisson_RHS>
__host__ void cuVEC_VC<VType>::IteratePoisson_aSOR(size_t arr_size, Class_Poisson_RHS& obj, bool start_iters, cuReal err_limit, cuReal& max_error, cuReal& max_val)
{
	IteratePoisson_aSOR_red_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, obj);

	IteratePoisson_aSOR_black_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, obj, max_error, max_val);

	adjust_aSOR_damping_kernel <<<1, CUDATHREADS >>> (*this, start_iters, err_limit, max_error, max_val);
}

//---------------------------------------------------- POISSON with non-homogeneous Neumann boundaries

//-------------------- GLOBAL RED

template <typename VType, typename Class_Poisson_NNeu>
__global__ void IteratePoisson_NNeu_SOR_red_kernel(cuVEC_VC<VType>& vec, Class_Poisson_NNeu& obj, cuReal& damping)
{
	vec.IteratePoisson_NNeu_SOR_red(obj, damping);
}

//-------------------- GLOBAL BLACK

template <typename VType, typename Class_Poisson_NNeu>
__global__ void IteratePoisson_NNeu_SOR_black_kernel(cuVEC_VC<VType>& vec, Class_Poisson_NNeu& obj, cuReal& damping, cuReal& max_error, cuReal& max_val)
{
	vec.IteratePoisson_NNeu_SOR_black(obj, damping, max_error, max_val);
}

//-------------------- DEVICE RED

template <typename VType>
template <typename Class_Poisson_NNeu>
__device__ void cuVEC_VC<VType>::IteratePoisson_NNeu_SOR_red(Class_Poisson_NNeu& obj, cuReal damping)
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

	//VType old_value = VType();

	if (calculate_idx) {

		//get maximum cell side
		cuReal h_max = cu_maximum(h.x, h.y, h.z);

		//get weights
		cuReal w_x = (h_max / h.x) * (h_max / h.x);
		cuReal w_y = (h_max / h.y) * (h_max / h.y);
		cuReal w_z = (h_max / h.z) * (h_max / h.z);

		VType weighted_sum = VType();
		cuReal total_weight = 0;

		//x direction
		if (ngbrFlags[idx] & NF_BOTHX) {

			total_weight += 2 * w_x;
			weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETX) {

			total_weight += 6 * w_x;

			if (ngbrFlags[idx] & NF_DIRICHLETPX) weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
			else								 weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
		}
		else if (ngbrFlags[idx] & NF_NGBRX) {

			total_weight += w_x;

			if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * (quantity[idx + 1] - obj.bdiff(idx).x * h.x);
			else						 weighted_sum += w_x * (quantity[idx - 1] + obj.bdiff(idx).x * h.x);
		}

		//y direction
		if (ngbrFlags[idx] & NF_BOTHY) {

			total_weight += 2 * w_y;
			weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETY) {

			total_weight += 6 * w_y;

			if (ngbrFlags[idx] & NF_DIRICHLETPY) weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
			else								 weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
		}
		else if (ngbrFlags[idx] & NF_NGBRY) {

			total_weight += w_y;

			if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * (quantity[idx + n.x] - obj.bdiff(idx).y * h.y);
			else						 weighted_sum += w_y * (quantity[idx - n.x] + obj.bdiff(idx).y * h.y);
		}

		//z direction
		if (ngbrFlags[idx] & NF_BOTHZ) {

			total_weight += 2 * w_z;
			weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

			total_weight += 6 * w_z;

			if (ngbrFlags[idx] & NF_DIRICHLETPZ) weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
			else								 weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
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
__device__ void cuVEC_VC<VType>::IteratePoisson_NNeu_SOR_black(Class_Poisson_NNeu& obj, cuReal damping, cuReal& max_error, cuReal& max_val)
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
		cuReal h_max = cu_maximum(h.x, h.y, h.z);

		//get weights
		cuReal w_x = (h_max / h.x) * (h_max / h.x);
		cuReal w_y = (h_max / h.y) * (h_max / h.y);
		cuReal w_z = (h_max / h.z) * (h_max / h.z);

		VType weighted_sum = VType();
		cuReal total_weight = 0;

		//x direction
		if (ngbrFlags[idx] & NF_BOTHX) {

			total_weight += 2 * w_x;
			weighted_sum += w_x * (quantity[idx - 1] + quantity[idx + 1]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETX) {

			total_weight += 6 * w_x;

			if (ngbrFlags[idx] & NF_DIRICHLETPX) weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETPX, idx) + 2 * quantity[idx + 1]);
			else								 weighted_sum += w_x * (4 * get_dirichlet_value(NF_DIRICHLETNX, idx) + 2 * quantity[idx - 1]);
		}
		else if (ngbrFlags[idx] & NF_NGBRX) {

			total_weight += w_x;

			if (ngbrFlags[idx] & NF_NPX) weighted_sum += w_x * (quantity[idx + 1] - obj.bdiff(idx).x * h.x);
			else						 weighted_sum += w_x * (quantity[idx - 1] + obj.bdiff(idx).x * h.x);
		}

		//y direction
		if (ngbrFlags[idx] & NF_BOTHY) {

			total_weight += 2 * w_y;
			weighted_sum += w_y * (quantity[idx - n.x] + quantity[idx + n.x]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETY) {

			total_weight += 6 * w_y;

			if (ngbrFlags[idx] & NF_DIRICHLETPY) weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETPY, idx) + 2 * quantity[idx + n.x]);
			else								 weighted_sum += w_y * (4 * get_dirichlet_value(NF_DIRICHLETNY, idx) + 2 * quantity[idx - n.x]);
		}
		else if (ngbrFlags[idx] & NF_NGBRY) {

			total_weight += w_y;

			if (ngbrFlags[idx] & NF_NPY) weighted_sum += w_y * (quantity[idx + n.x] - obj.bdiff(idx).y * h.y);
			else						 weighted_sum += w_y * (quantity[idx - n.x] + obj.bdiff(idx).y * h.y);
		}

		//z direction
		if (ngbrFlags[idx] & NF_BOTHZ) {

			total_weight += 2 * w_z;
			weighted_sum += w_z * (quantity[idx - n.x*n.y] + quantity[idx + n.x*n.y]);
		}
		else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

			total_weight += 6 * w_z;

			if (ngbrFlags[idx] & NF_DIRICHLETPZ) weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETPZ, idx) + 2 * quantity[idx + n.x*n.y]);
			else								 weighted_sum += w_z * (4 * get_dirichlet_value(NF_DIRICHLETNZ, idx) + 2 * quantity[idx - n.x*n.y]);
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
__host__ void cuVEC_VC<VType>::IteratePoisson_NNeu_SOR(size_t arr_size, Class_Poisson_NNeu& obj, cuReal& damping, cuReal& max_error, cuReal& max_val)
{
	IteratePoisson_NNeu_SOR_red_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, obj, damping);

	IteratePoisson_NNeu_SOR_black_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, obj, damping, max_error, max_val);
}

//---------------------------------------------------- POISSON with non-homogeneous Neumann boundaries and adaptive SOR algorithm

//-------------------- GLOBAL RED

template <typename VType, typename Class_Poisson_NNeu>
__global__ void IteratePoisson_NNeu_aSOR_red_kernel(cuVEC_VC<VType>& vec, Class_Poisson_NNeu& obj)
{
	vec.IteratePoisson_NNeu_SOR_red(obj, vec.aSOR_damping);
}

//-------------------- GLOBAL BLACK

template <typename VType, typename Class_Poisson_NNeu>
__global__ void IteratePoisson_NNeu_aSOR_black_kernel(cuVEC_VC<VType>& vec, Class_Poisson_NNeu& obj, cuReal& max_error, cuReal& max_val)
{
	vec.IteratePoisson_NNeu_SOR_black(obj, vec.aSOR_damping, max_error, max_val);
}

//-------------------- LAUNCHER

template <typename VType>
template <typename Class_Poisson_NNeu>
__host__ void cuVEC_VC<VType>::IteratePoisson_NNeu_aSOR(size_t arr_size, Class_Poisson_NNeu& obj, bool start_iters, cuReal err_limit, cuReal& max_error, cuReal& max_val)
{
	IteratePoisson_NNeu_aSOR_red_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, obj);

	IteratePoisson_NNeu_aSOR_black_kernel <<< (arr_size / 2 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, obj, max_error, max_val);

	adjust_aSOR_damping_kernel <<<1, CUDATHREADS >>> (*this, start_iters, err_limit, max_error, max_val);
}