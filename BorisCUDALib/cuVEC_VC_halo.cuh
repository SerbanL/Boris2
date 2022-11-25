#pragma once

#include "cuVEC_VC.h"
#include "launchers.h"

//------------------------------------------------------------------- EXTRACT NGBRFLAGS

__global__ static void extract_ngbrFlags_kernel(cuBox cells_box, const cuSZ3& n, int*& ngbrFlags, int* linear_storage)
{
	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (box_idx < cells_box.size().dim()) {
		
		//i, j, k values inside cells_box
		int box_i = box_idx % (cells_box.size().x);
		int box_j = (box_idx / (cells_box.size().x)) % cells_box.size().y;
		int box_k = box_idx / (cells_box.size().x*cells_box.size().y);

		//now form index in ngbrFlags
		int ngbr_idx = (box_i + cells_box.s.i) + (box_j + cells_box.s.j) * n.x + (box_k + cells_box.s.k) * n.x*n.y;
		//linear_storage size matches cells_box.size().dim(), thus use box_idx
		linear_storage[box_idx] = ngbrFlags[ngbr_idx];
	}
}

template bool cuVEC_VC<float>::extract_ngbrFlags(cuBox cells_box, cu_arr<int>& linear_storage);
template bool cuVEC_VC<double>::extract_ngbrFlags(cuBox cells_box, cu_arr<int>& linear_storage);

template bool cuVEC_VC<cuFLT3>::extract_ngbrFlags(cuBox cells_box, cu_arr<int>& linear_storage);
template bool cuVEC_VC<cuDBL3>::extract_ngbrFlags(cuBox cells_box, cu_arr<int>& linear_storage);

//extract ngbrFlags values within the given cells_box (end point not included, so e.g. cuBox(n) is the entire mesh), and place them inside the provided linear_storage space
template <typename VType>
__host__ bool cuVEC_VC<VType>::extract_ngbrFlags(cuBox cells_box, cu_arr<int>& linear_storage)
{
	if (linear_storage.size() != cells_box.size().dim()) return false;

	extract_ngbrFlags_kernel <<< (cells_box.size().dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cells_box, cuVEC<VType>::n, ngbrFlags, linear_storage);
	return true;
}

//------------------------------------------------------------------- EXTRACT QUANTITY TO HALO

template <typename VType>
__global__ static void extract_halonx_kernel(const cuSZ3& n, VType*& quantity, VType*& halo, size_t& halo_size)
{
	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	int halo_depth = halo_size / (n.y * n.z);

	if (box_idx < halo_size) {

		//i, j, k values inside cells_box
		int box_i = box_idx % halo_depth;
		int box_j = (box_idx / halo_depth) % n.y;
		int box_k = box_idx / (halo_depth*n.y);

		//now form index in quantity
		int quant_idx = (box_i + n.x - halo_depth) + box_j* n.x + box_k * n.x*n.y;
		
		halo[box_idx] = quantity[quant_idx];
	}
}

template <typename VType>
__global__ static void extract_halopx_kernel(const cuSZ3& n, VType*& quantity, VType*& halo, size_t& halo_size)
{
	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	int halo_depth = halo_size / (n.y * n.z);

	if (box_idx < halo_size) {

		//i, j, k values inside cells_box
		int box_i = box_idx % halo_depth;
		int box_j = (box_idx / halo_depth) % n.y;
		int box_k = box_idx / (halo_depth*n.y);

		//now form index in quantity
		int quant_idx = box_i + box_j * n.x + box_k * n.x*n.y;

		halo[box_idx] = quantity[quant_idx];
	}
}

template <typename VType>
__global__ static void extract_halony_kernel(const cuSZ3& n, VType*& quantity, VType*& halo, size_t& halo_size)
{
	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	int halo_depth = halo_size / (n.x * n.z);

	if (box_idx < halo_size) {

		//i, j, k values inside cells_box
		int box_i = box_idx % n.x;
		int box_j = (box_idx / n.x) % halo_depth;
		int box_k = box_idx / (n.x*halo_depth);

		//now form index in quantity
		int quant_idx = box_i + (box_j + n.y - halo_depth) * n.x + box_k * n.x*halo_depth;

		halo[box_idx] = quantity[quant_idx];
	}
}

template <typename VType>
__global__ static void extract_halopy_kernel(const cuSZ3& n, VType*& quantity, VType*& halo, size_t& halo_size)
{
	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	int halo_depth = halo_size / (n.x * n.z);

	if (box_idx < halo_size) {

		//i, j, k values inside cells_box
		int box_i = box_idx % n.x;
		int box_j = (box_idx / n.x) % halo_depth;
		int box_k = box_idx / (n.x*halo_depth);

		//now form index in quantity
		int quant_idx = box_i + box_j * n.x + box_k * n.x*halo_depth;

		halo[box_idx] = quantity[quant_idx];
	}
}

template <typename VType>
__global__ static void extract_halonz_kernel(const cuSZ3& n, VType*& quantity, VType*& halo, size_t& halo_size)
{
	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	int halo_depth = halo_size / (n.x * n.y);

	if (box_idx < halo_size) {

		//i, j, k values inside cells_box
		int box_i = box_idx % n.x;
		int box_j = (box_idx / n.x) % n.y;
		int box_k = box_idx / (n.x*n.y);

		//now form index in quantity
		int quant_idx = box_i + box_j * n.x + (box_k + n.z - halo_depth) * n.x*n.y;

		halo[box_idx] = quantity[quant_idx];
	}
}

template <typename VType>
__global__ static void extract_halopz_kernel(const cuSZ3& n, VType*& quantity, VType*& halo, size_t& halo_size)
{
	int box_idx = blockIdx.x * blockDim.x + threadIdx.x;

	int halo_depth = halo_size / (n.x * n.y);

	if (box_idx < halo_size) {

		//i, j, k values inside cells_box
		int box_i = box_idx % n.x;
		int box_j = (box_idx / n.x) % n.y;
		int box_k = box_idx / (n.x*n.y);

		//now form index in quantity
		int quant_idx = box_i + box_j * n.x + box_k * n.x*n.y;

		halo[box_idx] = quantity[quant_idx];
	}
}

template void cuVEC_VC<float>::extract_halo(int halo_id, size_t size_transfer);
template void cuVEC_VC<double>::extract_halo(int halo_id, size_t size_transfer);

template void cuVEC_VC<cuFLT3>::extract_halo(int halo_id, size_t size_transfer);
template void cuVEC_VC<cuDBL3>::extract_halo(int halo_id, size_t size_transfer);

//extract values from quantity into halo arrays as specified by halo_id
//halo_id can be an individual flag, and identifies the required halos to extract values to, e.g. NF2_HALOPX, NF2_HALONX
//values from p side of quantity are extracted to n halo and values from n side are extracted to p halo
//e.g. if called with NF2_HALOPX, then values from n side are extracted to halo_px; if called with NF2_HALONX then values from p side are extracted to halo_nx
template <typename VType>
__host__ void cuVEC_VC<VType>::extract_halo(int halo_id, size_t size_transfer)
{
	switch (halo_id) {

	case NF2_HALONX:
		extract_halonx_kernel <<< (size_transfer + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, cuVEC<VType>::quantity, halo_nx, halo_x_size);
		break;
	case NF2_HALOPX:
		extract_halopx_kernel <<< (size_transfer + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, cuVEC<VType>::quantity, halo_px, halo_x_size);
		break;
	case NF2_HALONY:
		extract_halony_kernel <<< (size_transfer + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, cuVEC<VType>::quantity, halo_ny, halo_y_size);
		break;
	case NF2_HALOPY:
		extract_halopy_kernel <<< (size_transfer + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, cuVEC<VType>::quantity, halo_py, halo_y_size);
		break;
	case NF2_HALONZ:
		extract_halonz_kernel <<< (size_transfer + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, cuVEC<VType>::quantity, halo_nz, halo_z_size);
		break;
	case NF2_HALOPZ:
		extract_halopz_kernel <<< (size_transfer + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, cuVEC<VType>::quantity, halo_pz, halo_z_size);
		break;
	}
}