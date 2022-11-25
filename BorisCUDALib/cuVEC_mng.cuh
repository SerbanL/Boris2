#pragma once

#include "Funcs_Vectors.h"		//need malloc_vector from this BorisLib header

#include "cuVEC.h"
#include "launchers.h"
#include "cuVEC_mng.h"
#include "cuVEC_MeshTransfer.h"

//------------------------------------------------------------------- MAPMESH_NEWDIMS

template <typename VType>
__global__ void mapmesh_newdims_kernel(cuSZ3 new_n, cuSZ3& old_n, VType*& new_quantity, VType* old_quantity)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//now also transfer mesh values to new dimensions
	cuReal3 sourceIdx = (cuReal3)old_n / new_n;

	if (idx < new_n.dim()) {

		int _x = (int)floor((idx % new_n.x) * sourceIdx.x);
		int _y = (int)floor(((idx / new_n.x) % new_n.y) * sourceIdx.y);
		int _z = (int)floor((idx / (new_n.x*new_n.y)) * sourceIdx.z);

		new_quantity[idx] = old_quantity[_x + _y * int(old_n.x) + _z * (int(old_n.x*old_n.y))];
	}
}

template bool cuVEC<float>::mapmesh_newdims(const cuSZ3& new_n);
template bool cuVEC<double>::mapmesh_newdims(const cuSZ3& new_n);

template bool cuVEC<cuFLT3>::mapmesh_newdims(const cuSZ3& new_n);
template bool cuVEC<cuDBL3>::mapmesh_newdims(const cuSZ3& new_n);

template bool cuVEC<cuFLT4>::mapmesh_newdims(const cuSZ3& new_n);
template bool cuVEC<cuDBL4>::mapmesh_newdims(const cuSZ3& new_n);

template bool cuVEC<cuReIm3>::mapmesh_newdims(const cuSZ3& new_n);

template <typename VType>
__host__ bool cuVEC<VType>::mapmesh_newdims(const cuSZ3& new_n)
{
	cudaError_t error = cudaSuccess;

	//We need to allocate memory for new size new_n.dim() and map values from old size to new size.
	//1. The fastest way to do this is to allocate new space new_n and do the mapping using a gpu kernel, all in gpu memory : we need an extra new_n.dim() elements in gpu memory
	//2. If method 1 fails we can try to transfer array to cpu memory, do the mapping using the cpu then transfer it back to gpu : much slower but only need an extra new_n.dim() - n.dim() elements in gpu memory.

	//method 1 : need extra new_n.dim() elements in gpu memory

	bool enough_space_for_method1 = (new_n.dim() <= cudaMemGetFree() / sizeof(VType));
	bool enough_space_for_method2 = (new_n.dim() <= get_gpu_value(n).dim() || new_n.dim() - get_gpu_value(n).dim() <= cudaMemGetFree() / sizeof(VType));

	//allocate space for new quantity (call it old as we'll swap)
	VType* old_quantity = nullptr;
	
	if (enough_space_for_method1) {

		error = gpu_alloc(old_quantity, new_n.dim());
	}

	if (!enough_space_for_method1 || error != cudaSuccess) {

		//memory allocation has failed so try method 2. instead ...
		gpu_free(old_quantity);

		//... unless of course we know there isn't enough space even for this
		if (!enough_space_for_method2) return false;

		//could work so allocate cpu memory 
		std::vector<VType> cpu_quantity_old, cpu_quantity_new;
		
		if (!malloc_vector(cpu_quantity_old, get_gpu_value(n).dim())) return false;
		if (!malloc_vector(cpu_quantity_new, new_n.dim())) return false;

		//transfer values from gpu to cpu
		error = gpu_to_cpu_managed(cpu_quantity_old.data(), quantity, get_gpu_value(n).dim());
		if (error != cudaSuccess) return false;

		//resize gpu memory to new size
		error = gpu_alloc_managed(quantity, new_n.dim());
		if (error != cudaSuccess) {

			//failed : go back to original size and values
			error = gpu_alloc_managed(quantity, get_gpu_value(n).dim());
			if (error != cudaSuccess) return false;
			
			error = cpu_to_gpu_managed(quantity, cpu_quantity_old.data(), get_gpu_value(n).dim());
			if (error != cudaSuccess) return false;

			return true;
		}

		//now map values from old to new size using the cpu
		cuReal3 sourceIdx = (cuReal3)get_gpu_value(n) / new_n;

		#pragma omp parallel for
		for (int idx = 0; idx < new_n.dim(); idx++) {

			int _x = (int)floor((idx % new_n.x) * sourceIdx.x);
			int _y = (int)floor(((idx / new_n.x) % new_n.y) * sourceIdx.y);
			int _z = (int)floor((idx / (new_n.x*new_n.y)) * sourceIdx.z);

			cpu_quantity_new[idx] = cpu_quantity_old[_x + _y * int(get_gpu_value(n).x) + _z * (int(get_gpu_value(n).x*get_gpu_value(n).y))];
		}

		//transfer values to new size
		error = cpu_to_gpu_managed(quantity, cpu_quantity_new.data(), new_n.dim());
		if (error != cudaSuccess) return false;

		//set new size
		set_n(new_n);

		return true;
	}
	else {

		//enough memory for method 1 so continue all on gpu

		//swap new and old : quantity will now point to new memory space (of size new_n.dim()), old_quantity to old space (of size n.dim())
		gpu_swap_managed(quantity, old_quantity);

		//map values from old quantity to new quantity size using a kernel
		mapmesh_newdims_kernel <<< (new_n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (new_n, n, quantity, old_quantity);

		//free old quantity
		gpu_free(old_quantity);

		//set new size
		set_n(new_n);

		return true;
	}
}

//------------------------------------------------------------------- EXTRACT CUVEC

//Copy from quantity_in to quantity_out (considered as 3D arrays), values from cells containing the center of cells of quantity_out.
//Note : kernel launched at size n.dim() >= n_coarse.dim()
template <typename VType>
__global__ void strided_copy_3d(VType*& quantity_in, cuSZ3& n, VType*& quantity_out, cuSZ3& n_coarse)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < n_coarse.dim()) {

		//ijk_coarse for quantity_out (with size n_coarse)
		cuINT3 ijk_coarse = cuINT3(idx % n_coarse.x, (idx / n_coarse.x) % n_coarse.y, idx / (n_coarse.x*n_coarse.y));

		//we don't need the real cellsizes since we know the rectangle of the two meshes must be the same - take it as a normalized unit rectangle.
		cuReal3 h_coarse, h;
		h = cuReal3(1) / n;
		h_coarse = cuReal3(1) / n_coarse;

		//cell center in normalized rectangle
		cuReal3 cell_center = (ijk_coarse + cuReal3(0.5)) & h_coarse;

		//index of cell in quantity_in containing cell_center
		cuINT3 ijk = cu_floor(cell_center / h);

		//copy value
		quantity_out[idx] = quantity_in[ijk.i + ijk.j * n.x + ijk.k * n.x*n.y];
	}
}

template void cuVEC<float>::extract_cuvec(size_t size, cuVEC<float>& cuvec);
template void cuVEC<double>::extract_cuvec(size_t size, cuVEC<double>& cuvec);

template void cuVEC<cuFLT3>::extract_cuvec(size_t size, cuVEC<cuFLT3>& cuvec);
template void cuVEC<cuDBL3>::extract_cuvec(size_t size, cuVEC<cuDBL3>& cuvec);

template void cuVEC<cuFLT4>::extract_cuvec(size_t size, cuVEC<cuFLT4>& cuvec);
template void cuVEC<cuDBL4>::extract_cuvec(size_t size, cuVEC<cuDBL4>& cuvec);

template void cuVEC<cuReIm3>::extract_cuvec(size_t size, cuVEC<cuReIm3>& cuvec);

template <typename VType>
__host__ void cuVEC<VType>::extract_cuvec(size_t size, cuVEC<VType>& cuvec)
{
	strided_copy_3d <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (quantity, n, cuvec.quantity, cuvec.n);
}

//------------------------------------------------------------------- CUARR LOAD / STORE

template <typename VType>
__global__ void load_cuarr_kernel(VType*& quantity, VType* input, cuSZ3& n)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < n.dim()) {

		quantity[idx] = input[idx];
	}
}

template void cuVEC<float>::load_cuarr(size_t size, cu_arr<float>& input);
template void cuVEC<double>::load_cuarr(size_t size, cu_arr<double>& input);

template void cuVEC<cuFLT3>::load_cuarr(size_t size, cu_arr<cuFLT3>& input);
template void cuVEC<cuDBL3>::load_cuarr(size_t size, cu_arr<cuDBL3>& input);

template void cuVEC<cuFLT4>::load_cuarr(size_t size, cu_arr<cuFLT4>& input);
template void cuVEC<cuDBL4>::load_cuarr(size_t size, cu_arr<cuDBL4>& input);

template void cuVEC<cuReIm3>::load_cuarr(size_t size, cu_arr<cuReIm3>& input);

//copy values from a cu_arr of same type -> sizes must match
template <typename VType>
__host__ void cuVEC<VType>::load_cuarr(size_t size, cu_arr<VType>& input)
{
	load_cuarr_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (quantity, (VType*)input, n);
}

template <typename VType>
__global__ void store_cuarr_kernel(VType*& quantity, VType* output, cuSZ3& n)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < n.dim()) {

		output[idx] = quantity[idx];
	}
}

template void cuVEC<float>::store_cuarr(size_t size, cu_arr<float>& output);
template void cuVEC<double>::store_cuarr(size_t size, cu_arr<double>& output);

template void cuVEC<cuFLT3>::store_cuarr(size_t size, cu_arr<cuFLT3>& output);
template void cuVEC<cuDBL3>::store_cuarr(size_t size, cu_arr<cuDBL3>& output);

template void cuVEC<cuFLT4>::store_cuarr(size_t size, cu_arr<cuFLT4>& output);
template void cuVEC<cuDBL4>::store_cuarr(size_t size, cu_arr<cuDBL4>& output);

template void cuVEC<cuReIm3>::store_cuarr(size_t size, cu_arr<cuReIm3>& output);

//copy values to a cu_arr of same type -> sizes must match
template <typename VType>
__host__ void cuVEC<VType>::store_cuarr(size_t size, cu_arr<VType>& output)
{
	store_cuarr_kernel << < (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (quantity, (VType*)output, n);
}