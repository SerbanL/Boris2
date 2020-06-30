#pragma once

#include "cuVEC_VC.h"
#include "launchers.h"


//------------------------------------------------------------------- EXTRACT CUVEC

//Copy from quantity_in to quantity_out (considered as a 3D arrays), values from cells containing the center of cells of quantity_out.
//Note : kernel launched at size n.dim() >= n_coarse.dim()
template <typename VType>
__global__ void strided_copy_3d(VType*& quantity_in, int*& ngbrFlags_in, cuSZ3& n, VType*& quantity_out, int*& ngbrFlags_out, cuSZ3& n_coarse)
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
		ngbrFlags_out[idx] = ngbrFlags_in[ijk.i + ijk.j * n.x + ijk.k * n.x*n.y];
	}
}

template void cuVEC_VC<char>::extract_cuvec(size_t size, cuVEC_VC<char>& cuvec);
template void cuVEC_VC<int>::extract_cuvec(size_t size, cuVEC_VC<int>& cuvec);
template void cuVEC_VC<unsigned>::extract_cuvec(size_t size, cuVEC_VC<unsigned>& cuvec);
template void cuVEC_VC<long>::extract_cuvec(size_t size, cuVEC_VC<long>& cuvec);
template void cuVEC_VC<size_t>::extract_cuvec(size_t size, cuVEC_VC<size_t>& cuvec);
template void cuVEC_VC<float>::extract_cuvec(size_t size, cuVEC_VC<float>& cuvec);
template void cuVEC_VC<double>::extract_cuvec(size_t size, cuVEC_VC<double>& cuvec);

template void cuVEC_VC<cuINT3>::extract_cuvec(size_t size, cuVEC_VC<cuINT3>& cuvec);
template void cuVEC_VC<cuSZ3>::extract_cuvec(size_t size, cuVEC_VC<cuSZ3>& cuvec);
template void cuVEC_VC<cuFLT3>::extract_cuvec(size_t size, cuVEC_VC<cuFLT3>& cuvec);
template void cuVEC_VC<cuDBL3>::extract_cuvec(size_t size, cuVEC_VC<cuDBL3>& cuvec);

template void cuVEC_VC<cuINT33>::extract_cuvec(size_t size, cuVEC_VC<cuINT33>& cuvec);
template void cuVEC_VC<cuFLT33>::extract_cuvec(size_t size, cuVEC_VC<cuFLT33>& cuvec);
template void cuVEC_VC<cuDBL33>::extract_cuvec(size_t size, cuVEC_VC<cuDBL33>& cuvec);

template <typename VType>
__host__ void cuVEC_VC<VType>::extract_cuvec(size_t size, cuVEC_VC<VType>& cuvec)
{
	strided_copy_3d <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::quantity, ngbrFlags, cuVEC<VType>::n, cuvec.quantity, cuvec.ngbrFlags, cuvec.n);
}