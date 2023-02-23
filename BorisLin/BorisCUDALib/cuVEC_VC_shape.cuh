#pragma once

#include "cuVEC_VC.h"

//--------------------------------------------MULTIPLE ENTRIES SETTERS - SHAPE CHANGERS

//------------------------------------------------------------------- SETBOX (for cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_setbox_kernel(const cuSZ3& n, int*& ngbrFlags, VType*& quantity, cuBox box, VType value)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));
	
	if (idx < n.dim() && box.Contains(ijk)) {

		quantity[idx] = value;
		ngbrFlags[idx] |= NF_NOTEMPTY;
	}
}


template void cuVEC_VC<float>::setbox(cuBox box, float value);
template void cuVEC_VC<double>::setbox(cuBox box, double value);

template void cuVEC_VC<cuFLT3>::setbox(cuBox box, cuFLT3 value);
template void cuVEC_VC<cuDBL3>::setbox(cuBox box, cuDBL3 value);

//set value in box (i.e. in cells entirely included in box) - all cells become non-empty cells irrespective of value set
template <typename VType>
__host__ void cuVEC_VC<VType>::setbox(cuBox box, VType value)
{
	cuvec_vc_setbox_kernel <<< (get_gpu_value(cuVEC<VType>::n).dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, ngbrFlags, cuVEC<VType>::quantity, box, value);

	set_ngbrFlags();
}

//------------------------------------------------------------------- DELRECT

template <typename VType>
__global__ void delrect_kernel(const cuSZ3& n, int*& ngbrFlags, VType*& quantity, cuBox box)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	if (idx < n.dim() && box.Contains(ijk)) {

		quantity[idx] = VType();
		ngbrFlags[idx] &= ~NF_NOTEMPTY;
	}
}

template void cuVEC_VC<float>::delrect(cuRect rectangle);
template void cuVEC_VC<double>::delrect(cuRect rectangle);

template void cuVEC_VC<cuFLT3>::delrect(cuRect rectangle);
template void cuVEC_VC<cuDBL3>::delrect(cuRect rectangle);

//delete rectangle, where the rectangle is relative to this VEC's rectangle, by setting empty cell values - all cells become empty cells
template <typename VType>
__host__ void cuVEC_VC<VType>::delrect(cuRect rectangle)
{
	cuBox box = cuVEC<VType>::box_from_rect_max_cpu(rectangle + get_gpu_value(cuVEC<VType>::rect).s);

	delrect_kernel <<< (get_gpu_value(cuVEC<VType>::n).dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, ngbrFlags, cuVEC<VType>::quantity, box);

	set_ngbrFlags();
}

//--------------------------------------------MULTIPLE ENTRIES SETTERS - SHAPE CHANGERS

template <typename VType>
__global__ void apply_bitmap_mask_kernel(const cuSZ3& n, int*& ngbrFlags, VType*& quantity, unsigned char* bitmap, int zDepth)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	if (idx < n.dim()) {

		//In the file the 0, 0 point is in the left-top corner, but in the mesh the 0, 0 point is in the left-bottom corner
		int indexBitMap = ijk.i + (n.y - 1 - ijk.j)*n.x;

		float B = (float)bitmap[indexBitMap * 4] / 255;
		//float G = (float)bitmap[indexBitMap * 4 + 1] / 255;
		//float R = (float)bitmap[indexBitMap * 4 + 2] / 255;

		int depthCut = 0;

		if (zDepth) {

			//divide grayscale in zDepth equal intervals of length intervalLength
			float intervalLength = (float)1.0 / abs(zDepth);

			//find out in which interval the grayscale pixel lies (assume B=G=R)
			depthCut = (int)floor(B / intervalLength);
		}
		else {

			if (cuIsNZ(B)) depthCut = n.z;
		}

		if (zDepth > 0) {

			//void to given depth
			if (ijk.k > n.z - 1 - depthCut) {

				quantity[idx] = VType();
				ngbrFlags[idx] &= ~NF_NOTEMPTY;
			}
		}

		if (zDepth <= 0) {

			//void to given height
			if (ijk.k < depthCut) {

				quantity[idx] = VType();
				ngbrFlags[idx] &= ~NF_NOTEMPTY;
			}
		}
	}
}

template bool cuVEC_VC<float>::apply_bitmap_mask(std::vector<unsigned char>& bitmap, int zDepth);
template bool cuVEC_VC<double>::apply_bitmap_mask(std::vector<unsigned char>& bitmap, int zDepth);

template bool cuVEC_VC<cuFLT3>::apply_bitmap_mask(std::vector<unsigned char>& bitmap, int zDepth);
template bool cuVEC_VC<cuDBL3>::apply_bitmap_mask(std::vector<unsigned char>& bitmap, int zDepth);

template <typename VType>
bool cuVEC_VC<VType>::apply_bitmap_mask(std::vector<unsigned char>& bitmap, int zDepth)
{
	cuSZ3 n_ = get_gpu_value(cuVEC<VType>::n);

	//bitmap must have the right size (i.e. have n.x * n.y pixels, remembering each pixel has 4 bytes as B-G-R-A)
	if (bitmap.size() != n_.x*n_.y * 4) return false;

	unsigned char* bitmap_gpu = nullptr;
	if (gpu_alloc(bitmap_gpu, bitmap.size()) != cudaSuccess) return false;

	cpu_to_gpu(bitmap_gpu, bitmap.data(), bitmap.size());

	apply_bitmap_mask_kernel <<< (n_.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, ngbrFlags, cuVEC<VType>::quantity, bitmap_gpu, zDepth);

	gpu_free(bitmap_gpu);

	set_ngbrFlags();

	return true;
}