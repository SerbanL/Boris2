#pragma once

#include "cuVEC_VC.h"
#include "cuFuncs_Math.h"
#include "launchers.h"
#include "Reduction.cuh"

//------------------------------------------------------------------- AVERAGE NONEMPTY (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_average_nonempty_kernel(cuSZ3& n, int*& ngbrFlags, VType*& quantity, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//need the idx < n.dim() check before ngbrFlags[idx] to avoid bad memory access
	reduction_avg(idx, n.dim(), quantity, average_value, points_count, idx < n.dim() && ngbrFlags[idx] & NF_NOTEMPTY);
}

template <typename VType>
__global__ void cuvec_vc_average_nonempty_kernel(cuSZ3& n, cuBox box, int*& ngbrFlags, VType*& quantity, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	//need the idx < n.dim() check before ngbrFlags[idx] to avoid bad memory access
	reduction_avg(idx, n.dim(), quantity, average_value, points_count, idx < n.dim() && box.Contains(ijk) && (ngbrFlags[idx] & NF_NOTEMPTY));
}

template float cuVEC_VC<float>::average_nonempty(size_t arr_size, cuBox box);
template double cuVEC_VC<double>::average_nonempty(size_t arr_size, cuBox box);

template cuFLT3 cuVEC_VC<cuFLT3>::average_nonempty(size_t arr_size, cuBox box);
template cuDBL3 cuVEC_VC<cuDBL3>::average_nonempty(size_t arr_size, cuBox box);

template <typename VType>
__host__ VType cuVEC_VC<VType>::average_nonempty(size_t arr_size, cuBox box)
{
	zero_aux_values << <1, CUDATHREADS >> > (aux_value, aux_value2, aux_value3, aux_real, aux_integer);

	if (box.IsNull()) {

		cuvec_vc_average_nonempty_kernel << < (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (n, ngbrFlags, quantity, aux_value, aux_integer);
	}
	else cuvec_vc_average_nonempty_kernel << < (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (n, box, ngbrFlags, quantity, aux_value, aux_integer);

	VType av = get_gpu_value(aux_value);
	size_t points = get_gpu_value(aux_integer);

	if (points) return av / points;
	else return VType();
}

//------------------------------------------------------------------- AVERAGE NONEMPTY using a passed Rect

template <typename VType>
__global__ void average_nonempty_kernel(cuSZ3& n, cuRect rectangle, cuVEC_VC<VType>& cuvec, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	reduction_avg(idx, n.dim(), cuvec.data(), average_value, points_count, cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(ijk));
}

template float cuVEC_VC<float>::average_nonempty(size_t arr_size, cuRect rectangle);
template double cuVEC_VC<double>::average_nonempty(size_t arr_size, cuRect rectangle);

template cuFLT3 cuVEC_VC<cuFLT3>::average_nonempty(size_t arr_size, cuRect rectangle);
template cuDBL3 cuVEC_VC<cuDBL3>::average_nonempty(size_t arr_size, cuRect rectangle);

template <typename VType>
__host__ VType cuVEC_VC<VType>::average_nonempty(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average_nonempty(arr_size, cuBox());

	zero_aux_values << <1, CUDATHREADS >> > (aux_value, aux_value2, aux_value3, aux_real, aux_integer);

	average_nonempty_kernel << < (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (n, rectangle, *this, aux_value, aux_integer);

	VType av = get_gpu_value(aux_value);
	size_t points = get_gpu_value(aux_integer);

	if (points) return av / points;
	else return VType();
}