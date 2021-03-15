#pragma once

#include "cuVEC.h"
#include "cuFuncs_Math.h"
#include "launchers.h"
#include "Reduction.cuh"
#include "cuVEC_avg.h"

//------------------------------------------------------------------- AVERAGE

template <typename VType>
__global__ void average_kernel(cuSZ3& n, VType*& quantity, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	reduction_avg(idx, n.dim(), quantity, average_value, points_count, true);
}

template <typename VType>
__global__ void average_kernel(cuSZ3& n, cuBox box, VType*& quantity, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	reduction_avg(idx, n.dim(), quantity, average_value, points_count, box.Contains(ijk));
}

template float cuVEC<float>::average(size_t arr_size, cuBox box);
template double cuVEC<double>::average(size_t arr_size, cuBox box);

template cuFLT3 cuVEC<cuFLT3>::average(size_t arr_size, cuBox box);
template cuDBL3 cuVEC<cuDBL3>::average(size_t arr_size, cuBox box);

template <typename VType>
__host__ VType cuVEC<VType>::average(size_t arr_size, cuBox box)
{
	zero_aux_values <<<1, CUDATHREADS >>> (aux_value, aux_value2, aux_value3, aux_real, aux_real2, aux_integer);

	if (box.IsNull()) {

		average_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, quantity, aux_value, aux_integer);
	}
	else average_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, box, quantity, aux_value, aux_integer);

	VType av = get_gpu_value(aux_value);
	size_t points = get_gpu_value(aux_integer);

	if (points) return av / points;
	else return VType();
}

//------------------------------------------------------------------- AVERAGE using a passed cuRect

template <typename VType>
__global__ void average_kernel(cuSZ3& n, cuRect rectangle, cuVEC<VType>& cuvec, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	reduction_avg(idx, n.dim(), cuvec.data(), average_value, points_count, cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk));
}

template float cuVEC<float>::average(size_t arr_size, cuRect rectangle);
template double cuVEC<double>::average(size_t arr_size, cuRect rectangle);

template cuFLT3 cuVEC<cuFLT3>::average(size_t arr_size, cuRect rectangle);
template cuDBL3 cuVEC<cuDBL3>::average(size_t arr_size, cuRect rectangle);

//average over given rectangle (relative to this cuVEC's rect)
template <typename VType>
__host__ VType cuVEC<VType>::average(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average(arr_size, cuBox());

	zero_aux_values <<<1, CUDATHREADS >>> (aux_value, aux_value2, aux_value3, aux_real, aux_real2, aux_integer);

	average_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, rectangle, *this, aux_value, aux_integer);

	VType av = get_gpu_value(aux_value);
	size_t points = get_gpu_value(aux_integer);

	if (points) return av / points;
	else return VType();
}

//------------------------------------------------------------------- AVERAGE NONEMPTY

template <typename VType>
__global__ void average_nonempty_kernel(cuSZ3& n, VType*& quantity, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//need the idx < n.dim() check before quantity[idx] to avoid bad memory access
	reduction_avg(idx, n.dim(), quantity, average_value, points_count, idx < n.dim() && cuIsNZ(cu_GetMagnitude(quantity[idx])));
}

template <typename VType>
__global__ void average_nonempty_kernel(cuSZ3& n, cuBox box, VType*& quantity, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	//need the idx < n.dim() check before quantity[idx] to avoid bad memory access
	reduction_avg(idx, n.dim(), quantity, average_value, points_count, idx < n.dim() && box.Contains(ijk) && cuIsNZ(cu_GetMagnitude(quantity[idx])));
}

template float cuVEC<float>::average_nonempty(size_t arr_size, cuBox box);
template double cuVEC<double>::average_nonempty(size_t arr_size, cuBox box);

template cuFLT3 cuVEC<cuFLT3>::average_nonempty(size_t arr_size, cuBox box);
template cuDBL3 cuVEC<cuDBL3>::average_nonempty(size_t arr_size, cuBox box);

template <typename VType>
__host__ VType cuVEC<VType>::average_nonempty(size_t arr_size, cuBox box)
{
	zero_aux_values << <1, CUDATHREADS >> > (aux_value, aux_value2, aux_value3, aux_real, aux_real2, aux_integer);

	if (box.IsNull()) {

		average_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, quantity, aux_value, aux_integer);
	}
	else average_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, box, quantity, aux_value, aux_integer);

	VType av = get_gpu_value(aux_value);
	size_t points = get_gpu_value(aux_integer);

	if (points) return av / points;
	else return VType();
}

//------------------------------------------------------------------- AVERAGE NONEMPTY using a passed Rect

template <typename VType>
__global__ void average_nonempty_kernel(cuSZ3& n, cuRect rectangle, cuVEC<VType>& cuvec, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	//need the idx < n.dim() check before cuvec.is_not_empty(ijk) to avoid bad memory access
	reduction_avg(idx, n.dim(), cuvec.data(), average_value, points_count, idx < n.dim() && cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(ijk));
}

template float cuVEC<float>::average_nonempty(size_t arr_size, cuRect rectangle);
template double cuVEC<double>::average_nonempty(size_t arr_size, cuRect rectangle);

template cuFLT3 cuVEC<cuFLT3>::average_nonempty(size_t arr_size, cuRect rectangle);
template cuDBL3 cuVEC<cuDBL3>::average_nonempty(size_t arr_size, cuRect rectangle);

template <typename VType>
__host__ VType cuVEC<VType>::average_nonempty(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average_nonempty(arr_size, cuBox());

	zero_aux_values <<<1, CUDATHREADS >>> (aux_value, aux_value2, aux_value3, aux_real, aux_real2, aux_integer);

	average_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, rectangle, *this, aux_value, aux_integer);

	VType av = get_gpu_value(aux_value);
	size_t points = get_gpu_value(aux_integer);

	if (points) return av / points;
	else return VType();
}
