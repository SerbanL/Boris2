#pragma once

#include "cuVEC.h"
#include "cuFuncs_Math.h"
#include "launchers.h"
#include "Reduction.cuh"

//------------------------------------------------------------------- SETBOX

template <typename VType>
__global__ void setbox_kernel(cuSZ3& n, cuBox box, VType value, VType*& quantity)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	if (idx < n.dim()) {

		if (box.Contains(ijk)) quantity[idx] = value;
	}
}

template void cuVEC<char>::setbox(cuBox box, char value);
template void cuVEC<int>::setbox(cuBox box, int value);
template void cuVEC<unsigned>::setbox(cuBox box, unsigned value);
template void cuVEC<long>::setbox(cuBox box, long value);
template void cuVEC<size_t>::setbox(cuBox box, size_t value);
template void cuVEC<float>::setbox(cuBox box, float value);
template void cuVEC<double>::setbox(cuBox box, double value);

template void cuVEC<cuINT3>::setbox(cuBox box, cuINT3 value);
template void cuVEC<cuSZ3>::setbox(cuBox box, cuSZ3 value);
template void cuVEC<cuFLT3>::setbox(cuBox box, cuFLT3 value);
template void cuVEC<cuDBL3>::setbox(cuBox box, cuDBL3 value);

template void cuVEC<cuINT33>::setbox(cuBox box, cuINT33 value);
template void cuVEC<cuFLT33>::setbox(cuBox box, cuFLT33 value);
template void cuVEC<cuDBL33>::setbox(cuBox box, cuDBL33 value);

template void cuVEC<cuReIm>::setbox(cuBox box, cuReIm value);
template void cuVEC<cuReIm3>::setbox(cuBox box, cuReIm3 value);

template <typename VType>
__host__ void cuVEC<VType>::setbox(cuBox box, VType value)
{
	setbox_kernel <<< (get_gpu_value(n).dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, box, value, quantity);
}

//------------------------------------------------------------------- SET

template <typename VType>
__global__ void set_kernel(size_t size, VType value, VType*& quantity)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size) {

		quantity[idx] = value;
	}
}

template void cuVEC<char>::set(size_t size, char value);
template void cuVEC<int>::set(size_t size, int value);
template void cuVEC<unsigned>::set(size_t size, unsigned value);
template void cuVEC<long>::set(size_t size, long value);
template void cuVEC<size_t>::set(size_t size, size_t value);
template void cuVEC<float>::set(size_t size, float value);
template void cuVEC<double>::set(size_t size, double value);

template void cuVEC<cuINT3>::set(size_t size, cuINT3 value);
template void cuVEC<cuSZ3>::set(size_t size, cuSZ3 value);
template void cuVEC<cuFLT3>::set(size_t size, cuFLT3 value);
template void cuVEC<cuDBL3>::set(size_t size, cuDBL3 value);

template void cuVEC<cuINT33>::set(size_t size, cuINT33 value);
template void cuVEC<cuFLT33>::set(size_t size, cuFLT33 value);
template void cuVEC<cuDBL33>::set(size_t size, cuDBL33 value);

template void cuVEC<cuReIm>::set(size_t size, cuReIm value);
template void cuVEC<cuReIm3>::set(size_t size, cuReIm3 value);

template <typename VType>
__host__ void cuVEC<VType>::set(size_t size, VType value)
{
	set_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (size, value, quantity);
}

//------------------------------------------------------------------- RENORMALIZE

template <typename VType, typename PType>
__global__ void renormalize_kernel(cuSZ3& n, VType*& quantity, PType new_norm)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < n.dim()) {

		PType curr_norm = cu_GetMagnitude(quantity[idx]);

		if (cuIsNZ(curr_norm)) quantity[idx] *= new_norm / curr_norm;
	}
}

template void cuVEC<float>::renormalize(size_t arr_size, float new_norm);
template void cuVEC<double>::renormalize(size_t arr_size, double new_norm);

template void cuVEC<cuFLT3>::renormalize(size_t arr_size, float new_norm);
template void cuVEC<cuDBL3>::renormalize(size_t arr_size, double new_norm);

template <typename VType>
template <typename PType>
__host__ void cuVEC<VType>::renormalize(size_t arr_size, PType new_norm)
{
	renormalize_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, quantity, new_norm);
}

//------------------------------------------------------------------- AVERAGE

template <typename VType>
__global__ void average_kernel(cuSZ3& n, VType*& quantity, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	
	reduction_avg(idx, n.dim(), quantity, average_value, points_count);
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
	zero_aux_values <<<1, CUDATHREADS >>> (aux_value, aux_value2, aux_value3, aux_real, aux_integer);

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

	zero_aux_values <<<1, CUDATHREADS >>> (aux_value, aux_value2, aux_value3, aux_real, aux_integer);

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

	reduction_avg(idx, n.dim(), quantity, average_value, points_count, cuIsNZ(cu_GetMagnitude(quantity[idx])));
}

template <typename VType>
__global__ void average_nonempty_kernel(cuSZ3& n, cuBox box, VType*& quantity, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	reduction_avg(idx, n.dim(), quantity, average_value, points_count, box.Contains(ijk) && cuIsNZ(cu_GetMagnitude(quantity[idx])));
}

template float cuVEC<float>::average_nonempty(size_t arr_size, cuBox box);
template double cuVEC<double>::average_nonempty(size_t arr_size, cuBox box);

template cuFLT3 cuVEC<cuFLT3>::average_nonempty(size_t arr_size, cuBox box);
template cuDBL3 cuVEC<cuDBL3>::average_nonempty(size_t arr_size, cuBox box);

template <typename VType>
__host__ VType cuVEC<VType>::average_nonempty(size_t arr_size, cuBox box)
{
	zero_aux_values << <1, CUDATHREADS >> > (aux_value, aux_value2, aux_value3, aux_real, aux_integer);

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

	reduction_avg(idx, n.dim(), cuvec.data(), average_value, points_count, cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(ijk));
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

	zero_aux_values <<<1, CUDATHREADS >>> (aux_value, aux_value2, aux_value3, aux_real, aux_integer);

	average_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (n, rectangle, *this, aux_value, aux_integer);

	VType av = get_gpu_value(aux_value);
	size_t points = get_gpu_value(aux_integer);

	if (points) return av / points;
	else return VType();
}