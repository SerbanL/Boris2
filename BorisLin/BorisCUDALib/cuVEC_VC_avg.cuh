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
	reduction_avg(idx, n.dim(), quantity, average_value, points_count, idx < n.dim() && (ngbrFlags[idx] & NF_NOTEMPTY));
}

template <typename VType>
__global__ void cuvec_vc_average_nonempty_kernel(cuSZ3& n, cuBox box, int*& ngbrFlags, VType*& quantity, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	//need the idx < n.dim() check before ngbrFlags[idx] to avoid bad memory access
	reduction_avg(idx, n.dim(), quantity, average_value, points_count, idx < n.dim() && box.Contains(ijk) && (ngbrFlags[idx] & NF_NOTEMPTY));
}

template float cuVEC_VC<float>::average_nonempty(size_t arr_size, cuBox box, bool transfer_to_cpu);
template double cuVEC_VC<double>::average_nonempty(size_t arr_size, cuBox box, bool transfer_to_cpu);

template cuFLT3 cuVEC_VC<cuFLT3>::average_nonempty(size_t arr_size, cuBox box, bool transfer_to_cpu);
template cuDBL3 cuVEC_VC<cuDBL3>::average_nonempty(size_t arr_size, cuBox box, bool transfer_to_cpu);

template <typename VType>
__host__ VType cuVEC_VC<VType>::average_nonempty(size_t arr_size, cuBox box, bool transfer_to_cpu)
{
	zero_aux_values <<<1, CUDATHREADS >>> (cuVEC<VType>::aux_value, cuVEC<VType>::aux_value2, cuVEC<VType>::aux_value3, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2, cuVEC<VType>::aux_integer);

	if (box.IsNull()) {

		cuvec_vc_average_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, ngbrFlags, cuVEC<VType>::quantity, cuVEC<VType>::aux_value, cuVEC<VType>::aux_integer);
	}
	else cuvec_vc_average_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, box, ngbrFlags, cuVEC<VType>::quantity, cuVEC<VType>::aux_value, cuVEC<VType>::aux_integer);

	if (transfer_to_cpu) {

		VType av = get_gpu_value(cuVEC<VType>::aux_value);
		size_t points = get_gpu_value(cuVEC<VType>::aux_integer);

		if (points) return av / points;
		else return VType();
	}
	else {

		//not transferring average value to cpu, keep it in gpu memory in auxiliary value
		return VType();
	}
}

//------------------------------------------------------------------- AVERAGE NONEMPTY using a passed Rect

template <typename VType>
__global__ void average_nonempty_kernel(cuSZ3& n, cuRect rectangle, cuVEC_VC<VType>& cuvec, VType& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	//need the idx < n.dim() check before cuvec.is_not_empty(ijk) to avoid bad memory access
	reduction_avg(idx, n.dim(), cuvec.data(), average_value, points_count, idx < n.dim() && cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(ijk));
}

template float cuVEC_VC<float>::average_nonempty(size_t arr_size, cuRect rectangle, bool transfer_to_cpu);
template double cuVEC_VC<double>::average_nonempty(size_t arr_size, cuRect rectangle, bool transfer_to_cpu);

template cuFLT3 cuVEC_VC<cuFLT3>::average_nonempty(size_t arr_size, cuRect rectangle, bool transfer_to_cpu);
template cuDBL3 cuVEC_VC<cuDBL3>::average_nonempty(size_t arr_size, cuRect rectangle, bool transfer_to_cpu);

template <typename VType>
__host__ VType cuVEC_VC<VType>::average_nonempty(size_t arr_size, cuRect rectangle, bool transfer_to_cpu)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average_nonempty(arr_size, cuBox(), transfer_to_cpu);

	//if rect start and end point are the same, then just read single value
	if (rectangle.s == rectangle.e) {
		
		VType value = VType();
		int index = cuVEC<VType>::position_to_cellidx_cpu(rectangle.s);
		if (index < arr_size) gpu_to_cpu_managed(&value, cuVEC<VType>::quantity, 1, index);
		return value;
	}

	zero_aux_values <<<1, CUDATHREADS >>> (cuVEC<VType>::aux_value, cuVEC<VType>::aux_value2, cuVEC<VType>::aux_value3, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2, cuVEC<VType>::aux_integer);

	average_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, rectangle, *this, cuVEC<VType>::aux_value, cuVEC<VType>::aux_integer);

	if (transfer_to_cpu) {

		VType av = get_gpu_value(cuVEC<VType>::aux_value);
		size_t points = get_gpu_value(cuVEC<VType>::aux_integer);

		if (points) return av / points;
		else return VType();
	}
	else {

		//not transferring average value to cpu, keep it in gpu memory in auxiliary value
		return VType();
	}
}

//------------------------------------------------------------------- SUM NONEMPTY (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_sum_nonempty_kernel(cuSZ3& n, int*& ngbrFlags, VType*& quantity, VType& sum_value)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//need the idx < n.dim() check before ngbrFlags[idx] to avoid bad memory access
	reduction_sum(idx, n.dim(), quantity, sum_value, idx < n.dim() && (ngbrFlags[idx] & NF_NOTEMPTY));
}

template <typename VType>
__global__ void cuvec_vc_sum_nonempty_kernel(cuSZ3& n, cuBox box, int*& ngbrFlags, VType*& quantity, VType& sum_value)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	//need the idx < n.dim() check before ngbrFlags[idx] to avoid bad memory access
	reduction_sum(idx, n.dim(), quantity, sum_value, idx < n.dim() && box.Contains(ijk) && (ngbrFlags[idx] & NF_NOTEMPTY));
}

template float cuVEC_VC<float>::sum_nonempty(size_t arr_size, cuBox box, bool transfer_to_cpu);
template double cuVEC_VC<double>::sum_nonempty(size_t arr_size, cuBox box, bool transfer_to_cpu);

template cuFLT3 cuVEC_VC<cuFLT3>::sum_nonempty(size_t arr_size, cuBox box, bool transfer_to_cpu);
template cuDBL3 cuVEC_VC<cuDBL3>::sum_nonempty(size_t arr_size, cuBox box, bool transfer_to_cpu);

template <typename VType>
__host__ VType cuVEC_VC<VType>::sum_nonempty(size_t arr_size, cuBox box, bool transfer_to_cpu)
{
	zero_aux_values <<<1, CUDATHREADS >>> (cuVEC<VType>::aux_value, cuVEC<VType>::aux_value2, cuVEC<VType>::aux_value3, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2, cuVEC<VType>::aux_integer);

	if (box.IsNull()) {

		cuvec_vc_sum_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, ngbrFlags, cuVEC<VType>::quantity, cuVEC<VType>::aux_value);
	}
	else cuvec_vc_sum_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, box, ngbrFlags, cuVEC<VType>::quantity, cuVEC<VType>::aux_value);

	if (transfer_to_cpu) {

		return get_gpu_value(cuVEC<VType>::aux_value);
	}
	else {

		//not transferring average value to cpu, keep it in gpu memory in auxiliary value
		return VType();
	}
}

//------------------------------------------------------------------- AVERAGE NONEMPTY using a passed Rect

template <typename VType>
__global__ void sum_nonempty_kernel(cuSZ3& n, cuRect rectangle, cuVEC_VC<VType>& cuvec, VType& sum_value)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	//need the idx < n.dim() check before cuvec.is_not_empty(ijk) to avoid bad memory access
	reduction_sum(idx, n.dim(), cuvec.data(), sum_value, idx < n.dim() && cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(ijk));
}

template float cuVEC_VC<float>::sum_nonempty(size_t arr_size, cuRect rectangle, bool transfer_to_cpu);
template double cuVEC_VC<double>::sum_nonempty(size_t arr_size, cuRect rectangle, bool transfer_to_cpu);

template cuFLT3 cuVEC_VC<cuFLT3>::sum_nonempty(size_t arr_size, cuRect rectangle, bool transfer_to_cpu);
template cuDBL3 cuVEC_VC<cuDBL3>::sum_nonempty(size_t arr_size, cuRect rectangle, bool transfer_to_cpu);

template <typename VType>
__host__ VType cuVEC_VC<VType>::sum_nonempty(size_t arr_size, cuRect rectangle, bool transfer_to_cpu)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return sum_nonempty(arr_size, cuBox(), transfer_to_cpu);

	zero_aux_values <<<1, CUDATHREADS >>> (cuVEC<VType>::aux_value, cuVEC<VType>::aux_value2, cuVEC<VType>::aux_value3, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2, cuVEC<VType>::aux_integer);

	sum_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, rectangle, *this, cuVEC<VType>::aux_value);

	if (transfer_to_cpu) {

		return get_gpu_value(cuVEC<VType>::aux_value);
	}
	else {

		//not transferring average value to cpu, keep it in gpu memory in auxiliary value
		return VType();
	}
}