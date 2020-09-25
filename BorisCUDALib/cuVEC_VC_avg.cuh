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

template float cuVEC_VC<float>::average_nonempty(size_t arr_size, cuBox box);
template double cuVEC_VC<double>::average_nonempty(size_t arr_size, cuBox box);

template cuFLT3 cuVEC_VC<cuFLT3>::average_nonempty(size_t arr_size, cuBox box);
template cuDBL3 cuVEC_VC<cuDBL3>::average_nonempty(size_t arr_size, cuBox box);

template <typename VType>
__host__ VType cuVEC_VC<VType>::average_nonempty(size_t arr_size, cuBox box)
{
	zero_aux_values <<<1, CUDATHREADS >>> (cuVEC<VType>::aux_value, cuVEC<VType>::aux_value2, cuVEC<VType>::aux_value3, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2, cuVEC<VType>::aux_integer);

	if (box.IsNull()) {

		cuvec_vc_average_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, ngbrFlags, cuVEC<VType>::quantity, cuVEC<VType>::aux_value, cuVEC<VType>::aux_integer);
	}
	else cuvec_vc_average_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, box, ngbrFlags, cuVEC<VType>::quantity, cuVEC<VType>::aux_value, cuVEC<VType>::aux_integer);

	VType av = get_gpu_value(cuVEC<VType>::aux_value);
	size_t points = get_gpu_value(cuVEC<VType>::aux_integer);

	if (points) return av / points;
	else return VType();
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

template float cuVEC_VC<float>::average_nonempty(size_t arr_size, cuRect rectangle);
template double cuVEC_VC<double>::average_nonempty(size_t arr_size, cuRect rectangle);

template cuFLT3 cuVEC_VC<cuFLT3>::average_nonempty(size_t arr_size, cuRect rectangle);
template cuDBL3 cuVEC_VC<cuDBL3>::average_nonempty(size_t arr_size, cuRect rectangle);

template <typename VType>
__host__ VType cuVEC_VC<VType>::average_nonempty(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average_nonempty(arr_size, cuBox());

	zero_aux_values <<<1, CUDATHREADS >>> (cuVEC<VType>::aux_value, cuVEC<VType>::aux_value2, cuVEC<VType>::aux_value3, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2, cuVEC<VType>::aux_integer);

	average_nonempty_kernel << < (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cuVEC<VType>::n, rectangle, *this, cuVEC<VType>::aux_value, cuVEC<VType>::aux_integer);

	VType av = get_gpu_value(cuVEC<VType>::aux_value);
	size_t points = get_gpu_value(cuVEC<VType>::aux_integer);

	if (points) return av / points;
	else return VType();
}

//------------------------------------------------------------------- AVERAGE XSQ NONEMPTY using a passed Rect

template <typename VType>
__global__ void average_xsq_nonempty_kernel(cuSZ3& n, cuRect rectangle, cuVEC_VC<VType>& cuvec, cuBReal& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(ijk)) {

		value = cuvec[idx].x * cuvec[idx].x;
		include = true;
	}

	reduction_avg(0, 1, &value, average_value, points_count, include);
}

template float cuVEC_VC<cuFLT3>::average_xsq_nonempty(size_t arr_size, cuRect rectangle);
template double cuVEC_VC<cuDBL3>::average_xsq_nonempty(size_t arr_size, cuRect rectangle);

template <typename VType>
template <typename PType>
__host__ PType cuVEC_VC<VType>::average_xsq_nonempty(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average_xsq_nonempty(arr_size, cuBox());

	zero_aux_values <<<1, CUDATHREADS >>> (cuVEC<VType>::aux_value, cuVEC<VType>::aux_value2, cuVEC<VType>::aux_value3, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2, cuVEC<VType>::aux_integer);

	average_xsq_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, rectangle, *this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_integer);

	cuBReal av = get_gpu_value(cuVEC<VType>::aux_real);
	size_t points = get_gpu_value(cuVEC<VType>::aux_integer);

	if (points) return av / points;
	else return PType();
}

//------------------------------------------------------------------- AVERAGE XSQ NONEMPTY (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_average_xsq_nonempty_kernel(cuSZ3& n, cuVEC_VC<VType>& cuvec, cuBReal& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].x * cuvec[idx].x;
		include = true;
	}

	reduction_avg(0, 1, &value, average_value, points_count, include);
}

template <typename VType>
__global__ void cuvec_vc_average_xsq_nonempty_kernel(cuSZ3& n, cuBox box, cuVEC_VC<VType>& cuvec, cuBReal& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && box.Contains(ijk) && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].x * cuvec[idx].x;
		include = true;
	}

	reduction_avg(0, 1, &value, average_value, points_count, include);
}

template cuFLT3 cuVEC_VC<cuFLT3>::average_xsq_nonempty(size_t arr_size, cuBox box);
template cuDBL3 cuVEC_VC<cuDBL3>::average_xsq_nonempty(size_t arr_size, cuBox box);

template <typename VType>
template <typename PType>
__host__ PType cuVEC_VC<VType>::average_xsq_nonempty(size_t arr_size, cuBox box)
{
	zero_aux_values <<<1, CUDATHREADS >>> (cuVEC<VType>::aux_value, cuVEC<VType>::aux_value2, cuVEC<VType>::aux_value3, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2, cuVEC<VType>::aux_integer);

	if (box.IsNull()) {

		cuvec_vc_average_xsq_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, *this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_integer);
	}
	else cuvec_vc_average_xsq_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, box, *this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_integer);

	cuBReal av = get_gpu_value(cuVEC<VType>::aux_real);
	size_t points = get_gpu_value(cuVEC<VType>::aux_integer);

	if (points) return av / points;
	else return PType();
}

//------------------------------------------------------------------- AVERAGE YSQ NONEMPTY using a passed Rect

template <typename VType>
__global__ void average_ysq_nonempty_kernel(cuSZ3& n, cuRect rectangle, cuVEC_VC<VType>& cuvec, cuBReal& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(ijk)) {

		value = cuvec[idx].y * cuvec[idx].y;
		include = true;
	}

	reduction_avg(0, 1, &value, average_value, points_count, include);
}

template float cuVEC_VC<cuFLT3>::average_ysq_nonempty(size_t arr_size, cuRect rectangle);
template double cuVEC_VC<cuDBL3>::average_ysq_nonempty(size_t arr_size, cuRect rectangle);

template <typename VType>
template <typename PType>
__host__ PType cuVEC_VC<VType>::average_ysq_nonempty(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average_ysq_nonempty(arr_size, cuBox());

	zero_aux_values <<<1, CUDATHREADS >>> (cuVEC<VType>::aux_value, cuVEC<VType>::aux_value2, cuVEC<VType>::aux_value3, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2, cuVEC<VType>::aux_integer);

	average_ysq_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, rectangle, *this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_integer);

	cuBReal av = get_gpu_value(cuVEC<VType>::aux_real);
	size_t points = get_gpu_value(cuVEC<VType>::aux_integer);

	if (points) return av / points;
	else return PType();
}

//------------------------------------------------------------------- AVERAGE YSQ NONEMPTY (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_average_ysq_nonempty_kernel(cuSZ3& n, cuVEC_VC<VType>& cuvec, cuBReal& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].y * cuvec[idx].y;
		include = true;
	}

	reduction_avg(0, 1, &value, average_value, points_count, include);
}

template <typename VType>
__global__ void cuvec_vc_average_ysq_nonempty_kernel(cuSZ3& n, cuBox box, cuVEC_VC<VType>& cuvec, cuBReal& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && box.Contains(ijk) && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].y * cuvec[idx].y;
		include = true;
	}

	reduction_avg(0, 1, &value, average_value, points_count, include);
}

template cuFLT3 cuVEC_VC<cuFLT3>::average_ysq_nonempty(size_t arr_size, cuBox box);
template cuDBL3 cuVEC_VC<cuDBL3>::average_ysq_nonempty(size_t arr_size, cuBox box);

template <typename VType>
template <typename PType>
__host__ PType cuVEC_VC<VType>::average_ysq_nonempty(size_t arr_size, cuBox box)
{
	zero_aux_values <<<1, CUDATHREADS >>> (cuVEC<VType>::aux_value, cuVEC<VType>::aux_value2, cuVEC<VType>::aux_value3, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2, cuVEC<VType>::aux_integer);

	if (box.IsNull()) {

		cuvec_vc_average_ysq_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, *this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_integer);
	}
	else cuvec_vc_average_ysq_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, box, *this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_integer);

	cuBReal av = get_gpu_value(cuVEC<VType>::aux_real);
	size_t points = get_gpu_value(cuVEC<VType>::aux_integer);

	if (points) return av / points;
	else return PType();
}

//------------------------------------------------------------------- AVERAGE ZSQ NONEMPTY using a passed Rect

template <typename VType>
__global__ void average_zsq_nonempty_kernel(cuSZ3& n, cuRect rectangle, cuVEC_VC<VType>& cuvec, cuBReal& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(ijk)) {

		value = cuvec[idx].z * cuvec[idx].z;
		include = true;
	}

	reduction_avg(0, 1, &value, average_value, points_count, include);
}

template float cuVEC_VC<cuFLT3>::average_zsq_nonempty(size_t arr_size, cuRect rectangle);
template double cuVEC_VC<cuDBL3>::average_zsq_nonempty(size_t arr_size, cuRect rectangle);

template <typename VType>
template <typename PType>
__host__ PType cuVEC_VC<VType>::average_zsq_nonempty(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return average_zsq_nonempty(arr_size, cuBox());

	zero_aux_values <<<1, CUDATHREADS >>> (cuVEC<VType>::aux_value, cuVEC<VType>::aux_value2, cuVEC<VType>::aux_value3, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2, cuVEC<VType>::aux_integer);

	average_zsq_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, rectangle, *this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_integer);

	cuBReal av = get_gpu_value(cuVEC<VType>::aux_real);
	size_t points = get_gpu_value(cuVEC<VType>::aux_integer);

	if (points) return av / points;
	else return PType();
}

//------------------------------------------------------------------- AVERAGE ZSQ NONEMPTY (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_average_zsq_nonempty_kernel(cuSZ3& n, cuVEC_VC<VType>& cuvec, cuBReal& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].z * cuvec[idx].z;
		include = true;
	}

	reduction_avg(0, 1, &value, average_value, points_count, include);
}

template <typename VType>
__global__ void cuvec_vc_average_zsq_nonempty_kernel(cuSZ3& n, cuBox box, cuVEC_VC<VType>& cuvec, cuBReal& average_value, size_t& points_count)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && box.Contains(ijk) && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].z * cuvec[idx].z;
		include = true;
	}

	reduction_avg(0, 1, &value, average_value, points_count, include);
}

template cuFLT3 cuVEC_VC<cuFLT3>::average_zsq_nonempty(size_t arr_size, cuBox box);
template cuDBL3 cuVEC_VC<cuDBL3>::average_zsq_nonempty(size_t arr_size, cuBox box);

template <typename VType>
template <typename PType>
__host__ PType cuVEC_VC<VType>::average_zsq_nonempty(size_t arr_size, cuBox box)
{
	zero_aux_values <<<1, CUDATHREADS >>> (cuVEC<VType>::aux_value, cuVEC<VType>::aux_value2, cuVEC<VType>::aux_value3, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2, cuVEC<VType>::aux_integer);

	if (box.IsNull()) {

		cuvec_vc_average_zsq_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, *this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_integer);
	}
	else cuvec_vc_average_zsq_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, box, *this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_integer);

	cuBReal av = get_gpu_value(cuVEC<VType>::aux_real);
	size_t points = get_gpu_value(cuVEC<VType>::aux_integer);

	if (points) return av / points;
	else return PType();
}