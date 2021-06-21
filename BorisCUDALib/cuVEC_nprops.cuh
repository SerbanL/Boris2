#pragma once

#include "cuVEC.h"
#include "cuFuncs_Math.h"
#include "launchers.h"
#include "Reduction.cuh"

inline __global__ void set_minmaxreduction_aux_values(cuBReal& min, cuBReal& max)
{
	if (threadIdx.x == 0) min = 1e38;
	if (threadIdx.x == 1) max = -1e38;
}

//------------------------------------------------------------------- MINMAX by MAGNITUDE (cuVEC)

template <typename VType>
__global__ void cuvec_minmax_nonempty_kernel(cuVEC<VType>& cuvec, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = 0.0;
	bool include = false;

	if (idx < cuvec.linear_size() && cuvec.is_not_empty(idx)) {

		value = cu_GetMagnitude(cuvec[idx]);
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template <typename VType>
__global__ void cuvec_minmax_nonempty_kernel(cuVEC<VType>& cuvec, cuBox box, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < cuvec.linear_size() && box.Contains(ijk) && cuvec.is_not_empty(idx)) {

		value = cu_GetMagnitude(cuvec[idx]);
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC<float>::get_minmax(size_t arr_size, cuBox box);
template cuDBL2 cuVEC<double>::get_minmax(size_t arr_size, cuBox box);

template cuFLT2 cuVEC<cuFLT3>::get_minmax(size_t arr_size, cuBox box);
template cuDBL2 cuVEC<cuDBL3>::get_minmax(size_t arr_size, cuBox box);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC<VType>::get_minmax(size_t arr_size, cuBox box)
{
	set_minmaxreduction_aux_values <<<1, CUDATHREADS >>> (aux_real, aux_real2);

	if (box.IsNull()) {

		cuvec_minmax_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, aux_real, aux_real2);
	}
	else cuvec_minmax_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, box, aux_real, aux_real2);

	cuBReal min = get_gpu_value(aux_real);
	cuBReal max = get_gpu_value(aux_real2);

	return cuVAL2<PType>(min, max);
}

//------------------------------------------------------------------- MINMAX by MAGNITUDE using a passed Rect (cuVEC)

template <typename VType>
__global__ void cuvec_minmax_nonempty_kernel(cuVEC<VType>& cuvec, cuRect rectangle, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(idx)) {

		value = cu_GetMagnitude(cuvec[idx]);
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC<float>::get_minmax(size_t arr_size, cuRect rectangle);
template cuDBL2 cuVEC<double>::get_minmax(size_t arr_size, cuRect rectangle);

template cuFLT2 cuVEC<cuFLT3>::get_minmax(size_t arr_size, cuRect rectangle);
template cuDBL2 cuVEC<cuDBL3>::get_minmax(size_t arr_size, cuRect rectangle);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC<VType>::get_minmax(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return get_minmax(arr_size, cuBox());

	set_minmaxreduction_aux_values <<<1, CUDATHREADS >>> (aux_real, aux_real2);

	cuvec_minmax_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, rectangle, aux_real, aux_real2);

	cuBReal min = get_gpu_value(aux_real);
	cuBReal max = get_gpu_value(aux_real2);

	return cuVAL2<PType>(min, max);
}

//--------------------------------------------GET MIN-MAX COMPONENT X (cuVEC)

template <typename VType>
__global__ void cuvec_minmax_component_x_nonempty_kernel(cuVEC<VType>& cuvec, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].x;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template <typename VType>
__global__ void cuvec_minmax_component_x_nonempty_kernel(cuVEC<VType>& cuvec, cuBox box, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && box.Contains(ijk) && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].x;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC<cuFLT3>::get_minmax_component_x(size_t arr_size, cuBox box);
template cuDBL2 cuVEC<cuDBL3>::get_minmax_component_x(size_t arr_size, cuBox box);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC<VType>::get_minmax_component_x(size_t arr_size, cuBox box)
{
	set_minmaxreduction_aux_values <<<1, CUDATHREADS >>> (aux_real, aux_real2);

	if (box.IsNull()) {

		cuvec_minmax_component_x_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, aux_real, aux_real2);
	}
	else cuvec_minmax_component_x_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, box, aux_real, aux_real2);

	cuBReal min = get_gpu_value(aux_real);
	cuBReal max = get_gpu_value(aux_real2);

	return cuVAL2<PType>(min, max);
}

//------------------------------------------------------------------- MINMAX COMPONENT X using a passed Rect (cuVEC)

template <typename VType>
__global__ void cuvec_minmax_component_x_nonempty_kernel(cuVEC<VType>& cuvec, cuRect rectangle, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].x;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC<cuFLT3>::get_minmax_component_x(size_t arr_size, cuRect rectangle);
template cuDBL2 cuVEC<cuDBL3>::get_minmax_component_x(size_t arr_size, cuRect rectangle);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC<VType>::get_minmax_component_x(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return get_minmax_component_x(arr_size, cuBox());

	set_minmaxreduction_aux_values <<<1, CUDATHREADS >>> (aux_real, aux_real2);

	cuvec_minmax_component_x_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, rectangle, aux_real, aux_real2);

	cuBReal min = get_gpu_value(aux_real);
	cuBReal max = get_gpu_value(aux_real2);

	return cuVAL2<PType>(min, max);
}

//--------------------------------------------GET MIN-MAX COMPONENT Y (cuVEC)

template <typename VType>
__global__ void cuvec_minmax_component_y_nonempty_kernel(cuVEC<VType>& cuvec, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].y;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template <typename VType>
__global__ void cuvec_minmax_component_y_nonempty_kernel(cuVEC<VType>& cuvec, cuBox box, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && box.Contains(ijk) && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].y;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC<cuFLT3>::get_minmax_component_y(size_t arr_size, cuBox box);
template cuDBL2 cuVEC<cuDBL3>::get_minmax_component_y(size_t arr_size, cuBox box);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC<VType>::get_minmax_component_y(size_t arr_size, cuBox box)
{
	set_minmaxreduction_aux_values <<<1, CUDATHREADS >>> (aux_real, aux_real2);

	if (box.IsNull()) {

		cuvec_minmax_component_y_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, aux_real, aux_real2);
	}
	else cuvec_minmax_component_y_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, box, aux_real, aux_real2);

	cuBReal min = get_gpu_value(aux_real);
	cuBReal max = get_gpu_value(aux_real2);

	return cuVAL2<PType>(min, max);
}

//------------------------------------------------------------------- MINMAX COMPONENT Y using a passed Rect (cuVEC)

template <typename VType>
__global__ void cuvec_minmax_component_y_nonempty_kernel(cuVEC<VType>& cuvec, cuRect rectangle, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].y;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC<cuFLT3>::get_minmax_component_y(size_t arr_size, cuRect rectangle);
template cuDBL2 cuVEC<cuDBL3>::get_minmax_component_y(size_t arr_size, cuRect rectangle);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC<VType>::get_minmax_component_y(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return get_minmax_component_y(arr_size, cuBox());

	set_minmaxreduction_aux_values <<<1, CUDATHREADS >>> (aux_real, aux_real2);

	cuvec_minmax_component_y_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, rectangle, aux_real, aux_real2);

	cuBReal min = get_gpu_value(aux_real);
	cuBReal max = get_gpu_value(aux_real2);

	return cuVAL2<PType>(min, max);
}

//--------------------------------------------GET MIN-MAX COMPONENT Z (cuVEC)

template <typename VType>
__global__ void cuvec_minmax_component_z_nonempty_kernel(cuVEC<VType>& cuvec, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].z;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template <typename VType>
__global__ void cuvec_minmax_component_z_nonempty_kernel(cuVEC<VType>& cuvec, cuBox box, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && box.Contains(ijk) && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].z;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC<cuFLT3>::get_minmax_component_z(size_t arr_size, cuBox box);
template cuDBL2 cuVEC<cuDBL3>::get_minmax_component_z(size_t arr_size, cuBox box);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC<VType>::get_minmax_component_z(size_t arr_size, cuBox box)
{
	set_minmaxreduction_aux_values <<<1, CUDATHREADS >>> (aux_real, aux_real2);

	if (box.IsNull()) {

		cuvec_minmax_component_z_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, aux_real, aux_real2);
	}
	else cuvec_minmax_component_z_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, box, aux_real, aux_real2);

	cuBReal min = get_gpu_value(aux_real);
	cuBReal max = get_gpu_value(aux_real2);

	return cuVAL2<PType>(min, max);
}

//------------------------------------------------------------------- MINMAX COMPONENT Z using a passed Rect (cuVEC)

template <typename VType>
__global__ void cuvec_minmax_component_z_nonempty_kernel(cuVEC<VType>& cuvec, cuRect rectangle, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(idx)) {

		value = cuvec[idx].z;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC<cuFLT3>::get_minmax_component_z(size_t arr_size, cuRect rectangle);
template cuDBL2 cuVEC<cuDBL3>::get_minmax_component_z(size_t arr_size, cuRect rectangle);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC<VType>::get_minmax_component_z(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return get_minmax_component_z(arr_size, cuBox());

	set_minmaxreduction_aux_values <<<1, CUDATHREADS >>> (aux_real, aux_real2);

	cuvec_minmax_component_z_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, rectangle, aux_real, aux_real2);

	cuBReal min = get_gpu_value(aux_real);
	cuBReal max = get_gpu_value(aux_real2);

	return cuVAL2<PType>(min, max);
}