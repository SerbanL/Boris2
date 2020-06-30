#pragma once

#include "cuVEC_VC.h"
#include "cuFuncs_Math.h"
#include "launchers.h"
#include "Reduction.cuh"

//------------------------------------------------------------------- ZERO/SET AUX VALUES

template <typename VType>
__global__ void set_aux_values_vc_x(VType*& quantity, cuBReal& aux_real, cuBReal& aux_real2)
{
	if (threadIdx.x == 0) aux_real = quantity[0].x;
	if (threadIdx.x == 1) aux_real2 = quantity[0].x;
}

template <typename VType>
__global__ void set_aux_values_vc_y(VType*& quantity, cuBReal& aux_real, cuBReal& aux_real2)
{
	if (threadIdx.x == 0) aux_real = quantity[0].y;
	if (threadIdx.x == 1) aux_real2 = quantity[0].y;
}

template <typename VType>
__global__ void set_aux_values_vc_z(VType*& quantity, cuBReal& aux_real, cuBReal& aux_real2)
{
	if (threadIdx.x == 0) aux_real = quantity[0].z;
	if (threadIdx.x == 1) aux_real2 = quantity[0].z;
}

//------------------------------------------------------------------- MINMAX by MAGNITUDE (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_nonempty_kernel(cuSZ3& n, int*& ngbrFlags, VType*& quantity, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = cu_GetMagnitude(quantity[0]);
	bool include = false;

	if (idx < n.dim() && (ngbrFlags[idx] & NF_NOTEMPTY)) {

		value = cu_GetMagnitude(quantity[idx]);
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template <typename VType>
__global__ void cuvec_vc_minmax_nonempty_kernel(cuSZ3& n, cuBox box, int*& ngbrFlags, VType*& quantity, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = cu_GetMagnitude(quantity[0]);
	bool include = false;

	if (idx < n.dim() && box.Contains(ijk) && (ngbrFlags[idx] & NF_NOTEMPTY)) {

		value = cu_GetMagnitude(quantity[idx]);
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC_VC<float>::get_minmax(size_t arr_size, cuBox box);
template cuDBL2 cuVEC_VC<double>::get_minmax(size_t arr_size, cuBox box);

template cuFLT2 cuVEC_VC<cuFLT3>::get_minmax(size_t arr_size, cuBox box);
template cuDBL2 cuVEC_VC<cuDBL3>::get_minmax(size_t arr_size, cuBox box);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC_VC<VType>::get_minmax(size_t arr_size, cuBox box)
{
	zero_aux_values <<<1, CUDATHREADS >>> (cuVEC<VType>::aux_value, cuVEC<VType>::aux_value2, cuVEC<VType>::aux_value3, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2, cuVEC<VType>::aux_integer);

	if (box.IsNull()) {

		cuvec_vc_minmax_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, ngbrFlags, cuVEC<VType>::quantity, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);
	}
	else cuvec_vc_minmax_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, box, ngbrFlags, cuVEC<VType>::quantity, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//------------------------------------------------------------------- MINMAX by MAGNITUDE using a passed Rect (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_nonempty_kernel(cuSZ3& n, cuRect rectangle, cuVEC_VC<VType>& cuvec, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = cu_GetMagnitude(cuvec[0]);
	bool include = false;

	if (idx < n.dim() && cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(ijk)) {

		value = cu_GetMagnitude(cuvec[idx]);
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC_VC<float>::get_minmax(size_t arr_size, cuRect rectangle);
template cuDBL2 cuVEC_VC<double>::get_minmax(size_t arr_size, cuRect rectangle);

template cuFLT2 cuVEC_VC<cuFLT3>::get_minmax(size_t arr_size, cuRect rectangle);
template cuDBL2 cuVEC_VC<cuDBL3>::get_minmax(size_t arr_size, cuRect rectangle);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC_VC<VType>::get_minmax(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return get_minmax(arr_size, cuBox());

	zero_aux_values << <1, CUDATHREADS >> > (cuVEC<VType>::aux_value, cuVEC<VType>::aux_value2, cuVEC<VType>::aux_value3, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2, cuVEC<VType>::aux_integer);

	cuvec_vc_minmax_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, rectangle, *this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//--------------------------------------------GET MIN-MAX COMPONENT X (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_component_x_nonempty_kernel(cuSZ3& n, int*& ngbrFlags, VType*& quantity, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = quantity[0].x;
	bool include = false;

	if (idx < n.dim() && (ngbrFlags[idx] & NF_NOTEMPTY)) {

		value = quantity[idx].x;
		include = true;
	}
	
	reduction_minmax(0, 1, &value, min, max, include);
}

template <typename VType>
__global__ void cuvec_vc_minmax_component_x_nonempty_kernel(cuSZ3& n, cuBox box, int*& ngbrFlags, VType*& quantity, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = quantity[0].x;
	bool include = false;

	if (idx < n.dim() && box.Contains(ijk) && (ngbrFlags[idx] & NF_NOTEMPTY)) {

		value = quantity[idx].x;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC_VC<cuFLT3>::get_minmax_component_x(size_t arr_size, cuBox box);
template cuDBL2 cuVEC_VC<cuDBL3>::get_minmax_component_x(size_t arr_size, cuBox box);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC_VC<VType>::get_minmax_component_x(size_t arr_size, cuBox box)
{
	set_aux_values_vc_x <<<1, CUDATHREADS >>> (cuVEC<VType>::quantity, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	if (box.IsNull()) {

		cuvec_vc_minmax_component_x_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, ngbrFlags, cuVEC<VType>::quantity, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);
	}
	else cuvec_vc_minmax_component_x_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, box, ngbrFlags, cuVEC<VType>::quantity, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//------------------------------------------------------------------- MINMAX COMPONENT X using a passed Rect (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_component_x_nonempty_kernel(cuSZ3& n, cuRect rectangle, cuVEC_VC<VType>& cuvec, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = cuvec[0].x;
	bool include = false;

	if (idx < n.dim() && cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(ijk)) {

		value = cuvec[idx].x;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC_VC<cuFLT3>::get_minmax_component_x(size_t arr_size, cuRect rectangle);
template cuDBL2 cuVEC_VC<cuDBL3>::get_minmax_component_x(size_t arr_size, cuRect rectangle);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC_VC<VType>::get_minmax_component_x(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return get_minmax_component_x(arr_size, cuBox());

	set_aux_values_vc_x <<<1, CUDATHREADS >>> (cuVEC<VType>::quantity, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuvec_vc_minmax_component_x_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, rectangle, *this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//--------------------------------------------GET MIN-MAX COMPONENT Y (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_component_y_nonempty_kernel(cuSZ3& n, int*& ngbrFlags, VType*& quantity, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = quantity[0].y;
	bool include = false;

	if (idx < n.dim() && (ngbrFlags[idx] & NF_NOTEMPTY)) {

		value = quantity[idx].y;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template <typename VType>
__global__ void cuvec_vc_minmax_component_y_nonempty_kernel(cuSZ3& n, cuBox box, int*& ngbrFlags, VType*& quantity, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = quantity[0].y;
	bool include = false;

	if (idx < n.dim() && box.Contains(ijk) && (ngbrFlags[idx] & NF_NOTEMPTY)) {

		value = quantity[idx].y;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC_VC<cuFLT3>::get_minmax_component_y(size_t arr_size, cuBox box);
template cuDBL2 cuVEC_VC<cuDBL3>::get_minmax_component_y(size_t arr_size, cuBox box);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC_VC<VType>::get_minmax_component_y(size_t arr_size, cuBox box)
{
	set_aux_values_vc_y <<<1, CUDATHREADS >>> (cuVEC<VType>::quantity, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	if (box.IsNull()) {

		cuvec_vc_minmax_component_y_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, ngbrFlags, cuVEC<VType>::quantity, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);
	}
	else cuvec_vc_minmax_component_y_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, box, ngbrFlags, cuVEC<VType>::quantity, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//------------------------------------------------------------------- MINMAX COMPONENT Y using a passed Rect (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_component_y_nonempty_kernel(cuSZ3& n, cuRect rectangle, cuVEC_VC<VType>& cuvec, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = cuvec[0].y;
	bool include = false;

	if (idx < n.dim() && cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(ijk)) {
		
		value = cuvec[idx].y;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC_VC<cuFLT3>::get_minmax_component_y(size_t arr_size, cuRect rectangle);
template cuDBL2 cuVEC_VC<cuDBL3>::get_minmax_component_y(size_t arr_size, cuRect rectangle);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC_VC<VType>::get_minmax_component_y(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return get_minmax_component_y(arr_size, cuBox());

	set_aux_values_vc_y <<<1, CUDATHREADS >>> (cuVEC<VType>::quantity, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuvec_vc_minmax_component_y_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, rectangle, *this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//--------------------------------------------GET MIN-MAX COMPONENT Z (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_component_z_nonempty_kernel(cuSZ3& n, int*& ngbrFlags, VType*& quantity, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value = quantity[0].z;
	bool include = false;

	if (idx < n.dim() && (ngbrFlags[idx] & NF_NOTEMPTY)) {

		value = quantity[idx].z;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template <typename VType>
__global__ void cuvec_vc_minmax_component_z_nonempty_kernel(cuSZ3& n, cuBox box, int*& ngbrFlags, VType*& quantity, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = quantity[0].z;
	bool include = false;

	if (idx < n.dim() && box.Contains(ijk) && (ngbrFlags[idx] & NF_NOTEMPTY)) {

		value = quantity[idx].z;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC_VC<cuFLT3>::get_minmax_component_z(size_t arr_size, cuBox box);
template cuDBL2 cuVEC_VC<cuDBL3>::get_minmax_component_z(size_t arr_size, cuBox box);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC_VC<VType>::get_minmax_component_z(size_t arr_size, cuBox box)
{
	set_aux_values_vc_z <<<1, CUDATHREADS >>> (cuVEC<VType>::quantity, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	if (box.IsNull()) {

		cuvec_vc_minmax_component_z_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, ngbrFlags, cuVEC<VType>::quantity, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);
	}
	else cuvec_vc_minmax_component_z_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, box, ngbrFlags, cuVEC<VType>::quantity, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//------------------------------------------------------------------- MINMAX COMPONENT Z using a passed Rect (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_component_z_nonempty_kernel(cuSZ3& n, cuRect rectangle, cuVEC_VC<VType>& cuvec, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = cuvec[0].z;
	bool include = false;

	if (idx < n.dim() && cuvec.box_from_rect_max(rectangle + cuvec.rect.s).Contains(ijk) && cuvec.is_not_empty(ijk)) {

		value = cuvec[idx].z;
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template cuFLT2 cuVEC_VC<cuFLT3>::get_minmax_component_z(size_t arr_size, cuRect rectangle);
template cuDBL2 cuVEC_VC<cuDBL3>::get_minmax_component_z(size_t arr_size, cuRect rectangle);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC_VC<VType>::get_minmax_component_z(size_t arr_size, cuRect rectangle)
{
	//if empty rectangle then average ove the entire mesh
	if (rectangle.IsNull()) return get_minmax_component_z(arr_size, cuBox());

	set_aux_values_vc_z <<<1, CUDATHREADS >>> (cuVEC<VType>::quantity, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuvec_vc_minmax_component_z_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuVEC<VType>::n, rectangle, *this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}