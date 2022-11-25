#pragma once

#include "cuVEC_VC.h"
#include "cuFuncs_Math.h"
#include "launchers.h"
#include "Reduction.cuh"

//------------------------------------------------------------------- ZERO/SET AUX VALUES

inline __global__ void set_minmaxreduction_aux_values_vc(cuBReal& min, cuBReal& max)
{
	if (threadIdx.x == 0) min = 1e38;
	if (threadIdx.x == 1) max = -1e38;
}

//------------------------------------------------------------------- MINMAX by MAGNITUDE (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_nonempty_kernel(cuVEC_VC<VType>& cuvec, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;
	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && cuvec.is_not_empty(idx)) {

		value = cu_GetMagnitude(cuvec[idx]);
		include = true;
	}

	reduction_minmax(0, 1, &value, min, max, include);
}

template <typename VType>
__global__ void cuvec_vc_minmax_nonempty_kernel(cuVEC_VC<VType>& cuvec, cuBox box, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
	bool include = false;

	if (idx < n.dim() && box.Contains(ijk) && cuvec.is_not_empty(idx)) {

		value = cu_GetMagnitude(cuvec[idx]);
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
	set_minmaxreduction_aux_values_vc <<< 1, CUDATHREADS >>> (cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	if (box.IsNull()) {

		cuvec_vc_minmax_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);
	}
	else cuvec_vc_minmax_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, box, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//------------------------------------------------------------------- MINMAX by MAGNITUDE using a passed Rect (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_nonempty_kernel(cuVEC_VC<VType>& cuvec, cuRect rectangle, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
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

	set_minmaxreduction_aux_values_vc <<< 1, CUDATHREADS >>> (cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuvec_vc_minmax_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, rectangle, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//--------------------------------------------GET MIN-MAX COMPONENT X (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_component_x_nonempty_kernel(cuVEC_VC<VType>& cuvec, cuBReal& min, cuBReal& max)
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
__global__ void cuvec_vc_minmax_component_x_nonempty_kernel(cuVEC_VC<VType>& cuvec, cuBox box, cuBReal& min, cuBReal& max)
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

template cuFLT2 cuVEC_VC<cuFLT3>::get_minmax_component_x(size_t arr_size, cuBox box);
template cuDBL2 cuVEC_VC<cuDBL3>::get_minmax_component_x(size_t arr_size, cuBox box);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC_VC<VType>::get_minmax_component_x(size_t arr_size, cuBox box)
{
	set_minmaxreduction_aux_values_vc <<< 1, CUDATHREADS >>> (cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	if (box.IsNull()) {

		cuvec_vc_minmax_component_x_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);
	}
	else cuvec_vc_minmax_component_x_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, box, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//------------------------------------------------------------------- MINMAX COMPONENT X using a passed Rect (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_component_x_nonempty_kernel(cuVEC_VC<VType>& cuvec, cuRect rectangle, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
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

	set_minmaxreduction_aux_values_vc <<< 1, CUDATHREADS >>> (cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuvec_vc_minmax_component_x_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, rectangle, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//--------------------------------------------GET MIN-MAX COMPONENT Y (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_component_y_nonempty_kernel(cuVEC_VC<VType>& cuvec, cuBReal& min, cuBReal& max)
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
__global__ void cuvec_vc_minmax_component_y_nonempty_kernel(cuVEC_VC<VType>& cuvec, cuBox box, cuBReal& min, cuBReal& max)
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

template cuFLT2 cuVEC_VC<cuFLT3>::get_minmax_component_y(size_t arr_size, cuBox box);
template cuDBL2 cuVEC_VC<cuDBL3>::get_minmax_component_y(size_t arr_size, cuBox box);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC_VC<VType>::get_minmax_component_y(size_t arr_size, cuBox box)
{
	set_minmaxreduction_aux_values_vc <<< 1, CUDATHREADS >>> (cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	if (box.IsNull()) {

		cuvec_vc_minmax_component_y_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);
	}
	else cuvec_vc_minmax_component_y_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, box, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//------------------------------------------------------------------- MINMAX COMPONENT Y using a passed Rect (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_component_y_nonempty_kernel(cuVEC_VC<VType>& cuvec, cuRect rectangle, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
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

	set_minmaxreduction_aux_values_vc <<< 1, CUDATHREADS >>> (cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuvec_vc_minmax_component_y_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, rectangle, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//--------------------------------------------GET MIN-MAX COMPONENT Z (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_component_z_nonempty_kernel(cuVEC_VC<VType>& cuvec, cuBReal& min, cuBReal& max)
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
__global__ void cuvec_vc_minmax_component_z_nonempty_kernel(cuVEC_VC<VType>& cuvec, cuBox box, cuBReal& min, cuBReal& max)
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

template cuFLT2 cuVEC_VC<cuFLT3>::get_minmax_component_z(size_t arr_size, cuBox box);
template cuDBL2 cuVEC_VC<cuDBL3>::get_minmax_component_z(size_t arr_size, cuBox box);

template <typename VType>
template <typename PType>
__host__ cuVAL2<PType> cuVEC_VC<VType>::get_minmax_component_z(size_t arr_size, cuBox box)
{
	set_minmaxreduction_aux_values_vc <<< 1, CUDATHREADS >>> (cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	if (box.IsNull()) {

		cuvec_vc_minmax_component_z_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);
	}
	else cuvec_vc_minmax_component_z_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, box, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//------------------------------------------------------------------- MINMAX COMPONENT Z using a passed Rect (cuVEC_VC)

template <typename VType>
__global__ void cuvec_vc_minmax_component_z_nonempty_kernel(cuVEC_VC<VType>& cuvec, cuRect rectangle, cuBReal& min, cuBReal& max)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3& n = cuvec.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuBReal value = 0.0;
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

	set_minmaxreduction_aux_values_vc <<< 1, CUDATHREADS >>> (cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuvec_vc_minmax_component_z_nonempty_kernel <<< (arr_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (*this, rectangle, cuVEC<VType>::aux_real, cuVEC<VType>::aux_real2);

	cuBReal min = get_gpu_value(cuVEC<VType>::aux_real);
	cuBReal max = get_gpu_value(cuVEC<VType>::aux_real2);

	return cuVAL2<PType>(min, max);
}

//--------------------------------------------SPECIAL NUMERICAL PROPERTIES : cuVEC_VC_nprops.cuh

//Robin value is of the form alpha * (Tb - Ta). alpha and Ta are known from values set robin_nx, robin_px, ...
//Tb is quantity value at boundary. Here we will return Robin value for x, y, z axes for which any shift component is nonzero (otherwise zero for that component)
//e.g. if shift.x is non-zero then Tb value is obtained at rel_pos + (shift.x, 0, 0) using extrapolation from values at rel_pos and rel_pos - (shift.x, 0, 0) -> both these values should be inside the mesh, else return zero.
template <typename VType>
__device__ cuReal3 cuVEC_VC<VType>::get_robin_value(const cuReal3& rel_pos, const cuReal3& shift)
{
	cuReal3 robin_values;

	if (cuIsNZ(shift.x)) {

		VType val1 = (*this)[rel_pos];
		VType val2 = (*this)[rel_pos - cuReal3(shift.x, 0, 0)];
		VType val_bnd = (1.5 * val1 - 0.5 * val2);

		if (shift.x > 0) robin_values.x = robin_px.i * (val_bnd - robin_px.j);
		else robin_values.x = robin_nx.i * (val_bnd - robin_nx.j);
	}

	if (cuIsNZ(shift.y)) {

		VType val1 = (*this)[rel_pos];
		VType val2 = (*this)[rel_pos - cuReal3(0, shift.y, 0)];
		VType val_bnd = (1.5 * val1 - 0.5 * val2);

		if (shift.y > 0) robin_values.y = robin_py.i * (val_bnd - robin_py.j);
		else robin_values.y = robin_ny.i * (val_bnd - robin_ny.j);
	}

	if (cuIsNZ(shift.z)) {

		VType val1 = (*this)[rel_pos];
		VType val2 = (*this)[rel_pos - cuReal3(0, 0, shift.z)];
		VType val_bnd = (1.5 * val1 - 0.5 * val2);

		if (shift.z > 0) robin_values.z = robin_pz.i * (val_bnd - robin_pz.j);
		else robin_values.z = robin_nz.i * (val_bnd - robin_nz.j);
	}

	return robin_values;
}

//for given cell index, find if any neighboring cells are empty and get distance (shift) valuer to them along each axis
//if any shift is zero this means both cells are present either side, or both are missing
//NOTE : this is intended to be used with get_robin_value method to determine the shift value, and rel_pos will be position corresponding to idx
template <typename VType>
__device__ cuReal3 cuVEC_VC<VType>::get_shift_to_emptycell(int idx)
{
	cuReal3 shift;

	//x
	if ((ngbrFlags[idx] & NF_BOTHX) != NF_BOTHX) {

		if (ngbrFlags[idx] & NF_NPX) shift.x = -cuVEC<VType>::h.x / 2;
		else if (ngbrFlags[idx] & NF_NNX) shift.x = +cuVEC<VType>::h.x / 2;
	}

	//y
	if ((ngbrFlags[idx] & NF_BOTHY) != NF_BOTHY) {

		if (ngbrFlags[idx] & NF_NPY) shift.y = -cuVEC<VType>::h.y / 2;
		else if (ngbrFlags[idx] & NF_NNY) shift.y = +cuVEC<VType>::h.y / 2;
	}

	//z
	if ((ngbrFlags[idx] & NF_BOTHZ) != NF_BOTHZ) {

		if (ngbrFlags[idx] & NF_NPZ) shift.z = -cuVEC<VType>::h.z / 2;
		else if (ngbrFlags[idx] & NF_NNZ) shift.z = +cuVEC<VType>::h.z / 2;
	}

	return shift;
}