#pragma once

#include "VEC_VC.h"

//--------------------------------------------GET MIN-MAX by MAGNITUDE

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax(const Box& box) const
{
	VEC<VType>::magnitude_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * VEC<VType>::n.x + (k + box.s.k) * VEC<VType>::n.x * VEC<VType>::n.y;

		if (idx < 0 || idx >= VEC<VType>::n.dim() || !(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

		VEC<VType>::magnitude_reduction.reduce_minmax(GetMagnitude(VEC<VType>::quantity[idx]));
	}

	return VEC<VType>::magnitude_reduction.minmax();
}

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax(const Rect& rectangle) const
{
	//if empty rectangle then use entire mesh
	if (rectangle.IsNull()) return get_minmax(Box(VEC<VType>::n));

	//... otherwise rectangle must intersect with this mesh
	if (!VEC<VType>::rect.intersects(rectangle + VEC<VType>::rect.s)) return VAL2<PType>();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return get_minmax(VEC<VType>::box_from_rect_max(rectangle + VEC<VType>::rect.s));
}

//--------------------------------------------GET MIN-MAX COMPONENT X

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax_component_x(const Box& box) const
{
	VEC<VType>::magnitude_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * VEC<VType>::n.x + (k + box.s.k) * VEC<VType>::n.x * VEC<VType>::n.y;

		if (idx < 0 || idx >= VEC<VType>::n.dim() || !(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

		VEC<VType>::magnitude_reduction.reduce_minmax(VEC<VType>::quantity[idx].x);
	}

	return VEC<VType>::magnitude_reduction.minmax();
}

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax_component_x(const Rect& rectangle) const
{
	//if empty rectangle then use entire mesh
	if (rectangle.IsNull()) return get_minmax_component_x(Box(VEC<VType>::n));

	//... otherwise rectangle must intersect with this mesh
	if (!VEC<VType>::rect.intersects(rectangle + VEC<VType>::rect.s)) return VAL2<PType>();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return get_minmax_component_x(VEC<VType>::box_from_rect_max(rectangle + VEC<VType>::rect.s));
}

//--------------------------------------------GET MIN-MAX COMPONENT Y

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax_component_y(const Box& box) const
{
	VEC<VType>::magnitude_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * VEC<VType>::n.x + (k + box.s.k) * VEC<VType>::n.x * VEC<VType>::n.y;

		if (idx < 0 || idx >= VEC<VType>::n.dim() || !(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

		VEC<VType>::magnitude_reduction.reduce_minmax(VEC<VType>::quantity[idx].y);
	}

	return VEC<VType>::magnitude_reduction.minmax();
}

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax_component_y(const Rect& rectangle) const
{
	//if empty rectangle then use entire mesh
	if (rectangle.IsNull()) return get_minmax_component_y(Box(VEC<VType>::n));

	//... otherwise rectangle must intersect with this mesh
	if (!VEC<VType>::rect.intersects(rectangle + VEC<VType>::rect.s)) return VAL2<PType>();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return get_minmax_component_y(VEC<VType>::box_from_rect_max(rectangle + VEC<VType>::rect.s));
}

//--------------------------------------------GET MIN-MAX COMPONENT Z

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax_component_z(const Box& box) const
{
	VEC<VType>::magnitude_reduction.new_minmax_reduction();

#pragma omp parallel for
	for (int idx_box = 0; idx_box < box.size().dim(); idx_box++) {

		//i, j, k values inside the box only calculated from the box cell index
		int i = (idx_box % box.size().x);
		int j = ((idx_box / box.size().x) % box.size().y);
		int k = (idx_box / (box.size().x * box.size().y));

		//index inside the mesh for this box cell index
		int idx = (i + box.s.i) + (j + box.s.j) * VEC<VType>::n.x + (k + box.s.k) * VEC<VType>::n.x * VEC<VType>::n.y;

		if (idx < 0 || idx >= VEC<VType>::n.dim() || !(ngbrFlags[idx] & NF_NOTEMPTY)) continue;

		VEC<VType>::magnitude_reduction.reduce_minmax(VEC<VType>::quantity[idx].z);
	}

	return VEC<VType>::magnitude_reduction.minmax();
}

template <typename VType>
template <typename PType>
VAL2<PType> VEC_VC<VType>::get_minmax_component_z(const Rect& rectangle) const
{
	//if empty rectangle then use entire mesh
	if (rectangle.IsNull()) return get_minmax_component_z(Box(VEC<VType>::n));

	//... otherwise rectangle must intersect with this mesh
	if (!VEC<VType>::rect.intersects(rectangle + VEC<VType>::rect.s)) return VAL2<PType>();

	//convert rectangle to box (include all cells intersecting with the rectangle)
	return get_minmax_component_z(VEC<VType>::box_from_rect_max(rectangle + VEC<VType>::rect.s));
}

//--------------------------------------------SPECIAL NUMERICAL PROPERTIES : VEC_VC_nprops.h

//Robin value is of the form alpha * (Tb - Ta). alpha and Ta are known from values set robin_nx, robin_px, ...
//Tb is quantity value at boundary. Here we will return Robin value for x, y, z axes for which any shift component is nonzero (otherwise zero for that component)
//e.g. if shift.x is non-zero then Tb value is obtained at rel_pos + (shift.x, 0, 0) using extrapolation from values at rel_pos and rel_pos - (shift.x, 0, 0) -> both these values should be inside the mesh, else return zero.
template <typename VType>
DBL3 VEC_VC<VType>::get_robin_value(const DBL3& rel_pos, const DBL3& shift)
{
	DBL3 robin_values;

	if (IsNZ(shift.x)) {

		VType val1 = (*this)[rel_pos];
		VType val2 = (*this)[rel_pos - DBL3(shift.x, 0, 0)];
		VType val_bnd = (1.5 * val1 - 0.5 * val2);

		if (shift.x > 0) robin_values.x = robin_px.i * (val_bnd - robin_px.j);
		else robin_values.x = robin_nx.i * (val_bnd - robin_nx.j);
	}

	if (IsNZ(shift.y)) {

		VType val1 = (*this)[rel_pos];
		VType val2 = (*this)[rel_pos - DBL3(0, shift.y, 0)];
		VType val_bnd = (1.5 * val1 - 0.5 * val2);

		if (shift.y > 0) robin_values.y = robin_py.i * (val_bnd - robin_py.j);
		else robin_values.y = robin_ny.i * (val_bnd - robin_ny.j);
	}

	if (IsNZ(shift.z)) {

		VType val1 = (*this)[rel_pos];
		VType val2 = (*this)[rel_pos - DBL3(0, 0, shift.z)];
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
DBL3 VEC_VC<VType>::get_shift_to_emptycell(int idx)
{
	DBL3 shift;

	//x
	if ((ngbrFlags[idx] & NF_BOTHX) != NF_BOTHX) {

		if (ngbrFlags[idx] & NF_NPX) shift.x = -VEC<VType>::h.x / 2;
		else if (ngbrFlags[idx] & NF_NNX) shift.x = +VEC<VType>::h.x / 2;
	}

	//y
	if ((ngbrFlags[idx] & NF_BOTHY) != NF_BOTHY) {

		if (ngbrFlags[idx] & NF_NPY) shift.y = -VEC<VType>::h.y / 2;
		else if (ngbrFlags[idx] & NF_NNY) shift.y = +VEC<VType>::h.y / 2;
	}

	//z
	if ((ngbrFlags[idx] & NF_BOTHZ) != NF_BOTHZ) {

		if (ngbrFlags[idx] & NF_NPZ) shift.z = -VEC<VType>::h.z / 2;
		else if (ngbrFlags[idx] & NF_NNZ) shift.z = +VEC<VType>::h.z / 2;
	}

	return shift;
}