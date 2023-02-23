#pragma once

#include "VEC.h"
//--------------------------------------------EXTRACT A LINE PROFILE

//extract profile in profile_storage temporary vector, returned through reference: extract starting at start in the direction end - step, with given step; use average to extract profile with given stencil, excluding zero points (assumed empty)
//all coordinates are relative positions
template <typename VType>
std::vector<VType>& VEC<VType>::extract_profile(DBL3 start, DBL3 end, double step, DBL3 stencil)
{
	if (step > 0) {

		size_t size = round((end - start).norm() / step) + 1;
		if (line_profile_storage.size() != size && !malloc_vector(line_profile_storage, size)) {

			line_profile_storage.clear();
			return line_profile_storage;
		}

		DBL3 meshDim = rect.size();

#pragma omp parallel for
		for (int idx = 0; idx < size; idx++) {

			//position wrapped-around
			DBL3 position = (start + (double)idx * step * (end - start).normalized()) % meshDim;

			VType value = average_nonempty(Rect(position - stencil / 2, position + stencil / 2));
			line_profile_storage[idx] = value;
		}
	}
	else line_profile_storage.clear();

	return line_profile_storage;
}

template std::vector<DBL2>& VEC<FLT3>::extract_profile_component_x(DBL3 start, DBL3 end, double step, DBL3 stencil);
template std::vector<DBL2>& VEC<DBL3>::extract_profile_component_x(DBL3 start, DBL3 end, double step, DBL3 stencil);

//as above but only component x (for VAL3 floating types only). Return xy data as line profile position and extracted component at that position.
template <typename VType>
std::vector<DBL2>& VEC<VType>::extract_profile_component_x(DBL3 start, DBL3 end, double step, DBL3 stencil)
{
	if (step > 0) {

		size_t size = round((end - start).norm() / step) + 1;
		if (line_profile_component.size() != size && !malloc_vector(line_profile_component, size)) {

			line_profile_component.clear();
			return line_profile_component;
		}

		DBL3 meshDim = rect.size();

#pragma omp parallel for
		for (int idx = 0; idx < size; idx++) {

			//position wrapped-around
			DBL3 position = (start + (double)idx * step * (end - start).normalized()) % meshDim;
			line_profile_component[idx].i = (double)idx * step;

			VType value = average_nonempty(Rect(position - stencil / 2, position + stencil / 2));
			line_profile_component[idx].j = value.x;
		}
	}
	else line_profile_component.clear();

	return line_profile_component;
}

template std::vector<DBL2>& VEC<FLT3>::extract_profile_component_y(DBL3 start, DBL3 end, double step, DBL3 stencil);
template std::vector<DBL2>& VEC<DBL3>::extract_profile_component_y(DBL3 start, DBL3 end, double step, DBL3 stencil);

//as above but only component y (for VAL3 floating types only). Return xy data as line profile position and extracted component at that position.
template <typename VType>
std::vector<DBL2>& VEC<VType>::extract_profile_component_y(DBL3 start, DBL3 end, double step, DBL3 stencil)
{
	if (step > 0) {

		size_t size = round((end - start).norm() / step) + 1;
		if (line_profile_component.size() != size && !malloc_vector(line_profile_component, size)) {

			line_profile_component.clear();
			return line_profile_component;
		}

		DBL3 meshDim = rect.size();

#pragma omp parallel for
		for (int idx = 0; idx < size; idx++) {

			//position wrapped-around
			DBL3 position = (start + (double)idx * step * (end - start).normalized()) % meshDim;
			line_profile_component[idx].i = (double)idx * step;

			VType value = average_nonempty(Rect(position - stencil / 2, position + stencil / 2));
			line_profile_component[idx].j = value.y;
		}
	}
	else line_profile_component.clear();

	return line_profile_component;
}

template std::vector<DBL2>& VEC<FLT3>::extract_profile_component_z(DBL3 start, DBL3 end, double step, DBL3 stencil);
template std::vector<DBL2>& VEC<DBL3>::extract_profile_component_z(DBL3 start, DBL3 end, double step, DBL3 stencil);

//as above but only component z (for VAL3 floating types only). Return xy data as line profile position and extracted component at that position.
template <typename VType>
std::vector<DBL2>& VEC<VType>::extract_profile_component_z(DBL3 start, DBL3 end, double step, DBL3 stencil)
{
	if (step > 0) {

		size_t size = round((end - start).norm() / step) + 1;
		if (line_profile_component.size() != size && !malloc_vector(line_profile_component, size)) {

			line_profile_component.clear();
			return line_profile_component;
		}

		DBL3 meshDim = rect.size();

#pragma omp parallel for
		for (int idx = 0; idx < size; idx++) {

			//position wrapped-around
			DBL3 position = (start + (double)idx * step * (end - start).normalized()) % meshDim;
			line_profile_component[idx].i = (double)idx * step;

			VType value = average_nonempty(Rect(position - stencil / 2, position + stencil / 2));
			line_profile_component[idx].j = value.z;
		}
	}
	else line_profile_component.clear();

	return line_profile_component;
}

template std::vector<DBL2>& VEC<FLT3>::extract_profile_component_max(DBL3 start, DBL3 end, double step, DBL3 stencil);
template std::vector<DBL2>& VEC<DBL3>::extract_profile_component_max(DBL3 start, DBL3 end, double step, DBL3 stencil);

//as above but only component which has largest value for the first point (after stencil averaging) (for VAL3 floating types only). Return xy data as line profile position and extracted component at that position.
template <typename VType>
std::vector<DBL2>& VEC<VType>::extract_profile_component_max(DBL3 start, DBL3 end, double step, DBL3 stencil)
{
	VType value = average_nonempty(Rect(start - stencil / 2, start + stencil / 2));

	if (abs(value.z) > abs(value.y) && abs(value.z) > abs(value.x)) {

		//Get Z component
		return extract_profile_component_z(start, end, step, stencil);
	}
	else if (abs(value.y) > abs(value.z) && abs(value.y) > abs(value.x)) {

		//Get Y component
		return extract_profile_component_y(start, end, step, stencil);
	}
	else {

		//Get X component
		return extract_profile_component_x(start, end, step, stencil);
	}
}