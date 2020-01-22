#pragma once

#include "VEC.h"
//--------------------------------------------EXTRACT A LINE PROFILE

//extract profile to a vector : extract size points starting at (start + step * 0.5) in the direction step; use weighted average to extract profile with stencil given by h
//e.g. if you have a start and end point with given step, then setting size = |end - start| / |step| means the profile must be extracted between (start + 0.5*step) and (end - 0.5*step). e.g.: |.|.|.|.|
template <typename VType>
void VEC<VType>::extract_profile(size_t size, std::vector<VType>& profile, DBL3 start, DBL3 step)
{
#pragma omp parallel for
	for (int idx = 0; idx < size; idx++) {

		DBL3 position = start + ((double)idx + 0.5) * step;

		VType value = weighted_average(position, h);
		profile[idx] = value;
	}
}

//these specifically apply for VType == cuReal3, allowing extraction of the x, y, z components separately
template <typename VType>
template <typename PType>
void VEC<VType>::extract_profile_component_x(size_t size, std::vector<PType>& profile, DBL3 start, DBL3 step)
{
#pragma omp parallel for
	for (int idx = 0; idx < size; idx++) {

		DBL3 position = start + ((double)idx + 0.5) * step;

		VType value = weighted_average(position, h);
		profile[idx] = value.x;
	}
}

template <typename VType>
template <typename PType>
void VEC<VType>::extract_profile_component_y(size_t size, std::vector<PType>& profile, DBL3 start, DBL3 step)
{
#pragma omp parallel for
	for (int idx = 0; idx < size; idx++) {

		DBL3 position = start + ((double)idx + 0.5) * step;

		VType value = weighted_average(position, h);
		profile[idx] = value.y;
	}
}

template <typename VType>
template <typename PType>
void VEC<VType>::extract_profile_component_z(size_t size, std::vector<PType>& profile, DBL3 start, DBL3 step)
{
#pragma omp parallel for
	for (int idx = 0; idx < size; idx++) {

		DBL3 position = start + ((double)idx + 0.5) * step;

		VType value = weighted_average(position, h);
		profile[idx] = value.z;
	}
}