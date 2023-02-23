#pragma once

#include "mcuVEC.h"

//linear : use interpolation to set values in this VEC based on projected distance between position1 and position2 and given fixed end values.
template <typename VType, typename MType>
bool mcuVEC<VType, MType>::generate_linear(cuReal3 new_h, cuRect new_rect, cuRect contact1, VType value1, cuRect contact2, VType value2)
{
	if (!resize(new_h, new_rect)) return false;

	set_linear(contact1, value1, contact2, value2);

	return true;
}

//similar to generate_linear except new dimensions not set
//also allow 'degeneracy' : multiple linear generators may be superimposed with use of degeneracy : degeneracy.first is index, degeneracy.second is number of geenerators
//if using degeneracy make sure these are called in order, or at least index 0 goes first
template <typename VType, typename MType>
void mcuVEC<VType, MType>::set_linear(cuRect contact1, VType value1, cuRect contact2, VType value2, cuReal2 degeneracy)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->set_linear(contact1, value1, contact2, value2, degeneracy);
	}
}