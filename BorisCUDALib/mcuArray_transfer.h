#pragma once

#include "mcuArray.h"


//------------------------------------------- GPU <-> CPU TRANSFER : mcuArray_transfer.h

//copy values from a std::vector into gpu memory. Doesn't set size, but copies up to currently allocated size. Copy to all configured devices.
template <typename VType>
template <typename Type>
void mcu_arr<VType>::copy_from_vector(std::vector<Type>& cpuvec)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		pcuarr[mGPU]->copy_from_vector(cpuvec);
	}
}

//copy values to a std::vector into cpu memory. Doesn't set size, but copies up to currently allocated size.
//Copy by putting all data from devices head to head in order of devices.
template <typename VType>
template <typename Type>
void mcu_arr<VType>::copy_to_vector(std::vector<Type>& cpuvec)
{
	size_t size_accum = 0;

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//copy to vector starting at the corect offset
		pcuarr[mGPU]->copy_to_vector(cpuvec, size_accum);

		//next
		size_accum += pcuarr[mGPU]->size();
	}
}