#pragma once

#include "mcuArray.h"

//------------------------------------------- SET VALUE : mcuArray_aux.h

//set all entries to given value, for all devices
template <typename VType>
void mcu_arr<VType>::set(VType value)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		pcuarr[mGPU]->set(value);
	}
}

//set all entries to given value, for indexed device
template <typename VType>
void mcu_arr<VType>::set(int idx, VType value)
{
	mGPU.select_device(idx);
	pcuarr[idx]->set(value);
}

//set single value from cpu memory at given index for all devices configured
template <typename VType>
void mcu_arr<VType>::setvalue(int index, VType value)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		pcuarr[mGPU]->setvalue(index, value);
	}
}

//------------------------------------------- GET SIZE : mcuArray_aux.h

//get size of array on indexed device
template <typename VType>
size_t mcu_arr<VType>::size(int idx)
{
	//setting cuda device not needed here, just read size value, which is stored in cpu memory
	return pcuarr[idx]->size();
}

//get total size (i.e. sum of all allocated device sizes)
//get size of array on indexed device
template <typename VType>
size_t mcu_arr<VType>::size(void)
{
	size_t total_size = 0;

	for (int idx = 0; idx < mGPU.get_num_devices(); idx++) {

		total_size += pcuarr[idx]->size();
	}
	
	return total_size;
}