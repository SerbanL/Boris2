#pragma once

#include "mcuArray.h"

//------------------------------------------- RESIZING : mcuArray_sizing.h

//resize all devices to given size
template <typename VType>
bool mcu_arr<VType>::resize(size_t size)
{
	bool success = true;

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		success &= pcuarr[mGPU]->resize(size);
	}

	return success;
}

//resize indexed device to given size
template <typename VType>
bool mcu_arr<VType>::resize(int idx, size_t size)
{
	mGPU.select_device(idx);
	return pcuarr[idx]->resize(size);
}

//clear for all devices
template <typename VType>
void mcu_arr<VType>::clear(void)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		pcuarr[mGPU]->clear();
	}
}

//clear for indexed device only
template <typename VType>
void mcu_arr<VType>::clear(int idx)
{
	mGPU.select_device(idx);
	pcuarr[idx]->clear();
}

//------------------------------------------- STORE ENTRIES : mcuArray_sizing.h

//new_entry is a pointer in cpu memory to an object in gpu memory. add it to indexed cu_arr
template <typename VType>
void mcu_arr<VType>::push_back(int idx, VType*& new_entry)
{
	mGPU.select_device(idx);
	pcuarr[idx]->push_back(new_entry);
}