#pragma once

#include "mcuVEC.h"

//------------------------------------------------------------------- ADD

//add to this vec the values in add_this : must have same size : size
template <typename VType, typename MType>
void mcuVEC<VType, MType>::add_values(mcu_obj<MType, mcuVEC<VType, MType>>& add_this)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->add_values(pn_d[mGPU].dim(), add_this.get_deviceobject(mGPU));
	}
}

//------------------------------------------------------------------- SUB

//subtract from this vec the values in sub_this : must have same size : size
template <typename VType, typename MType>
void mcuVEC<VType, MType>::sub_values(mcu_obj<MType, mcuVEC<VType, MType>>& sub_this)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->sub_values(pn_d[mGPU].dim(), sub_this.get_deviceobject(mGPU));
	}
}

//------------------------------------------------------------------- SCALE

//scale all stored values by the given constant. Applicable to cuVEC_VC only
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::scale_values(cuBReal constant)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->scale_values(pn_d[mGPU].dim(), constant);
	}
}