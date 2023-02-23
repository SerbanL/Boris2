#pragma once

#include "mcuVEC.h"

//------------------------------------------------------------------- AVERAGE

//average in a box (which should be contained in the cuVEC dimensions)
template <typename VType, typename MType>
VType mcuVEC<VType, MType>::average(cuBox box)
{
	VType aux = VType();
	size_t points = 0;

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//average on each device but keep reduction values in auxiliary values on each gpu
		mng(mGPU)->average(pn_d[mGPU].dim(), (box.IsNull() ? cuBox() : box.get_intersection(pbox_d[mGPU]) - pbox_d[mGPU].s), false);

		aux += get_gpu_value(mng(mGPU)->aux_value_ref());
		points += get_gpu_value(mng(mGPU)->aux_integer_ref());
	}

	if (points) return aux / points;
	else return VType();
}

//average over given rectangle (relative to this cuVEC's rect)
template <typename VType, typename MType>
VType mcuVEC<VType, MType>::average(cuRect rectangle)
{
	VType aux = VType();
	size_t points = 0;

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//average on each device but keep reduction values in auxiliary values on each gpu
		mng(mGPU)->average(pn_d[mGPU].dim(), (rectangle.IsNull() ? cuRect() : rectangle.get_intersection(prect_d[mGPU]) - prect_d[mGPU].s), false);

		aux += get_gpu_value(mng(mGPU)->aux_value_ref());
		points += get_gpu_value(mng(mGPU)->aux_integer_ref());
	}

	if (points) return aux / points;
	else return VType();
}

//------------------------------------------------------------------- AVERAGE NONEMPTY

//as above but exclude empty points from averaging
template <typename VType, typename MType>
VType mcuVEC<VType, MType>::average_nonempty(cuBox box)
{
	VType aux = VType();
	size_t points = 0;

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//average on each device but keep reduction values in auxiliary values on each gpu
		mng(mGPU)->average_nonempty(pn_d[mGPU].dim(), (box.IsNull() ? cuBox() : box.get_intersection(pbox_d[mGPU]) - pbox_d[mGPU].s), false);

		aux += get_gpu_value(mng(mGPU)->aux_value_ref());
		points += get_gpu_value(mng(mGPU)->aux_integer_ref());
	}

	if (points) return aux / points;
	else return VType();
}

template <typename VType, typename MType>
VType mcuVEC<VType, MType>::average_nonempty(cuRect rectangle)
{
	VType aux = VType();
	size_t points = 0;

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//average on each device but keep reduction values in auxiliary values on each gpu
		mng(mGPU)->average_nonempty(pn_d[mGPU].dim(), (rectangle.IsNull() ? cuRect() : rectangle.get_intersection(prect_d[mGPU]) - prect_d[mGPU].s), false);

		aux += get_gpu_value(mng(mGPU)->aux_value_ref());
		points += get_gpu_value(mng(mGPU)->aux_integer_ref());
	}

	if (points) return aux / points;
	else return VType();
}

//------------------------------------------------------------------- SUM NONEMPTY

template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
VType mcuVEC<VType, MType>::sum_nonempty(cuBox box)
{
	VType aux = VType();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//average on each device but keep reduction values in auxiliary values on each gpu
		mng(mGPU)->sum_nonempty(pn_d[mGPU].dim(), (box.IsNull() ? cuBox() : box.get_intersection(pbox_d[mGPU]) - pbox_d[mGPU].s), false);

		aux += get_gpu_value(mng(mGPU)->aux_value_ref());
	}

	return aux;
}

template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
VType mcuVEC<VType, MType>::sum_nonempty(cuRect rectangle)
{
	VType aux = VType();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//average on each device but keep reduction values in auxiliary values on each gpu
		mng(mGPU)->sum_nonempty(pn_d[mGPU].dim(), (rectangle.IsNull() ? cuRect() : rectangle.get_intersection(prect_d[mGPU]) - prect_d[mGPU].s), false);

		aux += get_gpu_value(mng(mGPU)->aux_value_ref());
	}

	return aux;
}