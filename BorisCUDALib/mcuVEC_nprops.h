#pragma once

#include "mcuVEC.h"

//Find min and max values. rectangles are relative to this VEC.

//------------------------------------------------------------------- MIN-MAX

template <typename VType, typename MType>
template <typename PType>
cuVAL2<PType> mcuVEC<VType, MType>::get_minmax(cuBox box)
{
	cuVAL2<PType> minmax = cuVAL2<PType>(PType(1e38), PType(-1e38));

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		cuVAL2<PType> minmax_ = mng(mGPU)->get_minmax(pn_d[mGPU].dim(), (box.IsNull() ? cuBox() : box.get_intersection(pbox_d[mGPU]) - pbox_d[mGPU].s));

		if (minmax_.i < minmax.i) minmax.i = minmax_.i;
		if (minmax_.j > minmax.j) minmax.j = minmax_.j;
	}

	return minmax;
}

template <typename VType, typename MType>
template <typename PType>
cuVAL2<PType> mcuVEC<VType, MType>::get_minmax(cuRect rectangle)
{
	cuVAL2<PType> minmax = cuVAL2<PType>(PType(1e38), PType(-1e38));

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		cuVAL2<PType> minmax_ = mng(mGPU)->get_minmax(pn_d[mGPU].dim(), (rectangle.IsNull() ? cuRect() : rectangle.get_intersection(prect_d[mGPU]) - prect_d[mGPU].s));

		if (minmax_.i < minmax.i) minmax.i = minmax_.i;
		if (minmax_.j > minmax.j) minmax.j = minmax_.j;
	}

	return minmax;
}

//------------------------------------------------------------------- MIN-MAX COMPONENT X

template <typename VType, typename MType>
template <typename PType>
cuVAL2<PType> mcuVEC<VType, MType>::get_minmax_component_x(cuBox box)
{
	cuVAL2<PType> minmax = cuVAL2<PType>(PType(1e38), PType(-1e38));

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		cuVAL2<PType> minmax_ = mng(mGPU)->get_minmax_component_x(pn_d[mGPU].dim(), (box.IsNull() ? cuBox() : box.get_intersection(pbox_d[mGPU]) - pbox_d[mGPU].s));

		if (minmax_.i < minmax.i) minmax.i = minmax_.i;
		if (minmax_.j > minmax.j) minmax.j = minmax_.j;
	}

	return minmax;
}

template <typename VType, typename MType>
template <typename PType>
cuVAL2<PType> mcuVEC<VType, MType>::get_minmax_component_x(cuRect rectangle)
{
	cuVAL2<PType> minmax = cuVAL2<PType>(PType(1e38), PType(-1e38));

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		cuVAL2<PType> minmax_ = mng(mGPU)->get_minmax_component_x(pn_d[mGPU].dim(), (rectangle.IsNull() ? cuRect() : rectangle.get_intersection(prect_d[mGPU]) - prect_d[mGPU].s));

		if (minmax_.i < minmax.i) minmax.i = minmax_.i;
		if (minmax_.j > minmax.j) minmax.j = minmax_.j;
	}

	return minmax;
}

//------------------------------------------------------------------- MIN-MAX COMPONENT Y

template <typename VType, typename MType>
template <typename PType>
cuVAL2<PType> mcuVEC<VType, MType>::get_minmax_component_y(cuBox box)
{
	cuVAL2<PType> minmax = cuVAL2<PType>(PType(1e38), PType(-1e38));

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		cuVAL2<PType> minmax_ = mng(mGPU)->get_minmax_component_y(pn_d[mGPU].dim(), (box.IsNull() ? cuBox() : box.get_intersection(pbox_d[mGPU]) - pbox_d[mGPU].s));

		if (minmax_.i < minmax.i) minmax.i = minmax_.i;
		if (minmax_.j > minmax.j) minmax.j = minmax_.j;
	}

	return minmax;
}

template <typename VType, typename MType>
template <typename PType>
cuVAL2<PType> mcuVEC<VType, MType>::get_minmax_component_y(cuRect rectangle)
{
	cuVAL2<PType> minmax = cuVAL2<PType>(PType(1e38), PType(-1e38));

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		cuVAL2<PType> minmax_ = mng(mGPU)->get_minmax_component_y(pn_d[mGPU].dim(), (rectangle.IsNull() ? cuRect() : rectangle.get_intersection(prect_d[mGPU]) - prect_d[mGPU].s));

		if (minmax_.i < minmax.i) minmax.i = minmax_.i;
		if (minmax_.j > minmax.j) minmax.j = minmax_.j;
	}

	return minmax;
}

//------------------------------------------------------------------- MIN-MAX COMPONENT Z

template <typename VType, typename MType>
template <typename PType>
cuVAL2<PType> mcuVEC<VType, MType>::get_minmax_component_z(cuBox box)
{
	cuVAL2<PType> minmax = cuVAL2<PType>(PType(1e38), PType(-1e38));

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		cuVAL2<PType> minmax_ = mng(mGPU)->get_minmax_component_z(pn_d[mGPU].dim(), (box.IsNull() ? cuBox() : box.get_intersection(pbox_d[mGPU]) - pbox_d[mGPU].s));

		if (minmax_.i < minmax.i) minmax.i = minmax_.i;
		if (minmax_.j > minmax.j) minmax.j = minmax_.j;
	}

	return minmax;
}

template <typename VType, typename MType>
template <typename PType>
cuVAL2<PType> mcuVEC<VType, MType>::get_minmax_component_z(cuRect rectangle)
{
	cuVAL2<PType> minmax = cuVAL2<PType>(PType(1e38), PType(-1e38));

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		cuVAL2<PType> minmax_ = mng(mGPU)->get_minmax_component_z(pn_d[mGPU].dim(), (rectangle.IsNull() ? cuRect() : rectangle.get_intersection(prect_d[mGPU]) - prect_d[mGPU].s));

		if (minmax_.i < minmax.i) minmax.i = minmax_.i;
		if (minmax_.j > minmax.j) minmax.j = minmax_.j;
	}

	return minmax;
}