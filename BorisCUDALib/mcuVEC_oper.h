#pragma once

#include "mcuVEC.h"

//------------------------------------------------------------------- SETBOX

//set value in box
template <typename VType, typename MType>
void mcuVEC<VType, MType>::setbox(cuBox box, VType value)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->setbox(box.get_intersection(pbox_d[mGPU]) - pbox_d[mGPU].s, value);
	}
}

//------------------------------------------------------------------- SETRECT

//set value in rectangle (i.e. in cells intersecting the rectangle), where the rectangle is relative to this cuVEC's rectangle.
template <typename VType, typename MType>
void mcuVEC<VType, MType>::setrect(const cuRect& rectangle, VType value)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->setrect(rectangle.get_intersection(prect_d[mGPU]) - prect_d[mGPU].s, value);
	}
}

//------------------------------------------------------------------- DELRECT

//delete rectangle, where the rectangle is relative to this VEC's rectangle, by setting empty cell values - all cells become empty cells. Applicable to cuVEC_VC only
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::delrect(cuRect rectangle)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->delrect(rectangle.get_intersection(prect_d[mGPU]) - prect_d[mGPU].s);
	}
}

//------------------------------------------------------------------- SET

//set value in all cells
template <typename VType, typename MType>
void mcuVEC<VType, MType>::set(VType value)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->set(pn_d[mGPU].dim(), value);
	}
}

//------------------------------------------------------------------- RENORMALIZE

//re-normalize all non-zero values to have the new magnitude (multiply by new_norm and divide by current magnitude)
template <typename VType, typename MType>
template <typename PType>
void mcuVEC<VType, MType>::renormalize(PType new_norm)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->renormalize(pn_d[mGPU].dim(), new_norm);
	}
}

//------------------------------------------------------------------- BITMAP MASK

//mask values in cells using bitmap image : white -> empty cells. black -> keep values. Apply mask up to given z depth number of cells depending on grayscale value (zDepth, all if 0). Applicable to cuVEC_VC only
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
bool mcuVEC<VType, MType>::apply_bitmap_mask(std::vector<unsigned char>& bitmap, int zDepth)
{
	bool success = true;

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		cuBox box = pbox_d[mGPU];
		std::vector<unsigned char> sub_bitmap = subvec(bitmap, box.s.i, box.s.j, 0, box.e.i, box.e.j, 1, n.i, n.j, 1);

		success &= mng(mGPU)->apply_bitmap_mask(sub_bitmap, zDepth);
	}

	return success;
}

//------------------------------------------------------------------- SETNONEMPTY

//exactly the same as assign value - do not use assign as it is slow (sets flags). Applicable to cuVEC_VC only
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::setnonempty(VType value)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->setnonempty(value);
	}
}

//------------------------------------------------------------------- SET RECT NONEMPTY

//set value in non-empty cells only in given rectangle (relative coordinates). Applicable to cuVEC_VC only
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::setrectnonempty(const cuRect& rectangle, VType value)
{
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->setrectnonempty(rectangle.get_intersection(prect_d[mGPU]) - prect_d[mGPU].s, value);
	}
}

//------------------------------------------------------------------- SHIFT

//shift all the values in this cuVEC by the given delta (units same as cuVEC<VType>::h). Shift values in given shift_rect (absolute coordinates). Applicable to cuVEC_VC only
template <typename VType, typename MType>
template <typename VType_, typename MType_, std::enable_if_t<std::is_same<MType_, cuVEC_VC<VType_>>::value>*>
void mcuVEC<VType, MType>::shift_x(cuBReal delta, cuRect shift_rect)
{
	//before attempting shifts, must exchange halos as these are used for shifting values at boundaries for each gpu
	//for shift_x we only need to do this if sub-rectangles are arranged along x
	if (halo_flag == NF2_HALOX) exchange_halos();

	//now we can do the shifts
	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		mng(mGPU)->shift_x(pn_d[mGPU].dim(), delta, shift_rect);
	}
}