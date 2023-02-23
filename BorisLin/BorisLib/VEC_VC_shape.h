#pragma once

#include "VEC_VC.h"

//--------------------------------------------MULTIPLE ENTRIES SETTERS - SHAPE CHANGERS

//set value in box (i.e. in cells entirely included in box) - all cells become non-empty cells irrespective of value set
template <typename VType>
void VEC_VC<VType>::setbox(const Box& box, VType value)
{
#pragma omp parallel for
	for (int j = (box.s.y >= 0 ? box.s.y : 0); j < (box.e.y <= VEC<VType>::n.y ? box.e.y : VEC<VType>::n.y); j++) {
		for (int k = (box.s.z >= 0 ? box.s.z : 0); k < (box.e.z <= VEC<VType>::n.z ? box.e.z : VEC<VType>::n.z); k++) {
			for (int i = (box.s.x >= 0 ? box.s.x : 0); i < (box.e.x <= VEC<VType>::n.x ? box.e.x : VEC<VType>::n.x); i++) {

				int idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

				VEC<VType>::quantity[idx] = value;
				mark_not_empty(idx);
			}
		}
	}

	set_ngbrFlags(*this);
}

//set value in rectangle (i.e. in cells intersecting the rectangle), where the rectangle is relative to this VEC's rectangle - all cells become non-empty cells irrespective of value set
template <typename VType>
void VEC_VC<VType>::setrect(const Rect& rectangle, VType value)
{
	if (!VEC<VType>::rect.intersects(rectangle + VEC<VType>::rect.s)) return;

	Box cells_box = VEC<VType>::box_from_rect_max(rectangle + VEC<VType>::rect.s);

	setbox(cells_box, value);
}

//delete rectangle, where the rectangle is relative to this VEC's rectangle, by setting empty cell values - all cells become empty cells
template <typename VType>
void VEC_VC<VType>::delrect(const Rect& rectangle)
{
	Box box = VEC<VType>::box_from_rect_max(rectangle + VEC<VType>::rect.s);

#pragma omp parallel for
	for (int j = (box.s.y >= 0 ? box.s.y : 0); j < (box.e.y <= VEC<VType>::n.y ? box.e.y : VEC<VType>::n.y); j++) {
		for (int k = (box.s.z >= 0 ? box.s.z : 0); k < (box.e.z <= VEC<VType>::n.z ? box.e.z : VEC<VType>::n.z); k++) {
			for (int i = (box.s.x >= 0 ? box.s.x : 0); i < (box.e.x <= VEC<VType>::n.x ? box.e.x : VEC<VType>::n.x); i++) {

				int idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

				mark_empty(idx);
			}
		}
	}

	set_ngbrFlags(*this);
}

template <typename VType>
bool VEC_VC<VType>::apply_bitmap_mask(const std::vector<unsigned char>& bitmap, double zDepth)
{
	//bitmap must have the right size (i.e. have VEC<VType>::n.x * VEC<VType>::n.y pixels, remembering each pixel has 4 bytes as B-G-R-A)
	if (bitmap.size() != VEC<VType>::n.x*VEC<VType>::n.y * 4) return false;

#pragma omp parallel for
	for (int j = 0; j < VEC<VType>::n.y; j++) {
		for (int i = 0; i < VEC<VType>::n.x; i++) {

			//In the file the 0, 0 point is in the left-top corner, but in the mesh the 0, 0 point is in the left-bottom corner
			int indexBitMap = i + (VEC<VType>::n.y - 1 - j)*VEC<VType>::n.x;
			int indexMesh = i + j * VEC<VType>::n.x;

			double B = (double)bitmap[indexBitMap * 4] / 255;
			//double G = (double)bitmap[indexBitMap * 4 + 1] / 255;
			//double R = (double)bitmap[indexBitMap * 4 + 2] / 255;
			
			int depthCut = 0;
			
			int zCellsDepth = round(zDepth / VEC<VType>::h.z);
			
			if (zCellsDepth) {

				//divide grayscale in zDepth equal intervals of length intervalLength
				double intervalLength = (double)1.0 / (double)abs(zCellsDepth);

				//find out in which interval the grayscale pixel lies (assume B=G=R)
				depthCut = (int)ceil(B / intervalLength);
			}
			else {

				if (B > 0.5) depthCut = VEC<VType>::n.z;
			}

			if (zCellsDepth > 0) {
				//void magnetization to given depth
				for (int k = VEC<VType>::n.z - 1; k > (int)VEC<VType>::n.z - 1 - depthCut; k--) {

					if (k >= 0) {

						mark_empty(indexMesh + k * VEC<VType>::n.x*VEC<VType>::n.y);
					}
				}
			}
			
			if (zCellsDepth <= 0) {

				//void magnetization to given height
				for (int k = 0; k < depthCut; k++) {

					if (k < VEC<VType>::n.z) {

						mark_empty(indexMesh + k * VEC<VType>::n.x*VEC<VType>::n.y);
					}
				}
			}
		}
	}

	set_ngbrFlags(*this);

	return true;
}