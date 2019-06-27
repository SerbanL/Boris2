#pragma once

#include "VEC_VC.h"

//--------------------------------------------MULTIPLE ENTRIES SETTERS - SHAPE CHANGERS

//set value in box (i.e. in cells entirely included in box) - all cells become non-empty cells irrespective of value set
template <typename VType>
void VEC_VC<VType>::setbox(const Box& box, VType value)
{
#pragma omp parallel for
	for (int j = (box.s.y >= 0 ? box.s.y : 0); j < (box.e.y <= n.y ? box.e.y : n.y); j++) {
		for (int k = (box.s.z >= 0 ? box.s.z : 0); k < (box.e.z <= n.z ? box.e.z : n.z); k++) {
			for (int i = (box.s.x >= 0 ? box.s.x : 0); i < (box.e.x <= n.x ? box.e.x : n.x); i++) {

				int idx = i + j * n.x + k * n.x*n.y;

				quantity[idx] = value;
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
	if (!rect.intersects(rectangle + rect.s)) return;

	Box cells_box = box_from_rect_max(rectangle + rect.s);

	setbox(cells_box, value);
}

//delete rectangle, where the rectangle is relative to this VEC's rectangle, by setting empty cell values - all cells become empty cells
template <typename VType>
void VEC_VC<VType>::delrect(const Rect& rectangle)
{
	Box box = box_from_rect_max(rectangle + rect.s);

#pragma omp parallel for
	for (int j = (box.s.y >= 0 ? box.s.y : 0); j < (box.e.y <= n.y ? box.e.y : n.y); j++) {
		for (int k = (box.s.z >= 0 ? box.s.z : 0); k < (box.e.z <= n.z ? box.e.z : n.z); k++) {
			for (int i = (box.s.x >= 0 ? box.s.x : 0); i < (box.e.x <= n.x ? box.e.x : n.x); i++) {

				int idx = i + j * n.x + k * n.x*n.y;

				mark_empty(idx);
			}
		}
	}

	set_ngbrFlags(*this);
}

template <typename VType>
bool VEC_VC<VType>::apply_bitmap_mask(std::vector<BYTE>& bitmap, double zDepth)
{
	//bitmap must have the right size (i.e. have n.x * n.y pixels, remembering each pixel has 4 bytes as B-G-R-A)
	if (bitmap.size() != n.x*n.y * 4) return false;

#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {
		for (int i = 0; i < n.x; i++) {

			//In the file the 0, 0 point is in the left-top corner, but in the mesh the 0, 0 point is in the left-bottom corner
			int indexBitMap = i + (n.y - 1 - j)*n.x;
			int indexMesh = i + j * n.x;

			double B = (double)bitmap[indexBitMap * 4] / 255;
			//double G = (double)bitmap[indexBitMap * 4 + 1] / 255;
			//double R = (double)bitmap[indexBitMap * 4 + 2] / 255;
			
			int depthCut = 0;
			
			int zCellsDepth = round(zDepth / h.z);
			
			if (zCellsDepth) {

				//divide grayscale in zDepth equal intervals of length intervalLength
				double intervalLength = (double)1.0 / (double)abs(zCellsDepth);

				//find out in which interval the grayscale pixel lies (assume B=G=R)
				depthCut = (int)ceil(B / intervalLength);
			}
			else {

				if (B > 0.5) depthCut = n.z;
			}

			if (zCellsDepth > 0) {
				//void magnetization to given depth
				for (int k = n.z - 1; k > (int)n.z - 1 - depthCut; k--) {

					if (k >= 0) {

						mark_empty(indexMesh + k * n.x*n.y);
					}
				}
			}
			
			if (zCellsDepth <= 0) {

				//void magnetization to given height
				for (int k = 0; k < depthCut; k++) {

					if (k < n.z) {

						mark_empty(indexMesh + k * n.x*n.y);
					}
				}
			}
		}
	}

	set_ngbrFlags(*this);

	return true;
}