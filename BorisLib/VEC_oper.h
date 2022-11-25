#pragma once

#include "VEC.h"

//--------------------------------------------MULTIPLE ENTRIES SETTERS

template <typename VType>
void VEC<VType>::set(VType value)
{
#pragma omp parallel for
	for (int idx = 0; idx < quantity.size(); idx++) {

		quantity[idx] = value;
	}
}

//set value in box (i.e. in cells entirely included in box)
template <typename VType>
void VEC<VType>::setbox(const Box& box, VType value)
{
#pragma omp parallel for
	for (int j = (box.s.y >= 0 ? box.s.y : 0); j < (box.e.y <= n.y ? box.e.y : n.y); j++) {
		for (int k = (box.s.z >= 0 ? box.s.z : 0); k < (box.e.z <= n.z ? box.e.z : n.z); k++) {
			for (int i = (box.s.x >= 0 ? box.s.x : 0); i < (box.e.x <= n.x ? box.e.x : n.x); i++) {

				quantity[i + j * n.x + k * n.x*n.y] = value;
			}
		}
	}
}

//set value in rectangle (i.e. in cells intersecting the rectangle), where the rectangle is relative to this VEC's rectangle.
template <typename VType>
void VEC<VType>::setrect(const Rect& rectangle, VType value)
{
	if (!rect.intersects(rectangle + rect.s)) return;

	Box cells_box = box_from_rect_max(rectangle + rect.s);

	setbox(cells_box, value);
}

template <typename VType>
template <typename PType>
void VEC<VType>::renormalize(PType new_norm)
{
#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		PType curr_norm = GetMagnitude(quantity[idx]);

		if (IsNZ(curr_norm)) quantity[idx] *= new_norm / curr_norm;
	}
}

//copy values from copy_this but keep current dimensions - if necessary map values from copy_this to local dimensions
template <typename VType>
void VEC<VType>::copy_values(const VEC<VType>& copy_this, Rect dstRect, Rect srcRect)
{
	if (dstRect.IsNull()) dstRect = rect - rect.s;
	if (srcRect.IsNull()) srcRect = copy_this.rect - copy_this.rect.s;

	Box cells_box_dst = box_from_rect_max(dstRect + rect.s);
	Box cells_box_src = copy_this.box_from_rect_max(srcRect + copy_this.rect.s);

	SZ3 dst_n = cells_box_dst.size();
	SZ3 src_n = cells_box_src.size();

	DBL3 sourceIdx = (DBL3)src_n / dst_n;
	
#pragma omp parallel for
	for (int j = 0; j < dst_n.j; j++) {
		for (int k = 0; k < dst_n.k; k++) {
			for (int i = 0; i < dst_n.i; i++) {

				int idx_box_dst = i + j * dst_n.x + k * dst_n.x*dst_n.y;

				int _x = (int)floor(i * sourceIdx.x);
				int _y = (int)floor(j * sourceIdx.y);
				int _z = (int)floor(k * sourceIdx.z);

				int idx_out = (i + cells_box_dst.s.i) + (j + cells_box_dst.s.j) * n.x + (k + cells_box_dst.s.k) * n.x*n.y;
				int idx_in = (_x + cells_box_src.s.i) + (_y + cells_box_src.s.j) * copy_this.n.x + (_z + cells_box_src.s.k) * (copy_this.n.x*copy_this.n.y);

				if (idx_out < n.dim() && idx_in < copy_this.n.dim()) quantity[idx_out] = copy_this.quantity[idx_in];
			}
		}
	}
}
