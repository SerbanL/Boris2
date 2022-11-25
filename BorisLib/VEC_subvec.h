#pragma once

#include "VEC.h"

//--------------------------------------------HELPER METHODS

//get a copy from this VEC, as a sub-VEC defined by box; same cellsize maintained; any transfer data not copied.
template <typename VType>
VEC<VType> VEC<VType>::subvec(Box box)
{
	Rect newvec_rect = Rect(box.s&h, box.e&h);
	//make appropriately sized new VEC
	VEC<VType> newvec(h, newvec_rect);

	//copy data to subvec
	for (int k = box.s.k; k < box.e.k; k++) {
#pragma omp parallel for
		for (int j = box.s.j; j < box.e.j; j++) {
			for (int i = box.s.i; i < box.e.i; i++) {

				int idx = i + j*n.x + k*n.x*n.y;

				int sidx = (i - box.s.i) + (j - box.s.j)*newvec.n.x + (k - box.s.k)*newvec.n.x*newvec.n.y;

				newvec[sidx] = quantity[idx];
			}
		}
	}

	//NOTE : the VEC which captures the return will retain same memory address for quantity (stored in a std::vector)
	return newvec;
}
