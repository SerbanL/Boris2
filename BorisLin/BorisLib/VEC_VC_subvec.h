#pragma once

#include "VEC_VC.h"

//--------------------------------------------HELPER METHODS

//get a copy from this VEC_VC, as a sub-VEC_VC defined by box; same cellsize maintained; any transfer data not copied.
//dirichlet conditions not copied
template <typename VType>
VEC_VC<VType> VEC_VC<VType>::subvec(Box box)
{
	Rect newvec_rect = Rect(box.s&VEC<VType>::h, box.e&VEC<VType>::h);
	//make appropriately sized new VEC_VC
	VEC_VC<VType> newvec(VEC<VType>::h, newvec_rect);

	//copy data to subvec (values and flags)
	for (int k = box.s.k; k < box.e.k; k++) {
#pragma omp parallel for
		for (int j = box.s.j; j < box.e.j; j++) {
			for (int i = box.s.i; i < box.e.i; i++) {

				int idx = i + j * VEC<VType>::n.x + k * VEC<VType>::n.x*VEC<VType>::n.y;

				int sidx = (i - box.s.i) + (j - box.s.j)*newvec.n.x + (k - box.s.k)*newvec.n.x*newvec.n.y;

				newvec[sidx] = VEC<VType>::quantity[idx];
				newvec.ngbrFlags_ref()[sidx] = ngbrFlags[idx];
				//extended flags not copied here, but set after
			}
		}
	}

	//dirichlet vectors not copied - these should be set externally if needed

	//copy robin values
	newvec.robin_px = robin_px;
	newvec.robin_nx = robin_nx;
	newvec.robin_py = robin_py;
	newvec.robin_ny = robin_ny;
	newvec.robin_pz = robin_pz;
	newvec.robin_nz = robin_nz;
	newvec.robin_v = robin_v;

	//copy shift debt
	newvec.shift_debt = shift_debt;

	//copy pbc parameters
	newvec.pbc_x = pbc_x;
	newvec.pbc_y = pbc_y;
	newvec.pbc_z = pbc_z;

	//now make sure flags in new VEC_VC are self consistent
	newvec.set_ngbrFlags();

	//NOTE : the VEC which captures the return will retain same memory address for quantity (stored in a std::vector)
	return newvec;
}
