#pragma once

#include "VEC_VC.h"

//-------------------------------- CURL OPERATOR

//curl operator. Use Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
//can only be applied if VType is a VAL3
template <typename VType>
VType VEC_VC<VType>::curl_neu(int idx) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {
		
		curl.y -= (VEC<VType>::quantity[idx + 1].z - VEC<VType>::quantity[idx - 1].z) / (2 * VEC<VType>::h.x);
		curl.z += (VEC<VType>::quantity[idx + 1].y - VEC<VType>::quantity[idx - 1].y) / (2 * VEC<VType>::h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (VEC<VType>::quantity[idx + 1].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.x;
			curl.z += (VEC<VType>::quantity[idx + 1].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - 1].z) / VEC<VType>::h.x;
			curl.z += (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - 1].y) / VEC<VType>::h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (VEC<VType>::quantity[idx + 1].z - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].z) / (2 * VEC<VType>::h.x);
				curl.z += (VEC<VType>::quantity[idx + 1].y - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].y) / (2 * VEC<VType>::h.x);
			}
			else {

				curl.y -= 0.5 * (VEC<VType>::quantity[idx + 1].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.x;
				curl.z += 0.5 * (VEC<VType>::quantity[idx + 1].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.x;
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].z - VEC<VType>::quantity[idx - 1].z) / (2 * VEC<VType>::h.x);
				curl.z += (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].y - VEC<VType>::quantity[idx - 1].y) / (2 * VEC<VType>::h.x);
			}
			else {

				curl.y -= 0.5 * (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - 1].z) / VEC<VType>::h.x;
				curl.z += 0.5 * (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - 1].y) / VEC<VType>::h.x;
			}
		}
	}

	//y direction differentials
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
		curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.y;
			curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / VEC<VType>::h.y;
			curl.z -= (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / VEC<VType>::h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
				curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
			}
			else {

				curl.x += 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.y;
				curl.z -= 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.y;
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
				curl.z -= (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
			}
			else {

				curl.x += 0.5 * (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / VEC<VType>::h.y;
				curl.z -= 0.5 * (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / VEC<VType>::h.y;
			}
		}
	}
	
	//z direction differentials
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
		curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.z;
			curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / VEC<VType>::h.z;
			curl.y += (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / VEC<VType>::h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
				curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
			}
			else {

				curl.x -= 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.z;
				curl.y += 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.z;
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
				curl.y += (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
			}
			else {

				curl.x -= 0.5 * (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / VEC<VType>::h.z;
				curl.y += 0.5 * (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / VEC<VType>::h.z;
			}
		}
	}

	return curl;
}

//curl operator. Use non-homogeneous Neumann boundary conditions.
//Can be used at composite media boundaries where sided differentials will be used instead.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
//can only be applied if VType is a VAL3
template <typename VType>
VType VEC_VC<VType>::curl_nneu(int idx, const VAL3<VType>& bdiff) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		curl.y -= (VEC<VType>::quantity[idx + 1].z - VEC<VType>::quantity[idx - 1].z) / (2 * VEC<VType>::h.x);
		curl.z += (VEC<VType>::quantity[idx + 1].y - VEC<VType>::quantity[idx - 1].y) / (2 * VEC<VType>::h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (VEC<VType>::quantity[idx + 1].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.x;
			curl.z += (VEC<VType>::quantity[idx + 1].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - 1].z) / VEC<VType>::h.x;
			curl.z += (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - 1].y) / VEC<VType>::h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (VEC<VType>::quantity[idx + 1].z - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].z) / (2 * VEC<VType>::h.x);
				curl.z += (VEC<VType>::quantity[idx + 1].y - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].y) / (2 * VEC<VType>::h.x);
			}
			else {

				curl.y -= 0.5 * ((VEC<VType>::quantity[idx + 1].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.x + bdiff.x.z);
				curl.z += 0.5 * ((VEC<VType>::quantity[idx + 1].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.x + bdiff.x.y);
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].z - VEC<VType>::quantity[idx - 1].z) / (2 * VEC<VType>::h.x);
				curl.z += (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].y - VEC<VType>::quantity[idx - 1].y) / (2 * VEC<VType>::h.x);
			}
			else {

				curl.y -= 0.5 * ((VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - 1].z) / VEC<VType>::h.x + bdiff.x.z);
				curl.z += 0.5 * ((VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - 1].y) / VEC<VType>::h.x + bdiff.x.y);
			}
		}
	}

	//y direction differentials
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
		curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.y;
			curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / VEC<VType>::h.y;
			curl.z -= (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / VEC<VType>::h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
				curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
			}
			else {

				curl.x += 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.y + bdiff.y.z);
				curl.z -= 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.y + bdiff.y.x);
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
				curl.z -= (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
			}
			else {

				curl.x += 0.5 * ((VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / VEC<VType>::h.y + bdiff.y.z);
				curl.z -= 0.5 * ((VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / VEC<VType>::h.y + bdiff.y.x);
			}
		}
	}

	//z direction differentials
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
		curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.z;
			curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / VEC<VType>::h.z;
			curl.y += (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / VEC<VType>::h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
				curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
			}
			else {

				curl.x -= 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.z + bdiff.z.y);
				curl.y += 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.z + bdiff.z.x);
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
				curl.y += (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
			}
			else {

				curl.x -= 0.5 * ((VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / VEC<VType>::h.z + bdiff.z.y);
				curl.y += 0.5 * ((VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / VEC<VType>::h.z + bdiff.z.x);
			}
		}
	}

	return curl;
}

//curl operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
//can only be applied if VType is a VAL3
template <typename VType>
VType VEC_VC<VType>::curl_diri(int idx) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		curl.y -= (VEC<VType>::quantity[idx + 1].z - VEC<VType>::quantity[idx - 1].z) / (2 * VEC<VType>::h.x);
		curl.z += (VEC<VType>::quantity[idx + 1].y - VEC<VType>::quantity[idx - 1].y) / (2 * VEC<VType>::h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

			curl.y -= 0.5 * (VEC<VType>::quantity[idx + 1].z + VEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).z) / VEC<VType>::h.x;
			curl.z += 0.5 * (VEC<VType>::quantity[idx + 1].y + VEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).y) / VEC<VType>::h.x;
		}
		else {

			curl.y -= 0.5 * (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).z - VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - 1].z) / VEC<VType>::h.x;
			curl.z += 0.5 * (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).y - VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - 1].y) / VEC<VType>::h.x;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (VEC<VType>::quantity[idx + 1].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.x;
			curl.z += (VEC<VType>::quantity[idx + 1].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - 1].z) / VEC<VType>::h.x;
			curl.z += (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - 1].y) / VEC<VType>::h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (VEC<VType>::quantity[idx + 1].z - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].z) / (2 * VEC<VType>::h.x);
				curl.z += (VEC<VType>::quantity[idx + 1].y - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].y) / (2 * VEC<VType>::h.x);
			}
			else {

				curl.y -= 0.5 * (VEC<VType>::quantity[idx + 1].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.x;
				curl.z += 0.5 * (VEC<VType>::quantity[idx + 1].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.x;
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].z - VEC<VType>::quantity[idx - 1].z) / (2 * VEC<VType>::h.x);
				curl.z += (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].y - VEC<VType>::quantity[idx - 1].y) / (2 * VEC<VType>::h.x);
			}
			else {

				curl.y -= 0.5 * (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - 1].z) / VEC<VType>::h.x;
				curl.z += 0.5 * (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - 1].y) / VEC<VType>::h.x;
			}
		}
	}

	//y direction differentials
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
		curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

			curl.x += 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x].z + VEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).z) / VEC<VType>::h.y;
			curl.z -= 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x].x + VEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).x) / VEC<VType>::h.y;
		}
		else {

			curl.x += 0.5 * (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).z - VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / VEC<VType>::h.y;
			curl.z -= 0.5 * (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).x - VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / VEC<VType>::h.y;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.y;
			curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / VEC<VType>::h.y;
			curl.z -= (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / VEC<VType>::h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
				curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
			}
			else {

				curl.x += 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.y;
				curl.z -= 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.y;
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
				curl.z -= (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
			}
			else {

				curl.x += 0.5 * (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / VEC<VType>::h.y;
				curl.z -= 0.5 * (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / VEC<VType>::h.y;
			}
		}
	}

	//z direction differentials
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
		curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

			curl.x -= 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y + VEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).y) / VEC<VType>::h.z;
			curl.y += 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x + VEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).x) / VEC<VType>::h.z;
		}
		else {

			curl.x -= 0.5 * (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).y - VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / VEC<VType>::h.z;
			curl.y += 0.5 * (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).x - VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / VEC<VType>::h.z;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.z;
			curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / VEC<VType>::h.z;
			curl.y += (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / VEC<VType>::h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
				curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
			}
			else {

				curl.x -= 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.z;
				curl.y += 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.z;
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
				curl.y += (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
			}
			else {

				curl.x -= 0.5 * (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / VEC<VType>::h.z;
				curl.y += 0.5 * (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / VEC<VType>::h.z;
			}
		}
	}

	return curl;
}

//curl operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
//Can be used at composite media boundaries where sided differentials will be used instead.
//can only be applied if VType is a VAL3
template <typename VType>
VType VEC_VC<VType>::curl_diri_nneu(int idx, const VAL3<VType>& bdiff) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		curl.y -= (VEC<VType>::quantity[idx + 1].z - VEC<VType>::quantity[idx - 1].z) / (2 * VEC<VType>::h.x);
		curl.z += (VEC<VType>::quantity[idx + 1].y - VEC<VType>::quantity[idx - 1].y) / (2 * VEC<VType>::h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPX) {

			curl.y -= 0.5 * (VEC<VType>::quantity[idx + 1].z + VEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).z) / VEC<VType>::h.x;
			curl.z += 0.5 * (VEC<VType>::quantity[idx + 1].y + VEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).y) / VEC<VType>::h.x;
		}
		else {

			curl.y -= 0.5 * (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).z - VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - 1].z) / VEC<VType>::h.x;
			curl.z += 0.5 * (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).y - VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - 1].y) / VEC<VType>::h.x;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (VEC<VType>::quantity[idx + 1].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.x;
			curl.z += (VEC<VType>::quantity[idx + 1].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - 1].z) / VEC<VType>::h.x;
			curl.z += (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - 1].y) / VEC<VType>::h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (VEC<VType>::quantity[idx + 1].z - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].z) / (2 * VEC<VType>::h.x);
				curl.z += (VEC<VType>::quantity[idx + 1].y - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].y) / (2 * VEC<VType>::h.x);
			}
			else {

				curl.y -= 0.5 * ((VEC<VType>::quantity[idx + 1].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.x + bdiff.x.z);
				curl.z += 0.5 * ((VEC<VType>::quantity[idx + 1].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.x + bdiff.x.y);
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].z - VEC<VType>::quantity[idx - 1].z) / (2 * VEC<VType>::h.x);
				curl.z += (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].y - VEC<VType>::quantity[idx - 1].y) / (2 * VEC<VType>::h.x);
			}
			else {

				curl.y -= 0.5 * ((VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - 1].z) / VEC<VType>::h.x + bdiff.x.z);
				curl.z += 0.5 * ((VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - 1].y) / VEC<VType>::h.x + bdiff.x.y);
			}
		}
	}

	//y direction differentials
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
		curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPY) {

			curl.x += 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x].z + VEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).z) / VEC<VType>::h.y;
			curl.z -= 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x].x + VEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).x) / VEC<VType>::h.y;
		}
		else {

			curl.x += 0.5 * (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).z - VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / VEC<VType>::h.y;
			curl.z -= 0.5 * (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).x - VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / VEC<VType>::h.y;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.y;
			curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / VEC<VType>::h.y;
			curl.z -= (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / VEC<VType>::h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
				curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
			}
			else {

				curl.x += 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.y + bdiff.y.z);
				curl.z -= 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.y + bdiff.y.x);
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
				curl.z -= (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
			}
			else {

				curl.x += 0.5 * ((VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / VEC<VType>::h.y + bdiff.y.z);
				curl.z -= 0.5 * ((VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / VEC<VType>::h.y + bdiff.y.x);
			}
		}
	}

	//z direction differentials
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
		curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) {

			curl.x -= 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y + VEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).y) / VEC<VType>::h.z;
			curl.y += 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x + VEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).x) / VEC<VType>::h.z;
		}
		else {

			curl.x -= 0.5 * (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).y - VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / VEC<VType>::h.z;
			curl.y += 0.5 * (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).x - VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / VEC<VType>::h.z;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.z;
			curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / VEC<VType>::h.z;
			curl.y += (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / VEC<VType>::h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
				curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
			}
			else {

				curl.x -= 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.z + bdiff.z.y);
				curl.y += 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.z + bdiff.z.x);
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
				curl.y += (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
			}
			else {

				curl.x -= 0.5 * ((VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / VEC<VType>::h.z + bdiff.z.y);
				curl.y += 0.5 * ((VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / VEC<VType>::h.z + bdiff.z.x);
			}
		}
	}

	return curl;
}

//curl operator. Use sided differentials at boundaries (including at composite media boundaries)
//can only be applied if VType is a VAL3
template <typename VType>
VType VEC_VC<VType>::curl_sided(int idx) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		curl.y -= (VEC<VType>::quantity[idx + 1].z - VEC<VType>::quantity[idx - 1].z) / (2 * VEC<VType>::h.x);
		curl.z += (VEC<VType>::quantity[idx + 1].y - VEC<VType>::quantity[idx - 1].y) / (2 * VEC<VType>::h.x);
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (VEC<VType>::quantity[idx + 1].z - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].z) / (2 * VEC<VType>::h.x);
				curl.z += (VEC<VType>::quantity[idx + 1].y - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].y) / (2 * VEC<VType>::h.x);
			}
			else {

				curl.y -= (VEC<VType>::quantity[idx + 1].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.x;
				curl.z += (VEC<VType>::quantity[idx + 1].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.x;
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].z - VEC<VType>::quantity[idx - 1].z) / (2 * VEC<VType>::h.x);
				curl.z += (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].y - VEC<VType>::quantity[idx - 1].y) / (2 * VEC<VType>::h.x);
			}
			else {

				curl.y -= (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - 1].z) / VEC<VType>::h.x;
				curl.z += (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - 1].y) / VEC<VType>::h.x;
			}
		}
	}

	//y direction differentials
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
		curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
				curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1) * VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
			}
			else {

				curl.x += (VEC<VType>::quantity[idx + VEC<VType>::n.x].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.y;
				curl.z -= (VEC<VType>::quantity[idx + VEC<VType>::n.x].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.y;
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / (2 * VEC<VType>::h.y);
				curl.z -= (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1) * VEC<VType>::n.x].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / (2 * VEC<VType>::h.y);
			}
			else {

				curl.x += (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x].z) / VEC<VType>::h.y;
				curl.z -= (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x].x) / VEC<VType>::h.y;
			}
		}
	}

	//z direction differentials
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
		curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
				curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
			}
			else {

				curl.x -= (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.z;
				curl.y += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.z;
			}
		}
		//we know there is exactly one neighbor along this axis, so if it's not one must be the other
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				curl.x -= (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / (2 * VEC<VType>::h.z);
				curl.y += (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / (2 * VEC<VType>::h.z);
			}
			else {

				curl.x -= (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].y) / VEC<VType>::h.z;
				curl.y += (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].x) / VEC<VType>::h.z;
			}
		}
	}

	return curl;
}