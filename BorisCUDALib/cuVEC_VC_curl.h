#pragma once

#include "cuVEC_VC.h"

//-------------------------------- CURL OPERATOR

//curl operator. Use Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
//can only be applied if VType is a VAL3
template <typename VType>
__device__ VType cuVEC_VC<VType>::curl_neu(int idx) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if (ngbrFlags[idx] & NF_BOTHX) {

		curl.y -= (quantity[idx + 1].z - quantity[idx - 1].z) / (2 * h.x);
		curl.z += (quantity[idx + 1].y - quantity[idx - 1].y) / (2 * h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (quantity[idx + 1].z - quantity[idx].z) / h.x;
			curl.z += (quantity[idx + 1].y - quantity[idx].y) / h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (quantity[idx].z - quantity[idx - 1].z) / h.x;
			curl.z += (quantity[idx].y - quantity[idx - 1].y) / h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (quantity[idx + 1].z - quantity[idx + n.x - 1].z) / (2 * h.x);
				curl.z += (quantity[idx + 1].y - quantity[idx + n.x - 1].y) / (2 * h.x);
			}
			else {

				curl.y -= 0.5 * (quantity[idx + 1].z - quantity[idx].z) / h.x;
				curl.z += 0.5 * (quantity[idx + 1].y - quantity[idx].y) / h.x;
			}
		}

		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (quantity[idx - (n.x - 1)].z - quantity[idx - 1].z) / (2 * h.x);
				curl.z += (quantity[idx - (n.x - 1)].y - quantity[idx - 1].y) / (2 * h.x);
			}
			else {

				curl.y -= 0.5 * (quantity[idx].z - quantity[idx - 1].z) / h.x;
				curl.z += 0.5 * (quantity[idx].y - quantity[idx - 1].y) / h.x;
			}
		}
	}

	//y direction differentials
	if (ngbrFlags[idx] & NF_BOTHY) {

		curl.x += (quantity[idx + n.x].z - quantity[idx - n.x].z) / (2 * h.y);
		curl.z -= (quantity[idx + n.x].x - quantity[idx - n.x].x) / (2 * h.y);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (quantity[idx + n.x].z - quantity[idx].z) / h.y;
			curl.z -= (quantity[idx + n.x].x - quantity[idx].x) / h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (quantity[idx].z - quantity[idx - n.x].z) / h.y;
			curl.z -= (quantity[idx].x - quantity[idx - n.x].x) / h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (quantity[idx + n.x].z - quantity[idx + (n.y - 1) * n.x].z) / (2 * h.y);
				curl.z -= (quantity[idx + n.x].x - quantity[idx + (n.y - 1) * n.x].x) / (2 * h.y);
			}
			else {

				curl.x += 0.5 * (quantity[idx + n.x].z - quantity[idx].z) / h.y;
				curl.z -= 0.5 * (quantity[idx + n.x].x - quantity[idx].x) / h.y;
			}
		}

		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (quantity[idx - (n.y - 1) * n.x].z - quantity[idx - n.x].z) / (2 * h.y);
				curl.z -= (quantity[idx - (n.y - 1) * n.x].x - quantity[idx - n.x].x) / (2 * h.y);
			}
			else {

				curl.x += 0.5 * (quantity[idx].z - quantity[idx - n.x].z) / h.y;
				curl.z -= 0.5 * (quantity[idx].x - quantity[idx - n.x].x) / h.y;
			}
		}
	}

	//z direction differentials
	if (ngbrFlags[idx] & NF_BOTHZ) {

		curl.x -= (quantity[idx + n.x*n.y].y - quantity[idx - n.x*n.y].y) / (2 * h.z);
		curl.y += (quantity[idx + n.x*n.y].x - quantity[idx - n.x*n.y].x) / (2 * h.z);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (quantity[idx + n.x*n.y].y - quantity[idx].y) / h.z;
			curl.y += (quantity[idx + n.x*n.y].x - quantity[idx].x) / h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z;
			curl.y += (quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= 0.5 * (quantity[idx + n.x*n.y].y - quantity[idx].y) / h.z;
			curl.y += 0.5 * (quantity[idx + n.x*n.y].x - quantity[idx].x) / h.z;
		}

		else {

			curl.x -= 0.5 * (quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z;
			curl.y += 0.5 * (quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z;
		}
	}

	return curl;
}


//curl operator. Use non-homogeneous Neumann boundary conditions.
//Can be used at composite media boundaries where sided differentials will be used instead.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
//can only be applied if VType is a VAL3
template <typename VType>
template <typename Class_BDiff>
__device__ VType cuVEC_VC<VType>::curl_nneu(int idx, Class_BDiff& bdiff_class) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if (ngbrFlags[idx] & NF_BOTHX) {

		curl.y -= (quantity[idx + 1].z - quantity[idx - 1].z) / (2 * h.x);
		curl.z += (quantity[idx + 1].y - quantity[idx - 1].y) / (2 * h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (quantity[idx + 1].z - quantity[idx].z) / h.x;
			curl.z += (quantity[idx + 1].y - quantity[idx].y) / h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (quantity[idx].z - quantity[idx - 1].z) / h.x;
			curl.z += (quantity[idx].y - quantity[idx - 1].y) / h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (quantity[idx + 1].z - quantity[idx + n.x - 1].z) / (2 * h.x);
				curl.z += (quantity[idx + 1].y - quantity[idx + n.x - 1].y) / (2 * h.x);
			}
			else {

				curl.y -= 0.5 * ((quantity[idx + 1].z - quantity[idx].z) / h.x + bdiff_val.x.z);
				curl.z += 0.5 * ((quantity[idx + 1].y - quantity[idx].y) / h.x + bdiff_val.x.y);
			}
		}

		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (quantity[idx - (n.x - 1)].z - quantity[idx - 1].z) / (2 * h.x);
				curl.z += (quantity[idx - (n.x - 1)].y - quantity[idx - 1].y) / (2 * h.x);
			}
			else {

				curl.y -= 0.5 * ((quantity[idx].z - quantity[idx - 1].z) / h.x + bdiff_val.x.z);
				curl.z += 0.5 * ((quantity[idx].y - quantity[idx - 1].y) / h.x + bdiff_val.x.y);
			}
		}
	}

	//y direction differentials
	if (ngbrFlags[idx] & NF_BOTHY) {

		curl.x += (quantity[idx + n.x].z - quantity[idx - n.x].z) / (2 * h.y);
		curl.z -= (quantity[idx + n.x].x - quantity[idx - n.x].x) / (2 * h.y);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (quantity[idx + n.x].z - quantity[idx].z) / h.y;
			curl.z -= (quantity[idx + n.x].x - quantity[idx].x) / h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (quantity[idx].z - quantity[idx - n.x].z) / h.y;
			curl.z -= (quantity[idx].x - quantity[idx - n.x].x) / h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (quantity[idx + n.x].z - quantity[idx + (n.y - 1) * n.x].z) / (2 * h.y);
				curl.z -= (quantity[idx + n.x].x - quantity[idx + (n.y - 1) * n.x].x) / (2 * h.y);
			}
			else {

				curl.x += 0.5 * ((quantity[idx + n.x].z - quantity[idx].z) / h.y + bdiff_val.y.z);
				curl.z -= 0.5 * ((quantity[idx + n.x].x - quantity[idx].x) / h.y + bdiff_val.y.x);
			}
		}

		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (quantity[idx - (n.y - 1) * n.x].z - quantity[idx - n.x].z) / (2 * h.y);
				curl.z -= (quantity[idx - (n.y - 1) * n.x].x - quantity[idx - n.x].x) / (2 * h.y);
			}
			else {

				curl.x += 0.5 * ((quantity[idx].z - quantity[idx - n.x].z) / h.y + bdiff_val.y.z);
				curl.z -= 0.5 * ((quantity[idx].x - quantity[idx - n.x].x) / h.y + bdiff_val.y.x);
			}
		}
	}

	//z direction differentials
	if (ngbrFlags[idx] & NF_BOTHZ) {

		curl.x -= (quantity[idx + n.x*n.y].y - quantity[idx - n.x*n.y].y) / (2 * h.z);
		curl.y += (quantity[idx + n.x*n.y].x - quantity[idx - n.x*n.y].x) / (2 * h.z);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (quantity[idx + n.x*n.y].y - quantity[idx].y) / h.z;
			curl.y += (quantity[idx + n.x*n.y].x - quantity[idx].x) / h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z;
			curl.y += (quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= 0.5 * ((quantity[idx + n.x*n.y].y - quantity[idx].y) / h.z + bdiff_val.z.y);
			curl.y += 0.5 * ((quantity[idx + n.x*n.y].x - quantity[idx].x) / h.z + bdiff_val.z.x);
		}

		else {

			curl.x -= 0.5 * ((quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z + bdiff_val.z.y);
			curl.y += 0.5 * ((quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z + bdiff_val.z.x);
		}
	}

	return curl;
}

//Same as above but boundary conditions specified using a constant
template <typename VType>
__device__ VType cuVEC_VC<VType>::curl_nneu(int idx, cuVAL3<VType>& bdiff) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if (ngbrFlags[idx] & NF_BOTHX) {

		curl.y -= (quantity[idx + 1].z - quantity[idx - 1].z) / (2 * h.x);
		curl.z += (quantity[idx + 1].y - quantity[idx - 1].y) / (2 * h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (quantity[idx + 1].z - quantity[idx].z) / h.x;
			curl.z += (quantity[idx + 1].y - quantity[idx].y) / h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (quantity[idx].z - quantity[idx - 1].z) / h.x;
			curl.z += (quantity[idx].y - quantity[idx - 1].y) / h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (quantity[idx + 1].z - quantity[idx + n.x - 1].z) / (2 * h.x);
				curl.z += (quantity[idx + 1].y - quantity[idx + n.x - 1].y) / (2 * h.x);
			}
			else {

				curl.y -= 0.5 * ((quantity[idx + 1].z - quantity[idx].z) / h.x + bdiff.x.z);
				curl.z += 0.5 * ((quantity[idx + 1].y - quantity[idx].y) / h.x + bdiff.x.y);
			}
		}

		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (quantity[idx - (n.x - 1)].z - quantity[idx - 1].z) / (2 * h.x);
				curl.z += (quantity[idx - (n.x - 1)].y - quantity[idx - 1].y) / (2 * h.x);
			}
			else {

				curl.y -= 0.5 * ((quantity[idx].z - quantity[idx - 1].z) / h.x + bdiff.x.z);
				curl.z += 0.5 * ((quantity[idx].y - quantity[idx - 1].y) / h.x + bdiff.x.y);
			}
		}
	}

	//y direction differentials
	if (ngbrFlags[idx] & NF_BOTHY) {

		curl.x += (quantity[idx + n.x].z - quantity[idx - n.x].z) / (2 * h.y);
		curl.z -= (quantity[idx + n.x].x - quantity[idx - n.x].x) / (2 * h.y);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (quantity[idx + n.x].z - quantity[idx].z) / h.y;
			curl.z -= (quantity[idx + n.x].x - quantity[idx].x) / h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (quantity[idx].z - quantity[idx - n.x].z) / h.y;
			curl.z -= (quantity[idx].x - quantity[idx - n.x].x) / h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (quantity[idx + n.x].z - quantity[idx + (n.y - 1) * n.x].z) / (2 * h.y);
				curl.z -= (quantity[idx + n.x].x - quantity[idx + (n.y - 1) * n.x].x) / (2 * h.y);
			}
			else {

				curl.x += 0.5 * ((quantity[idx + n.x].z - quantity[idx].z) / h.y + bdiff.y.z);
				curl.z -= 0.5 * ((quantity[idx + n.x].x - quantity[idx].x) / h.y + bdiff.y.x);
			}
		}

		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (quantity[idx - (n.y - 1) * n.x].z - quantity[idx - n.x].z) / (2 * h.y);
				curl.z -= (quantity[idx - (n.y - 1) * n.x].x - quantity[idx - n.x].x) / (2 * h.y);
			}
			else {

				curl.x += 0.5 * ((quantity[idx].z - quantity[idx - n.x].z) / h.y + bdiff.y.z);
				curl.z -= 0.5 * ((quantity[idx].x - quantity[idx - n.x].x) / h.y + bdiff.y.x);
			}
		}
	}

	//z direction differentials
	if (ngbrFlags[idx] & NF_BOTHZ) {

		curl.x -= (quantity[idx + n.x*n.y].y - quantity[idx - n.x*n.y].y) / (2 * h.z);
		curl.y += (quantity[idx + n.x*n.y].x - quantity[idx - n.x*n.y].x) / (2 * h.z);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (quantity[idx + n.x*n.y].y - quantity[idx].y) / h.z;
			curl.y += (quantity[idx + n.x*n.y].x - quantity[idx].x) / h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z;
			curl.y += (quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= 0.5 * ((quantity[idx + n.x*n.y].y - quantity[idx].y) / h.z + bdiff.z.y);
			curl.y += 0.5 * ((quantity[idx + n.x*n.y].x - quantity[idx].x) / h.z + bdiff.z.x);
		}

		else {

			curl.x -= 0.5 * ((quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z + bdiff.z.y);
			curl.y += 0.5 * ((quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z + bdiff.z.x);
		}
	}

	return curl;
}

//curl operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
//can only be applied if VType is a VAL3
template <typename VType>
__device__ VType cuVEC_VC<VType>::curl_diri(int idx) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if (ngbrFlags[idx] & NF_BOTHX) {

		curl.y -= (quantity[idx + 1].z - quantity[idx - 1].z) / (2 * h.x);
		curl.z += (quantity[idx + 1].y - quantity[idx - 1].y) / (2 * h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags[idx] & NF_DIRICHLETX) {

		if (ngbrFlags[idx] & NF_DIRICHLETPX) {

			curl.y -= 0.5 * (quantity[idx + 1].z + quantity[idx].z - 2 * get_dirichlet_value(NF_DIRICHLETPX, idx).z) / h.x;
			curl.z += 0.5 * (quantity[idx + 1].y + quantity[idx].y - 2 * get_dirichlet_value(NF_DIRICHLETPX, idx).y) / h.x;
		}
		else {

			curl.y -= 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNX, idx).z - quantity[idx].z - quantity[idx - 1].z) / h.x;
			curl.z += 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNX, idx).y - quantity[idx].y - quantity[idx - 1].y) / h.x;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (quantity[idx + 1].z - quantity[idx].z) / h.x;
			curl.z += (quantity[idx + 1].y - quantity[idx].y) / h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (quantity[idx].z - quantity[idx - 1].z) / h.x;
			curl.z += (quantity[idx].y - quantity[idx - 1].y) / h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (quantity[idx + 1].z - quantity[idx + n.x - 1].z) / (2 * h.x);
				curl.z += (quantity[idx + 1].y - quantity[idx + n.x - 1].y) / (2 * h.x);
			}
			else {

				curl.y -= 0.5 * (quantity[idx + 1].z - quantity[idx].z) / h.x;
				curl.z += 0.5 * (quantity[idx + 1].y - quantity[idx].y) / h.x;
			}
		}

		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (quantity[idx - (n.x - 1)].z - quantity[idx - 1].z) / (2 * h.x);
				curl.z += (quantity[idx - (n.x - 1)].y - quantity[idx - 1].y) / (2 * h.x);
			}
			else {

				curl.y -= 0.5 * (quantity[idx].z - quantity[idx - 1].z) / h.x;
				curl.z += 0.5 * (quantity[idx].y - quantity[idx - 1].y) / h.x;
			}
		}
	}

	//y direction differentials
	if (ngbrFlags[idx] & NF_BOTHY) {

		curl.x += (quantity[idx + n.x].z - quantity[idx - n.x].z) / (2 * h.y);
		curl.z -= (quantity[idx + n.x].x - quantity[idx - n.x].x) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETY) {

		if (ngbrFlags[idx] & NF_DIRICHLETPY) {

			curl.x += 0.5 * (quantity[idx + n.x].z + quantity[idx].z - 2 * get_dirichlet_value(NF_DIRICHLETPY, idx).z) / h.y;
			curl.z -= 0.5 * (quantity[idx + n.x].x + quantity[idx].x - 2 * get_dirichlet_value(NF_DIRICHLETPY, idx).x) / h.y;
		}
		else {

			curl.x += 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNY, idx).z - quantity[idx].z - quantity[idx - n.x].z) / h.y;
			curl.z -= 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNY, idx).x - quantity[idx].x - quantity[idx - n.x].x) / h.y;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (quantity[idx + n.x].z - quantity[idx].z) / h.y;
			curl.z -= (quantity[idx + n.x].x - quantity[idx].x) / h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (quantity[idx].z - quantity[idx - n.x].z) / h.y;
			curl.z -= (quantity[idx].x - quantity[idx - n.x].x) / h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (quantity[idx + n.x].z - quantity[idx + (n.y - 1) * n.x].z) / (2 * h.y);
				curl.z -= (quantity[idx + n.x].x - quantity[idx + (n.y - 1) * n.x].x) / (2 * h.y);
			}
			else {

				curl.x += 0.5 * (quantity[idx + n.x].z - quantity[idx].z) / h.y;
				curl.z -= 0.5 * (quantity[idx + n.x].x - quantity[idx].x) / h.y;
			}
		}

		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (quantity[idx - (n.y - 1) * n.x].z - quantity[idx - n.x].z) / (2 * h.y);
				curl.z -= (quantity[idx - (n.y - 1) * n.x].x - quantity[idx - n.x].x) / (2 * h.y);
			}
			else {

				curl.x += 0.5 * (quantity[idx].z - quantity[idx - n.x].z) / h.y;
				curl.z -= 0.5 * (quantity[idx].x - quantity[idx - n.x].x) / h.y;
			}
		}
	}

	//z direction differentials
	if (ngbrFlags[idx] & NF_BOTHZ) {

		curl.x -= (quantity[idx + n.x*n.y].y - quantity[idx - n.x*n.y].y) / (2 * h.z);
		curl.y += (quantity[idx + n.x*n.y].x - quantity[idx - n.x*n.y].x) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

		if (ngbrFlags[idx] & NF_DIRICHLETPZ) {

			curl.x -= 0.5 * (quantity[idx + n.x*n.y].y + quantity[idx].y - 2 * get_dirichlet_value(NF_DIRICHLETPZ, idx).y) / h.z;
			curl.y += 0.5 * (quantity[idx + n.x*n.y].x + quantity[idx].x - 2 * get_dirichlet_value(NF_DIRICHLETPZ, idx).x) / h.z;
		}
		else {

			curl.x -= 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNZ, idx).y - quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z;
			curl.y += 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNZ, idx).x - quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (quantity[idx + n.x*n.y].y - quantity[idx].y) / h.z;
			curl.y += (quantity[idx + n.x*n.y].x - quantity[idx].x) / h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z;
			curl.y += (quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= 0.5 * (quantity[idx + n.x*n.y].y - quantity[idx].y) / h.z;
			curl.y += 0.5 * (quantity[idx + n.x*n.y].x - quantity[idx].x) / h.z;
		}

		else {

			curl.x -= 0.5 * (quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z;
			curl.y += 0.5 * (quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z;
		}
	}

	return curl;
}

//curl operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
//Can be used at composite media boundaries where sided differentials will be used instead.
//can only be applied if VType is a VAL3
template <typename VType>
template <typename Class_BDiff>
__device__ VType cuVEC_VC<VType>::curl_diri_nneu(int idx, Class_BDiff& bdiff_class) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if (ngbrFlags[idx] & NF_BOTHX) {

		curl.y -= (quantity[idx + 1].z - quantity[idx - 1].z) / (2 * h.x);
		curl.z += (quantity[idx + 1].y - quantity[idx - 1].y) / (2 * h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags[idx] & NF_DIRICHLETX) {

		if (ngbrFlags[idx] & NF_DIRICHLETPX) {

			curl.y -= 0.5 * (quantity[idx + 1].z + quantity[idx].z - 2 * get_dirichlet_value(NF_DIRICHLETPX, idx).z) / h.x;
			curl.z += 0.5 * (quantity[idx + 1].y + quantity[idx].y - 2 * get_dirichlet_value(NF_DIRICHLETPX, idx).y) / h.x;
		}
		else {

			curl.y -= 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNX, idx).z - quantity[idx].z - quantity[idx - 1].z) / h.x;
			curl.z += 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNX, idx).y - quantity[idx].y - quantity[idx - 1].y) / h.x;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (quantity[idx + 1].z - quantity[idx].z) / h.x;
			curl.z += (quantity[idx + 1].y - quantity[idx].y) / h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (quantity[idx].z - quantity[idx - 1].z) / h.x;
			curl.z += (quantity[idx].y - quantity[idx - 1].y) / h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (quantity[idx + 1].z - quantity[idx + n.x - 1].z) / (2 * h.x);
				curl.z += (quantity[idx + 1].y - quantity[idx + n.x - 1].y) / (2 * h.x);
			}
			else {

				curl.y -= 0.5 * ((quantity[idx + 1].z - quantity[idx].z) / h.x + bdiff_val.x.z);
				curl.z += 0.5 * ((quantity[idx + 1].y - quantity[idx].y) / h.x + bdiff_val.x.y);
			}
		}

		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (quantity[idx - (n.x - 1)].z - quantity[idx - 1].z) / (2 * h.x);
				curl.z += (quantity[idx - (n.x - 1)].y - quantity[idx - 1].y) / (2 * h.x);
			}
			else {

				curl.y -= 0.5 * ((quantity[idx].z - quantity[idx - 1].z) / h.x + bdiff_val.x.z);
				curl.z += 0.5 * ((quantity[idx].y - quantity[idx - 1].y) / h.x + bdiff_val.x.y);
			}
		}
	}

	//y direction differentials
	if (ngbrFlags[idx] & NF_BOTHY) {

		curl.x += (quantity[idx + n.x].z - quantity[idx - n.x].z) / (2 * h.y);
		curl.z -= (quantity[idx + n.x].x - quantity[idx - n.x].x) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETY) {

		if (ngbrFlags[idx] & NF_DIRICHLETPY) {

			curl.x += 0.5 * (quantity[idx + n.x].z + quantity[idx].z - 2 * get_dirichlet_value(NF_DIRICHLETPY, idx).z) / h.y;
			curl.z -= 0.5 * (quantity[idx + n.x].x + quantity[idx].x - 2 * get_dirichlet_value(NF_DIRICHLETPY, idx).x) / h.y;
		}
		else {

			curl.x += 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNY, idx).z - quantity[idx].z - quantity[idx - n.x].z) / h.y;
			curl.z -= 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNY, idx).x - quantity[idx].x - quantity[idx - n.x].x) / h.y;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (quantity[idx + n.x].z - quantity[idx].z) / h.y;
			curl.z -= (quantity[idx + n.x].x - quantity[idx].x) / h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (quantity[idx].z - quantity[idx - n.x].z) / h.y;
			curl.z -= (quantity[idx].x - quantity[idx - n.x].x) / h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (quantity[idx + n.x].z - quantity[idx + (n.y - 1) * n.x].z) / (2 * h.y);
				curl.z -= (quantity[idx + n.x].x - quantity[idx + (n.y - 1) * n.x].x) / (2 * h.y);
			}
			else {

				curl.x += 0.5 * ((quantity[idx + n.x].z - quantity[idx].z) / h.y + bdiff_val.y.z);
				curl.z -= 0.5 * ((quantity[idx + n.x].x - quantity[idx].x) / h.y + bdiff_val.y.x);
			}
		}

		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (quantity[idx - (n.y - 1) * n.x].z - quantity[idx - n.x].z) / (2 * h.y);
				curl.z -= (quantity[idx - (n.y - 1) * n.x].x - quantity[idx - n.x].x) / (2 * h.y);
			}
			else {

				curl.x += 0.5 * ((quantity[idx].z - quantity[idx - n.x].z) / h.y + bdiff_val.y.z);
				curl.z -= 0.5 * ((quantity[idx].x - quantity[idx - n.x].x) / h.y + bdiff_val.y.x);
			}
		}
	}

	//z direction differentials
	if (ngbrFlags[idx] & NF_BOTHZ) {

		curl.x -= (quantity[idx + n.x*n.y].y - quantity[idx - n.x*n.y].y) / (2 * h.z);
		curl.y += (quantity[idx + n.x*n.y].x - quantity[idx - n.x*n.y].x) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

		if (ngbrFlags[idx] & NF_DIRICHLETPZ) {

			curl.x -= 0.5 * (quantity[idx + n.x*n.y].y + quantity[idx].y - 2 * get_dirichlet_value(NF_DIRICHLETPZ, idx).y) / h.z;
			curl.y += 0.5 * (quantity[idx + n.x*n.y].x + quantity[idx].x - 2 * get_dirichlet_value(NF_DIRICHLETPZ, idx).x) / h.z;
		}
		else {

			curl.x -= 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNZ, idx).y - quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z;
			curl.y += 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNZ, idx).x - quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (quantity[idx + n.x*n.y].y - quantity[idx].y) / h.z;
			curl.y += (quantity[idx + n.x*n.y].x - quantity[idx].x) / h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z;
			curl.y += (quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= 0.5 * ((quantity[idx + n.x*n.y].y - quantity[idx].y) / h.z + bdiff_val.z.y);
			curl.y += 0.5 * ((quantity[idx + n.x*n.y].x - quantity[idx].x) / h.z + bdiff_val.z.x);
		}

		else {

			curl.x -= 0.5 * ((quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z + bdiff_val.z.y);
			curl.y += 0.5 * ((quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z + bdiff_val.z.x);
		}
	}

	return curl;
}

//Same as above but boundary conditions specified using a constant
template <typename VType>
__device__ VType cuVEC_VC<VType>::curl_diri_nneu(int idx, cuVAL3<VType>& bdiff) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if (ngbrFlags[idx] & NF_BOTHX) {

		curl.y -= (quantity[idx + 1].z - quantity[idx - 1].z) / (2 * h.x);
		curl.z += (quantity[idx + 1].y - quantity[idx - 1].y) / (2 * h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags[idx] & NF_DIRICHLETX) {

		if (ngbrFlags[idx] & NF_DIRICHLETPX) {

			curl.y -= 0.5 * (quantity[idx + 1].z + quantity[idx].z - 2 * get_dirichlet_value(NF_DIRICHLETPX, idx).z) / h.x;
			curl.z += 0.5 * (quantity[idx + 1].y + quantity[idx].y - 2 * get_dirichlet_value(NF_DIRICHLETPX, idx).y) / h.x;
		}
		else {

			curl.y -= 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNX, idx).z - quantity[idx].z - quantity[idx - 1].z) / h.x;
			curl.z += 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNX, idx).y - quantity[idx].y - quantity[idx - 1].y) / h.x;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) {

			curl.y -= (quantity[idx + 1].z - quantity[idx].z) / h.x;
			curl.z += (quantity[idx + 1].y - quantity[idx].y) / h.x;
		}

		if (ngbrFlags[idx] & NF_NNX) {

			curl.y -= (quantity[idx].z - quantity[idx - 1].z) / h.x;
			curl.z += (quantity[idx].y - quantity[idx - 1].y) / h.x;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (quantity[idx + 1].z - quantity[idx + n.x - 1].z) / (2 * h.x);
				curl.z += (quantity[idx + 1].y - quantity[idx + n.x - 1].y) / (2 * h.x);
			}
			else {

				curl.y -= 0.5 * ((quantity[idx + 1].z - quantity[idx].z) / h.x + bdiff.x.z);
				curl.z += 0.5 * ((quantity[idx + 1].y - quantity[idx].y) / h.x + bdiff.x.y);
			}
		}

		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (quantity[idx - (n.x - 1)].z - quantity[idx - 1].z) / (2 * h.x);
				curl.z += (quantity[idx - (n.x - 1)].y - quantity[idx - 1].y) / (2 * h.x);
			}
			else {

				curl.y -= 0.5 * ((quantity[idx].z - quantity[idx - 1].z) / h.x + bdiff.x.z);
				curl.z += 0.5 * ((quantity[idx].y - quantity[idx - 1].y) / h.x + bdiff.x.y);
			}
		}
	}

	//y direction differentials
	if (ngbrFlags[idx] & NF_BOTHY) {

		curl.x += (quantity[idx + n.x].z - quantity[idx - n.x].z) / (2 * h.y);
		curl.z -= (quantity[idx + n.x].x - quantity[idx - n.x].x) / (2 * h.y);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETY) {

		if (ngbrFlags[idx] & NF_DIRICHLETPY) {

			curl.x += 0.5 * (quantity[idx + n.x].z + quantity[idx].z - 2 * get_dirichlet_value(NF_DIRICHLETPY, idx).z) / h.y;
			curl.z -= 0.5 * (quantity[idx + n.x].x + quantity[idx].x - 2 * get_dirichlet_value(NF_DIRICHLETPY, idx).x) / h.y;
		}
		else {

			curl.x += 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNY, idx).z - quantity[idx].z - quantity[idx - n.x].z) / h.y;
			curl.z -= 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNY, idx).x - quantity[idx].x - quantity[idx - n.x].x) / h.y;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) {

			curl.x += (quantity[idx + n.x].z - quantity[idx].z) / h.y;
			curl.z -= (quantity[idx + n.x].x - quantity[idx].x) / h.y;
		}

		if (ngbrFlags[idx] & NF_NNY) {

			curl.x += (quantity[idx].z - quantity[idx - n.x].z) / h.y;
			curl.z -= (quantity[idx].x - quantity[idx - n.x].x) / h.y;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (quantity[idx + n.x].z - quantity[idx + (n.y - 1) * n.x].z) / (2 * h.y);
				curl.z -= (quantity[idx + n.x].x - quantity[idx + (n.y - 1) * n.x].x) / (2 * h.y);
			}
			else {

				curl.x += 0.5 * ((quantity[idx + n.x].z - quantity[idx].z) / h.y + bdiff.y.z);
				curl.z -= 0.5 * ((quantity[idx + n.x].x - quantity[idx].x) / h.y + bdiff.y.x);
			}
		}

		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (quantity[idx - (n.y - 1) * n.x].z - quantity[idx - n.x].z) / (2 * h.y);
				curl.z -= (quantity[idx - (n.y - 1) * n.x].x - quantity[idx - n.x].x) / (2 * h.y);
			}
			else {

				curl.x += 0.5 * ((quantity[idx].z - quantity[idx - n.x].z) / h.y + bdiff.y.z);
				curl.z -= 0.5 * ((quantity[idx].x - quantity[idx - n.x].x) / h.y + bdiff.y.x);
			}
		}
	}

	//z direction differentials
	if (ngbrFlags[idx] & NF_BOTHZ) {

		curl.x -= (quantity[idx + n.x*n.y].y - quantity[idx - n.x*n.y].y) / (2 * h.z);
		curl.y += (quantity[idx + n.x*n.y].x - quantity[idx - n.x*n.y].x) / (2 * h.z);
	}
	else if (ngbrFlags[idx] & NF_DIRICHLETZ) {

		if (ngbrFlags[idx] & NF_DIRICHLETPZ) {

			curl.x -= 0.5 * (quantity[idx + n.x*n.y].y + quantity[idx].y - 2 * get_dirichlet_value(NF_DIRICHLETPZ, idx).y) / h.z;
			curl.y += 0.5 * (quantity[idx + n.x*n.y].x + quantity[idx].x - 2 * get_dirichlet_value(NF_DIRICHLETPZ, idx).x) / h.z;
		}
		else {

			curl.x -= 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNZ, idx).y - quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z;
			curl.y += 0.5 * (2 * get_dirichlet_value(NF_DIRICHLETNZ, idx).x - quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z;
		}
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (quantity[idx + n.x*n.y].y - quantity[idx].y) / h.z;
			curl.y += (quantity[idx + n.x*n.y].x - quantity[idx].x) / h.z;
		}

		if (ngbrFlags[idx] & NF_NNZ) {

			curl.x -= (quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z;
			curl.y += (quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z;
		}
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= 0.5 * ((quantity[idx + n.x*n.y].y - quantity[idx].y) / h.z + bdiff.z.y);
			curl.y += 0.5 * ((quantity[idx + n.x*n.y].x - quantity[idx].x) / h.z + bdiff.z.x);
		}

		else {

			curl.x -= 0.5 * ((quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z + bdiff.z.y);
			curl.y += 0.5 * ((quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z + bdiff.z.x);
		}
	}

	return curl;
}

//curl operator. Use sided differentials at boundaries (including at composite media boundaries)
//can only be applied if VType is a VAL3
template <typename VType>
__device__ VType cuVEC_VC<VType>::curl_sided(int idx) const
{
	VType curl = VType();

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return curl;

	//x direction differentials
	if (ngbrFlags[idx] & NF_BOTHX) {

		curl.y -= (quantity[idx + 1].z - quantity[idx - 1].z) / (2 * h.x);
		curl.z += (quantity[idx + 1].y - quantity[idx - 1].y) / (2 * h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (quantity[idx + 1].z - quantity[idx + n.x - 1].z) / (2 * h.x);
				curl.z += (quantity[idx + 1].y - quantity[idx + n.x - 1].y) / (2 * h.x);
			}
			else {

				curl.y -= (quantity[idx + 1].z - quantity[idx].z) / h.x;
				curl.z += (quantity[idx + 1].y - quantity[idx].y) / h.x;
			}
		}

		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				curl.y -= (quantity[idx - (n.x - 1)].z - quantity[idx - 1].z) / (2 * h.x);
				curl.z += (quantity[idx - (n.x - 1)].y - quantity[idx - 1].y) / (2 * h.x);
			}
			else {

				curl.y -= (quantity[idx].z - quantity[idx - 1].z) / h.x;
				curl.z += (quantity[idx].y - quantity[idx - 1].y) / h.x;
			}
		}
	}

	//y direction differentials
	if (ngbrFlags[idx] & NF_BOTHY) {

		curl.x += (quantity[idx + n.x].z - quantity[idx - n.x].z) / (2 * h.y);
		curl.z -= (quantity[idx + n.x].x - quantity[idx - n.x].x) / (2 * h.y);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (quantity[idx + n.x].z - quantity[idx + (n.y - 1) * n.x].z) / (2 * h.y);
				curl.z -= (quantity[idx + n.x].x - quantity[idx + (n.y - 1) * n.x].x) / (2 * h.y);
			}
			else {

				curl.x += (quantity[idx + n.x].z - quantity[idx].z) / h.y;
				curl.z -= (quantity[idx + n.x].x - quantity[idx].x) / h.y;
			}
		}

		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				curl.x += (quantity[idx - (n.y - 1) * n.x].z - quantity[idx - n.x].z) / (2 * h.y);
				curl.z -= (quantity[idx - (n.y - 1) * n.x].x - quantity[idx - n.x].x) / (2 * h.y);
			}
			else {

				curl.x += (quantity[idx].z - quantity[idx - n.x].z) / h.y;
				curl.z -= (quantity[idx].x - quantity[idx - n.x].x) / h.y;
			}
		}
	}

	//z direction differentials
	if (ngbrFlags[idx] & NF_BOTHZ) {

		curl.x -= (quantity[idx + n.x*n.y].y - quantity[idx - n.x*n.y].y) / (2 * h.z);
		curl.y += (quantity[idx + n.x*n.y].x - quantity[idx - n.x*n.y].x) / (2 * h.z);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			curl.x -= (quantity[idx + n.x*n.y].y - quantity[idx].y) / h.z;
			curl.y += (quantity[idx + n.x*n.y].x - quantity[idx].x) / h.z;
		}

		else {

			curl.x -= (quantity[idx].y - quantity[idx - n.x*n.y].y) / h.z;
			curl.y += (quantity[idx].x - quantity[idx - n.x*n.y].x) / h.z;
		}
	}

	return curl;
}