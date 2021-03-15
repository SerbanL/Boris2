#pragma once

#include "cuVEC_VC.h"

//-------------------------------- DIVERGENCE OPERATOR

//divergence operator. Use Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
//div operator can be applied if VType is a VAL3<Type>, returning Type
template __device__ cuBReal cuVEC_VC<cuReal3>::div_neu(int idx) const;

template <typename VType>
__device__ cuBReal cuVEC_VC<VType>::div_neu(int idx) const
{
	cuBReal div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		div += (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) div += (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.x;
		if (ngbrFlags[idx] & NF_NNX) div += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / cuVEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (cuVEC<VType>::quantity[idx + 1].x - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].x) / (2 * cuVEC<VType>::h.x);
			}
			else {

				div += 0.5 * (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
			}
			else {

				div += 0.5 * (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / cuVEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.y;
		if (ngbrFlags[idx] & NF_NNY) div += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cuVEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
			}
			else {

				div += 0.5 * (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
			}
			else {

				div += 0.5 * (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cuVEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.z;
		if (ngbrFlags[idx] & NF_NNZ) div += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / cuVEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
			}
			else {

				div += 0.5 * (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
			}
			else {

				div += 0.5 * (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / cuVEC<VType>::h.z;
			}
		}
	}

	return div;
}

//divergence operator. Use non-homogeneous Neumann boundary conditions.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
//Can be used at composite media boundaries where sided differentials will be used instead.
//div operator can be applied if VType is a VAL3<Type>, returning Type
template <typename VType>
template <typename Class_BDiff>
__device__ cuBReal cuVEC_VC<VType>::div_nneu(int idx, const Class_BDiff& bdiff_class) const
{
	cuBReal div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		div += (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) div += (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.x;
		if (ngbrFlags[idx] & NF_NNX) div += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / cuVEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (cuVEC<VType>::quantity[idx + 1].x - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].x) / (2 * cuVEC<VType>::h.x);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.x + bdiff_val.x.x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / cuVEC<VType>::h.x + bdiff_val.x.x);
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.y;
		if (ngbrFlags[idx] & NF_NNY) div += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cuVEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.y + bdiff_val.y.y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cuVEC<VType>::h.y + bdiff_val.y.y);
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.z;
		if (ngbrFlags[idx] & NF_NNZ) div += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / cuVEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.z + bdiff_val.z.z);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / cuVEC<VType>::h.z + bdiff_val.z.z);
			}
		}
	}

	return div;
}

//Same as above but boundary conditions specified using a constant
template <typename VType>
__device__ cuBReal cuVEC_VC<VType>::div_nneu(int idx, const cuVAL3<VType>& bdiff) const
{
	cuBReal div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		div += (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) div += (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.x;
		if (ngbrFlags[idx] & NF_NNX) div += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / cuVEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (cuVEC<VType>::quantity[idx + 1].x - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].x) / (2 * cuVEC<VType>::h.x);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.x + bdiff.x.x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / cuVEC<VType>::h.x + bdiff.x.x);
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.y;
		if (ngbrFlags[idx] & NF_NNY) div += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cuVEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.y + bdiff.y.y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cuVEC<VType>::h.y + bdiff.y.y);
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.z;
		if (ngbrFlags[idx] & NF_NNZ) div += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / cuVEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.z + bdiff.z.z);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / cuVEC<VType>::h.z + bdiff.z.z);
			}
		}
	}

	return div;
}

//divergence operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
//div operator can be applied if VType is a VAL3<Type>, returning Type
template __device__ cuBReal cuVEC_VC<cuReal3>::div_diri(int idx) const;

template <typename VType>
__device__ cuBReal cuVEC_VC<VType>::div_diri(int idx) const
{
	cuBReal div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		div += (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPX) div += (cuVEC<VType>::quantity[idx + 1].x + cuVEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).x) / (2 * cuVEC<VType>::h.x);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).x - cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
	}
	//Not Dirichlet, is it a CMBND boundary? - if not this either then use homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) div += (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.x;
		if (ngbrFlags[idx] & NF_NNX) div += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / cuVEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (cuVEC<VType>::quantity[idx + 1].x - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].x) / (2 * cuVEC<VType>::h.x);
			}
			else {

				div += 0.5 * (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
			}
			else {

				div += 0.5 * (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / cuVEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
	}
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPY) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y + cuVEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).y) / (2 * cuVEC<VType>::h.y);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).y - cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.y;
		if (ngbrFlags[idx] & NF_NNY) div += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cuVEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
			}
			else {

				div += 0.5 * (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
			}
			else {

				div += 0.5 * (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cuVEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
	}
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z + cuVEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).z) / (2 * cuVEC<VType>::h.z);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).z - cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.z;
		if (ngbrFlags[idx] & NF_NNZ) div += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / cuVEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
			}
			else {

				div += 0.5 * (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
			}
			else {

				div += 0.5 * (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / cuVEC<VType>::h.z;
			}
		}
	}

	return div;
}

//divergence operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions - the class Class_BDiff must define a method bdiff returning a cuVAL3<VType> and taking an int (the cell index)
//Can be used at composite media boundaries where sided differentials will be used instead.
//div operator can be applied if VType is a VAL3<Type>, returning Type
template <typename VType>
template <typename Class_BDiff>
__device__ cuBReal cuVEC_VC<VType>::div_diri_nneu(int idx, const Class_BDiff& bdiff_class) const
{
	cuBReal div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		div += (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPX) div += (cuVEC<VType>::quantity[idx + 1].x + cuVEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).x) / (2 * cuVEC<VType>::h.x);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).x - cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
	}
	//Not Dirichlet, is it a CMBND boundary? - if not this either then use homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) div += (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.x;
		if (ngbrFlags[idx] & NF_NNX) div += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / cuVEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (cuVEC<VType>::quantity[idx + 1].x - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].x) / (2 * cuVEC<VType>::h.x);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.x + bdiff_val.x.x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / cuVEC<VType>::h.x + bdiff_val.x.x);
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
	}
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPY) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y + cuVEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).y) / (2 * cuVEC<VType>::h.y);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).y - cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.y;
		if (ngbrFlags[idx] & NF_NNY) div += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cuVEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.y + bdiff_val.y.y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cuVEC<VType>::h.y + bdiff_val.y.y);
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
	}
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z + cuVEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).z) / (2 * cuVEC<VType>::h.z);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).z - cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.z;
		if (ngbrFlags[idx] & NF_NNZ) div += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / cuVEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		cuVAL3<VType> bdiff_val = bdiff_class.bdiff(idx);

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.z + bdiff_val.z.z);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / cuVEC<VType>::h.z + bdiff_val.z.z);
			}
		}
	}

	return div;
}

//Same as above but boundary conditions specified using a constant
template <typename VType>
__device__ cuBReal cuVEC_VC<VType>::div_diri_nneu(int idx, const cuVAL3<VType>& bdiff) const
{
	cuBReal div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		div += (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPX) div += (cuVEC<VType>::quantity[idx + 1].x + cuVEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).x) / (2 * cuVEC<VType>::h.x);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).x - cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
	}
	//Not Dirichlet, is it a CMBND boundary? - if not this either then use homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) div += (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.x;
		if (ngbrFlags[idx] & NF_NNX) div += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / cuVEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (cuVEC<VType>::quantity[idx + 1].x - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].x) / (2 * cuVEC<VType>::h.x);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.x + bdiff.x.x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / cuVEC<VType>::h.x + bdiff.x.x);
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
	}
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPY) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y + cuVEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).y) / (2 * cuVEC<VType>::h.y);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).y - cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.y;
		if (ngbrFlags[idx] & NF_NNY) div += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cuVEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.y + bdiff.y.y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cuVEC<VType>::h.y + bdiff.y.y);
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
	}
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z + cuVEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).z) / (2 * cuVEC<VType>::h.z);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).z - cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.z;
		if (ngbrFlags[idx] & NF_NNZ) div += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / cuVEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.z + bdiff.z.z);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
			}
			else {

				div += 0.5 * ((cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / cuVEC<VType>::h.z + bdiff.z.z);
			}
		}
	}

	return div;
}

//divergence operator. Use sided differentials (also at composite media boundaries)
//div operator can be applied if VType is a VAL3<Type>, returning Type
template __device__ cuBReal cuVEC_VC<cuReal3>::div_sided(int idx) const;

template <typename VType>
__device__ cuBReal cuVEC_VC<VType>::div_sided(int idx) const
{
	cuBReal div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		div += (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
	}
	//use sided differentials if one of the neighbors is present
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (cuVEC<VType>::quantity[idx + 1].x - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1].x) / (2 * cuVEC<VType>::h.x);
			}
			else {

				div += (cuVEC<VType>::quantity[idx + 1].x - cuVEC<VType>::quantity[idx].x) / cuVEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)].x - cuVEC<VType>::quantity[idx - 1].x) / (2 * cuVEC<VType>::h.x);
			}
			else {

				div += (cuVEC<VType>::quantity[idx].x - cuVEC<VType>::quantity[idx - 1].x) / cuVEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
			}
			else {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx].y) / cuVEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / (2 * cuVEC<VType>::h.y);
			}
			else {

				div += (cuVEC<VType>::quantity[idx].y - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x].y) / cuVEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
			}
			else {

				div += (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx].z) / cuVEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / (2 * cuVEC<VType>::h.z);
			}
			else {

				div += (cuVEC<VType>::quantity[idx].z - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y].z) / cuVEC<VType>::h.z;
			}
		}
	}

	return div;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

//-------------------------------- DIVERGENCE OPERATOR applied after multiplying with unit antisymmetric tensor (epsilon3)

//divergence operator of epsilon3(cuVEC<VType>::quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
//Use Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
template __device__ cuReal3 cuVEC_VC<cuReal3>::diveps3_neu(int idx) const;

template <typename VType>
__device__ VType cuVEC_VC<VType>::diveps3_neu(int idx) const
{
	//differentials along x, y, z directions
	VType diffx, diffy, diffz;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		diffx = (cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx - 1]) / (2 * cuVEC<VType>::h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) diffx = (cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.x;
		if (ngbrFlags[idx] & NF_NNX) diffx = (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - 1]) / cuVEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx = (cuVEC<VType>::quantity[idx + 1] - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1]) / (2 * cuVEC<VType>::h.x);
			}
			else {

				diffx = 0.5 * (cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx = (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)] - cuVEC<VType>::quantity[idx - 1]) / (2 * cuVEC<VType>::h.x);
			}
			else {

				diffx = 0.5 * (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - 1]) / cuVEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diffy = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]) / (2 * cuVEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) diffy = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.y;
		if (ngbrFlags[idx] & NF_NNY) diffy = (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]) / cuVEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]) / (2 * cuVEC<VType>::h.y);
			}
			else {

				diffy = 0.5 * (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy = (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]) / (2 * cuVEC<VType>::h.y);
			}
			else {

				diffy = 0.5 * (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]) / cuVEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diffz = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / (2 * cuVEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) diffz = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.z;
		if (ngbrFlags[idx] & NF_NNZ) diffz = (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / cuVEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diffz = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / (2 * cuVEC<VType>::h.z);
			}
			else {

				diffz = 0.5 * (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diffz = (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / (2 * cuVEC<VType>::h.z);
			}
			else {

				diffz = 0.5 * (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / cuVEC<VType>::h.z;
			}
		}
	}

	return VType(
		-diffy.z + diffz.y,
		-diffz.x + diffx.z,
		-diffx.y + diffy.x
	);
}

//divergence operator of epsilon3(cuVEC<VType>::quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
//Use Dirichlet conditions if set, else Neumann boundary conditions(homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
template __device__ cuReal3 cuVEC_VC<cuReal3>::diveps3_diri(int idx) const;

template <typename VType>
__device__ VType cuVEC_VC<VType>::diveps3_diri(int idx) const
{
	//differentials along x, y, z directions
	VType diffx, diffy, diffz;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		diffx = (cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx - 1]) / (2 * cuVEC<VType>::h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPX) diffx = (cuVEC<VType>::quantity[idx + 1] + cuVEC<VType>::quantity[idx] - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx)) / (2 * cuVEC<VType>::h.x);
		else								 diffx = (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx) - cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - 1]) / (2 * cuVEC<VType>::h.x);
	}
	//Not Dirichlet, is it a CMBND boundary? - if not this either then use homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) diffx = (cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.x;
		if (ngbrFlags[idx] & NF_NNX) diffx = (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - 1]) / cuVEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx = (cuVEC<VType>::quantity[idx + 1] - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1]) / (2 * cuVEC<VType>::h.x);
			}
			else {

				diffx = 0.5 * (cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx = (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)] - cuVEC<VType>::quantity[idx - 1]) / (2 * cuVEC<VType>::h.x);
			}
			else {

				diffx = 0.5 * (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - 1]) / cuVEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diffy = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]) / (2 * cuVEC<VType>::h.y);
	}
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPY) diffy = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] + cuVEC<VType>::quantity[idx] - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx)) / (2 * cuVEC<VType>::h.y);
		else								 diffy = (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx) - cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]) / (2 * cuVEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) diffy = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.y;
		if (ngbrFlags[idx] & NF_NNY) diffy = (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]) / cuVEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]) / (2 * cuVEC<VType>::h.y);
			}
			else {

				diffy = 0.5 * (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy = (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]) / (2 * cuVEC<VType>::h.y);
			}
			else {

				diffy = 0.5 * (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]) / cuVEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diffz = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / (2 * cuVEC<VType>::h.z);
	}
	else if (using_extended_flags && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) diffz = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] + cuVEC<VType>::quantity[idx] - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx)) / (2 * cuVEC<VType>::h.z);
		else								 diffz = (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) - cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / (2 * cuVEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) diffz = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.z;
		if (ngbrFlags[idx] & NF_NNZ) diffz = (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / cuVEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diffz = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / (2 * cuVEC<VType>::h.z);
			}
			else {

				diffz = 0.5 * (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.z;
			}
		}
		else {
			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diffz = (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / (2 * cuVEC<VType>::h.z);
			}
			else {

				diffz = 0.5 * (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / cuVEC<VType>::h.z;
			}
		}
	}

	return VType(
		-diffy.z + diffz.y,
		-diffz.x + diffx.z,
		-diffx.y + diffy.x
	);
}

//divergence operator of epsilon3(cuVEC<VType>::quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
//Use sided differentials (also at composite media boundaries)
template __device__ cuReal3 cuVEC_VC<cuReal3>::diveps3_sided(int idx) const;

template <typename VType>
__device__ VType cuVEC_VC<VType>::diveps3_sided(int idx) const
{
	//differentials along x, y, z directions
	VType diffx, diffy, diffz;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		diffx = (cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx - 1]) / (2 * cuVEC<VType>::h.x);
	}
	//use sided differentials if one of the neighbors is present
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx = (cuVEC<VType>::quantity[idx + 1] - cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x - 1]) / (2 * cuVEC<VType>::h.x);
			}
			else {

				diffx = (cuVEC<VType>::quantity[idx + 1] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx = (cu_get_sign(pbc_x) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.x - 1)] - cuVEC<VType>::quantity[idx - 1]) / (2 * cuVEC<VType>::h.x);
			}
			else {

				diffx = (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - 1]) / cuVEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diffy = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]) / (2 * cuVEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x]) / (2 * cuVEC<VType>::h.y);
			}
			else {

				diffy = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy = (cu_get_sign(pbc_y) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.y - 1)*cuVEC<VType>::n.x] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]) / (2 * cuVEC<VType>::h.y);
			}
			else {

				diffy = (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x]) / cuVEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diffz = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / (2 * cuVEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diffz = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx + (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / (2 * cuVEC<VType>::h.z);
			}
			else {

				diffz = (cuVEC<VType>::quantity[idx + cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx]) / cuVEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diffz = (cu_get_sign(pbc_z) * cuVEC<VType>::quantity[idx - (cuVEC<VType>::n.z - 1) * cuVEC<VType>::n.x*cuVEC<VType>::n.y] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / (2 * cuVEC<VType>::h.z);
			}
			else {

				diffz = (cuVEC<VType>::quantity[idx] - cuVEC<VType>::quantity[idx - cuVEC<VType>::n.x*cuVEC<VType>::n.y]) / cuVEC<VType>::h.z;
			}
		}
	}

	return VType(
		-diffy.z + diffz.y,
		-diffz.x + diffx.z,
		-diffx.y + diffy.x
	);
}