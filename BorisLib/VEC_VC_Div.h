#pragma once

#include "VEC_VC.h"

//-------------------------------- DIVERGENCE OPERATOR

//divergence operator. Use Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
//div operator can be applied if VType is a VAL3<Type>, returning Type
template double VEC_VC<DBL3>::div_neu(int idx) const;

template <typename VType>
double VEC_VC<VType>::div_neu(int idx) const
{
	double div = 0.0;
	
	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		div += (VEC<VType>::quantity[idx + 1].x - VEC<VType>::quantity[idx - 1].x) / (2 * VEC<VType>::h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) div += (VEC<VType>::quantity[idx + 1].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.x;
		if (ngbrFlags[idx] & NF_NNX) div += (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - 1].x) / VEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (VEC<VType>::quantity[idx + 1].x - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].x) / (2 * VEC<VType>::h.x);
			}
			else {

				div += 0.5 * (VEC<VType>::quantity[idx + 1].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].x - VEC<VType>::quantity[idx - 1].x) / (2 * VEC<VType>::h.x);
			}
			else {

				div += 0.5 * (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - 1].x) / VEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.y;
		if (ngbrFlags[idx] & NF_NNY) div += (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / VEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
			}
			else {

				div += 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
			}
			else {

				div += 0.5 * (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / VEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.z;
		if (ngbrFlags[idx] & NF_NNZ) div += (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / VEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
			}
			else {

				div += 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
			}
			else {

				div += 0.5 * (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / VEC<VType>::h.z;
			}
		}
	}
	
	return div;
}

//divergence operator. Use non-homogeneous Neumann boundary conditions.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
//Can be used at composite media boundaries where sided differentials will be used instead.
//div operator can be applied if VType is a VAL3<Type>, returning Type
template double VEC_VC<DBL3>::div_nneu(int idx, const VAL3<DBL3>& bdiff) const;

template <typename VType>
double VEC_VC<VType>::div_nneu(int idx, const VAL3<VType>& bdiff) const
{
	double div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		div += (VEC<VType>::quantity[idx + 1].x - VEC<VType>::quantity[idx - 1].x) / (2 * VEC<VType>::h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) div += (VEC<VType>::quantity[idx + 1].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.x;
		if (ngbrFlags[idx] & NF_NNX) div += (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - 1].x) / VEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (VEC<VType>::quantity[idx + 1].x - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].x) / (2 * VEC<VType>::h.x);
			}
			else {

				div += 0.5 * ((VEC<VType>::quantity[idx + 1].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.x + bdiff.x.x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].x - VEC<VType>::quantity[idx - 1].x) / (2 * VEC<VType>::h.x);
			}
			else {

				div += 0.5 * ((VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - 1].x) / VEC<VType>::h.x + bdiff.x.x);
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.y;
		if (ngbrFlags[idx] & NF_NNY) div += (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / VEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
			}
			else {

				div += 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.y + bdiff.y.y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
			}
			else {

				div += 0.5 * ((VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / VEC<VType>::h.y + bdiff.y.y);
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.z;
		if (ngbrFlags[idx] & NF_NNZ) div += (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / VEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
			}
			else {

				div += 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.z + bdiff.z.z);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
			}
			else {

				div += 0.5 * ((VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / VEC<VType>::h.z + bdiff.z.z);
			}
		}
	}

	return div;
}

//divergence operator. Use Dirichlet conditions if set, else Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
//div operator can be applied if VType is a VAL3<Type>, returning Type
template double VEC_VC<DBL3>::div_diri(int idx) const;

template <typename VType>
double VEC_VC<VType>::div_diri(int idx) const
{
	double div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		div += (VEC<VType>::quantity[idx + 1].x - VEC<VType>::quantity[idx - 1].x) / (2 * VEC<VType>::h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPX) div += (VEC<VType>::quantity[idx + 1].x + VEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).x) / (2 * VEC<VType>::h.x);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).x - VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - 1].x) / (2 * VEC<VType>::h.x);
	}
	//Not Dirichlet, is it a CMBND boundary? - if not this either then use homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) div += (VEC<VType>::quantity[idx + 1].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.x;
		if (ngbrFlags[idx] & NF_NNX) div += (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - 1].x) / VEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (VEC<VType>::quantity[idx + 1].x - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].x) / (2 * VEC<VType>::h.x);
			}
			else {

				div += 0.5 * (VEC<VType>::quantity[idx + 1].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].x - VEC<VType>::quantity[idx - 1].x) / (2 * VEC<VType>::h.x);
			}
			else {

				div += 0.5 * (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - 1].x) / VEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPY) div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y + VEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).y) / (2 * VEC<VType>::h.y);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).y - VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.y;
		if (ngbrFlags[idx] & NF_NNY) div += (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / VEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
			}
			else {

				div += 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
			}
			else {

				div += 0.5 * (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / VEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z + VEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).z) / (2 * VEC<VType>::h.z);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).z - VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.z;
		if (ngbrFlags[idx] & NF_NNZ) div += (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / VEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
			}
			else {

				div += 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
			}
			else {

				div += 0.5 * (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / VEC<VType>::h.z;
			}
		}
	}

	return div;
}

//divergence operator. Use Dirichlet conditions if set, else non-homogeneous Neumann boundary conditions.
//NOTE : the boundary differential is specified with 3 components, one for each of +x, +y, +z surface normal directions
//Can be used at composite media boundaries where sided differentials will be used instead.
//div operator can be applied if VType is a VAL3<Type>, returning Type
template double VEC_VC<DBL3>::div_diri_nneu(int idx, const VAL3<DBL3>& bdiff) const;

template <typename VType>
double VEC_VC<VType>::div_diri_nneu(int idx, const VAL3<VType>& bdiff) const
{
	double div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		div += (VEC<VType>::quantity[idx + 1].x - VEC<VType>::quantity[idx - 1].x) / (2 * VEC<VType>::h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPX) div += (VEC<VType>::quantity[idx + 1].x + VEC<VType>::quantity[idx].x - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx).x) / (2 * VEC<VType>::h.x);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx).x - VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - 1].x) / (2 * VEC<VType>::h.x);
	}
	//Not Dirichlet, is it a CMBND boundary? - if not this either then use homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) div += (VEC<VType>::quantity[idx + 1].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.x;
		if (ngbrFlags[idx] & NF_NNX) div += (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - 1].x) / VEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (VEC<VType>::quantity[idx + 1].x - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].x) / (2 * VEC<VType>::h.x);
			}
			else {

				div += 0.5 * ((VEC<VType>::quantity[idx + 1].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.x + bdiff.x.x);
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].x - VEC<VType>::quantity[idx - 1].x) / (2 * VEC<VType>::h.x);
			}
			else {

				div += 0.5 * ((VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - 1].x) / VEC<VType>::h.x + bdiff.x.x);
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPY) div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y + VEC<VType>::quantity[idx].y - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx).y) / (2 * VEC<VType>::h.y);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx).y - VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.y;
		if (ngbrFlags[idx] & NF_NNY) div += (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / VEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
			}
			else {

				div += 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.y + bdiff.y.y);
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
			}
			else {

				div += 0.5 * ((VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / VEC<VType>::h.y + bdiff.y.y);
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z + VEC<VType>::quantity[idx].z - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx).z) / (2 * VEC<VType>::h.z);
		else								 div += (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx).z - VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.z;
		if (ngbrFlags[idx] & NF_NNZ) div += (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / VEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
			}
			else {

				div += 0.5 * ((VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.z + bdiff.z.z);
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
			}
			else {

				div += 0.5 * ((VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / VEC<VType>::h.z + bdiff.z.z);
			}
		}
	}

	return div;
}

//divergence operator. Use sided differentials (also at composite media boundaries)
//div operator can be applied if VType is a VAL3<Type>, returning Type
template double VEC_VC<DBL3>::div_sided(int idx) const;

template <typename VType>
double VEC_VC<VType>::div_sided(int idx) const
{
	double div = 0.0;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return div;

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		div += (VEC<VType>::quantity[idx + 1].x - VEC<VType>::quantity[idx - 1].x) / (2 * VEC<VType>::h.x);
	}
	//use sided differentials if one of the neighbors is present
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (VEC<VType>::quantity[idx + 1].x - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1].x) / (2 * VEC<VType>::h.x);
			}
			else {

				div += (VEC<VType>::quantity[idx + 1].x - VEC<VType>::quantity[idx].x) / VEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				div += (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)].x - VEC<VType>::quantity[idx - 1].x) / (2 * VEC<VType>::h.x);
			}
			else {

				div += (VEC<VType>::quantity[idx].x - VEC<VType>::quantity[idx - 1].x) / VEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
			}
			else {

				div += (VEC<VType>::quantity[idx + VEC<VType>::n.x].y - VEC<VType>::quantity[idx].y) / VEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				div += (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / (2 * VEC<VType>::h.y);
			}
			else {

				div += (VEC<VType>::quantity[idx].y - VEC<VType>::quantity[idx - VEC<VType>::n.x].y) / VEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
			}
			else {

				div += (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx].z) / VEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				div += (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / (2 * VEC<VType>::h.z);
			}
			else {

				div += (VEC<VType>::quantity[idx].z - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y].z) / VEC<VType>::h.z;
			}
		}
	}

	return div;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

//-------------------------------- DIVERGENCE OPERATOR applied after multiplying with unit antisymmetric tensor (epsilon3)

//divergence operator of epsilon3(VEC<VType>::quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
//Use Neumann boundary conditions (homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
template DBL3 VEC_VC<DBL3>::diveps3_neu(int idx) const;
template FLT3 VEC_VC<FLT3>::diveps3_neu(int idx) const;

template <typename VType>
VType VEC_VC<VType>::diveps3_neu(int idx) const
{
	//differentials along x, y, z directions
	VType diffx, diffy, diffz;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		diffx = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
	}
	//Is it a CMBND boundary? - if not then use homogeneous Neumann condition (differential zero at the boundary)
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) diffx = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / VEC<VType>::h.x;
		if (ngbrFlags[idx] & NF_NNX) diffx = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / VEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx = (VEC<VType>::quantity[idx + 1] - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diffx = 0.5 * (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / VEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diffx = 0.5 * (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / VEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diffy = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) diffy = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / VEC<VType>::h.y;
		if (ngbrFlags[idx] & NF_NNY) diffy = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / VEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diffy = 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / VEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diffy = 0.5 * (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / VEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diffz = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) diffz = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / VEC<VType>::h.z;
		if (ngbrFlags[idx] & NF_NNZ) diffz = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / VEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diffz = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diffz = 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / VEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diffz = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diffz = 0.5 * (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / VEC<VType>::h.z;
			}
		}
	}

	return VType(
		-diffy.z + diffz.y,
		-diffz.x + diffx.z,
		-diffx.y + diffy.x
		);
}

//divergence operator of epsilon3(VEC<VType>::quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
//Use Dirichlet conditions if set, else Neumann boundary conditions(homogeneous).
//Can be used at composite media boundaries where sided differentials will be used instead.
template DBL3 VEC_VC<DBL3>::diveps3_diri(int idx) const;
template FLT3 VEC_VC<FLT3>::diveps3_diri(int idx) const;

template <typename VType>
VType VEC_VC<VType>::diveps3_diri(int idx) const
{
	//differentials along x, y, z directions
	VType diffx, diffy, diffz;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		diffx = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
	}
	//not an inner point along this direction - Use Dirichlet?
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETX)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPX) diffx = (VEC<VType>::quantity[idx + 1] + VEC<VType>::quantity[idx] - 2 * get_dirichlet_value(NF2_DIRICHLETPX, idx)) / (2 * VEC<VType>::h.x);
		else								   diffx = (2 * get_dirichlet_value(NF2_DIRICHLETNX, idx) - VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
	}
	//Not Dirichlet, is it a CMBND boundary? - if not this either then use homogeneous Neumann condition
	else if (ngbrFlags[idx] & NF_CMBNDX) {

		if (ngbrFlags[idx] & NF_NPX) diffx = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / VEC<VType>::h.x;
		if (ngbrFlags[idx] & NF_NNX) diffx = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / VEC<VType>::h.x;
	}
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx = (VEC<VType>::quantity[idx + 1] - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diffx = 0.5 * (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / VEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diffx = 0.5 * (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / VEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diffy = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETY)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPY) diffy = (VEC<VType>::quantity[idx + VEC<VType>::n.x] + VEC<VType>::quantity[idx] - 2 * get_dirichlet_value(NF2_DIRICHLETPY, idx)) / (2 * VEC<VType>::h.y);
		else								 diffy = (2 * get_dirichlet_value(NF2_DIRICHLETNY, idx) - VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_CMBNDY) {

		if (ngbrFlags[idx] & NF_NPY) diffy = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / VEC<VType>::h.y;
		if (ngbrFlags[idx] & NF_NNY) diffy = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / VEC<VType>::h.y;
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diffy = 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / VEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diffy = 0.5 * (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / VEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diffz = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags2.size() && (ngbrFlags2[idx] & NF2_DIRICHLETZ)) {

		if (ngbrFlags2[idx] & NF2_DIRICHLETPZ) diffz = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] + VEC<VType>::quantity[idx] - 2 * get_dirichlet_value(NF2_DIRICHLETPZ, idx)) / (2 * VEC<VType>::h.z);
		else								 diffz = (2 * get_dirichlet_value(NF2_DIRICHLETNZ, idx) - VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_CMBNDZ) {

		if (ngbrFlags[idx] & NF_NPZ) diffz = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / VEC<VType>::h.z;
		if (ngbrFlags[idx] & NF_NNZ) diffz = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / VEC<VType>::h.z;
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diffz = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diffz = 0.5 * (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / VEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diffz = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diffz = 0.5 * (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / VEC<VType>::h.z;
			}
		}
	}

	return VType(
		-diffy.z + diffz.y,
		-diffz.x + diffx.z,
		-diffx.y + diffy.x
	);
}

//divergence operator of epsilon3(VEC<VType>::quantity[idx]), i.e. multiply this VEC_VC by the unit antisymmetric tensor of rank3 before taking divergence. 
//Use sided differentials (also at composite media boundaries)
template DBL3 VEC_VC<DBL3>::diveps3_sided(int idx) const;
template FLT3 VEC_VC<FLT3>::diveps3_sided(int idx) const;

template <typename VType>
VType VEC_VC<VType>::diveps3_sided(int idx) const
{
	//differentials along x, y, z directions
	VType diffx, diffy, diffz;

	if (!(ngbrFlags[idx] & NF_NOTEMPTY)) return VType();

	//x direction
	if ((ngbrFlags[idx] & NF_BOTHX) == NF_BOTHX) {

		//inner point along this direction
		diffx = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
	}
	//use sided differentials if one of the neighbors is present
	else if (ngbrFlags[idx] & NF_NGBRX) {

		if (ngbrFlags[idx] & NF_NPX) {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx = (VEC<VType>::quantity[idx + 1] - get_sign(pbc_x) * VEC<VType>::quantity[idx + VEC<VType>::n.x - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diffx = (VEC<VType>::quantity[idx + 1] - VEC<VType>::quantity[idx]) / VEC<VType>::h.x;
			}
		}
		else {

			//is it a pbc along x? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCX) {

				diffx = (get_sign(pbc_x) * VEC<VType>::quantity[idx - (VEC<VType>::n.x - 1)] - VEC<VType>::quantity[idx - 1]) / (2 * VEC<VType>::h.x);
			}
			else {

				diffx = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - 1]) / VEC<VType>::h.x;
			}
		}
	}

	//y direction
	if ((ngbrFlags[idx] & NF_BOTHY) == NF_BOTHY) {

		diffy = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
	}
	else if (ngbrFlags[idx] & NF_NGBRY) {

		if (ngbrFlags[idx] & NF_NPY) {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - get_sign(pbc_y) * VEC<VType>::quantity[idx + (VEC<VType>::n.y - 1)*VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diffy = (VEC<VType>::quantity[idx + VEC<VType>::n.x] - VEC<VType>::quantity[idx]) / VEC<VType>::h.y;
			}
		}
		else {

			//is it a pbc along y? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCY) {

				diffy = (get_sign(pbc_y) * VEC<VType>::quantity[idx - (VEC<VType>::n.y - 1)*VEC<VType>::n.x] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / (2 * VEC<VType>::h.y);
			}
			else {

				diffy = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x]) / VEC<VType>::h.y;
			}
		}
	}

	//z direction
	if ((ngbrFlags[idx] & NF_BOTHZ) == NF_BOTHZ) {

		diffz = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
	}
	else if (ngbrFlags[idx] & NF_NGBRZ) {

		if (ngbrFlags[idx] & NF_NPZ) {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diffz = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - get_sign(pbc_z) * VEC<VType>::quantity[idx + (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diffz = (VEC<VType>::quantity[idx + VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx]) / VEC<VType>::h.z;
			}
		}
		else {

			//is it a pbc along z? If yes, then we are guaranteed to have a "neighbor" on the other side, so use it; otherwise apply boundary condition.
			if (ngbrFlags[idx] & NF_PBCZ) {

				diffz = (get_sign(pbc_z) * VEC<VType>::quantity[idx - (VEC<VType>::n.z - 1) * VEC<VType>::n.x*VEC<VType>::n.y] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / (2 * VEC<VType>::h.z);
			}
			else {

				diffz = (VEC<VType>::quantity[idx] - VEC<VType>::quantity[idx - VEC<VType>::n.x*VEC<VType>::n.y]) / VEC<VType>::h.z;
			}
		}
	}

	return VType(
		-diffy.z + diffz.y,
		-diffz.x + diffx.z,
		-diffx.y + diffy.x
	);
}