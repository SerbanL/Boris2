#pragma once

#include "VEC_VC.h"
#include "VEC_VC_CGSolve.h"

//--------------------------------------------CONSTRUCTORS : VEC_VC_mng.h

template <typename VType>
VEC_VC<VType>::VEC_VC(void) :
	VEC(),
	ProgramStateNames(this, { VINFO(n), VINFO(h), VINFO(rect), VINFO(quantity), VINFO(ngbrFlags), VINFO(nonempty_cells),
		VINFO(dirichlet_nx), VINFO(dirichlet_px), VINFO(dirichlet_ny), VINFO(dirichlet_py), VINFO(dirichlet_nz), VINFO(dirichlet_pz),
		VINFO(robin_px), VINFO(robin_nx), VINFO(robin_py), VINFO(robin_ny), VINFO(robin_pz), VINFO(robin_nz), VINFO(robin_v), VINFO(shift_debt), VINFO(aSOR_damping) }, {})
{
}

template <typename VType>
VEC_VC<VType>::VEC_VC(const SZ3& n_) :
	VEC(n_),
	ProgramStateNames(this, { VINFO(n), VINFO(h), VINFO(rect), VINFO(quantity), VINFO(ngbrFlags), VINFO(nonempty_cells),
		VINFO(dirichlet_nx), VINFO(dirichlet_px), VINFO(dirichlet_ny), VINFO(dirichlet_py), VINFO(dirichlet_nz), VINFO(dirichlet_pz),
		VINFO(robin_px), VINFO(robin_nx), VINFO(robin_py), VINFO(robin_ny), VINFO(robin_pz), VINFO(robin_nz), VINFO(robin_v), VINFO(shift_debt), VINFO(aSOR_damping) }, {})
{
	if (!mreserve_vector(ngbrFlags, n_.dim()) || quantity.size() != n_.dim()) {

		clear();
	}
	else {

		resize_ngbrFlags(n);
	}
}

template <typename VType>
VEC_VC<VType>::VEC_VC(const DBL3& h_, const Rect& rect_) :
	VEC(h_, rect_),
	ProgramStateNames(this, { VINFO(n), VINFO(h), VINFO(rect), VINFO(quantity), VINFO(ngbrFlags), VINFO(nonempty_cells),
		VINFO(dirichlet_nx), VINFO(dirichlet_px), VINFO(dirichlet_ny), VINFO(dirichlet_py), VINFO(dirichlet_nz), VINFO(dirichlet_pz),
		VINFO(robin_px), VINFO(robin_nx), VINFO(robin_py), VINFO(robin_ny), VINFO(robin_pz), VINFO(robin_nz), VINFO(robin_v), VINFO(shift_debt), VINFO(aSOR_damping) }, {})
{
	//new_n : the size we should have if everything succeeds - n will take on this value if so
	SZ3 new_n = get_n_from_h_and_rect(h_, rect_);

	if (!mreserve_vector(ngbrFlags, new_n.dim()) || quantity.size() != new_n.dim()) {

		clear();
	}
	else {

		resize_ngbrFlags(n);
		set_ngbrFlags();
	}
}

template <typename VType>
VEC_VC<VType>::VEC_VC(const DBL3& h_, const Rect& rect_, VType value) :
	VEC(h_, rect_, value),
	ProgramStateNames(this, { VINFO(n), VINFO(h), VINFO(rect), VINFO(quantity), VINFO(ngbrFlags), VINFO(nonempty_cells),
		VINFO(dirichlet_nx), VINFO(dirichlet_px), VINFO(dirichlet_ny), VINFO(dirichlet_py), VINFO(dirichlet_nz), VINFO(dirichlet_pz),
		VINFO(robin_px), VINFO(robin_nx), VINFO(robin_py), VINFO(robin_ny), VINFO(robin_pz), VINFO(robin_nz), VINFO(robin_v), VINFO(shift_debt), VINFO(aSOR_damping) }, {})
{
	//new_n : the size we should have if everything succeeds - n will take on this value if so
	SZ3 new_n = get_n_from_h_and_rect(h_, rect_);

	if (!mreserve_vector(ngbrFlags, new_n.dim()) || quantity.size() != new_n.dim()) {

		clear();
	}
	else {

		resize_ngbrFlags(n);
		set_ngbrFlags();
	}
}

//--------------------------------------------SIZING

//resize and set shape using linked vec
template <typename VType>
template <typename LVType>
bool VEC_VC<VType>::resize(const SZ3& new_n, VEC_VC<LVType> &linked_vec)
{
	//reserve memory for ngbrFlags at the new size
	if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

	//memory reserved so map flags to new size before changing the n value
	resize_ngbrFlags(new_n);

	if (!VEC::resize(new_n)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec);

	return true;
}

//resize but keep shape
template <typename VType>
bool VEC_VC<VType>::resize(const SZ3& new_n)
{
	//reserve memory for ngbrFlags at the new size
	if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

	//memory reserved so map flags to new size before changing the n value
	resize_ngbrFlags(new_n);

	if (!VEC::resize(new_n)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags();

	return true;
}

//resize and set shape using linked vec
template <typename VType>
template <typename LVType>
bool VEC_VC<VType>::resize(const DBL3& new_h, const Rect& new_rect, VEC_VC<LVType> &linked_vec)
{
	SZ3 new_n = get_n_from_h_and_rect(new_h, new_rect);

	//reserve memory for ngbrFlags at the new size
	if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

	//memory reserved so map flags to new size before changing the n value
	resize_ngbrFlags(new_n);

	if (!VEC::resize(new_h, new_rect)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec);

	return true;
}

//resize but keep shape
template <typename VType>
bool VEC_VC<VType>::resize(const DBL3& new_h, const Rect& new_rect)
{
	SZ3 new_n = get_n_from_h_and_rect(new_h, new_rect);

	//reserve memory for ngbrFlags at the new size
	if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

	//memory reserved so map flags to new size before changing the n value
	resize_ngbrFlags(new_n);

	if (!VEC::resize(new_h, new_rect)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags();

	return true;
}

//set value and shape from linked vec
template <typename VType>
template <typename LVType>
bool VEC_VC<VType>::assign(const SZ3& new_n, VType value, VEC_VC<LVType> &linked_vec)
{
	//reserve memory for ngbrFlags at the new size
	if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

	//memory reserved so map flags to new size before changing the n value
	resize_ngbrFlags(new_n);

	if (!VEC::assign(new_n, value)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec);

	return true;
}

//set value but keep shape - empty cells will retain zero value : i.e. set value everywhere but in empty cells
template <typename VType>
bool VEC_VC<VType>::assign(const SZ3& new_n, VType value)
{
	//reserve memory for ngbrFlags at the new size
	if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

	//memory reserved so map flags to new size before changing the n value
	resize_ngbrFlags(new_n);

	if (!VEC::assign(new_n, value)) return false;

	//all good, finish off by setting flags (noting that empty cells are set back to zero)
	set_ngbrFlags();

	return true;
}

//set value and shape from linked vec
template <typename VType>
template <typename LVType>
bool VEC_VC<VType>::assign(const DBL3& new_h, const Rect& new_rect, VType value, VEC_VC<LVType> &linked_vec)
{
	SZ3 new_n = get_n_from_h_and_rect(new_h, new_rect);

	//reserve memory for ngbrFlags at the new size
	if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

	//memory reserved so map flags to new size before changing the n value
	resize_ngbrFlags(new_n);

	if (!VEC::assign(new_h, new_rect, value)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec);

	return true;
}

//set value but keep shape - empty cells will retain zero value : i.e. set value everywhere but in empty cells
template <typename VType>
bool VEC_VC<VType>::assign(const DBL3& new_h, const Rect& new_rect, VType value)
{
	SZ3 new_n = get_n_from_h_and_rect(new_h, new_rect);

	//reserve memory for ngbrFlags at the new size
	if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

	//memory reserved so map flags to new size before changing the n value
	resize_ngbrFlags(new_n);

	if (!VEC::assign(new_h, new_rect, value)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags();

	return true;
}

template <typename VType>
void VEC_VC<VType>::clear(void)
{
	VEC::clear();
	ngbrFlags.clear();
	ngbrFlags.shrink_to_fit();

	clear_dirichlet_flags();

	shift_debt = DBL3();

	aSOR_damping = 1.0;
}

template <typename VType>
void VEC_VC<VType>::shrink_to_fit(void)
{
	VEC::shrink_to_fit();
	ngbrFlags.shrink_to_fit();

	dirichlet_px.shrink_to_fit();
	dirichlet_nx.shrink_to_fit();
	dirichlet_py.shrink_to_fit();
	dirichlet_ny.shrink_to_fit();
	dirichlet_pz.shrink_to_fit();
	dirichlet_nz.shrink_to_fit();
}