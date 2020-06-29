#pragma once

#include "VEC_VC.h"
#include "VEC_VC_CGSolve.h"

//--------------------------------------------CONSTRUCTORS : VEC_VC_mng.h

template <typename VType>
VEC_VC<VType>::VEC_VC(void) :
	VEC<VType>(),
	ProgramState<VEC_VC<VType>, 
	std::tuple<SZ3, DBL3, Rect, std::vector<VType>, std::vector<int>, std::vector<int>, int,
	std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>,
	DBL2, DBL2, DBL2, DBL2, DBL2, DBL2, DBL2, DBL3, int, int, int>,
	std::tuple<>>
	(this, 
	{ VINFO(VEC<VType>::n), VINFO(VEC<VType>::h), VINFO(VEC<VType>::rect), VINFO(VEC<VType>::quantity), VINFO(ngbrFlags), VINFO(ngbrFlags2), VINFO(nonempty_cells),
	  VINFO(dirichlet_nx), VINFO(dirichlet_px), VINFO(dirichlet_ny), VINFO(dirichlet_py), VINFO(dirichlet_nz), VINFO(dirichlet_pz),
	  VINFO(robin_px), VINFO(robin_nx), VINFO(robin_py), VINFO(robin_ny), VINFO(robin_pz), VINFO(robin_nz), VINFO(robin_v), VINFO(shift_debt),
	  VINFO(pbc_x), VINFO(pbc_y), VINFO(pbc_z) }, {})
{
}

template <typename VType>
VEC_VC<VType>::VEC_VC(const SZ3& n_) :
	VEC<VType>(n_),
	ProgramState<VEC_VC<VType>,
	std::tuple<SZ3, DBL3, Rect, std::vector<VType>, std::vector<int>, std::vector<int>, int,
	std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>,
	DBL2, DBL2, DBL2, DBL2, DBL2, DBL2, DBL2, DBL3, int, int, int>,
	std::tuple<>>
	(this,
		{ VINFO(VEC<VType>::n), VINFO(VEC<VType>::h), VINFO(VEC<VType>::rect), VINFO(VEC<VType>::quantity), VINFO(ngbrFlags), VINFO(ngbrFlags2), VINFO(nonempty_cells),
		  VINFO(dirichlet_nx), VINFO(dirichlet_px), VINFO(dirichlet_ny), VINFO(dirichlet_py), VINFO(dirichlet_nz), VINFO(dirichlet_pz),
		  VINFO(robin_px), VINFO(robin_nx), VINFO(robin_py), VINFO(robin_ny), VINFO(robin_pz), VINFO(robin_nz), VINFO(robin_v), VINFO(shift_debt),
		  VINFO(pbc_x), VINFO(pbc_y), VINFO(pbc_z) }, {})
{
	if (!mreserve_vector(ngbrFlags, n_.dim()) || VEC<VType>::quantity.size() != n_.dim()) {

		clear();
	}
	else {

		resize_ngbrFlags(VEC<VType>::n);
	}
}

template <typename VType>
VEC_VC<VType>::VEC_VC(const DBL3& h_, const Rect& rect_) :
	VEC<VType>(h_, rect_),
	ProgramState<VEC_VC<VType>,
	std::tuple<SZ3, DBL3, Rect, std::vector<VType>, std::vector<int>, std::vector<int>, int,
	std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>,
	DBL2, DBL2, DBL2, DBL2, DBL2, DBL2, DBL2, DBL3, int, int, int>,
	std::tuple<>>
	(this,
		{ VINFO(VEC<VType>::n), VINFO(VEC<VType>::h), VINFO(VEC<VType>::rect), VINFO(VEC<VType>::quantity), VINFO(ngbrFlags), VINFO(ngbrFlags2), VINFO(nonempty_cells),
		  VINFO(dirichlet_nx), VINFO(dirichlet_px), VINFO(dirichlet_ny), VINFO(dirichlet_py), VINFO(dirichlet_nz), VINFO(dirichlet_pz),
		  VINFO(robin_px), VINFO(robin_nx), VINFO(robin_py), VINFO(robin_ny), VINFO(robin_pz), VINFO(robin_nz), VINFO(robin_v), VINFO(shift_debt),
		  VINFO(pbc_x), VINFO(pbc_y), VINFO(pbc_z) }, {})
{
	//new_n : the size we should have if everything succeeds - VEC<VType>::n will take on this value if so
	SZ3 new_n = VEC<VType>::get_n_from_h_and_rect(h_, rect_);

	if (!mreserve_vector(ngbrFlags, new_n.dim()) || VEC<VType>::quantity.size() != new_n.dim()) {

		clear();
	}
	else {

		resize_ngbrFlags(VEC<VType>::n);
		set_ngbrFlags();
	}
}

template <typename VType>
VEC_VC<VType>::VEC_VC(const DBL3& h_, const Rect& rect_, VType value) :
	VEC<VType>(h_, rect_, value),
	ProgramState<VEC_VC<VType>,
	std::tuple<SZ3, DBL3, Rect, std::vector<VType>, std::vector<int>, std::vector<int>, int,
	std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>, std::vector<VType>,
	DBL2, DBL2, DBL2, DBL2, DBL2, DBL2, DBL2, DBL3, int, int, int>,
	std::tuple<>>
	(this,
		{ VINFO(VEC<VType>::n), VINFO(VEC<VType>::h), VINFO(VEC<VType>::rect), VINFO(VEC<VType>::quantity), VINFO(ngbrFlags), VINFO(ngbrFlags2), VINFO(nonempty_cells),
		  VINFO(dirichlet_nx), VINFO(dirichlet_px), VINFO(dirichlet_ny), VINFO(dirichlet_py), VINFO(dirichlet_nz), VINFO(dirichlet_pz),
		  VINFO(robin_px), VINFO(robin_nx), VINFO(robin_py), VINFO(robin_ny), VINFO(robin_pz), VINFO(robin_nz), VINFO(robin_v), VINFO(shift_debt),
		  VINFO(pbc_x), VINFO(pbc_y), VINFO(pbc_z) }, {})
{
	//new_n : the size we should have if everything succeeds - VEC<VType>::n will take on this value if so
	SZ3 new_n = VEC<VType>::get_n_from_h_and_rect(h_, rect_);

	if (!mreserve_vector(ngbrFlags, new_n.dim()) || VEC<VType>::quantity.size() != new_n.dim()) {

		clear();
	}
	else {

		resize_ngbrFlags(VEC<VType>::n);
		set_ngbrFlags();
	}
}

//--------------------------------------------SIZING

//resize and set shape using linked vec
template <typename VType>
template <typename LVType>
bool VEC_VC<VType>::resize(const SZ3& new_n, const VEC_VC<LVType> &linked_vec)
{
	if (new_n != VEC<VType>::n) {

		//reserve memory for ngbrFlags at the new size
		if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

		//reserve memory for extended ngbrFlags if needed
		if (use_extended_flags()) {

			if (!mreserve_vector(ngbrFlags2, new_n.dim())) return false;
		}

		//memory reserved so map flags to new size before changing the VEC<VType>::n value
		resize_ngbrFlags(new_n);

		if (!VEC<VType>::resize(new_n)) return false;
	}

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec);

	return true;
}

//resize but keep shape
template <typename VType>
bool VEC_VC<VType>::resize(const SZ3& new_n)
{
	if (new_n != VEC<VType>::n) {

		//reserve memory for ngbrFlags at the new size
		if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

		//reserve memory for extended ngbrFlags if needed
		if (use_extended_flags()) {

			if (!mreserve_vector(ngbrFlags2, new_n.dim())) return false;
		}

		//memory reserved so map flags to new size before changing the VEC<VType>::n value
		resize_ngbrFlags(new_n);

		if (!VEC<VType>::resize(new_n)) return false;

		//all good, finish off by setting flags.
		set_ngbrFlags();
	}

	return true;
}

//resize and set shape using linked vec
template <typename VType>
template <typename LVType>
bool VEC_VC<VType>::resize(const DBL3& new_h, const Rect& new_rect, const VEC_VC<LVType> &linked_vec)
{
	if (new_h != VEC<VType>::h || new_rect != VEC<VType>::rect) {

		SZ3 new_n = VEC<VType>::get_n_from_h_and_rect(new_h, new_rect);

		//reserve memory for ngbrFlags at the new size
		if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

		//reserve memory for extended ngbrFlags if needed
		if (use_extended_flags()) {

			if (!mreserve_vector(ngbrFlags2, new_n.dim())) return false;
		}

		//memory reserved so map flags to new size before changing the VEC<VType>::n value
		resize_ngbrFlags(new_n);

		if (!VEC<VType>::resize(new_h, new_rect)) return false;
	}

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec);

	return true;
}

//resize but keep shape
template <typename VType>
bool VEC_VC<VType>::resize(const DBL3& new_h, const Rect& new_rect)
{
	if (new_h != VEC<VType>::h || new_rect != VEC<VType>::rect) {

		SZ3 new_n = VEC<VType>::get_n_from_h_and_rect(new_h, new_rect);

		//reserve memory for ngbrFlags at the new size
		if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

		//reserve memory for extended ngbrFlags if needed
		if (use_extended_flags()) {

			if (!mreserve_vector(ngbrFlags2, new_n.dim())) return false;
		}

		//memory reserved so map flags to new size before changing the VEC<VType>::n value
		resize_ngbrFlags(new_n);

		if (!VEC<VType>::resize(new_h, new_rect)) return false;

		//all good, finish off by setting flags.
		set_ngbrFlags();
	}

	return true;
}

//set value and shape from linked vec
template <typename VType>
template <typename LVType>
bool VEC_VC<VType>::assign(const SZ3& new_n, VType value, const VEC_VC<LVType> &linked_vec)
{
	if (new_n != VEC<VType>::n) {

		//reserve memory for ngbrFlags at the new size
		if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

		//reserve memory for extended ngbrFlags if needed
		if (use_extended_flags()) {

			if (!mreserve_vector(ngbrFlags2, new_n.dim())) return false;
		}

		//memory reserved so map flags to new size before changing the VEC<VType>::n value
		resize_ngbrFlags(new_n);
	}

	if (!VEC<VType>::assign(new_n, value)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec);

	return true;
}

//set value but keep shape - empty cells will retain zero value : i.e. set value everywhere but in empty cells
template <typename VType>
bool VEC_VC<VType>::assign(const SZ3& new_n, VType value)
{
	if (new_n != VEC<VType>::n) {

		//reserve memory for ngbrFlags at the new size
		if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

		//reserve memory for extended ngbrFlags if needed
		if (use_extended_flags()) {

			if (!mreserve_vector(ngbrFlags2, new_n.dim())) return false;
		}

		//memory reserved so map flags to new size before changing the VEC<VType>::n value
		resize_ngbrFlags(new_n);
	}

	if (!VEC<VType>::assign(new_n, value)) return false;

	//all good, finish off by setting flags (noting that empty cells are set back to zero)
	set_ngbrFlags();

	return true;
}

//set value and shape from linked vec
template <typename VType>
template <typename LVType>
bool VEC_VC<VType>::assign(const DBL3& new_h, const Rect& new_rect, VType value, const VEC_VC<LVType> &linked_vec)
{
	if (new_h != VEC<VType>::h || new_rect != VEC<VType>::rect) {

		SZ3 new_n = VEC<VType>::get_n_from_h_and_rect(new_h, new_rect);

		//reserve memory for ngbrFlags at the new size
		if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

		//reserve memory for extended ngbrFlags if needed
		if (use_extended_flags()) {

			if (!mreserve_vector(ngbrFlags2, new_n.dim())) return false;
		}

		//memory reserved so map flags to new size before changing the VEC<VType>::n value
		resize_ngbrFlags(new_n);
	}

	if (!VEC<VType>::assign(new_h, new_rect, value)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags(linked_vec);

	return true;
}

//set value but keep shape - empty cells will retain zero value : i.e. set value everywhere but in empty cells
template <typename VType>
bool VEC_VC<VType>::assign(const DBL3& new_h, const Rect& new_rect, VType value)
{
	if (new_h != VEC<VType>::h || new_rect != VEC<VType>::rect) {

		SZ3 new_n = VEC<VType>::get_n_from_h_and_rect(new_h, new_rect);

		//reserve memory for ngbrFlags at the new size
		if (!mreserve_vector(ngbrFlags, new_n.dim())) return false;

		//reserve memory for extended ngbrFlags if needed
		if (use_extended_flags()) {

			if (!mreserve_vector(ngbrFlags2, new_n.dim())) return false;
		}

		//memory reserved so map flags to new size before changing the VEC<VType>::n value
		resize_ngbrFlags(new_n);
	}

	if (!VEC<VType>::assign(new_h, new_rect, value)) return false;

	//all good, finish off by setting flags.
	set_ngbrFlags();

	return true;
}

template <typename VType>
void VEC_VC<VType>::clear(void)
{
	VEC<VType>::clear();

	ngbrFlags.clear();
	ngbrFlags.shrink_to_fit();

	ngbrFlags2.clear();
	ngbrFlags2.shrink_to_fit();

	clear_dirichlet_flags();

	clear_pbc();

	shift_debt = DBL3();
}

template <typename VType>
void VEC_VC<VType>::shrink_to_fit(void)
{
	VEC<VType>::shrink_to_fit();
	ngbrFlags.shrink_to_fit();
	ngbrFlags2.shrink_to_fit();

	dirichlet_px.shrink_to_fit();
	dirichlet_nx.shrink_to_fit();
	dirichlet_py.shrink_to_fit();
	dirichlet_ny.shrink_to_fit();
	dirichlet_pz.shrink_to_fit();
	dirichlet_nz.shrink_to_fit();
}