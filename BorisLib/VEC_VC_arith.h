#pragma once

#include "VEC_VC.h"

//scale all stored values by the given constant
template <typename VType>
void VEC_VC<VType>::scale_values(double constant)
{
#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		quantity[idx] *= constant;
	}
}

//add values from add_this but keep current dimensions - if necessary map values from add_this to local dimensions
template <typename VType>
void VEC_VC<VType>::add_values(const VEC_VC<VType>& copy_this)
{
	VEC<VType>::add_values(copy_this);

	//copy shape

	//first clear current flags
	ngbrFlags.assign(n.dim(), 0);

	//now map shape from copy_this.ngbrFlags to ngbrFlags
	SZ3 source_n = copy_this.n;
	DBL3 sourceIdx = (DBL3)source_n / n;

#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		int _x = (int)floor((idx % n.x) * sourceIdx.x);
		int _y = (int)floor(((idx / n.x) % n.y) * sourceIdx.y);
		int _z = (int)floor((idx / (n.x*n.y)) * sourceIdx.z);

		if (copy_this.ngbrFlags[_x + _y * source_n.x + _z * (source_n.x*source_n.y)] & NF_NOTEMPTY)
			ngbrFlags[idx] = NF_NOTEMPTY;
	}

	//recalculate neighbor flags
	set_ngbrFlags();
}

//add values from add_this but keep current dimensions - if necessary map values from add_this to local dimensions
template <typename VType>
void VEC_VC<VType>::add_values(const VEC<VType>& copy_this)
{
	VEC<VType>::add_values(copy_this);

	//recalculate neighbor flags : current shape is maintained.
	set_ngbrFlags();
}

//subtract values from sub_this but keep current dimensions - if necessary map values from sub_this to local dimensions
template <typename VType>
void VEC_VC<VType>::sub_values(const VEC_VC<VType>& copy_this)
{
	VEC<VType>::sub_values(copy_this);

	//copy shape

	//first clear current flags
	ngbrFlags.assign(n.dim(), 0);

	//now map shape from copy_this.ngbrFlags to ngbrFlags
	SZ3 source_n = copy_this.n;
	DBL3 sourceIdx = (DBL3)source_n / n;

#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		int _x = (int)floor((idx % n.x) * sourceIdx.x);
		int _y = (int)floor(((idx / n.x) % n.y) * sourceIdx.y);
		int _z = (int)floor((idx / (n.x*n.y)) * sourceIdx.z);

		if (copy_this.ngbrFlags[_x + _y * source_n.x + _z * (source_n.x*source_n.y)] & NF_NOTEMPTY)
			ngbrFlags[idx] = NF_NOTEMPTY;
	}

	//recalculate neighbor flags
	set_ngbrFlags();
}

//add values from add_this but keep current dimensions - if necessary map values from add_this to local dimensions
template <typename VType>
void VEC_VC<VType>::sub_values(const VEC<VType>& copy_this)
{
	VEC<VType>::sub_values(copy_this);

	//recalculate neighbor flags : current shape is maintained.
	set_ngbrFlags();
}

template <typename VType>
void VEC_VC<VType>::operator+=(double constant)
{
#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		if (ngbrFlags[idx] & NF_NOTEMPTY) quantity[idx] += constant;
	}
}

template <typename VType>
void VEC_VC<VType>::operator-=(double constant)
{
#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		if (ngbrFlags[idx] & NF_NOTEMPTY) quantity[idx] -= constant;
	}
}

template <typename VType>
void VEC_VC<VType>::operator*=(double constant)
{
#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		if (ngbrFlags[idx] & NF_NOTEMPTY) quantity[idx] *= constant;
	}
}

template <typename VType>
void VEC_VC<VType>::operator/=(double constant)
{
#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		if (ngbrFlags[idx] & NF_NOTEMPTY) quantity[idx] /= constant;
	}
}
