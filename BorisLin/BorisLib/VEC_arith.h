#pragma once

#include "VEC.h"

//add values from add_this but keep current dimensions - if necessary map values from add_this to local dimensions
template <typename VType>
void VEC<VType>::add_values(const VEC<VType>& add_this)
{
	SZ3 source_n = add_this.n;
	
	if (!source_n.dim()) return;

	DBL3 sourceIdx = (DBL3)source_n / n;

#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		int _x = (int)floor((idx % n.x) * sourceIdx.x);
		int _y = (int)floor(((idx / n.x) % n.y) * sourceIdx.y);
		int _z = (int)floor((idx / (n.x*n.y)) * sourceIdx.z);

		quantity[idx] += add_this[_x + _y * source_n.x + _z * (source_n.x*source_n.y)];
	}
}

//subtract values from sub_this but keep current dimensions - if necessary map values from sub_this to local dimensions
template <typename VType>
void VEC<VType>::sub_values(const VEC<VType>& sub_this)
{
	SZ3 source_n = sub_this.n;

	if (!source_n.dim()) return;

	DBL3 sourceIdx = (DBL3)source_n / n;

#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		int _x = (int)floor((idx % n.x) * sourceIdx.x);
		int _y = (int)floor(((idx / n.x) % n.y) * sourceIdx.y);
		int _z = (int)floor((idx / (n.x*n.y)) * sourceIdx.z);

		quantity[idx] -= sub_this[_x + _y * source_n.x + _z * (source_n.x*source_n.y)];
	}
}

template <typename VType>
void VEC<VType>::operator+=(double constant)
{
#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		quantity[idx] += constant;
	}
}

template <typename VType>
void VEC<VType>::operator-=(double constant)
{
#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		quantity[idx] -= constant;
	}
}

template <typename VType>
void VEC<VType>::operator*=(double constant)
{
#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		quantity[idx] *= constant;
	}
}

template <typename VType>
void VEC<VType>::operator/=(double constant)
{
#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		quantity[idx] /= constant;
	}
}