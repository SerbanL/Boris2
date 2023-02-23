#pragma once

#include "VEC.h"

//NOTE : very basic out-of-place transpose routines. Just loop tiling (blocking) done to improve cache use. When you have time you can look at better algorithms (in-place).

#define BLOCK	4

//out-of-place xy transpose
template <typename VType>
template <typename SType>
void VEC<VType>::transpose_xy(VEC<SType>& out)
{
	if (out.linear_size() != n.dim()) return;

	for (int k = 0; k < n.z; k++) {

		#pragma omp parallel for
		for (int jb = 0; jb < n.y; jb += BLOCK) {
			for (int ib = 0; ib < n.x; ib += BLOCK) {

				for (int j = jb; j < (jb + BLOCK < n.y ? jb + BLOCK : n.y); j++) {
					for (int i = ib; i < (ib + BLOCK < n.x ? ib + BLOCK : n.x); i++) {

						int idx_src = i + j * n.x + k * n.x*n.y;

						int idx_dst = j + i * n.y + k * n.y*n.x;

						out[idx_dst] = quantity[idx_src];
					}
				}

			}
		}
	}
}

//out-of-place xz transpose
template <typename VType>
template <typename SType>
void VEC<VType>::transpose_xz(VEC<SType>& out)
{
	if (out.linear_size() != n.dim()) return;
	
	for (int j = 0; j < n.y; j++) {

		#pragma omp parallel for
		for (int kb = 0; kb < n.z; kb += BLOCK) {
			for (int ib = 0; ib < n.x; ib += BLOCK) {

				for (int k = kb; k < (kb + BLOCK < n.z ? kb + BLOCK : n.z); k++) {
					for (int i = ib; i < (ib + BLOCK < n.x ? ib + BLOCK : n.x); i++) {

						int idx_src = i + j * n.x + k * n.x*n.y;

						int idx_dst = k + j * n.z + i * n.z*n.y;

						out[idx_dst] = quantity[idx_src];
					}
				}

			}
		}
	}
}

//out-of-place yz transpose
template <typename VType>
template <typename SType>
void VEC<VType>::transpose_yz(VEC<SType>& out)
{
	if (out.linear_size() != n.dim()) return;

	for (int k = 0; k < n.z; k++) {

		#pragma omp parallel for
		for (int jb = 0; jb < n.y; jb += BLOCK) {
			for (int ib = 0; ib < n.x; ib += BLOCK) {

				for (int j = jb; j < (jb + BLOCK < n.y ? jb + BLOCK : n.y); j++) {
					for (int i = ib; i < (ib + BLOCK < n.x ? ib + BLOCK : n.x); i++) {

						int idx_src = i + j * n.x + k * n.x*n.y;

						int idx_dst = i + k * n.x + j * n.x*n.z;

						out[idx_dst] = quantity[idx_src];
					}
				}

			}
		}
	}
}

//out-of-place cyclic transpose up : xyz to zxy
template <typename VType>
template <typename SType>
void VEC<VType>::transpose_cycleup(VEC<SType>& out)
{
	if (out.linear_size() != n.dim()) return;

	for (int j = 0; j < n.y; j++) {

		#pragma omp parallel for
		for (int kb = 0; kb < n.z; kb += BLOCK) {
			for (int ib = 0; ib < n.x; ib += BLOCK) {

				for (int k = kb; k < (kb + BLOCK < n.z ? kb + BLOCK : n.z); k++) {
					for (int i = ib; i < (ib + BLOCK < n.x ? ib + BLOCK : n.x); i++) {

						int idx_src = i + j * n.x + k * n.x*n.y;

						int idx_dst = k + i * n.z + j * n.z*n.x;

						out[idx_dst] = quantity[idx_src];
					}
				}

			}
		}
	}

}

//out-of-place cyclic transpose up : xyz to yzx
template <typename VType>
template <typename SType>
void VEC<VType>::transpose_cycledn(VEC<SType>& out)
{
	if (out.linear_size() != n.dim()) return;

	for (int k = 0; k < n.z; k++) {

		#pragma omp parallel for
		for (int jb = 0; jb < n.y; jb += BLOCK) {
			for (int ib = 0; ib < n.x; ib += BLOCK) {

				for (int j = jb; j < (jb + BLOCK < n.y ? jb + BLOCK : n.y); j++) {
					for (int i = ib; i < (ib + BLOCK < n.x ? ib + BLOCK : n.x); i++) {

						int idx_src = i + j * n.x + k * n.x*n.y;

						int idx_dst = j + k * n.y + i * n.y*n.z;

						out[idx_dst] = quantity[idx_src];
					}
				}

			}
		}
	}
}