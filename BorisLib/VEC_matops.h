#pragma once

#include "VEC.h"

//--------------------------------------------MULTIPLICATIONS

//matrix multiplication using naive algorithm
//multiply xy planes of lvec and rvec considered as 2D matrices and place result in this VEC. Require lvec.n.x = rvec.n.y and lvec.n.z = rvec.n.z. Output matrix (this) sized as required.
template <>
inline void VEC<double>::matrix_mul(const VEC<double>& lvec, const VEC<double>& rvec)
{
	//check if we have correct dimensions for matrix multiplication
	if (lvec.n.x != rvec.n.y || lvec.n.z != rvec.n.z) return;

	//resize matrix if needed and set output matrix to zero
	assign(SZ3(rvec.n.x, lvec.n.y, lvec.n.z), 0.0);

	int length = lvec.n.x;

	for (int k = 0; k < n.z; k++) {

#pragma omp parallel for
		for (int j = 0; j < n.y; j++) {
			for (int i = 0; i < n.x; i++) {

				double value = 0;

				for (int idx = 0; idx < length; idx++) {

					int lidx = idx + j * lvec.n.x + k * lvec.n.x * lvec.n.y;
					int ridx = i + idx * rvec.n.x + k * rvec.n.x * rvec.n.y;

					value += lvec[lidx] * rvec[ridx];
				}

				quantity[i + j * n.x + k * n.x*n.y] = value;
			}
		}
	}
}

//multiply matrix by floating point constant
template <typename VType>
void VEC<VType>::matrix_mul(double constant)
{
#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		quantity[idx] = quantity[idx] * constant;
	}
}

//multiply diagonal values in each xy plane of this VEC (considered as a matrix) by value.
//If not square in the xy plane the "diagonal" starts at (0,0) and has min(n.x, n.y) points
template <typename VType>
void VEC<VType>::matrix_muldiag(double value)
{
	int length = minimum(n.x, n.y);

	for (int k = 0; k < n.z; k++) {
		
		#pragma omp parallel for
		for (int idx = 0; idx < length; idx++) {

			quantity[idx + idx * n.x + k * n.x*n.y] *= value;
		}
	}
}

//--------------------------------------------ADDITION and SUBTRACTION

//add matadd into this matrix point by point - sizes must match
template <typename VType>
void VEC<VType>::matrix_add(const VEC<VType>& matadd)
{
	//sizes must match
	if (n != matadd.n) assign(matadd.n, VType());

#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		quantity[idx] += matadd[idx];
	}
}

//add lvec and rvec (sizes must match) point by point, setting output in this matrix
template <typename VType>
void VEC<VType>::matrix_add(const VEC<VType>& lvec, const VEC<VType>& rvec)
{
	//sizes must match
	if (lvec.n != rvec.n) return;

	if(n != lvec.n) assign(lvec.n, VType());

#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		quantity[idx] = lvec[idx] + rvec[idx];
	}
}

//subtract matadd from this matrix point by point - sizes must match
template <typename VType>
void VEC<VType>::matrix_sub(const VEC<VType>& matadd)
{
	//sizes must match
	if (n != matadd.n) assign(matadd.n, VType());

#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		quantity[idx] -= matadd[idx];
	}
}

//subtract rvec from lvec (sizes must match) point by point, setting output in this matrix
template <typename VType>
void VEC<VType>::matrix_sub(const VEC<VType>& lvec, const VEC<VType>& rvec)
{
	//sizes must match
	if (lvec.n != rvec.n) return;

	if (n != lvec.n) assign(lvec.n, VType());

#pragma omp parallel for
	for (int idx = 0; idx < n.dim(); idx++) {

		quantity[idx] = rvec[idx] - lvec[idx];
	}
}

//--------------------------------------------INVERSE

//Invert each plane of this VEC considered as a matrix (must be square in xy plane) and return determinant of first matrix (first xy plane) - using algorithm from : A. Farooq, K. Hamid, "An Efficient and Simple Algorithm for Matrix Inversion" IJTD, 1, 20 (2010)
template <>
inline double VEC<double>::matrix_inverse(void)
{
	//must be a square input matrix
	if (n.x != n.y) return 0.0;

	int size = n.x;

	//Step 1 : start - setup determinant start value and pointers to vectors
	double d = 1;

	//swap space
	VEC<double> swap(SZ3(size, size, n.z));

	//setup pointers
	VEC<double> *pInput, *pOutput, *pTemp;
	pInput = this;
	pOutput = &swap;

	//Step 2: diagonal elements loop
	for (int k = 0; k < n.z; k++) {

		for (int p = 0; p < size; p++) {

			//(p,p) diagonal element
			double a_pp = (*pInput)[p + size * p + k*size*size];

			//Step 3: if zero cannot find inverse: finish and return 0 to signal inverse has failed
			if (!a_pp) return 0.0;

			//Step 4: update determinant
			d = d * a_pp;

#pragma omp parallel for
			for (int j = 0; j < size; j++) {

				if (j == p) continue;

				//Step 5: pivot row
				(*pOutput)[p*size + j + k * size*size] = (*pInput)[p*size + j + k * size*size] / a_pp;

				//Step 6: pivot column
				(*pOutput)[p + j * size + k * size*size] = -(*pInput)[p + j * size + k * size*size] / a_pp;
			}

			//Step 7: calculate remaining elements
#pragma omp parallel for
			for (int i = 0; i < size; i++) {

				if (i == p) continue;

				for (int j = 0; j < size; j++) {

					if (j == p) continue;

					(*pOutput)[j + i * size + k * size*size] = (*pInput)[j + i * size + k * size*size] + (*pInput)[j + p * size + k * size*size] * (*pOutput)[p + i * size + k * size*size];
				}
			}

			(*pOutput)[p + p * size + k * size*size] = 1 / a_pp;

			//Rotate 
			pTemp = pInput;
			pInput = pOutput;
			pOutput = pTemp;
		}
	}

	//For odd dimension the inverted matrix is in the swap space
	if (size % 2) {

#pragma omp parallel for
		for (int idx = 0; idx < n.dim(); idx++) {

			quantity[idx] = swap[idx];
		}
	}

	return d;
}

//--------------------------------------------OTHERS

//extract values from given xy plane (plane ranges from 0 to n.z - 1) diagonal into a std::vector
//If not square in the xy plane the "diagonal" starts at (0,0) and has min(n.x, n.y) points
template <typename VType>
void VEC<VType>::matrix_getdiagonal(std::vector<VType>& diagonal, int plane)
{
	if (!GoodIdx(n.z, plane)) return;

	int length = minimum(n.x, n.y);

	diagonal.resize(length);

#pragma omp parallel for
	for (int idx = 0; idx < length; idx++) {

		diagonal[idx] = quantity[idx + idx * n.x + plane * n.x*n.y];
	}
}