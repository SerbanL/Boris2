#pragma once

//Functions :

#include <cuda_runtime.h>

#include <type_traits>

#include "cuBRealComplex.h"

//If a floating point value is between epsilon negative and positive then it is considered to be zero (used for comparisons when we need to know if a value is zero or not)
#define CUEPSILON_POSITIVE 1e-11F	
#define CUEPSILON_NEGATIVE -1e-11F

//Fixed epsilon value used for floor and ceil functions with fixed precision.
#define CUEPSILON_ROUNDING	1e-4

///////////////////////////////////////////////////////////////////////////////
//COMPARISONS - Should not compare floating point values for equality directly (==, >= or <=). Just use one of the functions below instead.

//Is Zero

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
__host__ __device__ bool cuIsZ(bfloat fval) { return (fval < (bfloat)CUEPSILON_POSITIVE && fval >(bfloat)CUEPSILON_NEGATIVE); }

template <typename bint, std::enable_if_t<std::is_integral<bint>::value>* = nullptr>
__host__ __device__ bool cuIsZ(bint ival) { return (ival == 0); }

//Is Not Zero

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
__host__ __device__ bool cuIsNZ(bfloat fval) { return (fval > (bfloat)CUEPSILON_POSITIVE || fval < (bfloat)CUEPSILON_NEGATIVE); }

template <typename bint, std::enable_if_t<std::is_integral<bint>::value>* = nullptr>
__host__ __device__ bool cuIsNZ(bint ival) { return (ival != 0); }

//Is Zero or Positive

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
__host__ __device__ bool cuIsZoP(bfloat fval) { return (fval > 0 || (fval < (bfloat)CUEPSILON_POSITIVE && fval >(bfloat)CUEPSILON_NEGATIVE)); }

template <typename bint, std::enable_if_t<std::is_integral<bint>::value>* = nullptr>
__host__ __device__ bool cuIsZoP(bint ival) { return (ival >= 0); }

//Is Zero or Negative

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
__host__ __device__ bool cuIsZoN(bfloat fval) { return (fval < 0 || (fval < (bfloat)CUEPSILON_POSITIVE && fval >(bfloat)CUEPSILON_NEGATIVE)); }

template <typename bint, std::enable_if_t<std::is_integral<bint>::value>* = nullptr>
__host__ __device__ bool cuIsZoN(bint ival) { return (ival <= 0); }

//Is Equal

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
__host__ __device__ bool cuIsE(bfloat fval1, bfloat fval2) { return (fval2 - fval1 < (bfloat)CUEPSILON_POSITIVE && fval2 - fval1 >(bfloat)CUEPSILON_NEGATIVE); }

template <typename bint, std::enable_if_t<std::is_integral<bint>::value>* = nullptr>
__host__ __device__ bool cuIsE(bint ival1, bint ival2) { return (ival1 == ival2); }

//Is Greater or Equal

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
__host__ __device__ bool cuIsGE(bfloat fval1, bfloat fval2) { return (fval2 - fval1 < (bfloat)CUEPSILON_POSITIVE); }

template <typename bint, std::enable_if_t<std::is_integral<bint>::value>* = nullptr>
__host__ __device__ bool cuIsGE(bint ival1, bint ival2) { return (ival1 >= ival2); }

//Is Smaller or Equal

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
__host__ __device__ bool cuIsSE(bfloat fval1, bfloat fval2) { return (fval1 - fval2 < (bfloat)CUEPSILON_POSITIVE); }

template <typename bint, std::enable_if_t<std::is_integral<bint>::value>* = nullptr>
__host__ __device__ bool cuIsSE(bint ival1, bint ival2) { return (ival1 <= ival2); }

//Is Smaller - sometimes just checking fval1 < fval2 is not good enough : fval1 must also be at least an EPSILON_POSITIVE lower than fval2.

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
__host__ __device__ bool cuIsS(bfloat fval1, bfloat fval2) { return (fval1 < fval2 - (bfloat)CUEPSILON_POSITIVE); }

//Is Greater - sometimes just checking fval1 > fval2 is not good enough : fval1 must also be at least an EPSILON_POSITIVE greater than fval2.

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
__host__ __device__ bool cuIsG(bfloat fval1, bfloat fval2) { return (fval1 > fval2 + (bfloat)CUEPSILON_POSITIVE); }