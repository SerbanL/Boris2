#pragma once

#include <cuda_runtime.h>
#include <cufft.h>

#include "cuFuncs_Aux.h"

#include "alloc_cpy.h"

#include "Types_ReIm.h"

////////////////////////////////////////////////////////////////////////////////////////////////// cuReIm (complex number)
//
// Complex number and associated operators

template <typename Type = void>
struct __cuReIm {

	//----------------------------- DATA

	cuBReal Re;
	cuBReal Im;

	//----------------------------- cu_obj MANAGED CONSTRUCTORS / DESTRUCTOR

	__host__ void construct_cu_obj(void)
	{
		set_gpu_value(Re, (cuBReal)0.0);
		set_gpu_value(Im, (cuBReal)0.0);
	}

	__host__ void construct_cu_obj(const __cuReIm& copyThis)
	{
		assign_cu_obj(copyThis);
	}

	__host__ void assign_cu_obj(const __cuReIm& copyThis)
	{
		gpu_to_gpu(Re, copyThis.Re);
		gpu_to_gpu(Im, copyThis.Im);
	}

	__host__ void destruct_cu_obj(void)
	{
	}

	//----------------------------- VALUE CONSTRUCTORS

	__host__ __device__ __cuReIm(void) { Re = 0; Im = 0; }
	__host__ __device__ __cuReIm(cuBReal _Re, cuBReal _Im) { Re = _Re; Im = _Im; }

	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	__host__ __device__ __cuReIm(const __cuReIm &copyThis) { Re = copyThis.Re; Im = copyThis.Im; }

	//assignment operator
	__host__ __device__ __cuReIm& operator=(const __cuReIm &rhs) { Re = rhs.Re; Im = rhs.Im; return *this; }

	//----------------------------- CONVERSION TO/FROM NON-CUDA VERSION
	
	__host__ operator ReIm() const
	{
		return ReIm((double)Re, (double)Im);
	}

	__host__ __cuReIm& operator=(const ReIm &rhs)
	{
		Re = (cuBReal)rhs.Re; 
		Im = (cuBReal)rhs.Im;
		return *this;
	}

	__host__ __cuReIm(const ReIm &rhs)
	{
		Re = (cuBReal)rhs.Re; 
		Im = (cuBReal)rhs.Im;
	}
	
	//----------------------------- CONVERSION TO/FROM cuBComplex

	__host__ __device__ operator cuBComplex() const
	{
		return { Re, Im };
	}

	__host__ __device__ __cuReIm& operator=(const cuBComplex &rhs)
	{
		Re = (cuBReal)rhs.x;
		Im = (cuBReal)rhs.y;
		return *this;
	}

	__host__ __device__ __cuReIm(const cuBComplex &rhs)
	{
		Re = (cuBReal)rhs.x;
		Im = (cuBReal)rhs.y;
	}

	//----------------------------- ARITHMETIC OPERATORS

	//+=, -=
	__host__ __device__ void operator+=(const __cuReIm &rhs) { Re += rhs.Re; Im += rhs.Im; }
	__host__ __device__ void operator-=(const __cuReIm &rhs) { Re -= rhs.Re; Im -= rhs.Im; }

	//complex conjugate
	__host__ __device__ __cuReIm operator~(void) const { return __cuReIm(Re, -Im); }

	//multiplication by i
	__host__ __device__ __cuReIm operator!(void) const { return __cuReIm(-Im, Re); }

	//multiplication of two complex numbers
	__host__ __device__ __cuReIm operator*(const __cuReIm &rhs) const { return __cuReIm(Re * rhs.Re - Im * rhs.Im, Re * rhs.Im + Im * rhs.Re); }

	//multiplication with a cuBComplex : cuBComplex on the RHS
	__host__ __device__ __cuReIm operator*(const cuBComplex &rhs) const { return __cuReIm(Re * rhs.x - Im * rhs.y, Re * rhs.y + Im * rhs.x); }

	//multiplication with a cuBComplex : cuBComplex on the LHS
	__host__ __device__ friend __cuReIm operator*(const cuBComplex &lhs, const __cuReIm &rhs) { return __cuReIm( lhs.x * rhs.Re - lhs.y * rhs.Im, lhs.x * rhs.Im + lhs.y * rhs.Re ); }

	//multiplication by a constant - constant on the right (must be fundamental)
	template <typename VType, std::enable_if_t<std::is_fundamental<VType>::value>* = nullptr>
	__host__ __device__ __cuReIm operator*(const VType &rhs) const { return __cuReIm(Re * (cuBReal)rhs, Im * (cuBReal)rhs); }

	//multiplication by a constant - constant on the left (must be fundamental)
	template <typename VType, std::enable_if_t<std::is_fundamental<VType>::value>* = nullptr>
	__host__ __device__ friend __cuReIm operator*(const VType &lhs, const __cuReIm &rhs) { return __cuReIm((cuBReal)lhs * rhs.Re, (cuBReal)lhs * rhs.Im); }

	//division by a constant (must be fundamental)
	template <typename VType, std::enable_if_t<std::is_fundamental<VType>::value>* = nullptr>
	__host__ __device__ __cuReIm operator/(const VType &divisor) const { return __cuReIm(Re / (cuBReal)divisor, Im / (cuBReal)divisor); }

	//addition
	__host__ __device__ __cuReIm operator+(const __cuReIm &rhs) const { return __cuReIm(Re + rhs.Re, Im + rhs.Im); }

	//subtraction
	__host__ __device__ __cuReIm operator-(const __cuReIm &rhs) const { return __cuReIm(Re - rhs.Re, Im - rhs.Im); }

	//----------------------------- COMPARISON OPERATORS

	//comparison
	__host__ __device__ bool operator==(const __cuReIm &rhs) const { return (IsZ(this->Re - rhs.Re) && IsZ(this->Im - rhs.Im)); }
	__host__ __device__ bool operator!=(const __cuReIm &rhs) const { return (IsNZ(this->Re - rhs.Re) || IsNZ(this->Im - rhs.Im)); }
};

typedef __cuReIm<void> cuReIm;

////////////////////////////////////////////////////////////////////////////////////////////////// __cuReIm3
//
//

//this is just a container of 3 cuReIm numbers (labelled conveniently x, y, z) with all operators extended so a __cuReIm3 object can be manipulated just as a cuReIm one
template <typename Type = void>
struct __cuReIm3 {

	//----------------------------- DATA

	cuReIm x, y, z;

	//----------------------------- cu_obj MANAGED CONSTRUCTORS / DESTRUCTOR

	__host__ void construct_cu_obj(void)
	{
		set_gpu_value(x, cuReIm());
		set_gpu_value(y, cuReIm());
		set_gpu_value(z, cuReIm());
	}

	__host__ void construct_cu_obj(const __cuReIm3& copyThis)
	{
		assign_cu_obj(copyThis);
	}

	__host__ void assign_cu_obj(const __cuReIm3& copyThis)
	{
		gpu_to_gpu(x, copyThis.x);
		gpu_to_gpu(y, copyThis.y);
		gpu_to_gpu(z, copyThis.z);
	}

	__host__ void destruct_cu_obj(void)
	{
	}

	//----------------------------- VALUE CONSTRUCTORS

	__host__ __device__ __cuReIm3(void) { x = cuReIm(); y = cuReIm(); z = cuReIm(); }

	__host__ __device__ __cuReIm3(cuReal3 Re, cuReal3 Im) { x = cuReIm(Re.x, Im.x); y = cuReIm(Re.y, Im.y); z = cuReIm(Re.z, Im.z); }

	__host__ __device__ __cuReIm3(cuReIm x, cuReIm y, cuReIm z) { this->x = x; this->y = y; this->z = z; }

	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	__host__ __device__ __cuReIm3(const __cuReIm3 &copyThis) { x = copyThis.x; y = copyThis.y; z = copyThis.z; }

	//assignment operator
	__host__ __device__ __cuReIm3& operator=(const __cuReIm3 &rhs) { x = rhs.x; y = rhs.y; z = rhs.z; return *this; }

	//----------------------------- CONVERSION TO/FROM NON-CUDA VERSION
	
	__host__ operator ReIm3() const
	{
		return ReIm3((ReIm)x, (ReIm)y, (ReIm)z);
	}

	__host__ __cuReIm3& operator=(const ReIm3 &rhs)
	{
		x = (cuReIm)rhs.x;
		y = (cuReIm)rhs.y;
		z = (cuReIm)rhs.z;
		return *this;
	}

	__host__ __cuReIm3(const ReIm3 &rhs)
	{
		x = (cuReIm)rhs.x;
		y = (cuReIm)rhs.y;
		z = (cuReIm)rhs.z;
	}
	
	//----------------------------- ARITHMETIC OPERATORS

	//complex conjugate
	__host__ __device__ __cuReIm3 operator~(void) const { return __cuReIm3(~x, ~y, ~z); }

	//multiplication by i
	__host__ __device__ __cuReIm3 operator!(void) const { return __cuReIm3(!x, !y, !z); }

	//multiplication of two complex numbers
	__host__ __device__ __cuReIm3 operator*(const __cuReIm3 &rhs) const { return __cuReIm3(x * rhs.x, y * rhs.y, z * rhs.z); }

	//multiplication by a ReIm (this means multiply each component by a ReIm) - ReIm on the right
	__host__ __device__ __cuReIm3 operator*(const cuReIm &rhs) const { return __cuReIm3(x * rhs, y * rhs, z * rhs); }

	//multiplication by a ReIm (this means multiply each component by a ReIm) - ReIm on the left
	__host__ __device__ friend __cuReIm3 operator*(const cuReIm &lhs, const __cuReIm3 &rhs) { return __cuReIm3(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z); }

	//multiplication by a constant - constant on the right (must be fundamental)
	template <typename VType, std::enable_if_t<std::is_fundamental<VType>::value>* = nullptr>
	__host__ __device__ __cuReIm3 operator*(const VType &rhs) const { return __cuReIm3(x * (cuBReal)rhs, y * (cuBReal)rhs, z * (cuBReal)rhs); }

	//multiplication by a constant - constant on the left (must be fundamental)
	template <typename VType, std::enable_if_t<std::is_fundamental<VType>::value>* = nullptr>
	__host__ __device__ friend __cuReIm3 operator*(const VType &lhs, const __cuReIm3 &rhs) { return __cuReIm3((cuBReal)lhs * rhs.x, (cuBReal)lhs * rhs.y, (cuBReal)lhs * rhs.z); }

	//division by a constant (must be fundamental)
	template <typename VType, std::enable_if_t<std::is_fundamental<VType>::value>* = nullptr>
	__host__ __device__ __cuReIm3 operator/(const VType &divisor) const { return __cuReIm3(x / (cuBReal)divisor, y / (cuBReal)divisor, z / (cuBReal)divisor); }

	//addition
	__host__ __device__ __cuReIm3 operator+(const __cuReIm3 &rhs) const { return __cuReIm3(x + rhs.x, y + rhs.y, z + rhs.z); }

	//subtraction
	__host__ __device__ __cuReIm3 operator-(const __cuReIm3 &rhs) const { return __cuReIm3(x - rhs.x, y - rhs.y, z - rhs.z); }

	//----------------------------- COMPARISON OPERATORS

	//comparison
	__host__ __device__ bool operator==(const __cuReIm3 &rhs) const { return (x == rhs.x && y == rhs.y && z == rhs.z); }
	__host__ __device__ bool operator!=(const __cuReIm3 &rhs) const { return (x != rhs.x || y != rhs.y || z != rhs.z); }

	//----------------------------- OTHERS

	//extract real part as a cuReal3
	__host__ __device__ cuReal3 Re(void) const { return cuReal3(x.Re, y.Re, z.Re); }

	//extract imaginary part as a cuReal3
	__host__ __device__ cuReal3 Im(void) const { return cuReal3(x.Im, y.Im, z.Im); }
};

typedef __cuReIm3<void> cuReIm3;