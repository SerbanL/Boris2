#pragma once

#include <cuda_runtime.h>

#include "cuFuncs_Aux.h"

#include "Types_VAL.h"

#include "alloc_cpy.h"

////////////////////////////////////////////////////////////////////////////////////////////////// cuVAL2
//
//

template <typename VType> 
struct cuVAL2 {

	//----------------------------- DATA

	union {

		VType i, x, major, first;
	};

	union {

		VType j, y, minor, second;
	};

	//----------------------------- cu_obj MANAGED CONSTRUCTORS / DESTRUCTOR

	__host__ void construct_cu_obj(void)
	{
		set_gpu_value(x, VType());
		set_gpu_value(y, VType());
	}

	__host__ void construct_cu_obj(const cuVAL2& copyThis)
	{
		assign_cu_obj(copyThis);
	}

	__host__ void assign_cu_obj(const cuVAL2& copyThis)
	{
		gpu_to_gpu(x, copyThis.x);
		gpu_to_gpu(y, copyThis.y);
	}

	__host__ void destruct_cu_obj(void)
	{
	}

	//----------------------------- VALUE CONSTRUCTORS

	__host__ __device__ cuVAL2(void) { i = VType(); j = VType(); }

	__host__ __device__ cuVAL2(VType i_) : i(i_), j(i_) {}

	__host__ __device__ cuVAL2(VType i_, VType j_) : i(i_), j(j_) {}

	//----------------------------- CONVERTING CONSTRUCTORS

	//type conversion constructor
	template <typename CVType> 
	__host__ __device__ cuVAL2(const cuVAL2<CVType> &convThis)
	{ 
		x = (VType)convThis.x; 
		y = (VType)convThis.y; 
	}

	//copy constructor
	__host__ __device__ cuVAL2(const cuVAL2 &copyThis) { *this = copyThis; }

	//assignment operator
	__host__ __device__ cuVAL2& operator=(const cuVAL2 &rhs)
	{ 
		x = rhs.x; 
		y = rhs.y; 
		return *this;
	}

	//----------------------------- CONVERSION TO/FROM NON-CUDA VERSION

	template <typename SType>
	__host__ operator VAL2<SType>() const
	{
		return VAL2<SType>((SType)x, (SType)y);
	}

	template <typename SType>
	__host__ cuVAL2<VType>& operator=(const VAL2<SType> &rhs)
	{
		x = (VType)rhs.x; 
		y = (VType)rhs.y;
		return *this;
	}

	template <typename SType>
	__host__ cuVAL2(const VAL2<SType> &rhs)
	{
		x = (VType)rhs.x; 
		y = (VType)rhs.y;
	}

	//----------------------------- ARITHMETIC OPERATORS

	//For binary operators if the second operand uses a floating point and the first is integral type, the return type favours the second; otherwise the return type favours the first operand.

	//SUM
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	__host__ __device__ cuVAL2<RType> operator+(const cuVAL2<VType_> &rhs) const { return cuVAL2<RType>(x + rhs.x, y + rhs.y); }

	//DIFFERENCE
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	__host__ __device__ cuVAL2<RType> operator-(const cuVAL2<VType_> &rhs) const { return cuVAL2<RType>(x - rhs.x, y - rhs.y); }

	//SCALAR PRODUCT : (same fundamental types, e.g. DBL2 with a DBL2 -> double)
	template <
		typename VType_,
		std::enable_if_t<std::is_same<VType, VType_>::value>* = nullptr
	>
	__host__ __device__ cuBReal operator*(const cuVAL2<VType_> &rhs) const { return cuBReal(x * rhs.x + y * rhs.y); }

	//SIMPLE PRODUCT (component by component)
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	__host__ __device__ cuVAL2<RType> operator&(const cuVAL2<VType_> &rhs) const { return cuVAL2<RType>(x * rhs.x, y * rhs.y); }

	//DIVISION REMAINDER
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
		__host__ __device__ cuVAL2<RType> operator%(const cuVAL2<VType_> &rhs) const 
	{ 
		return cuVAL2<RType>(
			(x < rhs.x ? (x < 0.0 ? fmod(x, rhs.x) + rhs.x : x) : fmod(x, rhs.x)),
			(y < rhs.y ? (y < 0.0 ? fmod(y, rhs.y) + rhs.y : y) : fmod(y, rhs.y)));
	}

	//PRODUCTS WITH A CONSTANT

	//product with a constant (must be fundamental type) on the RHS
	template <
		typename MVType,
		typename RType = typename std::conditional<std::is_floating_point<MVType>::value && std::is_integral<VType>::value, MVType, VType>::type,
		std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr
	>
	__host__ __device__ cuVAL2<RType> operator*(const MVType &mult) const { return cuVAL2<RType>(x * mult, y * mult); }

	//product with a constant (must be fundamental type) on the LHS : floating point
	template <
		typename MVType,
		std::enable_if_t<std::is_floating_point<MVType>::value>* = nullptr
	>
	__host__ __device__ friend cuVAL2<MVType> operator*(const MVType &mult, const cuVAL2<VType> &rhs) { return cuVAL2<MVType>(rhs.x * mult, rhs.y * mult); }

	//product with a constant (must be fundamental type) on the LHS : integral
	template <
		typename MVType,
		std::enable_if_t<std::is_integral<MVType>::value>* = nullptr
	>
	__host__ __device__ friend cuVAL2<VType> operator*(const MVType &mult, const cuVAL2<VType> &rhs) { return cuVAL2<VType>(rhs.x * mult, rhs.y * mult); }

	//SIMPLE DIVISION (component by component)
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	__host__ __device__ cuVAL2<RType> operator/(const cuVAL2<VType_> &rhs) const { return cuVAL2<RType>(x / rhs.x, y / rhs.y); }

	//DIVISION BY A CONSTANT (must be fundamental type)
	template <
		class DVType,
		typename RType = typename std::conditional<std::is_floating_point<DVType>::value && std::is_integral<VType>::value, DVType, VType>::type,
		std::enable_if_t<std::is_fundamental<DVType>::value>* = nullptr
	>
	__host__ __device__ cuVAL2<RType> operator/(const DVType &divisor) const { return cuVAL2<RType>(x / divisor, y / divisor); }

	//OPERATION-ASSIGN : ADD
	template <typename VType_>
	__host__ __device__ void operator+=(const cuVAL2<VType_> &rhs) { x += rhs.x; y += rhs.y; }

	//OPERATION-ASSIGN : SUBTRACT
	template <typename VType_>
	__host__ __device__ void operator-=(const cuVAL2<VType_> &rhs) { x -= rhs.x; y -= rhs.y; }

	//OPERATION-ASSIGN : PRODUCT WITH CONSTANT (must be fundamental type)
	template <
		typename MVType,
		std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr
	>
	__host__ __device__ void operator*=(const MVType &mult) { x *= mult; y *= mult; }

	//OPERATION-ASSIGN : DIVISION BY CONSTANT (must be fundamental type)
	template <
		typename DVType,
		std::enable_if_t<std::is_fundamental<DVType>::value>* = nullptr
	>
	__host__ __device__ void operator/=(const DVType &divisor) { x /= divisor; y /= divisor; }

	//----------------------------- OTHER NUMERIC OPERATORS

	//the norm or magnitude
	__host__ __device__ VType norm(void) const { return sqrt(x * x + y * y); }

	//get this cuVAL3 in normalized form - not checking for zero magnitude here this must be ensured externally
	__host__ __device__ cuVAL2 normalized(void) const { return *this / norm(); }

	//set new magnitude (norm) - not checking for zero magnitude here this must be ensured externally
	__host__ __device__ void renormalize(VType new_norm) { *this *= new_norm / norm(); }

	//----------------------------- COMPARISON OPERATORS

	//comparison
	__host__ __device__ bool operator==(const cuVAL2 &rhs) const { return (cuIsZ(x - rhs.x) && cuIsZ(y - rhs.y)); }
	__host__ __device__ bool operator!=(const cuVAL2 &rhs) const { return (cuIsNZ(x - rhs.x) || cuIsNZ(y - rhs.y)); }

	//ordering
	__host__ __device__ bool operator>=(const cuVAL2 &rhs) const { return (cuIsZoP(x - rhs.x) && cuIsZoP(y - rhs.y)); }
	__host__ __device__ bool operator<=(const cuVAL2 &rhs) const { return (cuIsZoN(x - rhs.x) && cuIsZoN(y - rhs.y)); }
	__host__ __device__ bool operator>(const cuVAL2 &rhs) const { return ((x > rhs.x) && (y > rhs.y)); }
	__host__ __device__ bool operator<(const cuVAL2 &rhs) const { return ((x < rhs.x) && (y < rhs.y)); }

	__host__ __device__ bool IsNull(void) const { return (*this == cuVAL2()); }

	//----------------------------- OTHER

	__host__ __device__ VType dim(void) const { return x * y; }
};


typedef cuVAL2<int> cuINT2;
typedef cuVAL2<float> cuFLT2;
typedef cuVAL2<cuBReal> cuReal2;
typedef cuVAL2<double> cuDBL2;

////////////////////////////////////////////////////////////////////////////////////////////////// cuVAL3 and special cases INT3, FLT3, DBL3
//
//

template <typename VType> 
struct cuVAL3 {
	
	//----------------------------- DATA

	union {

		VType i, x;
	};

	union {

		VType j, y;
	};

	union {

		VType k, z;
	};

	//----------------------------- cu_obj MANAGED CONSTRUCTORS / DESTRUCTOR

	__host__ void construct_cu_obj(void) 
	{
		set_gpu_value(x, VType());
		set_gpu_value(y, VType());
		set_gpu_value(z, VType());
	}

	__host__ void construct_cu_obj(const cuVAL3& copyThis)
	{
		assign_cu_obj(copyThis);
	}

	__host__ void assign_cu_obj(const cuVAL3& copyThis)
	{
		gpu_to_gpu(x, copyThis.x);
		gpu_to_gpu(y, copyThis.y);
		gpu_to_gpu(z, copyThis.z);
	}

	__host__ void destruct_cu_obj(void)
	{
	}

	//----------------------------- VALUE CONSTRUCTORS

	__host__ __device__ cuVAL3(void) { x = VType(); y = VType(); z = VType(); }
	__host__ __device__ cuVAL3(VType val) { x = val; y = val; z = val; }
	__host__ __device__ cuVAL3(VType x, VType y, VType z) { this->x = x; this->y = y; this->z = z; }

	//----------------------------- CONVERTING CONSTRUCTORS

	//type conversion constructor
	template <typename CVType> 
	__host__ __device__ cuVAL3(const cuVAL3<CVType> &convThis) { x = (VType)convThis.x; y = (VType)convThis.y; z = (VType)convThis.z; }

	//copy constructor
	__host__ __device__ cuVAL3(const cuVAL3 &copyThis) { x = copyThis.x; y = copyThis.y; z = copyThis.z; }

	//assignment operator
	__host__ __device__ cuVAL3& operator=(const cuVAL3 &rhs) { x = rhs.x; y = rhs.y; z = rhs.z; return *this; }

	//----------------------------- CONVERSION TO/FROM NON-CUDA VERSION

	template <typename SType>
	__host__ operator VAL3<SType>() const
	{
		return VAL3<SType>((SType)x, (SType)y, (SType)z);
	}

	template <typename SType>
	__host__ cuVAL3<VType>& operator=(const VAL3<SType> &rhs)
	{
		x = (VType)rhs.x; 
		y = (VType)rhs.y; 
		z = (VType)rhs.z;
		return *this;
	}

	template <typename SType>
	__host__ cuVAL3(const VAL3<SType> &rhs)
	{
		x = (VType)rhs.x; 
		y = (VType)rhs.y; 
		z = (VType)rhs.z;
	}

	//----------------------------- ARITHMETIC OPERATORS

	//For binary operators if the second operand uses a floating point and the first is integral type, the return type favours the second; otherwise the return type favours the first operand.

	//SUM
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	__host__ __device__ cuVAL3<RType> operator+(const cuVAL3<VType_> &rhs) const { return cuVAL3<RType>(x + rhs.x, y + rhs.y, z + rhs.z); }

	//DIFFERENCE
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	__host__ __device__ cuVAL3<RType> operator-(const cuVAL3<VType_> &rhs) const { return cuVAL3<RType>(x - rhs.x, y - rhs.y, z - rhs.z); }

	//SCALAR PRODUCT : (same fundamental types, e.g. DBL3 with a DBL3 -> double)
	template <
		typename VType_,
		std::enable_if_t<std::is_same<VType, VType_>::value && std::is_fundamental<VType>::value>* = nullptr
	>
	__host__ __device__ cuBReal operator*(const cuVAL3<VType_> &rhs) const { return cuBReal(x * rhs.x + y * rhs.y + z * rhs.z); }

	//MATRIX PRODUCT : (same types but not fundamental, e.g. DBL33 with a DBL33 -> DBL33)
	template <
		typename VType_,
		std::enable_if_t<std::is_same<VType, VType_>::value && !std::is_fundamental<VType>::value>* = nullptr
	>
	__host__ __device__ cuVAL3<VType> operator*(const cuVAL3<VType_> &rhs) const
	{
		return cuVAL3<VType>(
			VType(i.i * rhs.i.i + i.j * rhs.j.i + i.k * rhs.k.i, i.i * rhs.i.j + i.j * rhs.j.j + i.k * rhs.k.j, i.i * rhs.i.k + i.j * rhs.j.k + i.k * rhs.k.k),
			VType(j.i * rhs.i.i + j.j * rhs.j.i + j.k * rhs.k.i, j.i * rhs.i.j + j.j * rhs.j.j + j.k * rhs.k.j, j.i * rhs.i.k + j.j * rhs.j.k + j.k * rhs.k.k),
			VType(k.i * rhs.i.i + k.j * rhs.j.i + k.k * rhs.k.i, k.i * rhs.i.j + k.j * rhs.j.j + k.k * rhs.k.j, k.i * rhs.i.k + k.j * rhs.j.k + k.k * rhs.k.k));
	}

	//MATRIX PRODUCT WITH VECTOR : (DBL33 * DBL3 -> DBL3)
	template <
		typename VType_,
		std::enable_if_t<!std::is_fundamental<VType>::value && std::is_fundamental<VType_>::value>* = nullptr
	>
	__host__ __device__ cuVAL3<VType_> operator*(const cuVAL3<VType_> &rhs) const { return cuVAL3<VType_>(x * rhs, y * rhs, z * rhs); }

	//VECTOR PRODUCT
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	__host__ __device__ cuVAL3<RType> operator^(const cuVAL3<VType_> &rhs) const { return cuVAL3<RType>(y*rhs.z - z * rhs.y, z*rhs.x - x * rhs.z, x*rhs.y - y * rhs.x); }

	//OUTER VECTOR PRODUCT
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type,
		std::enable_if_t<std::is_fundamental<VType_>::value && std::is_fundamental<VType>::value>* = nullptr
	>
	__host__ __device__ cuVAL3<cuVAL3<RType>> operator|(const cuVAL3<VType_> &rhs) const { return cuVAL3<cuVAL3<RType>>(cuVAL3<RType>(x * rhs.x, x * rhs.y, x * rhs.z), cuVAL3<RType>(y * rhs.x, y * rhs.y, y * rhs.z), cuVAL3<RType>(z * rhs.x, z * rhs.y, z * rhs.z)); }

	//MATRIX COLUMN PRODUCT WITH VECTOR : (not same types, but one is a fundamental type, e.g. DBL33 | DBL3 -> DBL3, where DBL33 | DBL3 = transpose(DBL33) * DBL3 )
	template <
		typename VType_,
		std::enable_if_t<!std::is_fundamental<VType>::value && std::is_fundamental<VType_>::value>* = nullptr
	>
	__host__ __device__ cuVAL3<VType_> operator|(const cuVAL3<VType_> &rhs) const { return cuVAL3<VType_>(x * rhs.x + y * rhs.y + z * rhs.z); }

	//SIMPLE PRODUCT (component by component)
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	__host__ __device__ cuVAL3<RType> operator&(const cuVAL3<VType_> &rhs) const { return cuVAL3<RType>(x * rhs.x, y * rhs.y, z * rhs.z); }

	//DIVISION REMAINDER
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
		__host__ __device__ cuVAL3<RType> operator%(const cuVAL3<VType_> &rhs) const 
	{ 
		return cuVAL3<RType>(
			(x < rhs.x ? (x < 0.0 ? fmod(x, rhs.x) + rhs.x : x) : fmod(x, rhs.x)),
			(y < rhs.y ? (y < 0.0 ? fmod(y, rhs.y) + rhs.y : y) : fmod(y, rhs.y)),
			(z < rhs.z ? (z < 0.0 ? fmod(z, rhs.z) + rhs.z : z) : fmod(z, rhs.z)));
	}

	//PRODUCTS WITH A CONSTANT

	//product with a constant (must be fundamental type) on the RHS
	template <
		typename MVType,
		std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr
	>
	__host__ __device__ cuVAL3<VType> operator*(const MVType &mult) const { return cuVAL3<VType>(x * mult, y * mult, z * mult); }

	//product with a constant (must be fundamental type) on the LHS
	template <
		typename MVType,
		std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr
	>
	__host__ __device__ friend cuVAL3<VType> operator*(const MVType &mult, const cuVAL3<VType> &rhs) { return cuVAL3<VType>(rhs.x * mult, rhs.y * mult, rhs.z * mult); }

	//SIMPLE DIVISION (component by component)
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	__host__ __device__ cuVAL3<RType> operator/(const cuVAL3<VType_> &rhs) const { return cuVAL3<RType>(x / rhs.x, y / rhs.y, z / rhs.z); }

	//DIVISION BY A CONSTANT (must be fundamental type)
	template <
		class DVType,
		std::enable_if_t<std::is_fundamental<DVType>::value>* = nullptr
	>
	__host__ __device__ cuVAL3<VType> operator/(const DVType &divisor) const { return cuVAL3<VType>(x / divisor, y / divisor, z / divisor); }

	//OPERATION-ASSIGN : ADD
	template <typename VType_>
	__host__ __device__ void operator+=(const cuVAL3<VType_> &rhs) { x += rhs.x; y += rhs.y; z += rhs.z; }

	//OPERATION-ASSIGN : SUBTRACT
	template <typename VType_>
	__host__ __device__ void operator-=(const cuVAL3<VType_> &rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; }

	//OPERATION-ASSIGN : PRODUCT WITH CONSTANT (must be fundamental type)
	template <
		typename MVType,
		std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr
	>
	__host__ __device__ void operator*=(const MVType &mult) { x *= mult; y *= mult; z *= mult; }

	//OPERATION-ASSIGN : DIVISION BY CONSTANT (must be fundamental type)
	template <
		typename DVType,
		std::enable_if_t<std::is_fundamental<DVType>::value>* = nullptr
	>
	__host__ __device__ void operator/=(const DVType &divisor) { x /= divisor; y /= divisor; z /= divisor; }

	//----------------------------- OTHER NUMERIC OPERATORS

	//the norm or magnitude
	__host__ __device__ VType norm(void) const { return sqrt(x * x + y * y + z * z); }

	//get this cuVAL3 in normalized form - not checking for zero magnitude here this must be ensured externally
	__host__ __device__ cuVAL3 normalized(void) const { return *this / norm(); }

	//set new magnitude (norm) - not checking for zero magnitude here this must be ensured externally
	__host__ __device__ void renormalize(VType new_norm) { *this *= new_norm / norm(); }

	//----------------------------- COMPARISON OPERATORS

	//comparison
	__host__ __device__ bool operator==(const cuVAL3 &rhs) const { return (cuIsZ(x - rhs.x) && cuIsZ(y - rhs.y) && cuIsZ(z - rhs.z)); }
	__host__ __device__ bool operator!=(const cuVAL3 &rhs) const { return (cuIsNZ(x - rhs.x) || cuIsNZ(y - rhs.y) || cuIsNZ(z - rhs.z)); }

	//ordering operators
	__host__ __device__ bool operator>=(const cuVAL3 &rhs) const { return (cuIsZoP(x - rhs.x) && cuIsZoP(y - rhs.y) && cuIsZoP(z - rhs.z)); }
	__host__ __device__ bool operator<=(const cuVAL3 &rhs) const { return (cuIsZoN(x - rhs.x) && cuIsZoN(y - rhs.y) && cuIsZoN(z - rhs.z)); }
	__host__ __device__ bool operator>(const cuVAL3 &rhs) const { return ((x > rhs.x) && (y > rhs.y) && (z > rhs.z)); }
	__host__ __device__ bool operator<(const cuVAL3 &rhs) const { return ((x < rhs.x) && (y < rhs.y) && (z < rhs.z)); }

	__host__ __device__ bool IsNull(void) const { return (*this == cuVAL3()); }

	//----------------------------- OTHER

	__host__ __device__ VType dim(void) const { return x * y * z; }
};

typedef cuVAL3<int> cuINT3;
typedef cuVAL3<size_t> cuSZ3;

typedef cuVAL3<float> cuFLT3;
typedef cuVAL3<cuBReal> cuReal3;
typedef cuVAL3<double> cuDBL3;

typedef cuVAL3<cuINT3> cuINT33;

typedef cuVAL3<cuFLT3> cuFLT33;
typedef cuVAL3<cuReal3> cuReal33;
typedef cuVAL3<cuDBL3> cuDBL33;

////////////////////////////////////////////////////////////////////////////////////////////////// cuVAL4 and special cases INT4, FLT4, DBL4
//
//

template <typename VType>
struct cuVAL4 {

	//----------------------------- DATA

	union {

		VType i, x;
	};

	union {

		VType j, y;
	};

	union {

		VType k, z;
	};

	union {

		VType l, t;
	};

	//----------------------------- cu_obj MANAGED CONSTRUCTORS / DESTRUCTOR

	__host__ void construct_cu_obj(void)
	{
		set_gpu_value(x, VType());
		set_gpu_value(y, VType());
		set_gpu_value(z, VType());
		set_gpu_value(t, VType());
	}

	__host__ void construct_cu_obj(const cuVAL4& copyThis)
	{
		assign_cu_obj(copyThis);
	}

	__host__ void assign_cu_obj(const cuVAL4& copyThis)
	{
		gpu_to_gpu(x, copyThis.x);
		gpu_to_gpu(y, copyThis.y);
		gpu_to_gpu(z, copyThis.z);
		gpu_to_gpu(t, copyThis.t);
	}

	__host__ void destruct_cu_obj(void)
	{
	}

	//----------------------------- VALUE CONSTRUCTORS

	__host__ __device__ cuVAL4(void) { x = VType(); y = VType(); z = VType(); t = VType(); }
	__host__ __device__ cuVAL4(VType val) { x = val; y = val; z = val; t = val; }
	__host__ __device__ cuVAL4(VType x, VType y, VType z, VType t) { this->x = x; this->y = y; this->z = z; this->t = t; }

	//----------------------------- CONVERTING CONSTRUCTORS

	//type conversion constructor
	template <typename CVType>
	__host__ __device__ cuVAL4(const cuVAL4<CVType> &convThis) { x = (VType)convThis.x; y = (VType)convThis.y; z = (VType)convThis.z; t = (VType)convThis.t; }

	//copy constructor
	__host__ __device__ cuVAL4(const cuVAL4 &copyThis) { x = copyThis.x; y = copyThis.y; z = copyThis.z; t = copyThis.t; }

	//assignment operator
	__host__ __device__ cuVAL4& operator=(const cuVAL4 &rhs) { x = rhs.x; y = rhs.y; z = rhs.z; t = rhs.t; return *this; }

	//----------------------------- CONVERSION TO/FROM NON-CUDA VERSION

	template <typename SType>
	__host__ operator VAL4<SType>() const
	{
		return VAL4<SType>((SType)x, (SType)y, (SType)z, (SType)t);
	}

	template <typename SType>
	__host__ cuVAL4<VType>& operator=(const VAL4<SType> &rhs)
	{
		x = (VType)rhs.x; y = (VType)rhs.y; z = (VType)rhs.z; t = (VType)rhs.t;
		return *this;
	}

	template <typename SType>
	__host__ cuVAL4(const VAL4<SType> &rhs)
	{
		x = (VType)rhs.x; y = (VType)rhs.y; z = (VType)rhs.z; t = (VType)rhs.t;
	}

	//----------------------------- ARITHMETIC OPERATORS

	//For binary operators if the second operand uses a floating point and the first is integral type, the return type favours the second; otherwise the return type favours the first operand.

	//SUM
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
		__host__ __device__ cuVAL4<RType> operator+(const cuVAL4<VType_> &rhs) const { return cuVAL4<RType>(x + rhs.x, y + rhs.y, z + rhs.z, t + rhs.t); }

	//DIFFERENCE
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
		__host__ __device__ cuVAL4<RType> operator-(const cuVAL4<VType_> &rhs) const { return cuVAL4<RType>(x - rhs.x, y - rhs.y, z - rhs.z, t - rhs.t); }

	//SCALAR PRODUCT : (same fundamental types, e.g. DBL4 with a DBL4 -> double)
	template <
		typename VType_,
		std::enable_if_t<std::is_same<VType, VType_>::value && std::is_fundamental<VType>::value>* = nullptr
	>
		__host__ __device__ cuBReal operator*(const cuVAL4<VType_> &rhs) const { return cuBReal(x * rhs.x + y * rhs.y + z * rhs.z + t * rhs.t); }

	//DIVISION REMAINDER
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
		__host__ __device__ cuVAL4<RType> operator%(const cuVAL4<VType_> &rhs) const 
	{ 
		return cuVAL4<RType>(
			(x < rhs.x ? (x < 0.0 ? fmod(x, rhs.x) + rhs.x : x) : fmod(x, rhs.x)),
			(y < rhs.y ? (y < 0.0 ? fmod(y, rhs.y) + rhs.y : y) : fmod(y, rhs.y)),
			(z < rhs.z ? (z < 0.0 ? fmod(z, rhs.z) + rhs.z : z) : fmod(z, rhs.z)),
			(t < rhs.t ? (t < 0.0 ? fmod(t, rhs.t) + rhs.t : t) : fmod(t, rhs.t)));
	}

	//PRODUCTS WITH A CONSTANT

	//product with a constant (must be fundamental type) on the RHS
	template <
		typename MVType,
		std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr
	>
		__host__ __device__ cuVAL4<VType> operator*(const MVType &mult) const { return cuVAL4<VType>(x * mult, y * mult, z * mult, t * mult); }

	//product with a constant (must be fundamental type) on the LHS
	template <
		typename MVType,
		std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr
	>
		__host__ __device__ friend cuVAL4<VType> operator*(const MVType &mult, const cuVAL4<VType> &rhs) { return cuVAL4<VType>(rhs.x * mult, rhs.y * mult, rhs.z * mult, rhs.t * mult); }

	//SIMPLE DIVISION (component by component)
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
		__host__ __device__ cuVAL4<RType> operator/(const cuVAL4<VType_> &rhs) const { return cuVAL4<RType>(x / rhs.x, y / rhs.y, z / rhs.z, t / rhs.t); }

	//DIVISION BY A CONSTANT (must be fundamental type)
	template <
		class DVType,
		std::enable_if_t<std::is_fundamental<DVType>::value>* = nullptr
	>
		__host__ __device__ cuVAL4<VType> operator/(const DVType &divisor) const { return cuVAL4<VType>(x / divisor, y / divisor, z / divisor, t / divisor); }

	//OPERATION-ASSIGN : ADD
	template <typename VType_>
	__host__ __device__ void operator+=(const cuVAL4<VType_> &rhs) { x += rhs.x; y += rhs.y; z += rhs.z; t += rhs.t; }

	//OPERATION-ASSIGN : SUBTRACT
	template <typename VType_>
	__host__ __device__ void operator-=(const cuVAL4<VType_> &rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; t -= rhs.t; }

	//OPERATION-ASSIGN : PRODUCT WITH CONSTANT (must be fundamental type)
	template <
		typename MVType,
		std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr
	>
		__host__ __device__ void operator*=(const MVType &mult) { x *= mult; y *= mult; z *= mult; t *= mult; }

	//OPERATION-ASSIGN : DIVISION BY CONSTANT (must be fundamental type)
	template <
		typename DVType,
		std::enable_if_t<std::is_fundamental<DVType>::value>* = nullptr
	>
		__host__ __device__ void operator/=(const DVType &divisor) { x /= divisor; y /= divisor; z /= divisor; t /= divisor; }

	//----------------------------- OTHER NUMERIC OPERATORS

	//the norm or magnitude
	__host__ __device__ VType norm(void) const { return sqrt(x * x + y * y + z * z + t * t); }

	//get this cuVAL4 in normalized form - not checking for zero magnitude here this must be ensured externally
	__host__ __device__ cuVAL4 normalized(void) const { return *this / norm(); }

	//set new magnitude (norm) - not checking for zero magnitude here this must be ensured externally
	__host__ __device__ void renormalize(VType new_norm) { *this *= new_norm / norm(); }

	//----------------------------- COMPARISON OPERATORS

	//comparison
	__host__ __device__ bool operator==(const cuVAL4 &rhs) const { return (cuIsZ(x - rhs.x) && cuIsZ(y - rhs.y) && cuIsZ(z - rhs.z) && cuIsZ(t - rhs.t)); }
	__host__ __device__ bool operator!=(const cuVAL4 &rhs) const { return (cuIsNZ(x - rhs.x) || cuIsNZ(y - rhs.y) || cuIsNZ(z - rhs.z) || cuIsNZ(t - rhs.t)); }

	//ordering operators
	__host__ __device__ bool operator>=(const cuVAL4 &rhs) const { return (cuIsZoP(x - rhs.x) && cuIsZoP(y - rhs.y) && cuIsZoP(z - rhs.z) && cuIsZoP(t - rhs.t)); }
	__host__ __device__ bool operator<=(const cuVAL4 &rhs) const { return (cuIsZoN(x - rhs.x) && cuIsZoN(y - rhs.y) && cuIsZoN(z - rhs.z) && cuIsZoN(t - rhs.t)); }
	__host__ __device__ bool operator>(const cuVAL4 &rhs) const { return ((x > rhs.x) && (y > rhs.y) && (z > rhs.z) && (t > rhs.t)); }
	__host__ __device__ bool operator<(const cuVAL4 &rhs) const { return ((x < rhs.x) && (y < rhs.y) && (z < rhs.z) && (t < rhs.t)); }

	__host__ __device__ bool IsNull(void) const { return (*this == cuVAL4()); }

	//----------------------------- OTHER

	__host__ __device__ VType dim(void) const { return x * y * z * t; }
};

typedef cuVAL4<int> cuINT4;
typedef cuVAL4<size_t> cuSZ4;

typedef cuVAL4<float> cuFLT4;
typedef cuVAL4<cuBReal> cuReal4;
typedef cuVAL4<double> cuDBL4;

typedef cuVAL4<cuINT4> cuINT44;

typedef cuVAL4<cuFLT4> cuFLT44;
typedef cuVAL4<cuReal4> cuReal44;
typedef cuVAL4<cuDBL4> cuDBL44;