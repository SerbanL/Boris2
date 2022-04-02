#pragma once

#include "Funcs_Conv.h"
#include "Funcs_Strings.h"
#include "Funcs_Aux_base.h"

////////////////////////////////////////////////////////////////////////////////////////////////// VAL2 and special cases INT2, PAIR
//
//

template <typename VType> 
struct VAL2 {

	//----------------------------- DATA

	union {

		VType i, x, major, first;
	};

	union {

		VType j, y, minor, second;
	};

	//----------------------------- VALUE CONSTRUCTORS

	VAL2(void) { i = VType(); j = VType(); }

	VAL2(VType i_) : i(i_), j(i_) {}

	VAL2(VType i_, VType j_) : i(i_), j(j_) {}

	//----------------------------- CONVERTING CONSTRUCTORS

	//type conversion constructor
	template <typename CVType> VAL2(const VAL2<CVType> &convThis) 
	{ 
		x = (VType)convThis.x; 
		y = (VType)convThis.y; 
	}

	//copy constructor
	VAL2(const VAL2 &copyThis) { *this = copyThis; }

	//assignment operator
	VAL2& operator=(const VAL2 &rhs) 
	{ 
		x = rhs.x; 
		y = rhs.y; 
		return *this;
	}

	//----------------------------- STREAM OPERATORS

	//allows conversion to std::string and other functionality (e.g. saving to file, or output to text console - stream types derived from std::ostream)
	friend std::ostream& operator<<(std::ostream &os, const VAL2 &rhs)
	{ 
		os << ToString(rhs.i) << ", " << ToString(rhs.j);
		return os; 
	}

	//this also does conversion to std::string, but allows further functionality through the stringconversion object, e.g. units.
	friend Conversion::tostringconversion& operator<<(Conversion::tostringconversion &lhs, const VAL2 &rhs) { lhs << rhs.x << std::string(", ") << rhs.y; return lhs; }

	//allows conversions from std::string
	friend VAL2& operator>>(const std::stringstream &ss, VAL2 &rhs) 
	{
		//normally std::string representation is given as: "x, y"
		std::vector<std::string> components = split(ss.str(), ",");
		//it could also be given as: "x y"
		if(components.size() == 1) components = split(ss.str(), " ");

		//now set values from strings as would be done in the available constructors if the values were passed directly
		switch (components.size()) {

		case 1:
			rhs.x = ToNum(trimspaces(components[0]), "");
			rhs.y = rhs.x;
			break;

		case 2:
			rhs.x = ToNum(trimspaces(components[0]), "");
			rhs.y = ToNum(trimspaces(components[1]), "");
			break;
		}

		return rhs;
	}

	//----------------------------- ARITHMETIC OPERATORS

	//For binary operators if the second operand uses a floating point and the first is integral type, the return type favours the second; otherwise the return type favours the first operand.

	//SUM
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	VAL2<RType> operator+(const VAL2<VType_> &rhs) const { return VAL2<RType>(x + rhs.x, y + rhs.y); }

	//DIFFERENCE
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	VAL2<RType> operator-(const VAL2<VType_> &rhs) const { return VAL2<RType>(x - rhs.x, y - rhs.y); }

	//SCALAR PRODUCT : (same fundamental types, e.g. DBL2 with a DBL2 -> double)
	template <
		typename VType_,
		std::enable_if_t<std::is_same<VType, VType_>::value>* = nullptr
	>
	double operator*(const VAL2<VType_> &rhs) const { return double(x * rhs.x + y * rhs.y); }

	//SIMPLE PRODUCT (component by component)
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	VAL2<RType> operator&(const VAL2<VType_> &rhs) const { return VAL2<RType>(x * rhs.x, y * rhs.y); }

	//DIVISION REMAINDER
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
		VAL2<RType> operator%(const VAL2<VType_> &rhs) const 
	{ 
		return VAL2<RType>(
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
	VAL2<RType> operator*(const MVType &mult) const { return VAL2<RType>(x * mult, y * mult); }

	//product with a constant (must be fundamental type) on the LHS : floating point
	template <
		typename MVType,
		std::enable_if_t<std::is_floating_point<MVType>::value>* = nullptr
	>
	friend VAL2<MVType> operator*(const MVType &mult, const VAL2<VType> &rhs) { return VAL2<MVType>(rhs.x * mult, rhs.y * mult); }

	//product with a constant (must be fundamental type) on the LHS : integral
	template <
		typename MVType,
		std::enable_if_t<std::is_integral<MVType>::value>* = nullptr
	>
	friend VAL2<VType> operator*(const MVType &mult, const VAL2<VType> &rhs) { return VAL2<VType>(rhs.x * mult, rhs.y * mult); }

	//SIMPLE DIVISION (component by component)
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	VAL2<RType> operator/(const VAL2<VType_> &rhs) const { return VAL2<RType>(x / rhs.x, y / rhs.y); }

	//DIVISION BY A CONSTANT (must be fundamental type)
	template <
		class DVType,
		typename RType = typename std::conditional<std::is_floating_point<DVType>::value && std::is_integral<VType>::value, DVType, VType>::type,
		std::enable_if_t<std::is_fundamental<DVType>::value>* = nullptr
	>
	VAL2<RType> operator/(const DVType &divisor) const { return VAL2<RType>(x / divisor, y / divisor); }

	//OPERATION-ASSIGN : ADD
	template <typename VType_>
	void operator+=(const VAL2<VType_> &rhs) { x += rhs.x; y += rhs.y; }

	//OPERATION-ASSIGN : SUBTRACT
	template <typename VType_>
	void operator-=(const VAL2<VType_> &rhs) { x -= rhs.x; y -= rhs.y; }

	//OPERATION-ASSIGN : PRODUCT WITH CONSTANT (must be fundamental type)
	template <
		typename MVType,
		std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr
	>
	void operator*=(const MVType &mult) { x *= mult; y *= mult; }

	//OPERATION-ASSIGN : DIVISION BY CONSTANT (must be fundamental type)
	template <
		typename DVType,
		std::enable_if_t<std::is_fundamental<DVType>::value>* = nullptr
	>
	void operator/=(const DVType &divisor) { x /= divisor; y /= divisor; }

	//----------------------------- OTHER NUMERIC OPERATORS

	//the norm or magnitude
	VType norm(void) const { return sqrt(x * x + y * y); }

	//get this VAL3 in normalized form - not checking for zero magnitude here this must be ensured externally
	VAL2 normalized(void) const { return *this / norm(); }

	//set new magnitude (norm) - not checking for zero magnitude here this must be ensured externally
	void renormalize(VType new_norm) { *this *= new_norm / norm(); }

	//----------------------------- OTHER

	VType dim(void) const { return x * y; }

	//get maximum dimension
	VType maxdim(void) const { return (x < y ? y : x); }

	//get minimum dimension
	VType mindim(void) const { return (x < y ? x : y); }

	//----------------------------- COMPARISON OPERATORS

	//comparison
	bool operator==(const VAL2 &rhs) const { return (IsZ(x - rhs.x) && IsZ(y - rhs.y)); }
	bool operator!=(const VAL2 &rhs) const { return (IsNZ(x - rhs.x) || IsNZ(y - rhs.y)); }

	//ordering
	bool operator>=(const VAL2 &rhs) const { return (IsZoP(x - rhs.x) && IsZoP(y - rhs.y)); }
	bool operator<=(const VAL2 &rhs) const { return (IsZoN(x - rhs.x) && IsZoN(y - rhs.y)); }
	bool operator>(const VAL2 &rhs) const { return ((x > rhs.x) && (y > rhs.y)); }
	bool operator<(const VAL2 &rhs) const { return ((x < rhs.x) && (y < rhs.y)); }

	bool IsNull(void) const { return (*this == VAL2()); }

	//----------------------------- PERMUTATIONS

	//return new VAL2 with swapped x and y values
	VAL2<VType> swap_xy(void) const
	{
		VAL2<VType> transp(y, x);
		return transp;
	}
};


typedef VAL2<int> INT2;
typedef VAL2<float> FLT2;
typedef VAL2<double> DBL2;
typedef std::pair<INT2, INT2> PAIR;

////////////////////////////////////////////////////////////////////////////////////////////////// VAL3 and special cases INT3, FLT3, DBL3
//
//

template <typename VType> 
struct VAL3 {
	
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

	//----------------------------- VALUE CONSTRUCTORS

	VAL3(void) { x = VType(); y = VType(); z = VType(); }
	VAL3(VType val) { x = val; y = val; z = val; }
	VAL3(VType x, VType y, VType z) { this->x = x; this->y = y; this->z = z; }

	//----------------------------- CONVERTING CONSTRUCTORS

	//type conversion constructor
	template <typename CVType> VAL3(const VAL3<CVType> &convThis) { x = (VType)convThis.x; y = (VType)convThis.y; z = (VType)convThis.z; }

	//copy constructor
	VAL3(const VAL3 &copyThis) { x = copyThis.x; y = copyThis.y; z = copyThis.z; }

	//assignment operator
	VAL3& operator=(const VAL3 &rhs) { x = rhs.x; y = rhs.y; z = rhs.z; return *this; }

	//----------------------------- STREAM OPERATORS

	//allows conversion to std::string and other functionality (e.g. saving to file, or output to text console - stream types derived from std::ostream)
	friend std::ostream& operator<<(std::ostream &os, const VAL3 &rhs) { os << ToString(rhs.x) << ", " << ToString(rhs.y) << ", " << ToString(rhs.z); return os; }

	//this also does conversion to std::string, but allows further functionality through the stringconversion object, e.g. units.
	friend Conversion::tostringconversion& operator<<(Conversion::tostringconversion &lhs, const VAL3 &rhs) { lhs << rhs.x << std::string(", ") << rhs.y << std::string(", ") << rhs.z; return lhs; }

	//allows conversions from std::string
	friend VAL3& operator>>(const std::stringstream &ss, VAL3 &rhs)
	{ 
		//normally std::string representation is given as: "x, y, z"
		std::vector<std::string> components = split(ss.str(), ",");
		
		//it could also be given as: "x y z"
		if (components.size() == 1) components = split(ss.str(), " ");

		//another possibility is to have it in polar form: 
		//the normal specification as x, y, z is in Cartesian form, but instead we can have r; t, p (r is the magnitude, t is theta - polar angle - p is phi - azimuthal angle)
		//in this case components.size() is two, and we should be able to split components[0] using ";"

		//now set values from strings as would be done in the available constructors if the values were passed directly
		switch (components.size()) {

		//all components are the same
		case 1:
			rhs.x = ToNum(trimspaces(components[0]), "");
			rhs.y = rhs.x;
			rhs.z = rhs.x;
			break;

		//std::string in Polar form - convert to Cartesian
		case 2:
		{
			std::vector<std::string> r_theta = split(components[0], ";");
			if (r_theta.size() != 2) return rhs;

			double r = ToNum(trimspaces(r_theta[0]), "");
			double pol = ToNum(trimspaces(r_theta[1]), "");
			double azim = ToNum(trimspaces(components[1]), "");

			//Convert to Cartesian, angles in degrees
			rhs.x = r * sin(pol*PI / 180) * cos(azim*PI / 180);
			rhs.y = r * sin(pol*PI / 180) * sin(azim*PI / 180);
			rhs.z = r * cos(pol*PI / 180);
		}
			break;

		//std::string in Cartesian form
		case 3:
			rhs.x = ToNum(trimspaces(components[0]), "");
			rhs.y = ToNum(trimspaces(components[1]), "");
			rhs.z = ToNum(trimspaces(components[2]), "");
			break;
		}

		return rhs;
	}

	//----------------------------- ARITHMETIC OPERATORS

	//For binary operators if the second operand uses a floating point and the first is integral type, the return type favours the second; otherwise the return type favours the first operand.
	
	//SUM
	template <
		typename VType_, 
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	VAL3<RType> operator+(const VAL3<VType_> &rhs) const { return VAL3<RType>(x + rhs.x, y + rhs.y, z + rhs.z); }
	
	//DIFFERENCE
	template <
		typename VType_, 
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	VAL3<RType> operator-(const VAL3<VType_> &rhs) const { return VAL3<RType>(x - rhs.x, y - rhs.y, z - rhs.z); }

	//SCALAR PRODUCT : (same fundamental types, e.g. DBL3 with a DBL3 -> double)
	template <
		typename VType_,
		std::enable_if_t<std::is_same<VType, VType_>::value && std::is_fundamental<VType>::value>* = nullptr
	>
	double operator*(const VAL3<VType_> &rhs) const { return double(x * rhs.x + y * rhs.y + z * rhs.z); }

	//MATRIX PRODUCT : (same types but not fundamental, e.g. DBL33 with a DBL33 -> DBL33)
	template <
		typename VType_,
		std::enable_if_t<std::is_same<VType, VType_>::value && !std::is_fundamental<VType>::value>* = nullptr
	>
		VAL3<VType> operator*(const VAL3<VType_> &rhs) const
	{
		return VAL3<VType>(
			VType(i.i * rhs.i.i + i.j * rhs.j.i + i.k * rhs.k.i, i.i * rhs.i.j + i.j * rhs.j.j + i.k * rhs.k.j, i.i * rhs.i.k + i.j * rhs.j.k + i.k * rhs.k.k),
			VType(j.i * rhs.i.i + j.j * rhs.j.i + j.k * rhs.k.i, j.i * rhs.i.j + j.j * rhs.j.j + j.k * rhs.k.j, j.i * rhs.i.k + j.j * rhs.j.k + j.k * rhs.k.k),
			VType(k.i * rhs.i.i + k.j * rhs.j.i + k.k * rhs.k.i, k.i * rhs.i.j + k.j * rhs.j.j + k.k * rhs.k.j, k.i * rhs.i.k + k.j * rhs.j.k + k.k * rhs.k.k));
	}

	//MATRIX PRODUCT WITH VECTOR : (DBL33 * DBL3 -> DBL3)
	template <
		typename VType_,
		std::enable_if_t<!std::is_fundamental<VType>::value && std::is_fundamental<VType_>::value>* = nullptr
	>
	VAL3<VType_> operator*(const VAL3<VType_> &rhs) const { return VAL3<VType_>(x * rhs, y * rhs, z * rhs); }

	//VECTOR PRODUCT (DBL3 ^ DBL3 -> DBL3)
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	VAL3<RType> operator^(const VAL3<VType_> &rhs) const { return VAL3<RType>(y*rhs.z - z * rhs.y, z*rhs.x - x * rhs.z, x*rhs.y - y * rhs.x); }

	//OUTER VECTOR PRODUCT : (both types are fundamental, e.g. DBL3 | DBL3 -> DBL33
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type,
		std::enable_if_t<std::is_fundamental<VType_>::value && std::is_fundamental<VType>::value>* = nullptr
	>
	VAL3<VAL3<RType>> operator|(const VAL3<VType_> &rhs) const { return VAL3<VAL3<RType>>(VAL3<RType>(x * rhs.x, x * rhs.y, x * rhs.z), VAL3<RType>(y * rhs.x, y * rhs.y, y * rhs.z), VAL3<RType>(z * rhs.x, z * rhs.y, z * rhs.z)); }

	//MATRIX COLUMN PRODUCT WITH VECTOR : (not same types, but one is a fundamental type, e.g. DBL33 | DBL3 -> DBL3, where DBL33 | DBL3 = transpose(DBL33) * DBL3 )
	template <
		typename VType_,
		std::enable_if_t<!std::is_fundamental<VType>::value && std::is_fundamental<VType_>::value>* = nullptr
	>
	VAL3<VType_> operator|(const VAL3<VType_> &rhs) const { return VAL3<VType_>(x * rhs.x + y * rhs.y + z * rhs.z); }

	//SIMPLE PRODUCT (component by component)
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	VAL3<RType> operator&(const VAL3<VType_> &rhs) const { return VAL3<RType>(x * rhs.x, y * rhs.y, z * rhs.z); }

	//DIVISION REMAINDER
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
		VAL3<RType> operator%(const VAL3<VType_> &rhs) const 
	{ 
		return VAL3<RType>(
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
	VAL3<VType> operator*(const MVType &mult) const { return VAL3<VType>(x * mult, y * mult, z * mult); }

	//product with a constant (must be fundamental type) on the LHS
	template <
		typename MVType,
		std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr
	>
	friend VAL3<VType> operator*(const MVType &mult, const VAL3<VType> &rhs) { return VAL3<VType>(rhs.x * mult, rhs.y * mult, rhs.z * mult); }

	//SIMPLE DIVISION (component by component)
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
	VAL3<RType> operator/(const VAL3<VType_> &rhs) const { return VAL3<RType>(x / rhs.x, y / rhs.y, z / rhs.z); }

	//DIVISION BY A CONSTANT (must be fundamental type)
	template <
		class DVType,
		std::enable_if_t<std::is_fundamental<DVType>::value>* = nullptr
	>
	VAL3<VType> operator/(const DVType &divisor) const { return VAL3<VType>(x / divisor, y / divisor, z / divisor); }

	//OPERATION-ASSIGN : ADD
	template <typename VType_>
	void operator+=(const VAL3<VType_> &rhs) { x += rhs.x; y += rhs.y; z += rhs.z; }
	
	//OPERATION-ASSIGN : SUBTRACT
	template <typename VType_>
	void operator-=(const VAL3<VType_> &rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; }

	//OPERATION-ASSIGN : PRODUCT WITH CONSTANT (must be fundamental type)
	template <
		typename MVType, 
		std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr
	>
	void operator*=(const MVType &mult) { x *= mult; y *= mult; z *= mult; }

	//OPERATION-ASSIGN : DIVISION BY CONSTANT (must be fundamental type)
	template <
		typename DVType, 
		std::enable_if_t<std::is_fundamental<DVType>::value>* = nullptr
	>
	void operator/=(const DVType &divisor) { x /= divisor; y /= divisor; z /= divisor; }

	//----------------------------- OTHER NUMERIC OPERATORS

	//the norm or magnitude
	VType norm(void) const { return sqrt(x * x + y * y + z * z); }

	//get this VAL3 in normalized form - not checking for zero magnitude here this must be ensured externally
	VAL3 normalized(void) const { return *this / norm(); }

	//set new magnitude (norm) - not checking for zero magnitude here this must be ensured externally
	void renormalize(VType new_norm) { *this *= new_norm / norm(); }

	//----------------------------- COMPARISON OPERATORS

	//comparison
	bool operator==(const VAL3 &rhs) const { return (IsZ(x - rhs.x) && IsZ(y - rhs.y) && IsZ(z - rhs.z)); }
	bool operator!=(const VAL3 &rhs) const { return (IsNZ(x - rhs.x) || IsNZ(y - rhs.y) || IsNZ(z - rhs.z)); }

	//ordering operators
	bool operator>=(const VAL3 &rhs) const { return (IsZoP(x - rhs.x) && IsZoP(y - rhs.y) && IsZoP(z - rhs.z)); }
	bool operator<=(const VAL3 &rhs) const { return (IsZoN(x - rhs.x) && IsZoN(y - rhs.y) && IsZoN(z - rhs.z)); }
	bool operator>(const VAL3 &rhs) const { return ((x > rhs.x) && (y > rhs.y) && (z > rhs.z)); }
	bool operator<(const VAL3 &rhs) const { return ((x < rhs.x) && (y < rhs.y) && (z < rhs.z)); }

	bool IsNull(void) const { return (*this == VAL3()); }

	//----------------------------- OTHER

	VType dim(void) const { return x * y * z; }

	//get maximum dimension
	VType maxdim(void) const { return maximum(x, y, z); }

	//get minimum dimension
	VType mindim(void) const { return minimum(x, y, z); }

	//----------------------------- PERMUTATIONS

	//return new VAL3 with swapped x and y values
	VAL3<VType> swap_xy(void) const
	{
		VAL3<VType> transp(y, x, z);
		return transp;
	}

	//return new VAL3 with swapped x and z values
	VAL3<VType> swap_xz(void) const
	{
		VAL3<VType> transp(z, y, x);
		return transp;
	}

	//return new VAL3 with swapped y and z values
	VAL3<VType> swap_yz(void) const
	{
		VAL3<VType> transp(x, z, y);
		return transp;
	}

	//return new VAL3 with down-cycled xyz values (xyz -> yzx)
	VAL3<VType> cycledn(void) const
	{
		VAL3<VType> transp(y, z, x);
		return transp;
	}

	//return new VAL3 with up-cycled xyz values (xyz -> zxy)
	VAL3<VType> cycleup(void) const
	{
		VAL3<VType> transp(z, x, y);
		return transp;
	}
};

typedef VAL3<int> INT3;
typedef VAL3<unsigned> uINT3;
typedef VAL3<size_t> SZ3;

typedef VAL3<float> FLT3;
typedef VAL3<double> DBL3;

typedef VAL3<INT3> INT33;

typedef VAL3<FLT3> FLT33;
typedef VAL3<DBL3> DBL33;

////////////////////////////////////////////////////////////////////////////////////////////////// VAL4 and special cases INT4, FLT4, DBL4
//
//

template <typename VType>
struct VAL4 {

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

	//----------------------------- VALUE CONSTRUCTORS

	VAL4(void) { x = VType(); y = VType(); z = VType(); t = VType(); }
	VAL4(VType val) { x = val; y = val; z = val; t = val; }
	VAL4(VType x, VType y, VType z, VType t) { this->x = x; this->y = y; this->z = z; this->t = t; }

	//----------------------------- CONVERTING CONSTRUCTORS

	//type conversion constructor
	template <typename CVType> VAL4(const VAL4<CVType> &convThis) { x = (VType)convThis.x; y = (VType)convThis.y; z = (VType)convThis.z; t = (VType)convThis.t; }

	//copy constructor
	VAL4(const VAL4 &copyThis) { x = copyThis.x; y = copyThis.y; z = copyThis.z; t = copyThis.t; }

	//assignment operator
	VAL4& operator=(const VAL4 &rhs) { x = rhs.x; y = rhs.y; z = rhs.z; t = rhs.t; return *this; }

	//----------------------------- STREAM OPERATORS

	//allows conversion to std::string and other functionality (e.g. saving to file, or output to text console - stream types derived from std::ostream)
	friend std::ostream& operator<<(std::ostream &os, const VAL4 &rhs) { os << ToString(rhs.x) << ", " << ToString(rhs.y) << ", " << ToString(rhs.z) << ", " << ToString(rhs.t); return os; }

	//this also does conversion to std::string, but allows further functionality through the stringconversion object, e.g. units.
	friend Conversion::tostringconversion& operator<<(Conversion::tostringconversion &lhs, const VAL4 &rhs) { lhs << rhs.x << std::string(", ") << rhs.y << std::string(", ") << rhs.z << std::string(", ") << rhs.t; return lhs; }

	//allows conversions from std::string
	friend VAL4& operator>>(const std::stringstream &ss, VAL4 &rhs)
	{
		//normally std::string representation is given as: "x, y, z, t"
		std::vector<std::string> components = split(ss.str(), ",");
		//it could also be given as: "x y z t"
		if (components.size() == 1) components = split(ss.str(), " ");

		//now set values from strings as would be done in the available constructors if the values were passed directly
		switch (components.size()) {

		case 1:
			rhs.x = ToNum(trimspaces(components[0]), "");
			rhs.y = rhs.x;
			rhs.z = rhs.x;
			rhs.t = rhs.x;
			break;

		case 4:
			rhs.x = ToNum(trimspaces(components[0]), "");
			rhs.y = ToNum(trimspaces(components[1]), "");
			rhs.z = ToNum(trimspaces(components[2]), "");
			rhs.t = ToNum(trimspaces(components[3]), "");
			break;
		}

		return rhs;
	}

	//----------------------------- ARITHMETIC OPERATORS

	//For binary operators if the second operand uses a floating point and the first is integral type, the return type favours the second; otherwise the return type favours the first operand.

	//SUM
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
		VAL4<RType> operator+(const VAL4<VType_> &rhs) const { return VAL4<RType>(x + rhs.x, y + rhs.y, z + rhs.z, t + rhs.t); }

	//DIFFERENCE
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
		VAL4<RType> operator-(const VAL4<VType_> &rhs) const { return VAL4<RType>(x - rhs.x, y - rhs.y, z - rhs.z, t - rhs.t); }

	//SCALAR PRODUCT : (same fundamental types, e.g. DBL4 with a DBL4 -> double)
	template <
		typename VType_,
		std::enable_if_t<std::is_same<VType, VType_>::value && std::is_fundamental<VType>::value>* = nullptr
	>
		double operator*(const VAL4<VType_> &rhs) const { return double(x * rhs.x + y * rhs.y + z * rhs.z + t * rhs.t); }

	//DIVISION REMAINDER
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
		VAL4<RType> operator%(const VAL4<VType_> &rhs) const 
	{ 
		return VAL4<RType>(
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
		VAL4<VType> operator*(const MVType &mult) const { return VAL4<VType>(x * mult, y * mult, z * mult, t * mult); }

	//product with a constant (must be fundamental type) on the LHS
	template <
		typename MVType,
		std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr
	>
		friend VAL4<VType> operator*(const MVType &mult, const VAL4<VType> &rhs) { return VAL4<VType>(rhs.x * mult, rhs.y * mult, rhs.z * mult, rhs.t * mult); }

	//SIMPLE DIVISION (component by component)
	template <
		typename VType_,
		typename RType = typename std::conditional<std::is_floating_point<VType_>::value && std::is_integral<VType>::value, VType_, VType>::type
	>
		VAL4<RType> operator/(const VAL4<VType_> &rhs) const { return VAL4<RType>(x / rhs.x, y / rhs.y, z / rhs.z, t / rhs.t); }

	//DIVISION BY A CONSTANT (must be fundamental type)
	template <
		class DVType,
		std::enable_if_t<std::is_fundamental<DVType>::value>* = nullptr
	>
		VAL4<VType> operator/(const DVType &divisor) const { return VAL4<VType>(x / divisor, y / divisor, z / divisor, t / divisor); }

	//OPERATION-ASSIGN : ADD
	template <typename VType_>
	void operator+=(const VAL4<VType_> &rhs) { x += rhs.x; y += rhs.y; z += rhs.z; t += rhs.t; }

	//OPERATION-ASSIGN : SUBTRACT
	template <typename VType_>
	void operator-=(const VAL4<VType_> &rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; t -= rhs.t; }

	//OPERATION-ASSIGN : PRODUCT WITH CONSTANT (must be fundamental type)
	template <
		typename MVType,
		std::enable_if_t<std::is_fundamental<MVType>::value>* = nullptr
	>
		void operator*=(const MVType &mult) { x *= mult; y *= mult; z *= mult; t *= mult; }

	//OPERATION-ASSIGN : DIVISION BY CONSTANT (must be fundamental type)
	template <
		typename DVType,
		std::enable_if_t<std::is_fundamental<DVType>::value>* = nullptr
	>
		void operator/=(const DVType &divisor) { x /= divisor; y /= divisor; z /= divisor; t /= divisor; }

	//----------------------------- OTHER NUMERIC OPERATORS

	//the norm or magnitude
	VType norm(void) const { return sqrt(x * x + y * y + z * z + t * t); }

	//get this VAL4 in normalized form - not checking for zero magnitude here this must be ensured externally
	VAL4 normalized(void) const { return *this / norm(); }

	//set new magnitude (norm) - not checking for zero magnitude here this must be ensured externally
	void renormalize(VType new_norm) { *this *= new_norm / norm(); }

	//----------------------------- COMPARISON OPERATORS

	//comparison
	bool operator==(const VAL4 &rhs) const { return (IsZ(x - rhs.x) && IsZ(y - rhs.y) && IsZ(z - rhs.z) && IsZ(t - rhs.t)); }
	bool operator!=(const VAL4 &rhs) const { return (IsNZ(x - rhs.x) || IsNZ(y - rhs.y) || IsNZ(z - rhs.z) || IsNZ(t - rhs.t)); }

	//ordering operators
	bool operator>=(const VAL4 &rhs) const { return (IsZoP(x - rhs.x) && IsZoP(y - rhs.y) && IsZoP(z - rhs.z) && IsZoP(t - rhs.t)); }
	bool operator<=(const VAL4 &rhs) const { return (IsZoN(x - rhs.x) && IsZoN(y - rhs.y) && IsZoN(z - rhs.z) && IsZoN(t - rhs.t)); }
	bool operator>(const VAL4 &rhs) const { return ((x > rhs.x) && (y > rhs.y) && (z > rhs.z) && (t > rhs.t)); }
	bool operator<(const VAL4 &rhs) const { return ((x < rhs.x) && (y < rhs.y) && (z < rhs.z) && (t < rhs.t)); }

	bool IsNull(void) const { return (*this == VAL4()); }

	//----------------------------- OTHER

	VType dim(void) const { return x * y * z * t; }

	//get maximum dimension
	VType maxdim(void) const { return maximum(x, y, z, t); }

	//get minimum dimension
	VType mindim(void) const { return minimum(x, y, z, t); }
};

typedef VAL4<int> INT4;
typedef VAL4<unsigned> uINT4;
typedef VAL4<size_t> SZ4;

typedef VAL4<float> FLT4;
typedef VAL4<double> DBL4;

typedef VAL4<INT4> INT44;

typedef VAL4<FLT4> FLT44;
typedef VAL4<DBL4> DBL44;