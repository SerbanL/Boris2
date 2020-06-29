#pragma once

#include "Funcs_Conv.h"

////////////////////////////////////////////////////////////////////////////////////////////////// ReIm (complex number)
//
// Complex number and associated operators

template <typename Type, std::enable_if_t<std::is_fundamental<Type>::value>* = nullptr>
struct __ReIm {

	//----------------------------- DATA

	Type Re;
	Type Im;
	
	//----------------------------- VALUE CONSTRUCTORS
	
	__ReIm(void) { Re = 0; Im = 0; }
	__ReIm(Type _Re, Type _Im) { Re = _Re; Im = _Im; }

	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	__ReIm(const __ReIm &copyThis) { Re = copyThis.Re; Im = copyThis.Im; }

	//assignment operator
	__ReIm& operator=(const __ReIm &rhs) { Re = rhs.Re; Im = rhs.Im; return *this; }

	//----------------------------- STREAM OPERATORS

	//allows conversion to std::string and other functionality (e.g. saving to file, or output to text console - stream types derived from std::ostream)
	friend std::ostream& operator<<(std::ostream &os, const __ReIm &rhs) { os << rhs.Re << ", " << rhs.Im; return os; }

	//this also does conversion to std::string, but allows further functionality through the stringconversion object, e.g. units.
	friend Conversion::tostringconversion& operator<<(Conversion::tostringconversion &lhs, const __ReIm &rhs) { lhs << rhs.Re << std::string(", ") << rhs.Im; return lhs; }

	//allows conversions from std::string to ReIm
	friend __ReIm& operator>>(const std::stringstream &ss, __ReIm &rhs)
	{ 
		//normally std::string representation is given as: "Re, Im"
		std::vector<std::string> components = split(ss.str(), ",");
		//it could also be given as: "Re Im"
		if(components.size() == 1) components = split(ss.str(), " ");

		//now set values from strings as would be done in the available constructors if the values were passed directly
		switch (components.size()) {

		case 2:
			rhs.Re = ToNum(trimspaces(components[0]), "");
			rhs.Im = ToNum(trimspaces(components[1]), "");
			break;
		}

		return rhs;
	}

	//----------------------------- ARITHMETIC OPERATORS

	//+=, -=
	void operator+=(const __ReIm &rhs) { Re += rhs.Re; Im += rhs.Im; }
	void operator-=(const __ReIm &rhs) { Re -= rhs.Re; Im -= rhs.Im; }

	//complex conjugate
	__ReIm operator~(void) const { return __ReIm(Re, -Im); }

	//multiplication by i
	__ReIm operator!(void) const { return __ReIm(-Im, Re); }

	//multiplication of two complex numbers
	__ReIm operator*(const __ReIm &rhs) const { return __ReIm(Re * rhs.Re - Im * rhs.Im, Re * rhs.Im + Im * rhs.Re); }

	//multiplication by a constant - constant on the right (must be fundamental)
	template <typename VType, std::enable_if_t<std::is_fundamental<VType>::value>* = nullptr>
	__ReIm operator*(const VType &rhs) const { return __ReIm(Re * (Type)rhs, Im * (Type)rhs); }

	//multiplication by a constant - constant on the left (must be fundamental)
	template <typename VType, std::enable_if_t<std::is_fundamental<VType>::value>* = nullptr>
	friend __ReIm operator*(const VType &lhs, const __ReIm &rhs) { return __ReIm((Type)lhs * rhs.Re, (Type)lhs * rhs.Im); }

	//division by a constant (must be fundamental)
	template <typename VType, std::enable_if_t<std::is_fundamental<VType>::value>* = nullptr>
	__ReIm operator/(const VType &divisor) const { return __ReIm(Re / (Type)divisor, Im / (Type)divisor); }

	//addition
	__ReIm operator+(const __ReIm &rhs) const { return __ReIm(Re + rhs.Re, Im + rhs.Im); }

	//subtraction
	__ReIm operator-(const __ReIm &rhs) const { return __ReIm(Re - rhs.Re, Im - rhs.Im); }

	//----------------------------- OTHER NUMERIC OPERATORS

	//the norm or magnitude
	Type norm(void) const { return sqrt(Re * Re + Im * Im); }

	//----------------------------- COMPARISON OPERATORS

	//comparison
	bool operator==(const __ReIm &rhs) const { return (IsZ(Re - rhs.Re) && IsZ(Im - rhs.Im)); }
	bool operator!=(const __ReIm &rhs) const { return (IsNZ(Re - rhs.Re) || IsNZ(Im - rhs.Im)); }
};

//this is the most common type to use : double precision.
typedef __ReIm<double> ReIm;

////////////////////////////////////////////////////////////////////////////////////////////////// __ReIm3
//
//

//this is just a container of 3 ReIm numbers (labelled conveniently x, y, z) with all operators extended so a __ReIm3 object can be manipulated just as a ReIm one
template <typename Type, std::enable_if_t<std::is_fundamental<Type>::value>* = nullptr>
struct __ReIm3 {

	//----------------------------- DATA

	__ReIm<Type> x, y, z;

	//----------------------------- VALUE CONSTRUCTORS

	__ReIm3(void) { x = __ReIm<Type>(); y = __ReIm<Type>(); z = __ReIm<Type>(); }
	
	__ReIm3(VAL3<Type> Re, VAL3<Type> Im) { x = __ReIm<Type>(Re.x, Im.x); y = __ReIm<Type>(Re.y, Im.y); z = __ReIm<Type>(Re.z, Im.z); }

	__ReIm3(__ReIm<Type> x, __ReIm<Type> y, __ReIm<Type> z) {this->x = x; this->y = y; this->z = z;}
	
	//----------------------------- CONVERTING CONSTRUCTORS

	//copy constructor
	__ReIm3(const __ReIm3 &copyThis) { x = copyThis.x; y = copyThis.y; z = copyThis.z; }

	//assignment operator
	__ReIm3& operator=(const __ReIm3 &rhs) { x = rhs.x; y = rhs.y; z = rhs.z; return *this; }

	//----------------------------- STREAM OPERATORS

	//allows conversion to std::string and other functionality (e.g. saving to file, or output to text console - stream types derived from std::ostream)
	friend std::ostream& operator<<(std::ostream &os, const __ReIm3 &rhs) 
	{
		os << rhs.x.Re << ", " << rhs.x.Im << "; " << rhs.y.Re << ", " << rhs.y.Im << "; " << rhs.z.Re << ", " << rhs.z.Im;

		return os;
	}

	//this also does conversion to std::string, but allows further functionality through the stringconversion object, e.g. units.
	friend Conversion::tostringconversion& operator<<(Conversion::tostringconversion &lhs, const __ReIm3 &rhs) 
	{
		lhs << rhs.x.Re << std::string(", ") << rhs.x.Im << std::string("; ")
			<< rhs.y.Re << std::string(", ") << rhs.y.Im << std::string("; ")
			<< rhs.z.Re << std::string(", ") << rhs.z.Im;

		return lhs;
	}

	//----------------------------- ARITHMETIC OPERATORS

	//complex conjugate
	__ReIm3 operator~(void) const { return __ReIm3(~x, ~y, ~z); }

	//multiplication by i
	__ReIm3 operator!(void) const { return __ReIm3(!x, !y, !z); }

	//multiplication of two complex numbers
	__ReIm3 operator*(const __ReIm3 &rhs) const { return __ReIm3(x * rhs.x, y * rhs.y, z * rhs.z); }

	//multiplication by a ReIm (this means multiply each component by a ReIm) - ReIm on the right
	__ReIm3 operator*(const ReIm &rhs) const { return __ReIm3(x * rhs, y * rhs, z * rhs); }

	//multiplication by a ReIm (this means multiply each component by a ReIm) - ReIm on the left
	friend __ReIm3 operator*(const ReIm &lhs, const __ReIm3 &rhs) { return __ReIm3(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z); }

	//multiplication by a constant - constant on the right (must be fundamental)
	template <typename VType, std::enable_if_t<std::is_fundamental<VType>::value>* = nullptr>
	__ReIm3 operator*(const VType &rhs) const { return __ReIm3(x * (double)rhs, y * (double)rhs, z * (double)rhs); }

	//multiplication by a constant - constant on the left (must be fundamental)
	template <typename VType, std::enable_if_t<std::is_fundamental<VType>::value>* = nullptr>
	friend __ReIm3 operator*(const VType &lhs, const __ReIm3 &rhs) { return __ReIm3((double)lhs * rhs.x, (double)lhs * rhs.y, (double)lhs * rhs.z); }

	//division by a constant (must be fundamental)
	template <typename VType, std::enable_if_t<std::is_fundamental<VType>::value>* = nullptr>
	__ReIm3 operator/(const VType &divisor) const { return __ReIm3(x / (double)divisor, y / (double)divisor, z / (double)divisor); }

	//addition
	__ReIm3 operator+(const __ReIm3 &rhs) const { return __ReIm3(x + rhs.x, y + rhs.y, z + rhs.z); }

	//subtraction
	__ReIm3 operator-(const __ReIm3 &rhs) const { return __ReIm3(x - rhs.x, y - rhs.y, z - rhs.z); }

	//----------------------------- OTHER NUMERIC OPERATORS

	//the norm or magnitude
	VAL3<Type> norm(void) const { return VAL3<Type>(x.norm(), y.norm(), z.norm()); }

	//----------------------------- COMPARISON OPERATORS

	//comparison
	bool operator==(const __ReIm3 &rhs) const { return (x == rhs.x && y == rhs.y && z == rhs.z); }
	bool operator!=(const __ReIm3 &rhs) const { return (x != rhs.x || y != rhs.y || z != rhs.z); }

	//----------------------------- OTHERS

	//extract real part as a VAL3
	VAL3<Type> Re(void) const { return VAL3<Type>(x.Re, y.Re, z.Re); }

	//extract imaginary part as a VAL3
	VAL3<Type> Im(void) const { return VAL3<Type>(x.Im, y.Im, z.Im); }
};

//this is the most common type to use : double precision.
typedef __ReIm3<double> ReIm3;