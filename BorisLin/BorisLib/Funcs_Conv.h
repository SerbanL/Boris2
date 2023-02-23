#pragma once

#include <string>
#include <vector>
#include <sstream>

#include "Types_Conversion.h"
#include "Introspection_base.h"

///////////////////////////////////////////////////////////////////////////////
// TO STRING

//convert numerical value to std::string.
template <typename Type> std::string ToString(const Type& value, const std::string& unit)
{ 
	Conversion::tostringconversion sc;
	sc.set_unit(unit);

	sc << value;

	return sc.str();
}

template <typename Type> 
std::string ToString(const Type& value)
{
	return Conversion::ToString_convertible(value, is_streamable_out<std::stringstream, Type>());
}

///////////////////////////////////////////////////////////////////////////////
// FROM STRING

//convert std::string to numerical value. Example call: float fval = ToNum(str);
//If unit specified, allow use of unit magnitudes in input std::string rather than e notation
//inline Conversion::tonumberconversion ToNum(const std::string& text) { return Conversion::tonumberconversion(text, ""); }
inline Conversion::tonumberconversion ToNum(const std::string& text, const std::string& unit = "") { return Conversion::tonumberconversion(text, unit); }
inline Conversion::tonumberconversion ToNum(std::string&& text, const std::string& unit = "") { return Conversion::tonumberconversion( move(text), unit ); }

///////////////////////////////////////////////////////////////////////////////
//

//OType is a std::string : use ToString
template <typename OType, typename IType>
OType Convert(const IType& value, std::true_type)
{
	return ToString(value);
}

//OType not a std::string : use ToNum
template <typename OType, typename IType>
OType Convert(const IType& value, std::false_type)
{
	return ToNum(value);
}

//Use this when either ToNum or ToString should be used depending on the output type
template <typename OType, typename IType> 
OType Convert(const IType& value)
{
	return Convert<OType, IType>(value, is_string<OType>{});
}
