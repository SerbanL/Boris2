#pragma once

#include <d2d1_1.h>
#include <string>
#include <vector>
#include <sstream>

#include "Types_Conversion.h"
#include "Introspection_base.h"

///////////////////////////////////////////////////////////////////////////////
// GENERAL CONVERSIONS

inline std::wstring StringtoWideString(const std::string& text) 
{
	//string to wstring conversion. To get LPCWSTR then attach .c_str() to returned wstring
	int len;
    int slength = (int)text.length() + 1;
    len = MultiByteToWideChar(CP_ACP, 0, text.c_str(), slength, 0, 0); 
    wchar_t* buf = new wchar_t[len];
    MultiByteToWideChar(CP_ACP, 0, text.c_str(), slength, buf, len);
	std::wstring wstring_text(buf);
    delete[] buf;

	return wstring_text;
}

inline std::string WideStringtoString(std::wstring wstr) 
{
	return std::string(wstr.begin(), wstr.end());
}

inline WCHAR* StringtoWCHARPointer(const std::string& text) 
{
	//string to wstring conversion. To get LPCWSTR then attach .c_str() to returned wstring
	int len;
    int slength = (int)text.length() + 1;
    len = MultiByteToWideChar(CP_ACP, 0, text.c_str(), slength, 0, 0);
    wchar_t* buf = new wchar_t[len];
	MultiByteToWideChar(CP_ACP, 0, text.c_str(), slength, buf, len);

	return buf;
}

inline void RECTtoD2D1RECT(RECT &rc, D2D1_RECT_F &d2d1_rc) { d2d1_rc.bottom = (FLOAT)rc.bottom; d2d1_rc.top = (FLOAT)rc.top; d2d1_rc.left = (FLOAT)rc.left; d2d1_rc.right = (FLOAT)rc.right; }

inline void D2D1RECTtoRECT(D2D1_RECT_F &d2d1_rc, RECT &rc) { rc.bottom = (LONG)d2d1_rc.bottom; rc.top = (LONG)d2d1_rc.top; rc.left = (LONG)d2d1_rc.left; rc.right = (LONG)d2d1_rc.right; }

///////////////////////////////////////////////////////////////////////////////
// TO STRING

//convert numerical value to string.
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

//convert string to numerical value. Example call: float fval = ToNum(str);
//If unit specified, allow use of unit magnitudes in input string rather than e notation
//inline Conversion::tonumberconversion ToNum(const string& text) { return Conversion::tonumberconversion(text, ""); }
inline Conversion::tonumberconversion ToNum(const std::string& text, const std::string& unit = "") { return Conversion::tonumberconversion(text, unit); }
inline Conversion::tonumberconversion ToNum(std::string&& text, const std::string& unit = "") { return Conversion::tonumberconversion( move(text), unit ); }

///////////////////////////////////////////////////////////////////////////////
//

//OType is a string : use ToString
template <typename OType, typename IType>
OType Convert(const IType& value, std::true_type)
{
	return ToString(value);
}

//OType not a string : use ToNum
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
