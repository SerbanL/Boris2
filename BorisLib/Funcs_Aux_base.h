#pragma once

#define _CRT_SECURE_NO_WARNINGS

#include <d2d1_1.h>
#include <string>
#include <vector>
#include <sstream>
#include <time.h>
#include <windows.h>

//If a floating point value is between epsilon negative and positive then it is considered to be zero (used for comparisons when we need to know if a value is zero or not)
#define EPSILON_POSITIVE 1e-12F						
#define EPSILON_NEGATIVE -1e-12F

//Used for floor and ceil functions to define a non-fixed epsilon : the epsilon depends on the value the function is used on.
#define FLOOR_CEIL_RATIO	1e-9

///////////////////////////////////////////////////////////////////////////////
//MEMORY

inline size_t MemGetFree(void)
{
	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof(statex);
	GlobalMemoryStatusEx(&statex);

	return statex.ullAvailPhys;
}

inline size_t MemGetTotal(void)
{
	MEMORYSTATUSEX statex;
	statex.dwLength = sizeof(statex);
	GlobalMemoryStatusEx(&statex);

	return statex.ullTotalPhys;
}

///////////////////////////////////////////////////////////////////////////////
//TIME

//Note, time and ctime cause deprecation warnings. Fine to use them, just remember they are not re-entrant.
inline std::string Get_Date_Time(void) { time_t rawtime; time(&rawtime); return std::string(ctime(&rawtime)); }

///////////////////////////////////////////////////////////////////////////////
//COMPARISONS - Should not compare floating point values for equality directly (==, >= or <=). Just use one of the functions below instead.

//Is Zero

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
bool IsZ(bfloat fval)  { return (fval < (bfloat)EPSILON_POSITIVE && fval > (bfloat)EPSILON_NEGATIVE); }

template <typename bint, std::enable_if_t<std::is_integral<bint>::value>* = nullptr>
bool IsZ(bint ival) { return (ival == 0); }

//Is Not Zero

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
bool IsNZ(bfloat fval) { return (fval > (bfloat)EPSILON_POSITIVE || fval < (bfloat)EPSILON_NEGATIVE); }

template <typename bint, std::enable_if_t<std::is_integral<bint>::value>* = nullptr>
bool IsNZ(bint ival) { return (ival != 0); }

//Is Zero or Positive

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
bool IsZoP(bfloat fval)  { return (fval > 0 || (fval < (bfloat)EPSILON_POSITIVE && fval > (bfloat)EPSILON_NEGATIVE)); }

template <typename bint, std::enable_if_t<std::is_integral<bint>::value>* = nullptr>
bool IsZoP(bint ival) { return (ival >= 0); }

//Is Zero or Negative

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
bool IsZoN(bfloat fval)  { return (fval < 0 || (fval < (bfloat)EPSILON_POSITIVE && fval > (bfloat)EPSILON_NEGATIVE)); }

template <typename bint, std::enable_if_t<std::is_integral<bint>::value>* = nullptr>
bool IsZoN(bint ival) { return (ival <= 0); }

//Is Equal

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
bool IsE(bfloat fval1, bfloat fval2) { return (fval2 - fval1 < (bfloat)EPSILON_POSITIVE && fval2 - fval1 > (bfloat)EPSILON_NEGATIVE); }

template <typename bint, std::enable_if_t<std::is_integral<bint>::value>* = nullptr>
bool IsE(bint ival1, bint ival2) { return (ival1 == ival2); }

//Is Greater or Equal

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
bool IsGE(bfloat fval1, bfloat fval2) { return (fval2 - fval1 < (bfloat)EPSILON_POSITIVE); }

template <typename bint, std::enable_if_t<std::is_integral<bint>::value>* = nullptr>
bool IsGE(bint ival1, bint ival2) { return (ival1 >= ival2); }

//Is Smaller or Equal

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
bool IsSE(bfloat fval1, bfloat fval2) { return (fval1 - fval2 < (bfloat)EPSILON_POSITIVE); }

template <typename bint, std::enable_if_t<std::is_integral<bint>::value>* = nullptr>
bool IsSE(bint ival1, bint ival2) { return (ival1 <= ival2); }

//Is Smaller - sometimes just checking fval1 < fval2 is not good enough : fval1 must also be at least an EPSILON_POSITIVE lower than fval2.

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
bool IsS(bfloat fval1, bfloat fval2) { return (fval1 < fval2 - (bfloat)EPSILON_POSITIVE); }

//Is Greater - sometimes just checking fval1 > fval2 is not good enough : fval1 must also be at least an EPSILON_POSITIVE greater than fval2.

template <typename bfloat, std::enable_if_t<std::is_floating_point<bfloat>::value>* = nullptr>
bool IsG(bfloat fval1, bfloat fval2) { return (fval1 > fval2 + (bfloat)EPSILON_POSITIVE); }

//not using variadic function here, since the number of indexes to be checked is usually small
inline bool GoodIdx(int upperLimit, int idx) { return (idx >= 0 && idx <= upperLimit); }
inline bool GoodIdx(int upperLimit, int idx1, int idx2) { return (idx1 >= 0 && idx1 <= upperLimit && idx2 >= 0 && idx2 <= upperLimit); }
inline bool GoodIdx(int upperLimit, int idx1, int idx2, int idx3) { return (idx1 >= 0 && idx1 <= upperLimit && idx2 >= 0 && idx2 <= upperLimit && idx3 >= 0 && idx3 <= upperLimit); }
inline bool GoodIdx(int upperLimit, int idx1, int idx2, int idx3, int idx4) { return (idx1 >= 0 && idx1 <= upperLimit && idx2 >= 0 && idx2 <= upperLimit && idx3 >= 0 && idx3 <= upperLimit && idx4 >= 0 && idx4 <= upperLimit); }

inline bool SameColor(D2D1_COLOR_F color1, D2D1_COLOR_F color2) { return (IsE(color1.r, color2.r) && IsE(color1.g, color2.g) && IsE(color1.b, color2.b) && IsE(color1.a, color2.a)); }