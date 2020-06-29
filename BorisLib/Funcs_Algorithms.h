#pragma once

#include <string>
#include <vector>
#include <functional>

#include "Types_VAL.h"
#include "Funcs_Math.h"
#include "Funcs_Algorithms_base.h"

///////////////////////////////////////////////////////////////////////////////
//SORTNG

//INCREASING ORDER BY MAGNITUDE

//Quick sort std::vector, which must be of a type which can be ordered (i.e. defines < operator) - note, this also works on strings. 
//Re-arrange all arr_dep in exactly the same way as arr (if included in the quicksort call): these are the dependent arrays, which can have any type. They must have dimensions equal to or greater than arr.
template <typename Type, typename ... PType>
void quicksortmag(std::vector<Type> &arr, std::vector<PType>& ... arr_dep)
{
	auto part = [&arr, &arr_dep...](int lo, int hi) -> int {

		int pi = (hi + lo) / 2;
		Type pv = arr[pi];

		int si = lo;
		swap(pi, hi, arr, arr_dep...);

		for (int i = lo; i <= hi - 1; i++) {

			if (GetMagnitude(arr[i]) < GetMagnitude(pv)) {

				swap(i, si, arr, arr_dep...);
				si++;
			}
		}

		swap(si, hi, arr, arr_dep...);

		return si;
	};

	std::function<void(int, int)> qs;
	qs = [&qs, &part](int lo, int hi) {

		if (lo < hi) {

			int pi = part(lo, hi);
			qs(lo, pi - 1);
			qs(pi + 1, hi);
		}
	};

	//now sort it
	qs(0, (int)arr.size() - 1);
}

//Quick sort array. Array must be of a type which can be ordered
template <typename Type>
void quicksortmag(std::vector<Type> &arr)
{
	auto part = [&arr](int lo, int hi) -> int {

		int pi = (hi + lo) / 2;
		Type pv = arr[pi];

		int si = lo;
		swap(pi, hi, arr);

		for (int i = lo; i <= hi - 1; i++) {

			if (GetMagnitude(arr[i]) < GetMagnitude(pv)) {

				swap(i, si, arr);
				si++;
			}
		}

		swap(si, hi, arr);

		return si;
	};

	std::function<void(int, int)> qs;
	qs = [&arr, &qs, &part](int lo, int hi) {

		if (lo < hi) {

			int pi = part(lo, hi);
			qs(lo, pi - 1);
			qs(pi + 1, hi);
		}
	};

	//now sort it
	qs(0, (int)arr.size() - 1);
}

//DECREASING ORDER BY MAGNITUDE

//Quick sort std::vector, which must be of a type which can be ordered (i.e. defines < operator) - note, this also works on strings. 
//Re-arrange all arr_dep in exactly the same way as arr (if included in the quicksort call): these are the dependent arrays, which can have any type. They must have dimensions equal to or greater than arr.
template <typename Type, typename ... PType>
void invquicksortmag(std::vector<Type> &arr, std::vector<PType>& ... arr_dep)
{
	auto part = [&arr, &arr_dep...](int lo, int hi) -> int {

		int pi = (hi + lo) / 2;
		Type pv = arr[pi];

		int si = lo;
		swap(pi, hi, arr, arr_dep...);

		for (int i = lo; i <= hi - 1; i++) {

			if (GetMagnitude(arr[i]) > GetMagnitude(pv)) {

				swap(i, si, arr, arr_dep...);
				si++;
			}
		}

		swap(si, hi, arr, arr_dep...);

		return si;
	};

	std::function<void(int, int)> qs;
	qs = [&qs, &part](int lo, int hi) {

		if (lo < hi) {

			int pi = part(lo, hi);
			qs(lo, pi - 1);
			qs(pi + 1, hi);
		}
	};

	//now sort it
	qs(0, (int)arr.size() - 1);
}

//Quick sort array. Array must be of a type which can be ordered
template <typename Type>
void invquicksortmag(std::vector<Type> &arr)
{
	auto part = [&arr](int lo, int hi) -> int {

		int pi = (hi + lo) / 2;
		Type pv = arr[pi];

		int si = lo;
		swap(pi, hi, arr);

		for (int i = lo; i <= hi - 1; i++) {

			if (GetMagnitude(arr[i]) > GetMagnitude(pv)) {

				swap(i, si, arr);
				si++;
			}
		}

		swap(si, hi, arr);

		return si;
	};

	std::function<void(int, int)> qs;
	qs = [&arr, &qs, &part](int lo, int hi) {

		if (lo < hi) {

			int pi = part(lo, hi);
			qs(lo, pi - 1);
			qs(pi + 1, hi);
		}
	};

	//now sort it
	qs(0, (int)arr.size() - 1);
}

///////////////////////////////////////////////////////////////////////////////
//ROOT FINDING

#define NEWTONRAPHSON_TIMEOUT	100

//Newton-Raphson single variable

//Solve the equation F(x) = 0 using the Newton-Raphson method 
//i.e. x_n+1 = x_n - F(x_n) / F'(x_n) as n -> inf, for a starting guess x0, to given accuracy
//e.g. if accuracy 1e-6 then solve for x s.t. |F(x)| <= 1e-6
inline double Root_NewtonRaphson(std::function<double(double)>& F, std::function<double(double)>& Fdiff, double x0, double accuracy, int iters_timeout = NEWTONRAPHSON_TIMEOUT)
{
	//starting values
	double xn = x0;
	double yn = F(x0);

	int iterations = 0;

	while (fabs(yn) > accuracy && iterations < iters_timeout) {

		//next values
		xn -= yn / Fdiff(xn);
		yn = F(xn);

		iterations++;
	}

	//return the root to given precision (unless timout reached)
	return xn;
}

//Newton-Raphson multiple variables : 2

//Simple version with off-diagonal Jacobian elements negligible

//Solve the equation system F1(x1, x2) = 0 and F2(x1, x2) = 0 using the Newton-Raphson method for x1, x2, if the off-diagonal Jacobian elements can be ignored.
//i.e. x1_n+1 = x1_n - F1(x1_n, x2_n) / F1'(x1_n, x2_n) and x2_n+1 = x2_n - F2(x1_n, x2_n) / F2'(x1_n, x2_n) as n -> inf, for a starting guess (x1_0, x2_0) to given accuracy
//e.g. if accuracy 1e-6 then solve for x s.t. |F1(x1, x2)| <= 1e-6 and |F2(x1, x2)| <= 1e-6
inline DBL2 Root_NewtonRaphson(
	std::function<double(DBL2)>& F1, std::function<double(DBL2)>& F2,
	std::function<double(DBL2)>& F1diff, std::function<double(DBL2)>& F2diff,
	DBL2 x0,
	double accuracy, int iters_timeout = NEWTONRAPHSON_TIMEOUT)
{
	//starting values
	DBL2 xn = x0;
	DBL2 yn = DBL2(F1(x0), F2(x0));

	int iterations = 0;

	while ((fabs(yn.i) > accuracy || fabs(yn.j) > accuracy) && iterations < iters_timeout) {

		//next values
		xn.i -= yn.i / F1diff(xn);
		xn.j -= yn.j / F2diff(xn);

		yn = DBL2(F1(xn), F2(xn));

		iterations++;
	}

	//return the root to given precision (unless timout reached)
	return xn;
}

//Full multi-variate Newton-Raphson root-finding for 2 variables (need off-diagonal Jacobian elements)

inline DBL2 Root_NewtonRaphson(
	std::function<double(DBL2)>& F1, std::function<double(DBL2)>& F2,
	std::function<double(DBL2)>& F11diff, std::function<double(DBL2)>& F12diff,
	std::function<double(DBL2)>& F21diff, std::function<double(DBL2)>& F22diff,
	DBL2 x0,
	double accuracy, int iters_timeout = NEWTONRAPHSON_TIMEOUT)
{
	//starting values
	DBL2 xn = x0;
	DBL2 yn = DBL2(F1(x0), F2(x0));

	int iterations = 0;

	while ((fabs(yn.i) > accuracy || fabs(yn.j) > accuracy) && iterations < iters_timeout) {

		double det = F11diff(xn) * F22diff(xn) - F12diff(xn) * F21diff(xn);

		//if (!det) break;

		//next values
		xn -= DBL2(F22diff(xn) * yn.i - F12diff(xn) * yn.j, F11diff(xn) * yn.j - F21diff(xn) * yn.i) / det;

		yn = DBL2(F1(xn), F2(xn));

		iterations++;
	}

	//return the root to given precision (unless timout reached)
	return xn;
}