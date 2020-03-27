#pragma once

#include <string>
#include <vector>
#include <functional>
#include "Types.h"

///////////////////////////////////////////////////////////////////////////////
//SORTNG

//INCREASING ORDER

//Quick sort std::vector, which must be of a type which can be ordered (i.e. defines < operator) - note, this also works on strings. 
//Re-arrange all arr_dep in exactly the same way as arr (if included in the quicksort call): these are the dependent arrays, which can have any type. They must have dimensions equal to or greater than arr.
template <typename Type, typename ... PType>
void quicksort(std::vector<Type> &arr, std::vector<PType>& ... arr_dep)
{
	auto part = [&arr, &arr_dep...](int lo, int hi) -> int {

		int pi = (hi + lo) / 2;
		Type pv = arr[pi];

		int si = lo;
		swap(pi, hi, arr, arr_dep...);

		for (int i = lo; i <= hi - 1; i++) {

			if (arr[i] < pv) {

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
void quicksort(std::vector<Type> &arr)
{
	auto part = [&arr](int lo, int hi) -> int {

		int pi = (hi + lo) / 2;
		Type pv = arr[pi];

		int si = lo;
		swap(pi, hi, arr);

		for (int i = lo; i <= hi - 1; i++) {

			if (arr[i] < pv) {

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

//DECREASING ORDER

//Quick sort std::vector, which must be of a type which can be ordered (i.e. defines < operator) - note, this also works on strings. 
//Re-arrange all arr_dep in exactly the same way as arr (if included in the quicksort call): these are the dependent arrays, which can have any type. They must have dimensions equal to or greater than arr.
template <typename Type, typename ... PType>
void invquicksort(std::vector<Type> &arr, std::vector<PType>& ... arr_dep)
{
	auto part = [&arr, &arr_dep...](int lo, int hi) -> int {

		int pi = (hi + lo) / 2;
		Type pv = arr[pi];

		int si = lo;
		swap(pi, hi, arr, arr_dep...);

		for (int i = lo; i <= hi - 1; i++) {

			if (arr[i] > pv) {

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
void invquicksort(std::vector<Type> &arr)
{
	auto part = [&arr](int lo, int hi) -> int {

		int pi = (hi + lo) / 2;
		Type pv = arr[pi];

		int si = lo;
		swap(pi, hi, arr);

		for (int i = lo; i <= hi - 1; i++) {

			if (arr[i] > pv) {

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

///////////////////////////////////////////////////////////////////////////////
//GREATEST COMMON DIVISOR - binary method for positive integers
//Very simple, no parallelization, could improve if needed
//This was taken from : https://en.wikipedia.org/wiki/Binary_GCD_algorithm

inline int gcd_pve(int u, int v)
{
	int shift;

	//GCD(0,v) == v; GCD(u,0) == u, GCD(0,0) == 0
	if (u == 0) return v;
	if (v == 0) return u;

	//Let shift := lg K, where K is the greatest power of 2
	//dividing both u and v.
	for (shift = 0; ((u | v) & 1) == 0; ++shift) {

		u >>= 1;
		v >>= 1;
	}

	while ((u & 1) == 0) {

		u >>= 1;
	}

	// From here on, u is always odd.
	do {
		// remove all factors of 2 in v -- they are not common
		//   note: v is not zero, so while will terminate
		// Loop X
		while ((v & 1) == 0) {

			v >>= 1;
		}

		// Now u and v are both odd. Swap if necessary so u <= v,
		//then set v = v - u (which is even). For bignums, the
		//swapping is just pointer movement, and the subtraction
		//can be done in-place.
		if (u > v) {

			int t = v; v = u; u = t; // Swap u and v.
		}

		v = v - u; // Here v >= u.
	} while (v != 0);

	// restore common factors of 2
	return u << shift;
}

template <typename ... Type>
int gcd_pve(int number1, int number2, Type ... numbers)
{
	return gcd_pve(gcd_pve(number1, number2), numbers...);
}

inline int gcd_pve(std::vector<int>& numbers)
{
	if (numbers.size() >= 2) {

		int gcd_value = gcd_pve(numbers[0], numbers[1]);

		for (int idx = 2; idx < numbers.size(); idx++) {

			gcd_value = gcd_pve(gcd_value, numbers[idx]);
		}

		return gcd_value;
	}
	else if (numbers.size()) {

		return numbers[0];
	}
	else return 0;
}