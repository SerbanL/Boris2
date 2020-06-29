#pragma once

#include <string>
#include <vector>
#include <functional>

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