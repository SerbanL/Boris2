#pragma once

#include "FFT.h"

//------------------- CONSTRUCTOR

template <typename fftReal>
FFTMethods_Cpp<fftReal>::FFTMethods_Cpp(void) 
{
	OmpThreads = omp_get_num_procs();

	//1. Calculate cos and sin values.

	cossin.resize(MAXFFT);
	sincos.resize(MAXFFT);

	for (int idx = 0; idx < MAXFFT; idx++) {

		cossin[idx] = fftComplex(cos(2 * PI*idx / (fftReal)MAXFFT), sin(2 * PI*idx / (fftReal)MAXFFT));
		sincos[idx] = fftComplex(sin(2 * PI*idx / (fftReal)MAXFFT), cos(2 * PI*idx / (fftReal)MAXFFT));
	}

	//2. Calculate shuffling tables

	//use this as a helper lambda
	auto ShuffleIndex = [](int n, int N) -> int {

		int nt = 0, powkm1 = 1;

		if (!n) return 0;

		for (int k = 1; n; k++, N /= 2) {

			int t = n / N;

			if (!t) {

				powkm1 *= 2;
				continue;
			}

			nt += t * powkm1;
			n -= N;
			powkm1 *= 2;
		}

		return nt;
	};

	//Calculate shuffling indices for all FFT lengths up to MAXFFT length - place these values contiguously in memory from shortest FFT length to highest
	//For an FFT of length N, the starting index for these values is N - 1.
	//For each fft of length N (which is a power of 2) we need to store N index values, thus we need 2 ^ (log2(MAXFFT) + 1) - 1 values, which is MAXFFT*2 - 1, total storage locations
	fftShuffle.resize(MAXFFT * 2 - 1);
	fftShuffleInv.resize(MAXFFT * 2 - 1);

	int m = 0;
	while (pow(2, m) < MAXFFT) m++;  //m = log2(MAXFFT)

	fftShuffle[0] = 0;
	fftShuffleInv[0] = 0;

	for (int idx = 1; idx <= m; idx++) {

		int fftLength = (1 << idx);

		//These are used to build shuffle inverse indices : fftShuffleCopy has a copy of the normal shuffle indices, fftShufInv has the indices in increasing order
		//After they are built, fftShufffleCopy is sorted in increasing order with ffftShufInv as a dependent array. The result is fftShufInv will contain the inverse shuffle indices .
		std::vector<int> fftShufInv, fftShufCopy;
		fftShufInv.resize(fftLength);
		fftShufCopy.resize(fftLength);

		fftShuffle[fftLength - 1] = 0;
		fftShuffle[(fftLength - 1) * 2] = fftLength - 1;

		fftShufInv[0] = 0;
		fftShufInv[fftLength - 1] = fftLength - 1;
		fftShufCopy[0] = 0;
		fftShufCopy[fftLength - 1] = fftLength - 1;

		for (int p = 1; p <= fftLength - 2; p++) {

			fftShuffle[fftLength - 1 + p] = ShuffleIndex(p, fftLength / 2);
			fftShufCopy[p] = fftShuffle[fftLength - 1 + p];
			fftShufInv[p] = p;
		}

		//make inverse shuffle indices
		quicksort(fftShufCopy, fftShufInv);

		//and copy them in place
		for (int p = 0; p <= fftLength - 1; p++) {

			fftShuffleInv[fftLength - 1 + p] = fftShufInv[p];
		}
	}
}