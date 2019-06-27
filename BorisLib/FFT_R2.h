#pragma once

#include "FFT.h"

//------------------- RADIX-2

//Radix-2 FFT, packed real and imaginary parts for inputs and outputs.
template <typename fftReal>
template <typename Type>
void FFTMethods_Cpp<fftReal>::FFT_Rad2DIT_Packed(Type *Fx, Type *FX, int ldn, int N)
{
	FX[0] = Fx[0];
	FX[(N - 1)] = Fx[(N - 1)];

	//Shuffle values, moving input data from input to output : input data is left unaffected.
	for (int p = 1; p <= N - 2; p++) {

		FX[p] = Fx[fftShuffle[N - 1 + p]];
	}

	Type u, v;
	fftComplex exp;

	//Simple multiplications by 1 + 0i are treated separately (starting part of outer loop)
	for (int r = 0; r <= N - 1; r += 2) {

		u = FX[r];
		v = FX[(r + 1)];

		FX[r] = u + v;

		FX[(r + 1)] = u - v;
	}

	int T = 2, Th;

	//Remaining part of outer loop
	for (int p = 2; p <= ldn; p++) {

		T *= 2;
		Th = T / 2;

		for (int j = 0; j < Th; j++) {

			exp = ~cossin[(MAXFFT / T) * j];

			for (int r = 0; r <= N - T; r += T) {

				u = FX[(r + j)];

				v = FX[(r + j + Th)] * exp;

				FX[(r + j)] = u + v;

				FX[(r + j + Th)] = u - v;
			}
		}
	}
}