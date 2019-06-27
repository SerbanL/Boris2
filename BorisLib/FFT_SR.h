#pragma once

#include "FFT.h"

template <typename fftReal>
template <typename Type>
void FFTMethods_Cpp<fftReal>::FFT_SRDIT_Packed(Type *Fx, Type *FX, int ldn, int N)
{
	int N2, N4;

	int ix, id, i0, i1, i2, i3;

	Type t0, t1;
	fftComplex exp1, exp3;

	//Shuffle values, moving input data from input to output : input data is left unaffected.
	FX[0] = Fx[0];
	FX[(N - 1)] = Fx[(N - 1)];

	for (int p = 1; p <= N - 2; p++) {

		FX[p] = Fx[fftShuffle[N - 1 + p]];
	}

	//Stage 1
	for (ix = 0, id = 4; ix < N; id *= 4) {

		for (i0 = ix; i0 < N; i0 += id) {

			t0 = FX[i0] - FX[(i0 + 1)];
			FX[i0] = FX[i0] + FX[(i0 + 1)];
			FX[(i0 + 1)] = t0;
		}

		ix = 2 * id - 2;
	}

	//Stage 2
	N2 = 2;

	for (int k = ldn - 2; k >= 0; k--) {

		N2 = N2 * 2;
		N4 = N2 / 4;

		//Treat case j = 0 separately due to simpler multiplications
		for (ix = 0, id = N2 * 2; ix < N; id *= 4) {

			for (i0 = ix; i0 < N; i0 += id) {

				i1 = i0 + N4;
				i2 = i1 + N4;
				i3 = i2 + N4;

				t0 = FX[i2] + FX[i3];
				t1 = !(FX[i3] - FX[i2]);  // ! is multiplication by i

				FX[i2] = FX[i0] - t0;
				FX[i0] = FX[i0] + t0;

				FX[i3] = FX[i1] - t1;
				FX[i1] = FX[i1] + t1;
			}

			ix = 2 * id - N2;
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			exp1 = ~cossin[(MAXFFT / N2) * j];
			exp3 = ~cossin[(MAXFFT / N2) * 3 * j];

			for (ix = j, id = N2 * 2; ix < N; id *= 4) {

				for (i0 = ix; i0 < N; i0 += id) {

					i1 = i0 + N4;
					i2 = i1 + N4;
					i3 = i2 + N4;

					t0 = FX[i3] * exp3 + FX[i2] * exp1;
					t1 = !(FX[i3] * exp3 - FX[i2] * exp1);  // ! is multiplication by i

					FX[i2] = FX[i0] - t0;
					FX[i0] = FX[i0] + t0;

					FX[i3] = FX[i1] - t1;
					FX[i1] = FX[i1] + t1;
				}

				ix = 2 * id - N2 + j;
			}
		}
	}
}

template <typename fftReal>
template <typename Type>
void FFTMethods_Cpp<fftReal>::IFFT_SRDIF_Packed(Type *Fx, Type *FX, int ldn, int N)
{
	int N2 = 2 * N, N4;

	int ix, id, i0, i1, i2, i3;

	Type t0, t1;
	fftComplex exp1, exp3;

	//Stage 1 (this is similar to stage 2 for DIT)
	for (int k = 0; k < ldn - 1; k++) {

		N2 = N2 / 2;
		N4 = N2 / 4;

		//Treat case j = 0 separately due to simpler multiplications
		for (ix = 0, id = N2 * 2; ix < N; id *= 4) {

			for (i0 = ix; i0 < N; i0 += id) {

				i1 = i0 + N4;
				i2 = i1 + N4;
				i3 = i2 + N4;

				t0 = Fx[i0] - Fx[i2];
				t1 = Fx[i1] - Fx[i3];

				Fx[i0] = Fx[i0] + Fx[i2];
				Fx[i1] = Fx[i1] + Fx[i3];

				Fx[i2] = t0 + !t1;
				Fx[i3] = t0 - !t1;
			}

			ix = 2 * id - N2;
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			//Note these are not conjugated (as you would think they should be since this is IFFT) since we are doing DIF, not DIT - butterflies require it
			exp1 = cossin[(MAXFFT / N2) * j];
			exp3 = cossin[(MAXFFT / N2) * 3 * j];

			for (ix = j, id = N2 * 2; ix < N; id *= 4) {

				for (i0 = ix; i0 < N; i0 += id) {

					i1 = i0 + N4;
					i2 = i1 + N4;
					i3 = i2 + N4;

					t0 = Fx[i0] - Fx[i2];
					t1 = Fx[i1] - Fx[i3];

					Fx[i0] = Fx[i0] + Fx[i2];
					Fx[i1] = Fx[i1] + Fx[i3];

					Fx[i2] = (t0 + !t1) * exp1;
					Fx[i3] = (t0 - !t1) * exp3;
				}

				ix = 2 * id - N2 + j;
			}
		}
	}

	//Stage 2 (same as stage 1 for DIT)
	for (ix = 0, id = 4; ix < N; id *= 4) {

		for (i0 = ix; i0 < N; i0 += id) {

			t0 = Fx[i0] - Fx[(i0 + 1)];

			Fx[i0] = Fx[i0] + Fx[(i0 + 1)];
			Fx[(i0 + 1)] = t0;
		}

		ix = 2 * id - 2;
	}

	//Unshuffle, moving from input data to output, and divide by N for IFFT
	for (int p = 1; p < N; p++) {

		FX[p] = Fx[fftShuffle[N - 1 + p]] / N;
	}

	FX[0] = Fx[0] / N;
	FX[(N - 1)] = Fx[(N - 1)] / N;
}