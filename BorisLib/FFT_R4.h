#pragma once

#include "FFT.h"

//------------------- RADIX-4

//In-place FFT, Radix 4 Decimation in Time. Input data must be shuffled prior to calling this.
template <typename fftReal>
template <typename Type>
void FFTMethods_Cpp<fftReal>::FFT_Radix4_DIT(Type *FX, int ldn, int N) 
{
	int i0, i1, i2, i3, p;
	int N2, N4;

	Type t0, t1, t2, t3;
	fftComplex exp, exp2, exp3;

	p = (ldn & 1);

	// n is not a power of 4, need a radix-2 step
	if (p != 0) {

		for (i0 = 0; i0 < N; i0 += 2) {

			t0 = FX[i0] - FX[(i0 + 1)];
			FX[i0] = FX[i0] + FX[(i0 + 1)];
			FX[(i0 + 1)] = t0;
		}
	}

	N2 = (1 << p);

	for (p = p + 2; p <= ldn; p += 2) {

		N4 = N2;
		N2 *= 4;

		//Special case j = 0 separated from loop : no multiplications.
		for (int r = 0; r < N; r += N2) {

			i1 = r + N4;
			i2 = i1 + N4;
			i3 = i2 + N4;

			t0 = (FX[r] + FX[i1]);
			t1 = (FX[r] - FX[i1]);
			t2 = (FX[i2] + FX[i3]);
			t3 = !(FX[i3] - FX[i2]);

			FX[r] = t0 + t2;
			FX[i1] = t1 + t3;
			FX[i2] = t0 - t2;
			FX[i3] = t1 - t3;
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			exp = ~cossin[(MAXFFT / N2) * j];
			exp2 = ~cossin[(MAXFFT / N2) * 2 * j];
			exp3 = ~cossin[(MAXFFT / N2) * 3 * j];

			for (int r = 0; r < N; r += N2) {

				i0 = j + r;
				i1 = i0 + N4;
				i2 = i1 + N4;
				i3 = i2 + N4;

				t0 = (FX[i0] + (FX[i1] * exp2));
				t1 = (FX[i0] - (FX[i1] * exp2));
				t2 = ((FX[i2] * exp) + (FX[i3] * exp3));
				t3 = !((FX[i3] * exp3) - (FX[i2] * exp));

				FX[i0] = t0 + t2;
				FX[i1] = t1 + t3;
				FX[i2] = t0 - t2;
				FX[i3] = t1 - t3;
			}
		}
	}
}

//In-place FFT, Radix 4 Decimation in Frequency. Output data is shuffled, so must unshuffle after calling this.
template <typename fftReal>
template <typename Type>
void FFTMethods_Cpp<fftReal>::FFT_Radix4_DIF(Type *FX, int ldn, int N)
{
	int i0, i1, i2, i3, p;
	int N2, N4;

	Type t0, t1, t2, t3;
	fftComplex exp, exp2, exp3;

	N2 = 4 * N;

	for (p = ldn; p >= 2; p -= 2) {

		N2 /= 4;
		N4 = N2 / 4;

		//Special case j = 0 separated from loop : no multiplications.
		for (int r = 0; r < N; r += N2) {

			i1 = r + N4;
			i2 = i1 + N4;
			i3 = i2 + N4;

			t0 = (FX[r] + FX[i2]);
			t1 = (FX[r] - FX[i2]);
			t2 = (FX[i1] + FX[i3]);
			t3 = !(FX[i3] - FX[i1]);

			FX[r] = (t0 + t2);
			FX[i1] = (t0 - t2);
			FX[i2] = (t1 + t3);
			FX[i3] = (t1 - t3);
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			exp = ~cossin[(MAXFFT / N2) * j];
			exp2 = ~cossin[(MAXFFT / N2) * 2 * j];
			exp3 = ~cossin[(MAXFFT / N2) * 3 * j];

			for (int r = 0; r < N; r += N2) {

				i0 = j + r;
				i1 = i0 + N4;
				i2 = i1 + N4;
				i3 = i2 + N4;

				t0 = (FX[i0] + FX[i2]);
				t1 = (FX[i0] - FX[i2]);
				t2 = (FX[i1] + FX[i3]);
				t3 = !(FX[i3] - FX[i1]);

				FX[i0] = (t0 + t2);
				FX[i1] = (t0 - t2) * exp2;
				FX[i2] = (t1 + t3) * exp;
				FX[i3] = (t1 - t3) * exp3;
			}
		}
	}

	p = (ldn & 1);

	// n is not a power of 4, need a radix-2 step
	if (p != 0) {

		for (i0 = 0; i0 < N; i0 += 2) {

			t0 = FX[i0] - FX[(i0 + 1)];
			FX[i0] = FX[i0] + FX[(i0 + 1)];
			FX[(i0 + 1)] = t0;
		}
	}
}

//In-place IFFT, Radix 4 Decimation in Time. Input data must be shuffled prior to calling this, output data not divided by N, so divide by N after calling this.
template <typename fftReal>
template <typename Type>
void FFTMethods_Cpp<fftReal>::IFFT_Radix4_DIT(Type *FX, int ldn, int N)
{
	int i0, i1, i2, i3, p;
	int N2, N4;

	Type t0, t1, t2, t3;
	fftComplex exp, exp2, exp3;

	p = (ldn & 1);

	// n is not a power of 4, need a radix-2 step
	if (p != 0) {

		for (i0 = 0; i0 < N; i0 += 2) {

			t0 = FX[i0] - FX[(i0 + 1)];
			FX[i0] = FX[i0] + FX[(i0 + 1)];
			FX[(i0 + 1)] = t0;
		}
	}

	N2 = (1 << p);

	for (p = p + 2; p <= ldn; p += 2) {

		N4 = N2;
		N2 *= 4;

		//Special case j = 0 separated from loop : no multiplications.
		for (int r = 0; r < N; r += N2) {

			i1 = r + N4;
			i2 = i1 + N4;
			i3 = i2 + N4;

			t0 = (FX[r] + FX[i1]);
			t1 = (FX[r] - FX[i1]);
			t2 = (FX[i2] + FX[i3]);
			t3 = !(FX[i3] - FX[i2]);

			FX[r] = t0 + t2;
			FX[i1] = t1 - t3;
			FX[i2] = t0 - t2;
			FX[i3] = t1 + t3;
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			exp = cossin[(MAXFFT / N2) * j];
			exp2 = cossin[(MAXFFT / N2) * 2 * j];
			exp3 = cossin[(MAXFFT / N2) * 3 * j];

			for (int r = 0; r < N; r += N2) {

				i0 = j + r;
				i1 = i0 + N4;
				i2 = i1 + N4;
				i3 = i2 + N4;

				t0 = (FX[i0] + FX[i1] * exp2);
				t1 = (FX[i0] - FX[i1] * exp2);
				t2 = (FX[i2] * exp + FX[i3] * exp3);
				t3 = !(FX[i3] * exp3 - FX[i2] * exp);

				FX[i0] = t0 + t2;
				FX[i1] = t1 - t3;
				FX[i2] = t0 - t2;
				FX[i3] = t1 + t3;
			}
		}
	}
}

//In-place IFFT, Radix 4 Decimation in Frequency. Output data is shuffled, so unshuffle after calling this, and not divided by N, so divide by N after calling this.
template <typename fftReal>
template <typename Type>
void FFTMethods_Cpp<fftReal>::IFFT_Radix4_DIF(Type *FX, int ldn, int N) 
{
	int i0, i1, i2, i3, p;
	int N2, N4;

	Type t0, t1, t2, t3;
	fftComplex exp, exp2, exp3;

	N2 = 4 * N;

	for (p = ldn; p >= 2; p -= 2) {

		N2 /= 4;
		N4 = N2 / 4;

		//Special case j = 0 separated from loop : no multiplications.
		for (int r = 0; r < N; r += N2) {

			i1 = r + N4;
			i2 = i1 + N4;
			i3 = i2 + N4;

			t0 = (FX[r] + FX[i2]);
			t1 = (FX[r] - FX[i2]);
			t2 = (FX[i1] + FX[i3]);
			t3 = !(FX[i3] - FX[i1]);

			FX[r] = (t0 + t2);
			FX[i1] = (t0 - t2);
			FX[i2] = (t1 - t3);
			FX[i3] = (t1 + t3);
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			exp = cossin[(MAXFFT / N2) * j];
			exp2 = cossin[(MAXFFT / N2) * 2 * j];
			exp3 = cossin[(MAXFFT / N2) * 3 * j];

			for (int r = 0; r < N; r += N2) {

				i0 = j + r;
				i1 = i0 + N4;
				i2 = i1 + N4;
				i3 = i2 + N4;

				t0 = (FX[i0] + FX[i2]);
				t1 = (FX[i0] - FX[i2]);
				t2 = (FX[i1] + FX[i3]);
				t3 = !(FX[i3] - FX[i1]);

				FX[i0] = (t0 + t2);
				FX[i1] = (t0 - t2) * exp2;
				FX[i2] = (t1 - t3) * exp;
				FX[i3] = (t1 + t3) * exp3;
			}
		}
	}

	p = (ldn & 1);

	// n is not a power of 4, need a radix-2 step
	if (p != 0) {

		for (i0 = 0; i0 < N; i0 += 2) {

			t0 = FX[i0] - FX[(i0 + 1)];
			FX[i0] = FX[i0] + FX[(i0 + 1)];
			FX[(i0 + 1)] = t0;
		}
	}
}

//------------------ - RADIX-4 - PACKED VERSIONS

//Split-Radix FFT with decimation-in-time. Packed means no pre-shuffling needed as it's done by this function : just call to do the FFT
template <typename fftReal>
template <typename Type>
void FFTMethods_Cpp<fftReal>::FFT_RAD4DIT_Packed(Type *Fx, Type *FX, int ldn, int N)
{
	int i0, i1, i2, i3, p;
	int N2, N4;

	Type t0, t1, t2, t3;
	fftComplex exp, exp2, exp3;

	//Shuffle values, moving input data from input to output : input data is left unaffected.
	FX[0] = Fx[0];
	FX[(N - 1)] = Fx[(N - 1)];

	for (p = 1; p <= N - 2; p++) {

		FX[p] = Fx[fftShuffle[N - 1 + p]];
	}

	p = (ldn & 1);

	// n is not a power of 4, need a radix-2 step
	if (p != 0) {

		for (i0 = 0; i0 < N; i0 += 2) {

			t0 = FX[i0] - FX[(i0 + 1)];
			FX[i0] = FX[i0] + FX[(i0 + 1)];
			FX[(i0 + 1)] = t0;
		}
	}

	N2 = (1 << p);

	for (p = p + 2; p <= ldn; p += 2) {

		N4 = N2;
		N2 *= 4;

		//Special case j = 0 separated from loop : no multiplications.
		for (int r = 0; r < N; r += N2) {

			i1 = r + N4;
			i2 = i1 + N4;
			i3 = i2 + N4;

			t0 = (FX[r] + FX[i1]);
			t1 = (FX[r] - FX[i1]);
			t2 = (FX[i2] + FX[i3]);
			t3 = !(FX[i3] - FX[i2]);

			FX[r] = t0 + t2;
			FX[i1] = t1 + t3;
			FX[i2] = t0 - t2;
			FX[i3] = t1 - t3;
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			exp = ~cossin[(MAXFFT / N2) * j];
			exp2 = ~cossin[(MAXFFT / N2) * 2 * j];
			exp3 = ~cossin[(MAXFFT / N2) * 3 * j];

			for (int r = 0; r < N; r += N2) {

				i0 = j + r;
				i1 = i0 + N4;
				i2 = i1 + N4;
				i3 = i2 + N4;

				t0 = (FX[i0] + FX[i1] * exp2);
				t1 = (FX[i0] - FX[i1] * exp2);
				t2 = (FX[i2] * exp + FX[i3] * exp3);
				t3 = !(FX[i3] * exp3 - FX[i2] * exp);

				FX[i0] = t0 + t2;
				FX[i1] = t1 + t3;
				FX[i2] = t0 - t2;
				FX[i3] = t1 - t3;
			}
		}
	}
}

//Split-Radix IFFT with decimation-in-frequency. Packed means no un-shuffling needed as it's done by this function : just call to do the IFFT
template <typename fftReal>
template <typename Type>
void FFTMethods_Cpp<fftReal>::IFFT_RAD4DIF_Packed(Type *Fx, Type *FX, int ldn, int N)
{
	int i0, i1, i2, i3, p;
	int N2, N4;

	Type t0, t1, t2, t3;
	fftComplex exp, exp2, exp3;

	N2 = 4 * N;

	for (p = ldn; p >= 2; p -= 2) {

		N2 /= 4;
		N4 = N2 / 4;

		//Special case j = 0 separated from loop : no multiplications.
		for (int r = 0; r < N; r += N2) {

			i1 = r + N4;
			i2 = i1 + N4;
			i3 = i2 + N4;

			t0 = (Fx[r] + Fx[i2]);
			t1 = (Fx[r] - Fx[i2]);
			t2 = (Fx[i1] + Fx[i3]);
			t3 = !(Fx[i3] - Fx[i1]);

			//Note change in sign for t3 for IFFT compared to FFT
			Fx[r] = (t0 + t2);
			Fx[i1] = (t0 - t2);
			Fx[i2] = (t1 - t3);
			Fx[i3] = (t1 + t3);
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			//Conjugated roots since IFFT
			exp = cossin[(MAXFFT / N2) * j];
			exp2 = cossin[(MAXFFT / N2) * 2 * j];
			exp3 = cossin[(MAXFFT / N2) * 3 * j];

			for (int r = 0; r < N; r += N2) {

				i0 = j + r;
				i1 = i0 + N4;
				i2 = i1 + N4;
				i3 = i2 + N4;

				t0 = (Fx[i0] + Fx[i2]);
				t1 = (Fx[i0] - Fx[i2]);
				t2 = (Fx[i1] + Fx[i3]);
				t3 = !(Fx[i3] - Fx[i1]);

				//Note change in sign for t3 for IFFT compared to FFT
				Fx[i0] = (t0 + t2);
				Fx[i1] = (t0 - t2) * exp2;
				Fx[i2] = (t1 - t3) * exp;
				Fx[i3] = (t1 + t3) * exp3;
			}
		}
	}

	p = (ldn & 1);

	// n is not a power of 4, need a radix-2 step
	if (p != 0) {

		for (i0 = 0; i0 < N; i0 += 2) {

			t0 = Fx[i0] - Fx[(i0 + 1)];
			Fx[i0] = Fx[i0] + Fx[(i0 + 1)];
			Fx[(i0 + 1)] = t0;
		}
	}

	//Shuffle values, moving input data from input to output. Also divide by N since this is IFFT
	FX[0] = Fx[0] / N;
	FX[(N - 1)] = Fx[(N - 1)] / N;

	for (p = 1; p <= N - 2; p++) {

		FX[p] = Fx[fftShuffle[N - 1 + p]] / N;
	}
}