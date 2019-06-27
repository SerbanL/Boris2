#pragma once

#include "FFT.h"

//------------------- REAL INPUT FFT / REAL OUTPUT IFFT (fftComplex and fftComplex3 versions)

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::RealfromComplexFFT(fftComplex *Fx, fftComplex *FX, int N) 
{
	//Use this after doing a half-length complex FFT on a real input (should have used CopyRealShuffle or CopyRealShuffleZeroPad before doing the FFT as well).
	//FFT result is in Fx (N points, N is the half-length complex FFT length (real input would have had 2N points)).
	//This method will write N + 1 values in FX. (Herimitian symmetric property will give the remaining N - 1 points. Also note n = 0 and n = N points have zero imaginary parts.)

	fftComplex Yn, YNmn;

	//n = 0
	FX[0] = fftComplex(Fx[0].Re + Fx[0].Im, 0);

	for (int n = 1; n < N; n++) {

		Yn = Fx[n];
		YNmn = Fx[N - n];

		FX[n] = (Yn + ~YNmn - sincos[(MAXFFT / (2 * N)) * n] * (Yn - ~YNmn)) * 0.5;
	}

	//n = N
	FX[N] = fftComplex(Fx[0].Re - Fx[0].Im, 0);
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::RealfromComplexFFT(fftComplex3 *Fx, fftComplex3 *FX, int N) 
{
	fftComplex3 Yn, YNmn;

	//n = 0
	FX[0] = fftComplex3(fftComplex(Fx[0].x.Re + Fx[0].x.Im, 0), fftComplex(Fx[0].y.Re + Fx[0].y.Im, 0), fftComplex(Fx[0].z.Re + Fx[0].z.Im, 0));

	for (int n = 1; n < N; n++) {

		Yn = Fx[n];
		YNmn = Fx[N - n];

		FX[n] = (Yn + ~YNmn - ((Yn - ~YNmn) * sincos[(MAXFFT / (2 * N)) * n])) * 0.5;
	}

	//n = N
	FX[N] = fftComplex3(fftComplex(Fx[0].x.Re - Fx[0].x.Im, 0), fftComplex(Fx[0].y.Re - Fx[0].y.Im, 0), fftComplex(Fx[0].z.Re - Fx[0].z.Im, 0));
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::RealfromComplexIFFT(fftComplex *Fx, fftComplex *FX, int N) 
{
	//Use this before doing a half-length complex IFFT on complex input to get real output (should also use UnShuffleRealTruncate or UnShuffleReal after the IFFT).
	//N is the half-length complex IFFT length (real output would have had 2N points)
	//This method will write N values in FX using N + 1 values in Fx.

	for (int n = 0; n < N; n++) {

		fftComplex Xn = Fx[n];
		fftComplex XNmn = Fx[N - n];

		fftComplex Yn = (Xn + ~XNmn - ~sincos[(MAXFFT / (2 * N)) * n] * (Xn - ~XNmn)) * 0.5;

		FX[n] = Yn;
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::RealfromComplexIFFT(fftComplex3 *Fx, fftComplex3 *FX, int N) 
{
	//Use this before doing a half-length complex IFFT on complex input to get real output (should also use UnShuffleRealTruncate or UnShuffleReal after the IFFT).
	//N is the half-length complex IFFT length (real output would have had 2N points)
	//This method will write N values in FX using N + 1 values in Fx.

	for (int n = 0; n < N; n++) {

		fftComplex3 Xn = Fx[n];
		fftComplex3 XNmn = Fx[N - n];

		fftComplex3 Yn = (Xn + ~XNmn - ((Xn - ~XNmn) * ~sincos[(MAXFFT / (2 * N)) * n])) * 0.5;

		FX[n] = Yn;
	}
}

//------------------- SHUFFLING (before DIT)

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::CopyShuffle(fftComplex *Fx, fftComplex *FX, int stride, int N) 
{
	//Shuffle values moving input data to output

	FX[0] = Fx[0];
	FX[N - 1] = Fx[(N - 1) * stride];

	for (int p = 1; p <= N - 2; p++) {

		FX[p] = Fx[fftShuffle[N - 1 + p] * stride];
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::CopyShuffle(fftComplex3 *Fx, fftComplex3 *FX, int stride, int N) 
{
	//Shuffle values moving input data to output

	FX[0] = Fx[0];
	FX[N - 1] = Fx[(N - 1) * stride];

	for (int p = 1; p <= N - 2; p++) {

		FX[p] = Fx[fftShuffle[N - 1 + p] * stride];
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::CopyRealShuffle(fftReal *Rex, fftComplex *FX, int N) 
{
	//Shuffle values moving input data to output
	//This routine is for real input to prepare half-length complex FFT : even indexed input values are set in Real part, odd indexed input values set in Complex part
	//we have 2N real points so N is the length of the complex fft

	//INVERSE METHOD
	for (int p = 0; p < N; p++) {

		int index = fftShuffleInv[N - 1 + p];

		FX[index] = fftComplex(Rex[2 * p], Rex[(2 * p + 1)]);
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::CopyRealShuffle(fftReal3 *Rex, fftComplex3 *FX, int N) 
{
	//Shuffle values moving input data to output
	//This routine is for real input to prepare half-length complex FFT : even indexed input values are set in Real part, odd indexed input values set in Complex part
	//we have 2N real points so N is the length of the complex fft

	//INVERSE METHOD
	for (int p = 0; p < N; p++) {

		int index = fftShuffleInv[N - 1 + p];

		FX[index] = fftComplex3(Rex[2 * p], Rex[(2 * p + 1)]);
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::CopyRealShuffle(fftReal *Rex, fftComplex *FX, int stride, int N)
{
	//Shuffle values moving input data to output
	//This routine is for real input to prepare half-length complex FFT : even indexed input values are set in Real part, odd indexed input values set in Complex part
	//we have 2N real points so N is the length of the complex fft

	//INVERSE METHOD
	for (int p = 0; p < N; p++) {

		int index = fftShuffleInv[N - 1 + p];

		FX[index] = fftComplex(Rex[2 * p * stride], Rex[(2 * p + 1) * stride]);
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::CopyRealShuffle(fftReal3 *Rex, fftComplex3 *FX, int stride, int N)
{
	//Shuffle values moving input data to output
	//This routine is for real input to prepare half-length complex FFT : even indexed input values are set in Real part, odd indexed input values set in Complex part
	//we have 2N real points so N is the length of the complex fft

	//INVERSE METHOD
	for (int p = 0; p < N; p++) {

		int index = fftShuffleInv[N - 1 + p];

		FX[index] = fftComplex3(Rex[2 * p * stride], Rex[(2 * p + 1) * stride]);
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::CopyShuffleZeroPad(fftComplex *Fx, fftComplex *FX, int stride, int N) 
{
	//Shuffle values moving input data to output
	//Zero padding on input data from N/2 upwards - results in odd points being zero in shuffled data


	//DIRECT METHOD
	//for(int p = 0; p <= N-2; p += 2) {
	//
	//	FX[p] = Fx[fftShuffle[N - 1 + p] * stride];
	//	FX[p + 1] = fftComplex(0, 0);
	//}
	
	//INVERSE METHOD
	for (int p = 0; p < N / 2; p++) {

		int index = fftShuffleInv[N - 1 + p];

		FX[index] = Fx[p * stride];

		//odd indexed shuffled inputs must be zero due to zero padding (note first half of non-shuffled inputs are always distributed to even-indexed shuffled values, thus index + 1 is odd)
		FX[index + 1] = fftComplex(0, 0);
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::CopyShuffleZeroPad(fftComplex3 *Fx, fftComplex3 *FX, int stride, int N) 
{
	//DIRECT METHOD
	//for(int p = 0; p <= N-2; p += 2) {
	//
	//	FX[p] = Fx[fftShuffle[N - 1 + p] * stride];
	//	FX[p + 1] = fftComplex3(0, 0);
	//}
	
	//INVERSE METHOD
	for (int p = 0; p < N / 2; p++) {

		int index = fftShuffleInv[N - 1 + p];

		FX[index] = Fx[p * stride];

		//odd indexed shuffled inputs must be zero due to zero padding (note first half of non-shuffled inputs are always distributed to even-indexed shuffled values, thus index + 1 is odd)
		FX[index + 1] = fftComplex3();
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::CopyRealShuffleZeroPad(fftReal *Rex, fftComplex *FX, int N) 
{
	//Shuffle values moving input data to output
	//Zero padding on input data from N/2 upwards - results in odd points being zero in shuffled data
	//This routine is for real input to prepare half-length complex FFT : even indexed input values are set in Real part, odd indexed input values set in Complex part
	//we have 2N real points (but the upper N are taken to be zero) so N is the length of the complex fft

	//INVERSE METHOD
	for (int p = 0; p < N / 2; p++) {

		int index = fftShuffleInv[N - 1 + p];

		FX[index] = fftComplex(Rex[2 * p], Rex[(2 * p + 1)]);

		//odd indexed shuffled inputs must be zero due to zero padding (note first half of non-shuffled inputs are always distributed to even-indexed shuffled values, thus index + 1 is odd)
		FX[index + 1] = fftComplex(0, 0);
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::CopyRealShuffleZeroPad(fftReal3 *Rex, fftComplex3 *FX, int n, int N) 
{
	//Shuffle values moving input data to output
	//Zero padding on input data from N/2 upwards - results in odd points being zero in shuffled data
	//This routine is for real input to prepare half-length complex FFT : even indexed input values are set in Real part, odd indexed input values set in Complex part
	//we have 2N real points (but the upper N are taken to be zero) so N is the length of the complex fft

	//n is the dimension of Rex along x. N is a number greater or equal to n, which is a power of 2.
	//n could be smaller than N, in which case fill the remaining space with zeroes (i.e. as if Rex has dimension N/2 along x, with zero values from n to N/2 - 1.

	//INVERSE METHOD
	for (int p = 0; p < N / 2; p++) {

		int index = fftShuffleInv[N - 1 + p];

		fftReal3 val1, val2; //these are initialized to zero by default in their constructor

		if (2 * p < n) val1 = Rex[2 * p];
		if (2 * p + 1 < n) val2 = Rex[2 * p + 1];

		FX[index] = fftComplex3(val1, val2);

		//odd indexed shuffled inputs must be zero due to zero padding (note first half of non-shuffled inputs are always distributed to even-indexed shuffled values, thus index + 1 is odd)
		FX[index + 1] = fftComplex3();
	}
}

//------------------- UN-SHUFFLING (after DIF)

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::UnShuffle(fftComplex *Fx, fftComplex *FX, int stride, int N) 
{
	//Unshuffle values, moving input data to output placed using stride.

	//DIRECT METHOD
	//for(int p = 0; p < N; p++) {
	//
	//	int index = fftShuffle[N - 1 + p];
	//	
	//	FX[p * stride] = Fx[index];
	//}

	//INVERSE METHOD
	for (int p = 0; p < N; p++) {

		int index = fftShuffleInv[N - 1 + p];

		FX[index * stride] = Fx[p];
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::UnShuffle(fftComplex3 *Fx, fftComplex3 *FX, int stride, int N) 
{
	//INVERSE METHOD
	//Even-indexed inputs are always distributed to lower half of output de-shuffled (or shuffled) values : we only want the lower half.
	for (int p = 0; p < N; p++) {

		int index = fftShuffleInv[N - 1 + p];

		FX[index * stride] = Fx[p];
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::UnShuffleTruncate(fftComplex *Fx, fftComplex *FX, int stride, int N) 
{
	//Unshuffle values, moving input data to output placed using stride. Only half the outputs are written.

	//DIRECT METHOD
	//for(int p = 0; p < N/2; p++) {
	//
	//	int index = fftShuffle[N - 1 + p];
	//	
	//	FX[p * stride] = Fx[index];
	//}

	//INVERSE METHOD
	//Even-indexed inputs are always distributed to lower half of output de-shuffled (or shuffled) values : we only want the lower half.
	for (int p = 0; p < N; p += 2) {

		int index = fftShuffleInv[N - 1 + p];

		FX[index * stride] = Fx[p];
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::UnShuffleTruncate(fftComplex3 *Fx, fftComplex3 *FX, int stride, int N) 
{
	//INVERSE METHOD
	//Even-indexed inputs are always distributed to lower half of output de-shuffled (or shuffled) values : we only want the lower half.
	for (int p = 0; p < N; p += 2) {

		int index = fftShuffleInv[N - 1 + p];

		FX[index * stride] = Fx[p];
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::UnShuffleRealTruncate(fftComplex *Fx, fftReal *ReX, int n, int N) 
{
	//INVERSE METHOD
	//Even-indexed inputs are always distributed to lower half of output de-shuffled (or shuffled) values : we only want the lower half.
	//This is routine for real output from half-length complex IFFT :real output values are set as even indexed real output, complex output values are set as odd indexed real output 
	//We have 2N real output points (so the IFFT which should have preceeded the call to this method should have been a N-point complex IFFT), but we only write the first N points out and ignore the upper half.

	for (int p = 0; p < N; p += 2) {

		int index = fftShuffleInv[N - 1 + p];

		if (2 * index < n)	ReX[2 * index] = Fx[p].Re;
		if (2 * index + 1 < n)	ReX[2 * index + 1] = Fx[p].Im;
	}
}

template <typename fftReal>
void FFTMethods_Cpp<fftReal>::UnShuffleRealTruncate(fftComplex3 *Fx, fftReal3 *ReX, int n, int N) 
{
	//INVERSE METHOD
	//Even-indexed inputs are always distributed to lower half of output de-shuffled (or shuffled) values : we only want the lower half.
	//This is routine for real output from half-length complex IFFT :real output values are set as even indexed real output, complex output values are set as odd indexed real output 
	//We have 2N real output points (so the IFFT which should have preceeded the call to this method should have been a N-point complex IFFT), but we only write the first N points out and ignore the upper half.

	//n is the dimension of ReX along x. N is a number greater or equal to n, which is a power of 2.
	//n could be smaller than N, in which case only extract values with indices up to n - 1.

	for (int p = 0; p < N; p += 2) {

		int index = fftShuffleInv[N - 1 + p];

		if (2 * index < n)	ReX[2 * index] = Fx[p].Re();
		if (2 * index + 1 < n)	ReX[2 * index + 1] = Fx[p].Im();
	}
}

template <typename fftReal>
fftReal FFTMethods_Cpp<fftReal>::UnShuffleRealTruncate_FinishConvolution_Set(fftComplex3 *Fx, fftReal3 *ReX, fftReal3 *ReIn, int n, int N, fftReal divisor)
{
	fftReal dot_prod = 0;

	for (int p = 0; p < N; p += 2) {

		int index = fftShuffleInv[N - 1 + p];

		if (2 * index < n) {

			fftReal3 value = Fx[p].Re() / divisor;
			//ReX could be same as ReIn so first find dot product before writing result to ReX
			dot_prod += value * ReIn[2 * index];
			ReX[2 * index] = value;
		}

		if (2 * index + 1 < n) {

			fftReal3 value = Fx[p].Im() / divisor;
			dot_prod += value * ReIn[2 * index + 1];
			ReX[2 * index + 1] = value;
		}
	}

	return dot_prod;
}

template <typename fftReal>
fftReal FFTMethods_Cpp<fftReal>::UnShuffleRealTruncate_FinishConvolution_Add(fftComplex3 *Fx, fftReal3 *ReX, fftReal3 *ReIn, int n, int N, fftReal divisor)
{
	fftReal dot_prod = 0;

	for (int p = 0; p < N; p += 2) {

		int index = fftShuffleInv[N - 1 + p];

		if (2 * index < n) {

			fftReal3 value = Fx[p].Re() / divisor;
			//ReX could be same as ReIn so first find dot product before writing result to ReX
			dot_prod += value * ReIn[2 * index];
			ReX[2 * index] += value;
		}

		if (2 * index + 1 < n) {

			fftReal3 value = Fx[p].Im() / divisor;
			dot_prod += value * ReIn[2 * index + 1];
			ReX[2 * index + 1] += value;
		}
	}

	return dot_prod;
}