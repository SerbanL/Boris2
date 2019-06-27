#pragma once

#include "BorisLib.h"

#include <omp.h>

#define MAXFFT 32768								//Maximum length of FFT - used to initialize the cos and sin look-up tables

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//FFT Methods for demag calculation using convolution -> uses cpp routines

class FFTMethods_Cpp {

private: //private data

	vector<ReIm> cossin, sincos;

	//look-up tables for fft shuffling. First index, m, is the power of 2 of the fft length. On that row, the entries (2^m of them) give the shuffling indices.
	//fftShuffle gives indices for input data when output data is parsed linearly, fftShuffleInv is the inverse (input data parsed linearly, gives indices for output data)
	vector<int> fftShuffle, fftShuffleInv;

	int OmpThreads;

private: //private methods

	int ShuffleIndex_Radix2(int n, int m);

protected: //protected methods

	//-----------------------------------------------------------------------------------------------
	// DBL3 and ReIm3 methods

	void CopyRealShuffleZeroPad(DBL3 *Rex, ReIm3 *FX, int n, int N);
	void CopyShuffleZeroPad(ReIm3 *Fx, ReIm3 *FX, int stride, int N);
	void CopyRealShuffle(DBL3 *Rex, ReIm3 *FX, int N);
	void CopyShuffle(ReIm3 *Fx, ReIm3 *FX, int stride, int N);

	void FFT_Radix4_DIT(ReIm3 *FX, int ldn, int N);

	void RealfromComplexFFT(ReIm3 *Fx, ReIm3 *FX, int N);

	void IFFT_Radix4_DIF(ReIm3 *FX, int ldn, int N);

	void RealfromComplexIFFT(ReIm3 *Fx, ReIm3 *FX, int N);

	void UnShuffleTruncate(ReIm3 *Fx, ReIm3 *FX, int stride, int N);
	void UnShuffleRealTruncate(ReIm3 *Fx, DBL3 *ReX, int n, int N);

	//-----------------------------------------------------------------------------------------------
	// Single component array versions

	///
	void CopyShuffle(double *Rex, double *Imx, ReIm *FX, int stride, int N);		//VERIFIED, CORRECT - use to bit-reverse inputs before calling a DIT fft
	void CopyShuffle(ReIm *Fx, ReIm *FX, int stride, int N);						//VERIFIED, CORRECT - use to bit-reverse inputs before calling a DIT fft
	void CopyRealShuffle(double *Rex, ReIm *FX, int N);								//VERIFIED, CORRECT - use to bit-reverse inputs before calling a half-length complex DIT fft: if we have 2N real points then place them alternately in real and complex parts and do an N-point complex FFT
	
	void CopyShuffleZeroPad(ReIm *Fx, ReIm *FX, int stride, int N);					//VERIFIED, CORRECT - use to bit-reverse inputs before calling a DIT fft. Upper half of inputs are taken as zero, so the input need only have N/2 points stored.
	void CopyShuffleZeroPad(double *Rex, double *Imx, ReIm *FX, int N);				//VERIFIED, CORRECT - use to bit-reverse inputs before calling a DIT fft. Upper half of inputs are taken as zero, so the input need only have N/2 points stored.
	void CopyRealShuffleZeroPad(double *Rex, ReIm *FX, int N);						//VERIFIED, CORRECT - use to bit-reverse inputs before calling a half-length complex DIT fft: if we have 2N real points then place them alternately in real and complex parts and do an N-point complex FFT. On top of this the upper half of inputs are taken as zero, so the input need only have N points stored.
	void CopyRealProdShuffleZeroPad(double *ReS, double *Rex, ReIm *FX, int N);		//VERIFIED, CORRECT - use to bit-reverse inputs before calling a half-length complex DIT fft: if we have 2N real points then place them alternately in real and complex parts and do an N-point complex FFT. On top of this the upper half of inputs are taken as zero, so the input need only have N points stored.
	
	void UnShuffleTruncate(ReIm *Fx, ReIm *FX, int stride, int N);					//VERIFIED, CORRECT - use to bit-reverse outputs after calling a DIF fft. This method just writes the first half of outputs and ignores the second half.
	void UnShuffleTruncate(ReIm *Fx, double *ReX, double *ImX, int N);				//VERIFIED, CORRECT - use to bit-reverse outputs after calling a DIF fft. This method just writes the first half of outputs and ignores the second half.
	void UnShuffleRealTruncate(ReIm *Fx, double *ReX, int N);						//VERIFIED, CORRECT - use to bit-reverse outputs after calling a DIF fft.

	///
	void RealfromComplexFFT(ReIm *Fx, ReIm *FX, int N);								//VERIFIED, CORRECT - use to obtain the fft of a real input after using the half-length complex fft method (i.e. after CopyRealShuffle to place the input in alternate real and imaginary parts, then doing a DIT fft of half-length on this)
	void RealfromComplexIFFT(ReIm *Fx, ReIm *FX, int N);							//VERIFIED, CORRECT - use before running a half-length ifft which will have only real output. Here N is the length of the half-length IFFT, and we need N + 1 input points (don't need to store 2N points due to Hermitian symmetric property). After this method, run a DFT IFFT followed by UnShuffleRealTruncate

	///

	//Radix-2 FFT, packed real and imaginary parts for inputs and/or outputs. 
	void FFT_Rad2DIT_Packed(ReIm *Fx, ReIm *FX, int ldn, int N, int mult);			//VERIFIED, CORRECT

	//Split-Radix FFT : lower multiplications count compared to Radix4 (~ 10%), same additions count, but much larger number of reads and writes (almost double!!!).
	//My tests showed these are slower compared to Radix4 method, including on odd powers of 2.
	void FFT_SRDIT_Packed(ReIm *Fx, ReIm *FX, int ldn, int N, int mult);			//VERIFIED, CORRECT

	void IFFT_SRDIF_Packed(ReIm *Fx, ReIm *FX, int ldn, int N, int mult);			//VERIFIED, CORRECT

	//Radix-4 FFT and IFFT, packed real and imaginary parts for inputs and/or outputs. Uses a Radix-2 step for odd powers of 2.
	//NOTE : these could be optimized further :
	//The special case j = 0, giving (exp, exp2 and exp3) factors of (1, 1, 1), I've already separated out, easy.
	//Other cases that could be optimized include multiplications by +/- 1, and +/- i. 
	//These can be separated out from the main loop by tracking special valus of N2 and j for which they occur and writing them out specifically (no extra if conditions used, just rewrite the loops)
	//This will result in a MUCH larger FFT function but should run faster (but probably not significantly faster).
	void FFT_RAD4DIT_Packed(ReIm *Fx, ReIm *FX, int ldn, int N, int mult);			//VERIFIED, CORRECT
	
	void IFFT_RAD4DIF_Packed(ReIm *Fx, ReIm *FX, int ldn, int N, int mult);			//VERIFIED, CORRECT

	//Methods without data pre- or post-processing : Need to shuffle (for DIT), unshuffle (after DIF), zero-pad, truncate or divide (for IFFT) data separately as necessary - see auxiliary methods above
	//THESE ARE USED FOR COMPUTATIONS
	//
	
	void FFT_Radix4_DIT(ReIm *FX, int ldn, int N);				//VERIFIED, CORRECT - must use a shuffle method before
	void FFT_Radix4_DIF(ReIm *FX, int ldn, int N);				//VERIFIED, CORRECT - must use an unshuffle method after

	void IFFT_Radix4_DIT(ReIm *FX, int ldn, int N);				//VERIFIED, CORRECT - must use a shuffle method before and divide all numbers by N after
	void IFFT_Radix4_DIF(ReIm *FX, int ldn, int N);				//VERIFIED, CORRECT - must use an unshuffle method after and divide all numbers by N

	//
	//
	///////////////////////////////////

	//Advanced routines combining 1D FFTs included for convenience.

	//3D FFT of real input (N*M*K) into complex output (N/2 + 1)*M*K. Output is transposed.
	void FFT3D_R2C(double *Src, ReIm *Dst, ReIm *FFTSpace, int N, int M, int K, int powN, int powM, int powK);
	//As above but the real input is of dimension (N/2) * (M/2) * (K/2) and must be zero padded to dimensions N*M*K.
	void FFT3D_RZP2C(double *Src, ReIm *Dst, ReIm *FFTSpace, int N, int M, int K, int powN, int powM, int powK);

	//3D Inverse FFT from (N/2 + 1)*M*K to truncated real output of dimension (N/2) * (M/2) * (K/2)
	void IFFT3D_C2RT(ReIm *Src, double *Dst, ReIm *FFTSpace, int N, int M, int K, int powN, int powM, int powK);

	//2D FFT of real input (N*M) into complex output (N/2 + 1)*M. Output is transposed.
	void FFT2D_R2C(double *Src, ReIm *Dst, ReIm *FFTSpace, int N, int M, int powN, int powM);
	//As above but the real input is of dimension (N/2) * (M/2) and must be zero padded to dimensions N*M.
	void FFT2D_RZP2C(double *Src, ReIm *Dst, ReIm *FFTSpace, int N, int M, int powN, int powM);

	//3D Inverse FFT from (N/2 + 1)*M to truncated real output of dimension (N/2) * (M/2)
	void IFFT2D_C2RT(ReIm *Src, double *Dst, ReIm *FFTSpace, int N, int M, int powN, int powM);

	void Product(ReIm *Src1, ReIm *Src2, ReIm *Dst, int elements);

public: //public methods

	FFTMethods_Cpp(void);
	virtual ~FFTMethods_Cpp();
};