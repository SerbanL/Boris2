#pragma once

#include <omp.h>

#include "Types.h"

#define MAXFFT 32768								//Maximum length of FFT - used to initialize the cos and sin look-up tables

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//FFT Methods for demag calculation using convolution -> uses cpp routines

template <typename fftReal>
class FFTMethods_Cpp {

	using fftReal3 = VAL3<fftReal>;
	using fftComplex = __ReIm<fftReal>;
	using fftComplex3 = __ReIm3<fftReal>;

private:

	std::vector<fftComplex> cossin, sincos;

	//look-up tables for fft shuffling. First index, m, is the power of 2 of the fft length. On that row, the entries (2^m of them) give the shuffling indices.
	//fftShuffle gives indices for input data when output data is parsed linearly, fftShuffleInv is the inverse (input data parsed linearly, gives indices for output data)
	std::vector<int> fftShuffle, fftShuffleInv;

	int OmpThreads;

public:

	//------------------- CONSTRUCTOR : FFT_mng.h

	FFTMethods_Cpp(void);
	virtual ~FFTMethods_Cpp() {}

	//------------------- RADIX-2 : FFT_R2.h

	//NOTE : for all FFTs, N is the size which must be a power of 2, ldn is the power of 2 : N = 2^ldn
	//The template Type should be either a fftComplex or a fftComplex3 (ReIm or ReIm3)

	//Radix-2 FFT with decimation-in-time. Packed means no pre-shuffling needed as it's done by this function : just call to do the FFT
	template <typename Type>
	void FFT_Rad2DIT_Packed(Type *Fx, Type *FX, int ldn, int N);

	//------------------- SPLIT-RADIX : FFT_SR.h

	//Split-Radix FFT : lower multiplications count compared to Radix4 (~ 10%), same additions count, but much larger number of reads and writes.
	//My tests showed these are slower compared to Radix4 method, including on odd powers of 2.
	
	//Split-Radix FFT with decimation-in-time. Packed means no pre-shuffling needed as it's done by this function : just call to do the FFT
	template <typename Type>
	void FFT_SRDIT_Packed(Type *Fx, Type *FX, int ldn, int N);

	//Split-Radix IFFT with decimation-in-frequency. Packed means no un-shuffling needed as it's done by this function : just call to do the IFFT
	template <typename Type>
	void IFFT_SRDIF_Packed(Type *Fx, Type *FX, int ldn, int N);

	//------------------- RADIX-4 (use these) : FFT_R4.h

	//Radix-4 FFT and IFFT. Uses a Radix-2 step for odd powers of 2.
	//NOTE : these could be optimized further :
	//The special case j = 0, giving (exp, exp2 and exp3) factors of (1, 1, 1), I've already separated out, easy.
	//Other cases that could be optimized include multiplications by +/- 1, and +/- i. 
	//These can be separated out from the main loop by tracking special valus of N2 and j for which they occur and writing them out specifically (no extra if conditions used, just rewrite the loops)
	//This will result in a MUCH larger FFT function but should run faster (but probably not significantly faster).

	//Packed versions for testing

	//Radix-4 FFT with decimation-in-time. Packed means no pre-shuffling needed as it's done by this function : just call to do the FFT
	template <typename Type>
	void FFT_RAD4DIT_Packed(Type *Fx, Type *FX, int ldn, int N);

	//Radix-4 IFFT with decimation-in-frequency. Packed means no pre-shuffling needed as it's done by this function : just call to do the FFT
	template <typename Type>
	void IFFT_RAD4DIF_Packed(Type *Fx, Type *FX, int ldn, int N);

	//Non-packed versions: 
	//1. Before a DIT FFT/IFFT use a shuffle method to prepare input, e.g. CopyShuffle or CopyRealShuffle etc.
	//1. After a DIF FFT/IFFT use an unshuffle method to prepare output.
	//3. If running an DIT-FFT on real input (used a CopyRealShuffle version) then after the FFT call RealfromComplexFFT
	//4. After an IFFT must divide all numbers by N (the IFFT size)
	//5. If using a half-length IFFT with only real output, call RealfromComplexIFFT before the IFFT

	//FFTs

	template <typename Type>
	void FFT_Radix4_DIT(Type *FX, int ldn, int N);

	template <typename Type>
	void FFT_Radix4_DIF(Type *FX, int ldn, int N);

	//IFFTs

	template <typename Type>
	void IFFT_Radix4_DIT(Type *FX, int ldn, int N);
	
	template <typename Type>
	void IFFT_Radix4_DIF(Type *FX, int ldn, int N);

	//------------------- REAL INPUT FFT / REAL OUTPUT IFFT : FFT_shuf.h

	//use to obtain the fft of a real input after using the half-length complex fft method (i.e. after CopyRealShuffle to place the input in alternate real and imaginary parts, then doing a DIT fft of half-length on this)
	void RealfromComplexFFT(fftComplex *Fx, fftComplex *FX, int N);
	void RealfromComplexFFT(fftComplex3 *Fx, fftComplex3 *FX, int N);

	//use before running a half-length ifft which will have only real output. 
	//Here N is the length of the half-length IFFT, and we need N + 1 input points (don't need to store 2N points due to Hermitian symmetric property). After this method, run a DFT IFFT followed by UnShuffleRealTruncate
	void RealfromComplexIFFT(fftComplex *Fx, fftComplex *FX, int N);
	void RealfromComplexIFFT(fftComplex3 *Fx, fftComplex3 *FX, int N);

	//------------------- SHUFFLING (before DIT) : FFT_shuf.h

	//NOTE : use stride to transpose data when performing higher order FFTs - FFT input data must always be contiguous to minimise cache misses
	//(e.g. if stride is the row length then columns will be transposed to rows - in shuffled order - etc.)

	//use to bit-reverse inputs before calling a DIT fft;  
	void CopyShuffle(fftComplex *Fx, fftComplex *FX, int stride, int N);						
	void CopyShuffle(fftComplex3 *Fx, fftComplex3 *FX, int stride, int N);

	//use to bit-reverse inputs before calling a half-length complex DIT fft: if we have 2N real points then place them alternately in real and complex parts and do an N-point complex FFT
	void CopyRealShuffle(fftReal *Rex, fftComplex *FX, int N);								
	void CopyRealShuffle(fftReal3 *Rex, fftComplex3 *FX, int N);

	//as above but strided versions
	void CopyRealShuffle(fftReal *Rex, fftComplex *FX, int stride, int N);
	void CopyRealShuffle(fftReal3 *Rex, fftComplex3 *FX, int stride, int N);

	//use to bit-reverse inputs before calling a DIT fft. Upper half of inputs are taken as zero, so the input need only have N/2 points stored.
	void CopyShuffleZeroPad(fftComplex *Fx, fftComplex *FX, int stride, int N);					
	void CopyShuffleZeroPad(fftComplex3 *Fx, fftComplex3 *FX, int stride, int N);
	
	//use to bit-reverse inputs before calling a half-length complex DIT fft: if we have 2N real points then place them alternately in real and complex parts and do an N-point complex FFT. 
	//Also the upper half of inputs are taken as zero, so the input need only have N points stored.
	void CopyRealShuffleZeroPad(fftReal *Rex, fftComplex *FX, int N);						
	void CopyRealShuffleZeroPad(fftReal3 *Rex, fftComplex3 *FX, int n, int N);

	//------------------- UN-SHUFFLING (after DIF) : FFT_shuf.h
	
	//NOTE : use stride to transpose data when performing higher order FFTs
	//int n : only write first nth outputs, where n < N

	//use to bit-reverse outputs after calling a DIF fft.
	void UnShuffle(fftComplex *Fx, fftComplex *FX, int stride, int N);
	void UnShuffle(fftComplex3 *Fx, fftComplex3 *FX, int stride, int N);

	//use to bit-reverse outputs after calling a DIF fft. This method just writes the first half of outputs and ignores the second half.
	void UnShuffleTruncate(fftComplex *Fx, fftComplex *FX, int stride, int N);
	void UnShuffleTruncate(fftComplex3 *Fx, fftComplex3 *FX, int stride, int N);

	//use to bit-reverse outputs after calling a DIF fft.
	void UnShuffleRealTruncate(fftComplex *Fx, fftReal *ReX, int n, int N);						
	void UnShuffleRealTruncate(fftComplex3 *Fx, fftReal3 *ReX, int n, int N);

	//similar to UnShuffleRealTruncate, but in addition each value set in ReX is divided by divisor. The return value is the sum of of ReX multipled by ReIn cell-by-cell.
	fftReal UnShuffleRealTruncate_FinishConvolution_Set(fftComplex3 *Fx, fftReal3 *ReX, fftReal3 *ReIn, int n, int N, fftReal divisor);
	fftReal UnShuffleRealTruncate_FinishConvolution_Add(fftComplex3 *Fx, fftReal3 *ReX, fftReal3 *ReIn, int n, int N, fftReal divisor);
};