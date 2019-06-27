#include "stdafx.h"
#include "FFT.h"

FFTMethods_Cpp::FFTMethods_Cpp(void) {

	//Calculate cos and sin values.
	
	cossin.resize(MAXFFT);
	sincos.resize(MAXFFT);
	
	for(int idx = 0; idx < MAXFFT; idx++) {

		cossin[idx] = ReIm( cos(2*PI*idx / (double)MAXFFT), sin(2*PI*idx / (double)MAXFFT) );
		sincos[idx] = ReIm( sin(2*PI*idx / (double)MAXFFT), cos(2*PI*idx / (double)MAXFFT) );
	}
	
	//Calculate Radix-2 shuffling indices for all FFT lengths up to MAXFFT length - place these values contiguously in memory from shortest FFT length to highest
	//For an FFT of length N, the starting index for these values is N - 1.
	//For each fft of length N (which is a power of 2) we need to store N index values, thus we need 2 ^ (log2(MAXFFT) + 1) - 1 values, which is MAXFFT*2 - 1, total storage locations
	fftShuffle.resize(MAXFFT * 2 - 1);
	fftShuffleInv.resize(MAXFFT * 2 - 1);
	
	int m = 0;
	while(pow(2, m) < MAXFFT) m++;  //m = log2(MAXFFT)

	fftShuffle[0] = 0;
	fftShuffleInv[0] = 0;

	for(int idx = 1; idx <= m; idx++) {

		int fftLength = (1<<idx);

		//These are used to build shuffle inverse indices : fftShuffleCopy has a copy of the normal shuffle indices, fftShufInv has the indices in increasing order
		//After they are built, fftShufffleCopy is sorted in increasing order with ffftShufInv as a dependent array. The result is fftShufInv will contain the inverse shuffle indices .
		vector<int> fftShufInv, fftShufCopy;
		fftShufInv.resize(fftLength);
		fftShufCopy.resize(fftLength);

		fftShuffle[fftLength - 1] = 0;
		fftShuffle[(fftLength - 1) * 2] = fftLength - 1;

		fftShufInv[0] = 0;
		fftShufInv[fftLength - 1] = fftLength - 1;
		fftShufCopy[0] = 0;
		fftShufCopy[fftLength - 1] = fftLength - 1;

		for(int p = 1; p <= fftLength - 2; p++) {

			 fftShuffle[fftLength - 1 + p] = ShuffleIndex_Radix2(p, fftLength/2);
			 fftShufCopy[p] = fftShuffle[fftLength - 1 + p];
			 fftShufInv[p] = p;
		}

		//make inverse shuffle indices
		quicksort(fftShufCopy, fftShufInv);

		//and copy them in place
		for(int p = 0; p <= fftLength - 1; p++) {

			fftShuffleInv[fftLength - 1 + p] = fftShufInv[p];
		}
	}

	OmpThreads = omp_get_num_procs();
}

FFTMethods_Cpp::~FFTMethods_Cpp() {
	
}

int FFTMethods_Cpp::ShuffleIndex_Radix2(int n, int N) {

	int nt = 0, powkm1 = 1;

	if(!n) return 0;

	for(int k = 1; n; k++, N /= 2) {

		int t = n/N;

		if(!t) {
		
			powkm1 *= 2;
			continue;
		}

		nt += t*powkm1;
		n -= N;
		powkm1 *= 2;
	}

	return nt;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FFTMethods_Cpp::FFT_Rad2DIT_Packed(ReIm *Fx, ReIm *FX, int ldn, int N, int mult) {

	FX[0] = Fx[0];
	FX[(N-1)*mult] = Fx[(N-1)*mult];

	for(int p = 1; p <= N-2; p++) {

		FX[p*mult] = Fx[fftShuffle[N - 1 + p]*mult];
	}

	ReIm u, v, exp;
	
	//Simple multiplications by 1 + 0i are treated separately (starting part of outer loop)
	for(int r = 0; r <= N - 1; r += 2) {

		u = FX[r*mult]; 
		v = FX[(r+1)*mult];

		FX[r*mult] = u + v;

		FX[(r+1)*mult] = u - v;
	}
	
	int T = 2, Th;

	//Remaining part of outer loop
	for(int p = 2; p <= ldn; p++) {

		T *= 2;
		Th = T/2;

		for(int j = 0; j < Th; j++) {

			exp = ~cossin[(MAXFFT / T) * j];

			for(int r = 0; r <= N - T; r+=T) {

				u = FX[(r + j)*mult];

				v = FX[(r + j + Th)*mult] * exp;

				FX[(r + j)*mult] = u + v;

				FX[(r + j + Th)*mult] = u - v;
			}
		}
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FFTMethods_Cpp::FFT_RAD4DIT_Packed(ReIm *Fx, ReIm *FX, int ldn, int N, int mult) {

	int i0, i1, i2, i3, p;
	int N2, N4;

	ReIm t0, t1, t2, t3;
	ReIm exp, exp2, exp3;

	//Shuffle values, moving input data from input to output : input data is left unaffected.
	FX[0] = Fx[0];
	FX[(N-1)*mult] = Fx[(N-1)*mult];

	for(p = 1; p <= N-2; p++) {

		FX[p*mult] = Fx[fftShuffle[N - 1 + p]*mult];
	}

    p = (ldn&1);

	// n is not a power of 4, need a radix-2 step
    if (p!=0)  {

		for(i0 = 0; i0 < N; i0 += 2) {

			t0 = FX[i0*mult] - FX[(i0 + 1)*mult];
			FX[i0*mult] = FX[i0*mult] + FX[(i0 + 1)*mult];
			FX[(i0 + 1)*mult] = t0;
		}
    }
	
	N2 = (1 << p);

    for (p = p + 2; p <= ldn; p+= 2) {
        
		N4 = N2;
		N2 *= 4;
		
		//Special case j = 0 separated from loop : no multiplications.
		for (int r = 0; r < N; r += N2) {

            i1 = r + N4;
            i2 = i1 + N4;
            i3 = i2 + N4;
               
			t0 = (FX[r*mult] + FX[i1*mult]);
			t1 = (FX[r*mult] - FX[i1*mult]);
			t2 = (FX[i2*mult] + FX[i3*mult]);
			t3 =!(FX[i3*mult] - FX[i2*mult]);

			FX[r*mult] = t0 + t2;
			FX[i1*mult] = t1 + t3;
			FX[i2*mult] = t0 - t2;
			FX[i3*mult] = t1 - t3;
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			exp  = ~cossin[(MAXFFT / N2) * j];
			exp2 = ~cossin[(MAXFFT / N2) * 2*j];
			exp3 = ~cossin[(MAXFFT / N2) * 3*j];

            for (int r = 0; r < N; r += N2) {

                i0 = j + r;
                i1 = i0 + N4;
                i2 = i1 + N4;
                i3 = i2 + N4;
               
				t0 = (FX[i0*mult]       + FX[i1*mult] * exp2);
				t1 = (FX[i0*mult]       - FX[i1*mult] * exp2);
				t2 = (FX[i2*mult] * exp + FX[i3*mult] * exp3);
				t3 =!(FX[i3*mult] * exp3 - FX[i2*mult] * exp);

				FX[i0*mult] = t0 + t2;
				FX[i1*mult] = t1 + t3;
				FX[i2*mult] = t0 - t2;
				FX[i3*mult] = t1 - t3;
            }
        }
    }
}

void FFTMethods_Cpp::IFFT_RAD4DIF_Packed(ReIm *Fx, ReIm *FX, int ldn, int N, int mult) {

	int i0, i1, i2, i3, p;
	int N2, N4;

	ReIm t0, t1, t2, t3;
	ReIm exp, exp2, exp3;

	N2 = 4*N;

    for (p = ldn; p >= 2; p-= 2) {
        
		N2 /= 4;
		N4 = N2 / 4;
		
		//Special case j = 0 separated from loop : no multiplications.
		for (int r = 0; r < N; r += N2) {

            i1 = r + N4;
            i2 = i1 + N4;
            i3 = i2 + N4;
               
			t0 = (Fx[r*mult] + Fx[i2*mult]);
			t1 = (Fx[r*mult] - Fx[i2*mult]);
			t2 = (Fx[i1*mult] + Fx[i3*mult]);
			t3 =!(Fx[i3*mult] - Fx[i1*mult]);

			//Note change in sign for t3 for IFFT compared to FFT
			Fx[r*mult] = (t0 + t2);
			Fx[i1*mult] = (t0 - t2);
			Fx[i2*mult] = (t1 - t3);
			Fx[i3*mult] = (t1 + t3);
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			//Conjugated roots since IFFT
			exp  = cossin[(MAXFFT / N2) * j];
			exp2 = cossin[(MAXFFT / N2) * 2*j];
			exp3 = cossin[(MAXFFT / N2) * 3*j];

            for (int r = 0; r < N; r += N2) {

                i0 = j + r;
                i1 = i0 + N4;
                i2 = i1 + N4;
                i3 = i2 + N4;
               
				t0 = (Fx[i0*mult] + Fx[i2*mult]);
				t1 = (Fx[i0*mult] - Fx[i2*mult]);
				t2 = (Fx[i1*mult] + Fx[i3*mult]);
				t3 =!(Fx[i3*mult] - Fx[i1*mult]);

				//Note change in sign for t3 for IFFT compared to FFT
				Fx[i0*mult] = (t0 + t2);
				Fx[i1*mult] = (t0 - t2) * exp2;
				Fx[i2*mult] = (t1 - t3) * exp;
				Fx[i3*mult] = (t1 + t3) * exp3;
            }
        }
    }

	p = (ldn&1);

	// n is not a power of 4, need a radix-2 step
    if (p!=0)  {

		for(i0 = 0; i0 < N; i0 += 2) {

			t0 = Fx[i0*mult] - Fx[(i0 + 1)*mult];
			Fx[i0*mult] = Fx[i0*mult] + Fx[(i0 + 1)*mult];
			Fx[(i0 + 1)*mult] = t0;
		}
    }

	//Shuffle values, moving input data from input to output. Also divide by N since this is IFFT
	FX[0] = Fx[0] / N;
	FX[(N-1)*mult] = Fx[(N-1)*mult] / N;

	for(p = 1; p <= N-2; p++) {

		FX[p*mult] = Fx[fftShuffle[N - 1 + p]*mult] / N;
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FFTMethods_Cpp::FFT_SRDIT_Packed(ReIm *Fx, ReIm *FX, int ldn, int N, int mult) {

	int N2, N4;

	int ix, id, i0, i1, i2, i3;

	ReIm exp1, exp3;
	ReIm t0, t1;

	//Shuffle values, moving input data from input to output : input data is left unaffected.
	FX[0] = Fx[0];
	FX[(N-1)*mult] = Fx[(N-1)*mult];

	for(int p = 1; p <= N-2; p++) {

		FX[p*mult] = Fx[fftShuffle[N - 1 + p]*mult];
	}

	//Stage 1
	for (ix = 0, id = 4; ix < N;  id *= 4) {
        
		for (i0 = ix; i0 < N; i0 += id)  {
			
			t0 = FX[i0*mult] - FX[(i0 + 1)*mult];
			FX[i0*mult] = FX[i0*mult] + FX[(i0 + 1)*mult];
			FX[(i0 + 1)*mult] = t0;
		}

        ix = 2 * id - 2;
    }
	
	//Stage 2
	N2 = 2;

	for(int k = ldn - 2; k >= 0; k--) {

		N2 = N2 * 2;
		N4 = N2 / 4;
		
		//Treat case j = 0 separately due to simpler multiplications
		for(ix = 0, id = N2 * 2; ix < N; id *= 4) {

			for (i0 = ix; i0 < N; i0 += id) {

				i1 = i0 + N4;
				i2 = i1 + N4;
				i3 = i2 + N4;

				t0 = FX[i2*mult] + FX[i3*mult];
				t1 = !(FX[i3*mult] - FX[i2*mult]);  // ! is multiplication by i

				FX[i2*mult] = FX[i0*mult] - t0;
				FX[i0*mult] = FX[i0*mult] + t0;
				
				FX[i3*mult] = FX[i1*mult] - t1;
				FX[i1*mult] = FX[i1*mult] + t1;
			}

			ix = 2 * id - N2;
		}

		//Remaining cases from j = 1 upwards
		for(int j = 1; j < N4; j++) {

			exp1 = ~cossin[(MAXFFT / N2) * j];
			exp3 = ~cossin[(MAXFFT / N2) * 3*j];

			for(ix = j, id = N2 * 2; ix < N; id *= 4) {

				for (i0 = ix; i0 < N; i0 += id) {

					i1 = i0 + N4;
					i2 = i1 + N4;
					i3 = i2 + N4;

					t0 = FX[i3*mult] * exp3 + FX[i2*mult] * exp1;
					t1 = !(FX[i3*mult] * exp3 - FX[i2*mult] * exp1);  // ! is multiplication by i

					FX[i2*mult] = FX[i0*mult] - t0;
					FX[i0*mult] = FX[i0*mult] + t0;
				
					FX[i3*mult] = FX[i1*mult] - t1;
					FX[i1*mult] = FX[i1*mult] + t1;
				}

				ix = 2 * id - N2 + j;
			}
		}
	}	
}

void FFTMethods_Cpp::IFFT_SRDIF_Packed(ReIm *Fx, ReIm *FX, int ldn, int N, int mult) {

	int N2 = 2*N, N4;

	int ix, id, i0, i1, i2, i3;

	ReIm exp1, exp3, t0, t1;

	//Stage 1 (this is similar to stage 2 for DIT)
	for(int k = 0; k < ldn - 1; k++) {

		N2 = N2 / 2;
		N4 = N2 / 4;
		
		//Treat case j = 0 separately due to simpler multiplications
		for(ix = 0, id = N2 * 2; ix < N; id *= 4) {

			for (i0 = ix; i0 < N; i0 += id) {

				i1 = i0 + N4;
				i2 = i1 + N4;
				i3 = i2 + N4;

				t0 = Fx[i0*mult] - Fx[i2*mult];
				t1 = Fx[i1*mult] - Fx[i3*mult];
				
				Fx[i0*mult] = Fx[i0*mult] + Fx[i2*mult];
				Fx[i1*mult] = Fx[i1*mult] + Fx[i3*mult];

				Fx[i2*mult] = t0 + !t1;
				Fx[i3*mult] = t0 - !t1;
			}

			ix = 2 * id - N2;
		}

		//Remaining cases from j = 1 upwards
		for(int j = 1; j < N4; j++) {

			//Note these are not conjugated (as you would think they should be since this is IFFT) since we are doing DIF, not DIT - butterflies require it
			exp1 = cossin[(MAXFFT / N2) * j];
			exp3 = cossin[(MAXFFT / N2) * 3*j];

			for(ix = j, id = N2 * 2; ix < N; id *= 4) {

				for (i0 = ix; i0 < N; i0 += id) {

					i1 = i0 + N4;
					i2 = i1 + N4;
					i3 = i2 + N4;
					
					t0 = Fx[i0*mult] - Fx[i2*mult];
					t1 = Fx[i1*mult] - Fx[i3*mult];

					Fx[i0*mult] = Fx[i0*mult] + Fx[i2*mult];
					Fx[i1*mult] = Fx[i1*mult] + Fx[i3*mult];

					Fx[i2*mult] = (t0 + !t1) * exp1;
					Fx[i3*mult] = (t0 - !t1) * exp3;
				}

				ix = 2 * id - N2 + j;
			}
		}
	}

	//Stage 2 (same as stage 1 for DIT)
	for (ix = 0, id = 4; ix < N;  id *= 4) {
        
		for (i0 = ix; i0 < N; i0 += id)  {
			
			t0 = Fx[i0*mult] - Fx[(i0 + 1)*mult];

			Fx[i0*mult] = Fx[i0*mult] + Fx[(i0 + 1)*mult];
			Fx[(i0 + 1)*mult] = t0;
		}

        ix = 2 * id - 2;
    }
		
	//Unshuffle, moving from input data to output, and divide by N for IFFT
	for(int p = 1; p < N; p++) {

		FX[p*mult] = Fx[fftShuffle[N - 1 + p]*mult] / N;
	}

	FX[0] = Fx[0] / N;
	FX[(N-1)*mult] = Fx[(N-1)*mult] / N;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//In-place FFT, Radix 4 Decimation in Time. Input data must be shuffled.
void FFTMethods_Cpp::FFT_Radix4_DIT(ReIm *FX, int ldn, int N) {

	int i0, i1, i2, i3, p;
	int N2, N4;

	ReIm t0, t1, t2, t3;
	ReIm exp, exp2, exp3;

    p = (ldn&1);

	// n is not a power of 4, need a radix-2 step
    if (p!=0)  {

		for(i0 = 0; i0 < N; i0 += 2) {

			t0 = FX[i0] - FX[(i0 + 1)];
			FX[i0] = FX[i0] + FX[(i0 + 1)];
			FX[(i0 + 1)] = t0;
		}
    }
	
	N2 = (1 << p);

    for (p = p + 2; p <= ldn; p+= 2) {
        
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
			t3 =!(FX[i3] - FX[i2]);

			FX[r]  = t0 + t2;
			FX[i1] = t1 + t3;
			FX[i2] = t0 - t2;
			FX[i3] = t1 - t3;
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			exp  = ~cossin[(MAXFFT / N2) * j];
			exp2 = ~cossin[(MAXFFT / N2) * 2*j];
			exp3 = ~cossin[(MAXFFT / N2) * 3*j];

            for (int r = 0; r < N; r += N2) {

                i0 = j + r;
                i1 = i0 + N4;
                i2 = i1 + N4;
                i3 = i2 + N4;
               
				t0 = (FX[i0]       + FX[i1] * exp2);
				t1 = (FX[i0]       - FX[i1] * exp2);
				t2 = (FX[i2] * exp + FX[i3] * exp3);
				t3 =!(FX[i3] * exp3 - FX[i2] * exp);

				FX[i0] = t0 + t2;
				FX[i1] = t1 + t3;
				FX[i2] = t0 - t2;
				FX[i3] = t1 - t3;
            }
        }
    }
}

void FFTMethods_Cpp::FFT_Radix4_DIT(ReIm3 *FX, int ldn, int N) {
	
	int i0, i1, i2, i3, p;
	int N2, N4;

	ReIm3 t0, t1, t2, t3;
	ReIm exp, exp2, exp3;

    p = (ldn&1);

	// n is not a power of 4, need a radix-2 step
    if (p!=0)  {

		for(i0 = 0; i0 < N; i0 += 2) {

			t0 = FX[i0] - FX[(i0 + 1)];
			FX[i0] = FX[i0] + FX[(i0 + 1)];
			FX[(i0 + 1)] = t0;
		}
    }
	
	N2 = (1 << p);

    for (p = p + 2; p <= ldn; p+= 2) {
        
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
			t3 =!(FX[i3] - FX[i2]);

			FX[r]  = t0 + t2;
			FX[i1] = t1 + t3;
			FX[i2] = t0 - t2;
			FX[i3] = t1 - t3;
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			exp  = ~cossin[(MAXFFT / N2) * j];
			exp2 = ~cossin[(MAXFFT / N2) * 2*j];
			exp3 = ~cossin[(MAXFFT / N2) * 3*j];

            for (int r = 0; r < N; r += N2) {

                i0 = j + r;
                i1 = i0 + N4;
                i2 = i1 + N4;
                i3 = i2 + N4;
               
				t0 = (FX[i0]          + (FX[i1] * exp2));
				t1 = (FX[i0]          - (FX[i1] * exp2));
				t2 = ((FX[i2] * exp)  + (FX[i3] * exp3));
				t3 =!((FX[i3] * exp3) - (FX[i2] * exp ));

				FX[i0] = t0 + t2;
				FX[i1] = t1 + t3;
				FX[i2] = t0 - t2;
				FX[i3] = t1 - t3;
            }
        }
    }
}

//In-place FFT, Radix 4 Decimation in Frequency. Output data is shuffled.
void FFTMethods_Cpp::FFT_Radix4_DIF(ReIm *FX, int ldn, int N) {

	int i0, i1, i2, i3, p;
	int N2, N4;

	ReIm t0, t1, t2, t3;
	ReIm exp, exp2, exp3;

	N2 = 4*N;

    for (p = ldn; p >= 2; p-= 2) {
        
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
			t3 =!(FX[i3] - FX[i1]);

			FX[r] = (t0 + t2);
			FX[i1] = (t0 - t2);
			FX[i2] = (t1 + t3);
			FX[i3] = (t1 - t3);
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			exp  = ~cossin[(MAXFFT / N2) * j];
			exp2 = ~cossin[(MAXFFT / N2) * 2*j];
			exp3 = ~cossin[(MAXFFT / N2) * 3*j];

            for (int r = 0; r < N; r += N2) {

                i0 = j + r;
                i1 = i0 + N4;
                i2 = i1 + N4;
                i3 = i2 + N4;
               
				t0 = (FX[i0] + FX[i2]);
				t1 = (FX[i0] - FX[i2]);
				t2 = (FX[i1] + FX[i3]);
				t3 =!(FX[i3] - FX[i1]);

				FX[i0] = (t0 + t2);
				FX[i1] = (t0 - t2) * exp2;
				FX[i2] = (t1 + t3) * exp;
				FX[i3] = (t1 - t3) * exp3;
            }
        }
    }

	p = (ldn&1);

	// n is not a power of 4, need a radix-2 step
    if (p!=0)  {

		for(i0 = 0; i0 < N; i0 += 2) {

			t0 = FX[i0] - FX[(i0 + 1)];
			FX[i0] = FX[i0] + FX[(i0 + 1)];
			FX[(i0 + 1)] = t0;
		}
    }
}

//In-place IFFT, Radix 4 Decimation in Time. Input data must be shuffled, output data not divided by N.
void FFTMethods_Cpp::IFFT_Radix4_DIT(ReIm *FX, int ldn, int N) {

	int i0, i1, i2, i3, p;
	int N2, N4;

	ReIm t0, t1, t2, t3;
	ReIm exp, exp2, exp3;

    p = (ldn&1);

	// n is not a power of 4, need a radix-2 step
    if (p!=0)  {

		for(i0 = 0; i0 < N; i0 += 2) {

			t0 = FX[i0] - FX[(i0 + 1)];
			FX[i0] = FX[i0] + FX[(i0 + 1)];
			FX[(i0 + 1)] = t0;
		}
    }
	
	N2 = (1 << p);

    for (p = p + 2; p <= ldn; p+= 2) {
        
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
			t3 =!(FX[i3] - FX[i2]);

			FX[r] = t0 + t2;
			FX[i1] = t1 - t3;
			FX[i2] = t0 - t2;
			FX[i3] = t1 + t3;
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			exp  = cossin[(MAXFFT / N2) * j];
			exp2 = cossin[(MAXFFT / N2) * 2*j];
			exp3 = cossin[(MAXFFT / N2) * 3*j];

            for (int r = 0; r < N; r += N2) {

                i0 = j + r;
                i1 = i0 + N4;
                i2 = i1 + N4;
                i3 = i2 + N4;
               
				t0 = (FX[i0]       + FX[i1] * exp2);
				t1 = (FX[i0]       - FX[i1] * exp2);
				t2 = (FX[i2] * exp + FX[i3] * exp3);
				t3 =!(FX[i3] * exp3 - FX[i2] * exp);

				FX[i0] = t0 + t2;
				FX[i1] = t1 - t3;
				FX[i2] = t0 - t2;
				FX[i3] = t1 + t3;
            }
        }
    }
}

//In-place IFFT, Radix 4 Decimation in Frequency. Output data is shuffled and not divided by N.
void FFTMethods_Cpp::IFFT_Radix4_DIF(ReIm *FX, int ldn, int N) {

	int i0, i1, i2, i3, p;
	int N2, N4;

	ReIm t0, t1, t2, t3;
	ReIm exp, exp2, exp3;

	N2 = 4*N;

    for (p = ldn; p >= 2; p-= 2) {
        
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
			t3 =!(FX[i3] - FX[i1]);

			FX[r] = (t0 + t2);
			FX[i1] = (t0 - t2);
			FX[i2] = (t1 - t3);
			FX[i3] = (t1 + t3);
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			exp  = cossin[(MAXFFT / N2) * j];
			exp2 = cossin[(MAXFFT / N2) * 2*j];
			exp3 = cossin[(MAXFFT / N2) * 3*j];

            for (int r = 0; r < N; r += N2) {

                i0 = j + r;
                i1 = i0 + N4;
                i2 = i1 + N4;
                i3 = i2 + N4;
               
				t0 = (FX[i0] + FX[i2]);
				t1 = (FX[i0] - FX[i2]);
				t2 = (FX[i1] + FX[i3]);
				t3 =!(FX[i3] - FX[i1]);

				FX[i0] = (t0 + t2);
				FX[i1] = (t0 - t2) * exp2;
				FX[i2] = (t1 - t3) * exp;
				FX[i3] = (t1 + t3) * exp3;
            }
        }
    }

	p = (ldn&1);

	// n is not a power of 4, need a radix-2 step
    if (p!=0)  {

		for(i0 = 0; i0 < N; i0 += 2) {

			t0 = FX[i0] - FX[(i0 + 1)];
			FX[i0] = FX[i0] + FX[(i0 + 1)];
			FX[(i0 + 1)] = t0;
		}
    }
}

void FFTMethods_Cpp::IFFT_Radix4_DIF(ReIm3 *FX, int ldn, int N) {

	int i0, i1, i2, i3, p;
	int N2, N4;

	ReIm3 t0, t1, t2, t3;
	ReIm exp, exp2, exp3;

	N2 = 4*N;

    for (p = ldn; p >= 2; p-= 2) {
        
		N2 /= 4;
		N4 = N2 / 4;
		
		//Special case j = 0 separated from loop : no multiplications.
		for (int r = 0; r < N; r += N2) {

            i1 = r + N4;
            i2 = i1 + N4;
            i3 = i2 + N4;
               
			t0 = (FX[r]  + FX[i2]);
			t1 = (FX[r]  - FX[i2]);
			t2 = (FX[i1] + FX[i3]);
			t3 =!(FX[i3] - FX[i1]);

			FX[r]  = (t0 + t2);
			FX[i1] = (t0 - t2);
			FX[i2] = (t1 - t3);
			FX[i3] = (t1 + t3);
		}

		//Remaining cases from j = 1 upwards
		for (int j = 1; j < N4; j++) {

			exp  = cossin[(MAXFFT / N2) * j];
			exp2 = cossin[(MAXFFT / N2) * 2*j];
			exp3 = cossin[(MAXFFT / N2) * 3*j];

            for (int r = 0; r < N; r += N2) {

                i0 = j + r;
                i1 = i0 + N4;
                i2 = i1 + N4;
                i3 = i2 + N4;
               
				t0 = (FX[i0] + FX[i2]);
				t1 = (FX[i0] - FX[i2]);
				t2 = (FX[i1] + FX[i3]);
				t3 =!(FX[i3] - FX[i1]);

				FX[i0] = (t0 + t2);
				FX[i1] = (t0 - t2) * exp2;
				FX[i2] = (t1 - t3) * exp;
				FX[i3] = (t1 + t3) * exp3;
            }
        }
    }

	p = (ldn&1);

	// n is not a power of 4, need a radix-2 step
    if (p!=0)  {

		for(i0 = 0; i0 < N; i0 += 2) {

			t0 = FX[i0] - FX[(i0 + 1)];
			FX[i0] = FX[i0] + FX[(i0 + 1)];
			FX[(i0 + 1)] = t0;
		}
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Auxiliary Methods

void FFTMethods_Cpp::CopyShuffle(double *Rex, double *Imx, ReIm *FX, int stride, int N) {
	
	//Shuffle values moving input data to output
	
	FX[0] = ReIm(Rex[0], Imx[0]);
	FX[N - 1] = ReIm(Rex[(N - 1) * stride], Imx[(N - 1) * stride]);

	for(int p = 1; p <= N-2; p++) {

		FX[p] = ReIm(Rex[fftShuffle[N - 1 + p] * stride], Imx[fftShuffle[N - 1 + p] * stride]);
	}
}

void FFTMethods_Cpp::CopyShuffle(ReIm *Fx, ReIm *FX, int stride, int N) {
	
	//Shuffle values moving input data to output
	
	FX[0] = Fx[0];
	FX[N - 1] = Fx[(N - 1) * stride];

	for(int p = 1; p <= N-2; p++) {

		FX[p] = Fx[fftShuffle[N - 1 + p] * stride];
	}
}

void FFTMethods_Cpp::CopyShuffle(ReIm3 *Fx, ReIm3 *FX, int stride, int N) {
	
	//Shuffle values moving input data to output
	
	FX[0] = Fx[0];
	FX[N - 1] = Fx[(N - 1) * stride];

	for(int p = 1; p <= N-2; p++) {

		FX[p] = Fx[fftShuffle[N - 1 + p] * stride];
	}
}

void FFTMethods_Cpp::CopyRealShuffle(double *Rex, ReIm *FX, int N) {

	//Shuffle values moving input data to output
	//This routine is for real input to prepare half-length complex FFT : even indexed input values are set in Real part, odd indexed input values set in Complex part
	//we have 2N real points so N is the length of the complex fft

	//INVERSE METHOD
	for(int p = 0; p < N; p++) {

		int index = fftShuffleInv[N - 1 + p];
		
		FX[index] = ReIm(Rex[2*p], Rex[(2*p + 1)]);
	}
}

void FFTMethods_Cpp::CopyRealShuffle(DBL3 *Rex, ReIm3 *FX, int N) {

	//Shuffle values moving input data to output
	//This routine is for real input to prepare half-length complex FFT : even indexed input values are set in Real part, odd indexed input values set in Complex part
	//we have 2N real points so N is the length of the complex fft

	//INVERSE METHOD
	for(int p = 0; p < N; p++) {

		int index = fftShuffleInv[N - 1 + p];
		
		FX[index] = ReIm3(Rex[2*p], Rex[(2*p + 1)]);
	}
}

void FFTMethods_Cpp::CopyShuffleZeroPad(ReIm *Fx, ReIm *FX, int stride, int N) {

	//Shuffle values moving input data to output
	//Zero padding on input data from N/2 upwards - results in odd points being zero in shuffled data
	
	/*
	//DIRECT METHOD
	for(int p = 0; p <= N-2; p += 2) {

		FX[p] = Fx[fftShuffle[N - 1 + p] * stride];
		FX[p + 1] = ReIm(0, 0);
	}
	*/
	
	//INVERSE METHOD
	for(int p = 0; p < N/2; p++) {

		int index = fftShuffleInv[N - 1 + p];
		
		FX[index] = Fx[p * stride];
		
		//odd indexed shuffled inputs must be zero due to zero padding (note first half of non-shuffled inputs are always distributed to even-indexed shuffled values, thus index + 1 is odd)
		FX[index + 1] = ReIm(0, 0);
	}
}

void FFTMethods_Cpp::CopyShuffleZeroPad(ReIm3 *Fx, ReIm3 *FX, int stride, int N) {

	//INVERSE METHOD
	for(int p = 0; p < N/2; p++) {

		int index = fftShuffleInv[N - 1 + p];
		
		FX[index] = Fx[p * stride];
		
		//odd indexed shuffled inputs must be zero due to zero padding (note first half of non-shuffled inputs are always distributed to even-indexed shuffled values, thus index + 1 is odd)
		FX[index + 1] = ReIm3();
	}
}

void FFTMethods_Cpp::CopyShuffleZeroPad(double *Rex, double *Imx, ReIm *FX, int N) {

	//Shuffle values moving input data to output
	//Zero padding on input data from N/2 upwards - results in odd points being zero in shuffled data

	/*
	//DIRECT METHOD
	for(int p = 0; p <= N-2; p += 2) {

		FX[p] = ReIm(Rex[fftShuffle[N - 1 + p] * stride], Imx[fftShuffle[N - 1 + p] * stride]);
		FX[p + 1] = ReIm(0, 0);
	}
	*/
	
	//INVERSE METHOD
	for(int p = 0; p < N/2; p++) {

		int index = fftShuffleInv[N - 1 + p];
		
		FX[index] = ReIm(Rex[p], Imx[p]);
		
		//odd indexed shuffled inputs must be zero due to zero padding (note first half of non-shuffled inputs are always distributed to even-indexed shuffled values, thus index + 1 is odd)
		FX[index + 1] = ReIm(0, 0);
	}
}

void FFTMethods_Cpp::CopyRealShuffleZeroPad(double *Rex, ReIm *FX, int N) {

	//Shuffle values moving input data to output
	//Zero padding on input data from N/2 upwards - results in odd points being zero in shuffled data
	//This routine is for real input to prepare half-length complex FFT : even indexed input values are set in Real part, odd indexed input values set in Complex part
	//we have 2N real points (but the upper N are taken to be zero) so N is the length of the complex fft

	//INVERSE METHOD
	for(int p = 0; p < N/2; p++) {

		int index = fftShuffleInv[N - 1 + p];
		
		FX[index] = ReIm(Rex[2*p], Rex[(2*p + 1)]);
		
		//odd indexed shuffled inputs must be zero due to zero padding (note first half of non-shuffled inputs are always distributed to even-indexed shuffled values, thus index + 1 is odd)
		FX[index + 1] = ReIm(0, 0);
	}
}

void FFTMethods_Cpp::CopyRealShuffleZeroPad(DBL3 *Rex, ReIm3 *FX, int n, int N) {

	//Shuffle values moving input data to output
	//Zero padding on input data from N/2 upwards - results in odd points being zero in shuffled data
	//This routine is for real input to prepare half-length complex FFT : even indexed input values are set in Real part, odd indexed input values set in Complex part
	//we have 2N real points (but the upper N are taken to be zero) so N is the length of the complex fft

	//n is the dimension of Rex along x. N is a number greater or equal to n, which is a power of 2.
	//n could be smaller than N, in which case fill the remaining space with zeroes (i.e. as if Rex has dimension N/2 along x, with zero values from n to N/2 - 1.

	//INVERSE METHOD
	for(int p = 0; p < N/2; p++) {

		int index = fftShuffleInv[N - 1 + p];
		
		DBL3 val1, val2; //these are initialized to zero by default in their constructor

		if(2*p < n) val1 = Rex[2*p];
		if(2*p + 1 < n) val2 = Rex[2*p + 1];

		FX[index] = ReIm3(val1, val2);
		
		//odd indexed shuffled inputs must be zero due to zero padding (note first half of non-shuffled inputs are always distributed to even-indexed shuffled values, thus index + 1 is odd)
		FX[index + 1] = ReIm3();
	}
}

void FFTMethods_Cpp::CopyRealProdShuffleZeroPad(double *ReS, double *Rex, ReIm *FX, int N) {

	//Shuffle values moving input data to output
	//Zero padding on input data from N/2 upwards - results in odd points being zero in shuffled data
	//This routine is for real input to prepare half-length complex FFT : even indexed input values are set in Real part, odd indexed input values set in Complex part
	//we have 2N real points (but the upper N are taken to be zero) so N is the length of the complex fft

	//INVERSE METHOD
	for(int p = 0; p < N/2; p++) {

		int index = fftShuffleInv[N - 1 + p];
		
		FX[index] = ReIm(Rex[2*p] * ReS[2*p], Rex[2*p + 1] * ReS[2*p + 1]);
		
		//odd indexed shuffled inputs must be zero due to zero padding (note first half of non-shuffled inputs are always distributed to even-indexed shuffled values, thus index + 1 is odd)
		FX[index + 1] = ReIm(0, 0);
	}
}

void FFTMethods_Cpp::RealfromComplexFFT(ReIm *Fx, ReIm *FX, int N) {

	//Use this after doing a half-length complex FFT on a real input (should have used CopyRealShuffle or CopyRealShuffleZeroPad before doing the FFT as well).
	//FFT result is in Fx (N points, N is the half-length complex FFT length (real input would have had 2N points)).
	//This method will write N + 1 values in FX. (Herimitian symmetric property will give the remaining N - 1 points. Also note n = 0 and n = N points have zero imaginary parts.)

	ReIm Yn, YNmn;

	//n = 0
	FX[0] = ReIm(Fx[0].Re + Fx[0].Im, 0);

	for(int n = 1; n < N; n++) {

		Yn   = Fx[n];
		YNmn = Fx[N - n];

		FX[n] = (Yn + ~YNmn - sincos[(MAXFFT/(2*N)) * n] * (Yn - ~YNmn)) * 0.5;
	}

	//n = N
	FX[N] = ReIm(Fx[0].Re - Fx[0].Im, 0);
}

void FFTMethods_Cpp::RealfromComplexFFT(ReIm3 *Fx, ReIm3 *FX, int N) {

	ReIm3 Yn, YNmn;

	//n = 0
	FX[0] = ReIm3( ReIm(Fx[0].x.Re + Fx[0].x.Im, 0), ReIm(Fx[0].y.Re + Fx[0].y.Im, 0), ReIm(Fx[0].z.Re + Fx[0].z.Im, 0) );

	for(int n = 1; n < N; n++) {

		Yn   = Fx[n];
		YNmn = Fx[N - n];

		FX[n] = (Yn + ~YNmn - ((Yn - ~YNmn) * sincos[(MAXFFT/(2*N)) * n])) * 0.5;
	}

	//n = N
	FX[N] = ReIm3( ReIm(Fx[0].x.Re - Fx[0].x.Im, 0), ReIm(Fx[0].y.Re - Fx[0].y.Im, 0), ReIm(Fx[0].z.Re - Fx[0].z.Im, 0) );
}

void FFTMethods_Cpp::RealfromComplexIFFT(ReIm *Fx, ReIm *FX, int N) {

	//Use this before doing a half-length complex IFFT on complex input to get real output (should also use UnShuffleRealTruncate or UnShuffleReal after the IFFT).
	//N is the half-length complex IFFT length (real output would have had 2N points)
	//This method will write N values in FX using N + 1 values in Fx.

	for(int n = 0; n < N; n++) {

		ReIm Xn   = Fx[n];
		ReIm XNmn = Fx[N - n];

		ReIm Yn = (Xn + ~XNmn - ~sincos[(MAXFFT/(2*N)) * n] * (Xn - ~XNmn)) * 0.5;

		FX[n] = Yn;
	}
}

void FFTMethods_Cpp::RealfromComplexIFFT(ReIm3 *Fx, ReIm3 *FX, int N) {
	
	//Use this before doing a half-length complex IFFT on complex input to get real output (should also use UnShuffleRealTruncate or UnShuffleReal after the IFFT).
	//N is the half-length complex IFFT length (real output would have had 2N points)
	//This method will write N values in FX using N + 1 values in Fx.

	for(int n = 0; n < N; n++) {

		ReIm3 Xn   = Fx[n];
		ReIm3 XNmn = Fx[N - n];

		ReIm3 Yn = (Xn + ~XNmn - ((Xn - ~XNmn) * ~sincos[(MAXFFT/(2*N)) * n])) * 0.5;

		FX[n] = Yn;
	}
}

void FFTMethods_Cpp::UnShuffleTruncate(ReIm *Fx, ReIm *FX, int stride, int N) {

	//Unshuffle values, moving input data to output placed using stride. Only half the outputs are written.

	/*
	//DIRECT METHOD
	for(int p = 0; p < N/2; p++) {

		FX[p * stride] = Fx[fftShuffle[N - 1 + p]];
	}
	*/
	
	//INVERSE METHOD
	//Even-indexed inputs are always distributed to lower half of output de-shuffled (or shuffled) values : we only want the lower half.
	for(int p = 0; p < N; p+=2) {

		int index = fftShuffleInv[N - 1 + p];
		
		FX[index * stride] = Fx[p];
	}
}

void FFTMethods_Cpp::UnShuffleTruncate(ReIm3 *Fx, ReIm3 *FX, int stride, int N) {
	
	//INVERSE METHOD
	//Even-indexed inputs are always distributed to lower half of output de-shuffled (or shuffled) values : we only want the lower half.
	for(int p = 0; p < N; p+=2) {

		int index = fftShuffleInv[N - 1 + p];
		
		FX[index * stride] = Fx[p];
	}
}

void FFTMethods_Cpp::UnShuffleTruncate(ReIm *Fx, double *ReX, double *ImX, int N) {

	//Unshuffle values, moving input data to output placed using stride. Only half the outputs are written.

	/*
	//DIRECT METHOD
	for(int p = 0; p < N/2; p++) {

		ReX[p * stride] = Fx[fftShuffle[N - 1 + p]].Re;
		ImX[p * stride] = Fx[fftShuffle[N - 1 + p]].Im;
	}
	*/
	
	//INVERSE METHOD
	//Even-indexed inputs are always distributed to lower half of output de-shuffled (or shuffled) values : we only want the lower half.
	for(int p = 0; p < N; p+=2) {

		int index = fftShuffleInv[N - 1 + p];
		
		ReX[index] = Fx[p].Re;
		ImX[index] = Fx[p].Im;
	}
}

void FFTMethods_Cpp::UnShuffleRealTruncate(ReIm *Fx, double *ReX, int N) {

	//INVERSE METHOD
	//Even-indexed inputs are always distributed to lower half of output de-shuffled (or shuffled) values : we only want the lower half.
	//This is routine for real output from half-length complex IFFT :real output values are set as even indexed real output, complex output values are set as odd indexed real output 
	//We have 2N real output points (so the IFFT which should have preceeded the call to this method should have been a N-point complex IFFT), but we only write the first N points out and ignore the upper half.

	for(int p = 0; p < N; p+=2) {

		int index = fftShuffleInv[N - 1 + p];
		
		ReX[2*index] = Fx[p].Re;
		ReX[2*index + 1] = Fx[p].Im;
	}
}

void FFTMethods_Cpp::UnShuffleRealTruncate(ReIm3 *Fx, DBL3 *ReX, int n, int N) {
	
	//INVERSE METHOD
	//Even-indexed inputs are always distributed to lower half of output de-shuffled (or shuffled) values : we only want the lower half.
	//This is routine for real output from half-length complex IFFT :real output values are set as even indexed real output, complex output values are set as odd indexed real output 
	//We have 2N real output points (so the IFFT which should have preceeded the call to this method should have been a N-point complex IFFT), but we only write the first N points out and ignore the upper half.

	//n is the dimension of ReX along x. N is a number greater or equal to n, which is a power of 2.
	//n could be smaller than N, in which case only extract values with indices up to n - 1.

	for(int p = 0; p < N; p+=2) {

		int index = fftShuffleInv[N - 1 + p];
		
		if(2*index < n)	ReX[2*index] = Fx[p].Re();
		if(2*index + 1 < n)	ReX[2*index + 1] = Fx[p].Im();
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Advanced routines combining 1D FFTs

void FFTMethods_Cpp::FFT3D_R2C(double *Src, ReIm *Dst, ReIm *FFTSpace, int N, int M, int K, int powN, int powM, int powK) {

	//3D FFT real to complex. Output is transposed.
	//

	for(int k = 0; k < K; k++) {		
		#pragma omp parallel for
		for(int m = 0; m < M; m++) {
		
			CopyRealShuffle(Src + m*N + k*N*M, FFTSpace + m*N/2 + k*N*M/2, N/2);
			FFT_Radix4_DIT(FFTSpace + m*N/2 + k*N*M/2, powN - 1, N/2);
			RealfromComplexFFT(FFTSpace + m*N/2 + k*N*M/2, Dst + m*(N/2 + 1) + k*(N/2 + 1)*M, N/2);
		}
	}

	for(int k = 0; k < K; k++) {
		#pragma omp parallel for
		for(int n = 0; n <= N/2; n++) {
		
			CopyShuffle(Dst + n + k*(N/2 + 1)*M, FFTSpace + n*M + k*(N/2 + 1)*M, N/2 + 1, M);
			FFT_Radix4_DIT(FFTSpace + n*M + k*(N/2 + 1)*M, powM, M);
		}
	}

	#pragma omp parallel for
	for(int n = 0; n <= N/2; n++) {
		for(int m = 0; m < M; m++) {
		
			CopyShuffle(FFTSpace + m + n*M, Dst + n*K + m*(N/2 + 1)*K, (N/2 + 1)*M, K);
			FFT_Radix4_DIT(Dst + n*K + m*(N/2 + 1)*K, powK, K);
		}
	}
}

void FFTMethods_Cpp::FFT3D_RZP2C(double *Src, ReIm *Dst, ReIm *FFTSpace, int N, int M, int K, int powN, int powM, int powK) {

	//3D FFT real to complex with (N/2) * (M/2) * (K/2) input zero padded to dimensions N*M*K. Output is transposed.
	//

	//FFT current input
	for(int k = 0; k < K/2; k++) {
		#pragma omp parallel for
		for(int m = 0; m < M/2; m++) {
		
			CopyRealShuffleZeroPad(Src + m*(N/2) + k*(N/2)*(M/2), FFTSpace + m*N/2 + k*N*M/2, N/2);
			FFT_Radix4_DIT(FFTSpace + m*N/2 + k*N*M/2, powN - 1, N/2);
			RealfromComplexFFT(FFTSpace + m*N/2 + k*N*M/2, Dst + m*(N/2 + 1) + k*(N/2 + 1)*M, N/2);
		}
	}

	for(int k = 0; k < K/2; k++) {
		#pragma omp parallel for
		for(int n = 0; n <= N/2; n++) {
		
			CopyShuffleZeroPad(Dst + n + k*(N/2 + 1)*M, FFTSpace + n*M + k*(N/2 + 1)*M, N/2 + 1, M);
			FFT_Radix4_DIT(FFTSpace + n*M + k*(N/2 + 1)*M, powM, M);
		}
	}

	#pragma omp parallel for
	for(int m = 0; m < M; m++) {
		for(int n = 0; n <= N/2; n++) {

			CopyShuffleZeroPad(FFTSpace + m + n*M, Dst + n*K + m*(N/2 + 1)*K, (N/2 + 1)*M, K);
			FFT_Radix4_DIT(Dst + n*K + m*(N/2 + 1)*K, powK, K);
		}
	}
}

void FFTMethods_Cpp::IFFT3D_C2RT(ReIm *Src, double *Dst, ReIm *FFTSpace, int N, int M, int K, int powN, int powM, int powK) {

	//Inverse 3D FFT from complex input (N/2 + 1)*M*K to truncated real output of dimension (N/2) * (M/2) * (K/2)

	#pragma omp parallel for
	for(int m = 0; m < M; m++) {
		for(int n = 0; n <= N/2; n++) {

			IFFT_Radix4_DIF(Src + n*K + m*(N/2 + 1)*K, powK, K);
			UnShuffleTruncate(Src + n*K + m*(N/2 + 1)*K, FFTSpace + n*M + m, (N/2 + 1)*M, K);
		}
	}

	for(int k = 0; k < K/2; k++) {
		#pragma omp parallel for
		for(int n = 0; n <= N/2; n++) {

			IFFT_Radix4_DIF(FFTSpace + n*M + k*(N/2 + 1)*M, powM, M);
			UnShuffleTruncate(FFTSpace + n*M + k*(N/2 + 1)*M, Src + n + k*(N/2 + 1)*M, N/2 + 1, M);
		}
	}

	for(int k = 0; k < K/2; k++) {
		#pragma omp parallel for
		for(int m = 0; m < M/2; m++) {

			RealfromComplexIFFT(Src + m*(N/2 + 1) + k*(N/2 + 1)*M, FFTSpace + m*N/2 + k*N*M/2, N/2);
			IFFT_Radix4_DIF(FFTSpace + m*N/2 + k*N*M/2, powN - 1, N/2);
			UnShuffleRealTruncate(FFTSpace + m*N/2 + k*N*M/2, Dst + m*(N/2) + k*(N/2)*(M/2), N/2);
		}
	}

	//Finish IFFT
	#pragma omp parallel for
	for(int n = 0; n < (N/2)*(M/2)*(K/2); n++) {

		Dst[n] = Dst[n] / ((N/2)*M*K);
	}
}

void FFTMethods_Cpp::FFT2D_R2C(double *Src, ReIm *Dst, ReIm *FFTSpace, int N, int M, int powN, int powM) {

	//2D FFT real to complex. Output is transposed.
	//

	#pragma omp parallel for
	for(int m = 0; m < M; m++) {
		
		CopyRealShuffle(Src + m*N, Dst + m*N/2, N/2);
		FFT_Radix4_DIT(Dst + m*N/2, powN - 1, N/2);
		RealfromComplexFFT(Dst + m*N/2, FFTSpace + m*(N/2 + 1), N/2);
	}

	#pragma omp parallel for
	for(int n = 0; n <= N/2; n++) {
		
		CopyShuffle(FFTSpace + n, Dst + n*M, N/2 + 1, M);
		FFT_Radix4_DIT(Dst + n*M, powM, M);
	}
}

void FFTMethods_Cpp::FFT2D_RZP2C(double *Src, ReIm *Dst, ReIm *FFTSpace, int N, int M, int powN, int powM) {

	//2D FFT real to complex with (N/2) * (M/2) input zero padded to dimensions N*M*K. Output is transposed.
	//

	//FFT current input

	#pragma omp parallel for
	for(int m = 0; m < M/2; m++) {
		
		CopyRealShuffleZeroPad(Src + m*(N/2), Dst + m*N/2, N/2);
		FFT_Radix4_DIT(Dst + m*N/2, powN - 1, N/2);
		RealfromComplexFFT(Dst + m*N/2, FFTSpace + m*(N/2 + 1), N/2);
	}



	#pragma omp parallel for
	for(int n = 0; n <= N/2; n++) {
		
		CopyShuffleZeroPad(FFTSpace + n, Dst + n*M, N/2 + 1, M);
		FFT_Radix4_DIT(Dst + n*M, powM, M);
	}
}

void FFTMethods_Cpp::IFFT2D_C2RT(ReIm *Src, double *Dst, ReIm *FFTSpace, int N, int M, int powN, int powM) {

	//Inverse 2D FFT from complex input (N/2 + 1)*M to truncated real output of dimension (N/2) * (M/2) * (K/2)


	#pragma omp parallel for
	for(int n = 0; n <= N/2; n++) {

		IFFT_Radix4_DIF(Src + n*M, powM, M);
		UnShuffleTruncate(Src + n*M, FFTSpace + n, N/2 + 1, M);
	}

	#pragma omp parallel for
	for(int m = 0; m < M/2; m++) {

		RealfromComplexIFFT(FFTSpace + m*(N/2 + 1), Src + m*N/2, N/2);
		IFFT_Radix4_DIF(Src + m*N/2, powN - 1, N/2);
		UnShuffleRealTruncate(Src + m*N/2, Dst + m*(N/2), N/2);
	}

	//Finish IFFT
	#pragma omp parallel for
	for(int n = 0; n < (N/2)*(M/2); n++) {

		Dst[n] = Dst[n] / ((N/2)*M);
	}
}

void FFTMethods_Cpp::Product(ReIm *Src1, ReIm *Src2, ReIm *Dst, int elements) {

	#pragma omp parallel for
	for(int i = 0; i < elements; i++) {

		Dst[i] = Src1[i] * Src2[i];
	}
}

