#pragma once

#include "cuBLib_Flags.h"

#include <cufft.h>

#if SINGLEPRECISION == 1

#define cuBReal cufftReal
#define cuBComplex cufftComplex

#else

#define cuBReal cufftDoubleReal
#define cuBComplex cufftDoubleComplex

#endif