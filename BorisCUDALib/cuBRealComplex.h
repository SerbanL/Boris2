#pragma once

#include "cuBLib_Flags.h"

#if SINGLEPRECISION == 1

#define cuReal cufftReal
#define cuComplex cufftComplex

#else

#define cuReal cufftDoubleReal
#define cuComplex cufftDoubleComplex

#endif