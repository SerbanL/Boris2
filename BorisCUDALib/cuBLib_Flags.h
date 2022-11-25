#pragma once

//make sure to compile with matching architecture, e.g. compute_50, sm_50 for __CUDA_ARCH__ 500, etc.
#define __CUDA_ARCH__ 800

//compile with cuda single precision (float types) or double precision (double types). Set SINGLEPRECISION to 1 for single, otherwise (0) for double.
#define SINGLEPRECISION 1
