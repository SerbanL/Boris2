#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ParametersDefs.h"

class MatPFormulaCUDA {

protected:

	//selector for pre-defined temperature scaling formula (if this is MATPFORM_NONE then t_scaling is used if set)
	int formula_selector;

	//coefficients
	cuBReal *pcoeff;

private:

	//--------- MATPFORM_NONE

	__device__ cuBReal __MATPFORM_NONE(cuBReal T) const { return 1.0; }

	//--------- MATPFORM_LINEAR

	// y = x0 * T + 1
	__device__ cuBReal __MATPFORM_LINEAR(cuBReal T) const { return (pcoeff[0] * T + 1); }

	//--------- MATPFORM_PARABOLIC

	// y = x0 * T^2 + x1 * T + 1
	__device__ cuBReal __MATPFORM_PARABOLIC(cuBReal T) const { return (pcoeff[0] * T*T + pcoeff[1] * T + 1); }

	//--------- MATPFORM_INVERSELINEAR

	// y = 1 / (x0 * T + 1). coeff[0] not allowed to be negative.
	__device__ cuBReal __MATPFORM_INVERSELINEAR(cuBReal T) const { return 1 / (pcoeff[0] * T + 1); }

protected:

	//---------Constructors / destructor (cu_obj managed)

	//void constructor
	__host__ void construct_cu_obj(void)
	{
		set_formula_selector(MATPFORM_NONE);

		nullgpuptr(pcoeff);
		gpu_alloc_managed(pcoeff, MAXFORMULACOEFFICIENTS);
		gpu_set_managed(pcoeff, (cuBReal)0.0, MAXFORMULACOEFFICIENTS);
	}

	//destructor
	__host__ void destruct_cu_obj(void)
	{
		gpu_free_managed(pcoeff);
	}

	//---------Helpers

	__host__ void set_formula_selector(int value)
	{
		set_gpu_value(formula_selector, value);
	}

	//---------Get scaling value

	//get value at given temperature for set formula
	__device__ cuBReal t_scaling_formula(cuBReal T)
	{
		switch (formula_selector) {

		case MATPFORM_NONE:
			return __MATPFORM_NONE(T);

		case MATPFORM_LINEAR:
			return __MATPFORM_LINEAR(T);

		case MATPFORM_PARABOLIC:
			return __MATPFORM_PARABOLIC(T);

		case MATPFORM_INVERSELINEAR:
			return __MATPFORM_INVERSELINEAR(T);
		}

		return 1.0;
	}
};

#endif