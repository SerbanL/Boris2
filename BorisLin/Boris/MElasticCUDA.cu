#include "MElasticCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_MELASTIC

#include "BorisCUDALib.cuh"

//------------------- Configuration

//set diagonal and shear strain text equations
BError MElasticCUDA::Set_Sd_Equation(std::vector<std::vector< std::vector<EqComp::FSPEC> >> fspec)
{
	BError error(CLASS_STR(MElasticCUDA));

	if (!Sd_equation.make_vector(fspec)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

BError MElasticCUDA::Set_Sod_Equation(std::vector<std::vector< std::vector<EqComp::FSPEC> >> fspec)
{
	BError error(CLASS_STR(MElasticCUDA));

	if (!Sod_equation.make_vector(fspec)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	return error;
}

#endif

#endif