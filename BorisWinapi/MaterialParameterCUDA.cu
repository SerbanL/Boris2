#include "MaterialParameterCUDA.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"

void MatPCUDA<cuBReal, cuBReal>::set_t_equation_from_cpu(TEquationCUDA<cuBReal>& Tscaling_CUDAeq, const std::vector< std::vector<EqComp::FSPEC> >& fspec)
{
	//make CUDA version of text equation for temperature dependence
	Tscaling_CUDAeq.make_scalar(fspec);
}

void MatPCUDA<cuReal2, cuBReal>::set_t_equation_from_cpu(TEquationCUDA<cuBReal>& Tscaling_CUDAeq, const std::vector< std::vector<EqComp::FSPEC> >& fspec)
{
	//make CUDA version of text equation for temperature dependence
	Tscaling_CUDAeq.make_scalar(fspec);
}

void MatPCUDA<cuReal3, cuBReal>::set_t_equation_from_cpu(TEquationCUDA<cuBReal>& Tscaling_CUDAeq, const std::vector< std::vector<EqComp::FSPEC> >& fspec)
{
	//make CUDA version of text equation for temperature dependence
	Tscaling_CUDAeq.make_scalar(fspec);
}

void MatPCUDA<cuReal3, cuReal3>::set_t_equation_from_cpu(TEquationCUDA<cuBReal>& Tscaling_CUDAeq, const std::vector< std::vector<EqComp::FSPEC> >& fspec)
{
	//make CUDA version of text equation for temperature dependence
	Tscaling_CUDAeq.make_scalar(fspec);
}


//set spatial variation equation from cpu version : scalar version
void MatPCUDA<cuBReal, cuBReal>::set_s_equation_from_cpu(TEquationCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sscaling_CUDAeq, const std::vector< std::vector<EqComp::FSPEC> >& fspec)
{
	//make CUDA version of text equation for spatial variation
	Sscaling_CUDAeq.make_scalar(fspec);
}

//set spatial variation equation from cpu version : scalar version
void MatPCUDA<cuReal2, cuBReal>::set_s_equation_from_cpu(TEquationCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sscaling_CUDAeq, const std::vector< std::vector<EqComp::FSPEC> >& fspec)
{
	//make CUDA version of text equation for spatial variation
	Sscaling_CUDAeq.make_scalar(fspec);
}

//set spatial variation equation from cpu version : scalar version
void MatPCUDA<cuReal3, cuBReal>::set_s_equation_from_cpu(TEquationCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sscaling_CUDAeq, const std::vector< std::vector<EqComp::FSPEC> >& fspec)
{
	//make CUDA version of text equation for spatial variation
	Sscaling_CUDAeq.make_scalar(fspec);
}

//set spatial variation equation from cpu version : vector version
void MatPCUDA<cuReal3, cuReal3>::set_s_equation_from_cpu(TEquationCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Sscaling_CUDAeq, const std::vector<std::vector< std::vector<EqComp::FSPEC> >>& fspec)
{
	//make CUDA version of text equation for spatial variation
	Sscaling_CUDAeq.make_vector(fspec);
}

#endif
