#include "stdafx.h"
#include "MaterialParameter.h"

#if COMPILECUDA == 1

//-------------- UPDATE CUDA VALUE ONLY

//update just the value at 0K and current value in the corresponding MatPCUDA
template <>
void MatP<float, double>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuBReal, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

template <>
void MatP<double, double>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuBReal, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

template <>
void MatP<FLT2, double>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal2, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

template <>
void MatP<DBL2, double>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal2, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

template <>
void MatP<FLT3, double>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

template <>
void MatP<DBL3, double>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

template <>
void MatP<FLT3, FLT3>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuReal3>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

template <>
void MatP<DBL3, DBL3>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuReal3>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

//-------------- UPDATE FULL CUDA OBJECT

//fully update the corresponding MatPCUDA
template <>
void MatP<float, double>::update_cuda_object(void)
{
	//temperature equation
	if (Tscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuBReal, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_t_equation_from_cpu(Tscaling_CUDAeq, Tscaling_eq.get_scalar_fspec());
	//clear CUDA version so when we call set_from_cpu below the ManagedFunctionCUDA pointer is nulled.
	else Tscaling_CUDAeq.clear();

	//spatial variation equation
	if (Sscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuBReal, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_s_equation_from_cpu(Sscaling_CUDAeq, Sscaling_eq.get_scalar_fspec());
	//clear CUDA version so when we call set_from_cpu below the ManagedFunctionCUDA pointer is nulled.
	else Sscaling_CUDAeq.clear();

	reinterpret_cast<cu_obj<MatPCUDA<cuBReal, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

template <>
void MatP<double, double>::update_cuda_object(void)
{
	//temperature equation
	if (Tscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuBReal, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_t_equation_from_cpu(Tscaling_CUDAeq, Tscaling_eq.get_scalar_fspec());
	else Tscaling_CUDAeq.clear();

	//spatial variation equation
	if (Sscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuBReal, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_s_equation_from_cpu(Sscaling_CUDAeq, Sscaling_eq.get_scalar_fspec());
	else Sscaling_CUDAeq.clear();

	reinterpret_cast<cu_obj<MatPCUDA<cuBReal, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

template <>
void MatP<FLT2, double>::update_cuda_object(void)
{
	//temperature equation
	if (Tscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuReal2, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_t_equation_from_cpu(Tscaling_CUDAeq, Tscaling_eq.get_dual_fspec());
	else Tscaling_CUDAeq.clear();

	//spatial variation equation
	if (Sscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuReal2, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_s_equation_from_cpu(Sscaling_CUDAeq, Sscaling_eq.get_scalar_fspec());
	else Sscaling_CUDAeq.clear();

	reinterpret_cast<cu_obj<MatPCUDA<cuReal2, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

template <>
void MatP<DBL2, double>::update_cuda_object(void)
{
	//temperature equation
	if (Tscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuReal2, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_t_equation_from_cpu(Tscaling_CUDAeq, Tscaling_eq.get_dual_fspec());
	else Tscaling_CUDAeq.clear();

	//spatial variation equation
	if (Sscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuReal2, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_s_equation_from_cpu(Sscaling_CUDAeq, Sscaling_eq.get_scalar_fspec());
	else Sscaling_CUDAeq.clear();

	reinterpret_cast<cu_obj<MatPCUDA<cuReal2, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

template <>
void MatP<FLT3, double>::update_cuda_object(void)
{
	//temperature equation
	if (Tscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_t_equation_from_cpu(Tscaling_CUDAeq, Tscaling_eq.get_vector_fspec());
	else Tscaling_CUDAeq.clear();

	//spatial variation equation
	if (Sscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_s_equation_from_cpu(Sscaling_CUDAeq, Sscaling_eq.get_scalar_fspec());
	else Sscaling_CUDAeq.clear();

	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

template <>
void MatP<DBL3, double>::update_cuda_object(void)
{
	//temperature equation
	if (Tscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_t_equation_from_cpu(Tscaling_CUDAeq, Tscaling_eq.get_vector_fspec());
	else Tscaling_CUDAeq.clear();

	//spatial variation equation
	if (Sscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_s_equation_from_cpu(Sscaling_CUDAeq, Sscaling_eq.get_scalar_fspec());
	else Sscaling_CUDAeq.clear();

	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

template <>
void MatP<FLT3, FLT3>::update_cuda_object(void)
{
	//temperature equation
	if (Tscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuReal3>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_t_equation_from_cpu(Tscaling_CUDAeq, Tscaling_eq.get_vector_fspec());
	else Tscaling_CUDAeq.clear();

	//spatial variation equation
	if (Sscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuReal3>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_s_equation_from_cpu(Sscaling_CUDAeq, Sscaling_eq.get_vector_fspec());
	else Sscaling_CUDAeq.clear();

	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuReal3>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

template <>
void MatP<DBL3, DBL3>::update_cuda_object(void)
{
	//temperature equation
	if (Tscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuReal3>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_t_equation_from_cpu(Tscaling_CUDAeq, Tscaling_eq.get_vector_fspec());
	else Tscaling_CUDAeq.clear();

	//spatial variation equation
	if (Sscaling_eq.is_set()) reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuReal3>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_s_equation_from_cpu(Sscaling_CUDAeq, Sscaling_eq.get_vector_fspec());
	else Sscaling_CUDAeq.clear();

	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuReal3>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

#endif