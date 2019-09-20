#include "stdafx.h"
#include "MaterialParameter.h"

#if COMPILECUDA == 1

//-------------- UPDATE CUDA VALUE ONLY

//update just the value at 0K and current value in the corresponding MatPCUDA
void MatP<float, double>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuBReal, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

void MatP<double, double>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuBReal, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

void MatP<FLT2, double>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal2, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

void MatP<DBL2, double>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal2, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

void MatP<FLT3, double>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

void MatP<DBL3, double>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

void MatP<FLT3, FLT3>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuReal3>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

void MatP<DBL3, DBL3>::update_cuda_value(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuReal3>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu_value(*this);
}

//-------------- UPDATE FULL CUDA OBJECT

//update just the value at 0K and current value in the corresponding MatPCUDA
void MatP<float, double>::update_cuda_object(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuBReal, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

void MatP<double, double>::update_cuda_object(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuBReal, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

void MatP<FLT2, double>::update_cuda_object(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal2, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

void MatP<DBL2, double>::update_cuda_object(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal2, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

void MatP<FLT3, double>::update_cuda_object(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

void MatP<DBL3, double>::update_cuda_object(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuBReal>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

void MatP<FLT3, FLT3>::update_cuda_object(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuReal3>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

void MatP<DBL3, DBL3>::update_cuda_object(void)
{
	reinterpret_cast<cu_obj<MatPCUDA<cuReal3, cuReal3>>*>(p_cu_obj_mpcuda)->get_managed_object()->set_from_cpu(*this);
}

#endif