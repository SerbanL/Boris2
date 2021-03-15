#include "stdafx.h"
#include "ManagedModulesCUDA.h"
#include "ModulesCUDA.h"

#if COMPILECUDA == 1

BError ManagedModulesCUDA::set_pointers(ModulesCUDA* pModulesCUDA)
{
	BError error(__FUNCTION__);

	//Mesh quantities

	if (set_gpu_value(penergy, pModulesCUDA->energy.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(ppoints_count, pModulesCUDA->points_count.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pModule_Heff, pModulesCUDA->Module_Heff.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pModule_Heff2, pModulesCUDA->Module_Heff2.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pModule_energy, pModulesCUDA->Module_energy.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pModule_energy2, pModulesCUDA->Module_energy2.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	return error;
}

#endif