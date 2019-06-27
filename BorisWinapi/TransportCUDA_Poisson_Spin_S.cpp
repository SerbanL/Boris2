#include "stdafx.h"
#include "TransportCUDA_Poisson_Spin_S.h"

#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include "MeshCUDA.h"
#include "DiffEqCUDA.h"

BError TransportCUDA_Spin_S_Funcs::set_pointers(MeshCUDA* pMeshCUDA, DifferentialEquationCUDA* pdiffEqCUDA)
{
	BError error(__FUNCTION__);

	if (set_gpu_value(pcuMesh, pMeshCUDA->cuMesh.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (pdiffEqCUDA) {

		if (set_gpu_value(pcuDiffEq, pdiffEqCUDA->Get_ManagedDiffEqCUDA().get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	}
	else {

		nullgpuptr(pcuDiffEq);
	}

	return error;
}

#endif

#endif