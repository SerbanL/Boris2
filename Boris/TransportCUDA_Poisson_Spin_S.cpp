#include "stdafx.h"
#include "TransportCUDA_Poisson_Spin_S.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "MeshCUDA.h"
#include "DiffEqFMCUDA.h"
#include "TransportCUDA.h"

BError TransportCUDA_Spin_S_Funcs::set_pointers(MeshCUDA* pMeshCUDA, DifferentialEquationFMCUDA* pdiffEqCUDA, TransportCUDA* pTransportCUDA)
{
	BError error(__FUNCTION__);

	if (set_gpu_value(pcuMesh, pMeshCUDA->cuMesh.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (pdiffEqCUDA) {

		if (set_gpu_value(pcuDiffEq, pdiffEqCUDA->Get_ManagedDiffEqCUDA().get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	}
	else {

		nullgpuptr(pcuDiffEq);
	}

	if (set_gpu_value(pPoisson_Spin_V, pTransportCUDA->poisson_Spin_V.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pdM_dt, pTransportCUDA->dM_dt.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pdelsq_S_fixed, pTransportCUDA->delsq_S_fixed.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(stsolve, pTransportCUDA->Get_STSolveType()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	return error;
}

#endif

#endif