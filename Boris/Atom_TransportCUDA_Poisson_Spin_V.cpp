#include "stdafx.h"
#include "TransportCUDA_Poisson_Spin_V.h"

#if COMPILECUDA == 1

#if defined(MODULE_COMPILATION_TRANSPORT) && ATOMISTIC == 1

#include "Atom_MeshCUDA.h"
#include "Atom_TransportCUDA.h"

BError Atom_TransportCUDA_Spin_V_Funcs::set_pointers(Atom_MeshCUDA* paMeshCUDA, Atom_TransportCUDA* pTransportCUDA)
{
	BError error(__FUNCTION__);

	//Mesh quantities

	if (set_gpu_value(pcuaMesh, paMeshCUDA->cuaMesh.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(pdM_dt, pTransportCUDA->dM_dt.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);
	if (set_gpu_value(pdelsq_V_fixed, pTransportCUDA->delsq_V_fixed.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	if (set_gpu_value(stsolve, pTransportCUDA->Get_STSolveType()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	return error;
}

#endif

#endif