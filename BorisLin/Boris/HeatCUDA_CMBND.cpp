#include "stdafx.h"
#include "HeatCUDA_CMBND.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_HEAT

#include "MeshCUDA.h"
#include "Atom_MeshCUDA.h"

BError HeatCUDA_CMBND::set_pointers(MeshCUDA* pMeshCUDA)
{
	BError error(__FUNCTION__);

	if (set_gpu_value(pcuMesh, pMeshCUDA->cuMesh.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	return error;
}

BError HeatCUDA_CMBND::set_pointers(Atom_MeshCUDA* paMeshCUDA)
{
	BError error(__FUNCTION__);

	if (set_gpu_value(pcuaMesh, paMeshCUDA->cuaMesh.get_managed_object()) != cudaSuccess) error(BERROR_GPUERROR_CRIT);

	return error;
}

//set pQ_equation as needed
BError HeatCUDA_CMBND::set_Q_equation(TEquationCUDA<cuBReal, cuBReal, cuBReal, cuBReal>& Q_equation_ref)
{
	BError error(__FUNCTION__);

	if (Q_equation_ref.is_set()) {
		
		if (set_gpu_value(pQ_equation, Q_equation_ref.get_pcu_obj_x()->get_managed_object()) != cudaSuccess) return error(BERROR_GPUERROR_CRIT);
	}
	else {

		nullgpuptr(pQ_equation);
	}

	return error;
}

#endif

#endif