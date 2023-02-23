#include "stdafx.h"
#include "TMR.h"

#ifdef MODULE_COMPILATION_TMR

#include "Mesh.h"
#include "SuperMesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "TMRCUDA.h"
#endif

//calculate charge current density over the mesh : applies to both charge-only transport and spin transport solvers (if check used inside this method)
//VERIFIED - CORRECT
VEC_VC<DBL3>& TMR::GetChargeCurrent(void)
{
	//make sure memory is allocated to the correct size
	if (!PrepareDisplayVEC_VC(pMesh->h_e)) return displayVEC_VC;

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		cu_obj<cuVEC_VC<cuReal3>>& cudisplayVEC_VC = GetChargeCurrentCUDA();
		cudisplayVEC_VC()->copy_to_cpuvec(displayVEC_VC);

		return displayVEC_VC;
	}
#endif

	//compute charge current and store result in displayVEC_VC

	if (!pSMesh->disabled_transport_solver) {

		//calculate current density using Jc = -sigma * grad V
#pragma omp parallel for
		for (int idx = 0; idx < displayVEC_VC.linear_size(); idx++) {

			//only calculate current on non-empty cells - empty cells have already been assigned 0 at UpdateConfiguration
			if (pMesh->V.is_not_empty(idx)) {

				displayVEC_VC[idx] = -pMesh->elC[idx] * pMesh->V.grad_diri(idx);
			}
			else displayVEC_VC[idx] = DBL3(0);
		}
	}
	else {

		//if transport solver disabled we need to set displayVEC_VC directly from E and elC as Jc = elC * E
#pragma omp parallel for
		for (int idx = 0; idx < pMesh->E.linear_size(); idx++) {

			if (pMesh->elC.is_not_empty(idx)) {

				displayVEC_VC[idx] = pMesh->E[idx] * pMesh->elC[idx];
			}
			else {

				displayVEC_VC[idx] = DBL3(0.0);
			}
		}
	}

	return displayVEC_VC;
}

#endif