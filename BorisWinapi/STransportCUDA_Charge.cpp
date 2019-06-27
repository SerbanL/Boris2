#include "stdafx.h"
#include "STransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include "STransport.h"
#include "SuperMesh.h"

void STransportCUDA::solve_charge_transport_sor(void)
{
	cuReal2 normalized_max_error = cuReal2();

	pSTrans->iters_to_conv = 0;

	bool start_iters = true;

	do {

		pSTrans->aSOR_damping = DBL2();

		Zero_Errors();
		normalized_max_error = cuReal2();

		//1. solve V in each mesh separately (1 iteration each) - but do not set CMBND cells yet
		for (int idx = 0; idx < (int)pTransport.size(); idx++) {

			if (pSTrans->fixed_SOR_damping) pTransport[idx]->IterateChargeSolver_SOR(SOR_damping_V, max_error, max_value);
			else pTransport[idx]->IterateChargeSolver_aSOR(start_iters, pSTrans->s_errorMax, max_error, max_value);

			//store minimum and maximum damping values
			double damping = (*pV[idx])()->aSOR_get_damping_cpu();
			if (pSTrans->aSOR_damping == DBL2()) pSTrans->aSOR_damping = DBL2(damping, damping);
			else pSTrans->aSOR_damping = DBL2(min(pSTrans->aSOR_damping.i, damping), max(pSTrans->aSOR_damping.j, damping));
		}

		//normalize error to maximum change in cpu memory
		normalized_max_error = cuReal2(max_error.to_cpu(), max_value.to_cpu());
		normalized_max_error.first = (normalized_max_error.second > 0 ? normalized_max_error.first / normalized_max_error.second : normalized_max_error.first);

		start_iters = false;

		//2. now set CMBND cells
		set_cmbnd_charge_transport();

		pSTrans->iters_to_conv++;

	} while (normalized_max_error.first > pSTrans->errorMaxLaplace && pSTrans->iters_to_conv < pSTrans->maxLaplaceIterations);

	//continue next iteration if iterations timeout reached - with this timeout built in the program doesn't block if errorMaxLaplace cannot be reached. 
	if (pSTrans->iters_to_conv == pSTrans->maxLaplaceIterations) pSTrans->recalculate_transport = true;

	//3. update Jc in all meshes if a significant change occured
	for (int idx = 0; idx < (int)pTransport.size(); idx++) {

		pTransport[idx]->CalculateCurrentDensity();
	}

	//store the current max error in the energy term so it can be read if requested
	energy.from_cpu(normalized_max_error.first);
}

#endif

#endif