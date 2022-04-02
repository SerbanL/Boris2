#include "stdafx.h"
#include "STransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "STransport.h"
#include "SuperMesh.h"

void STransportCUDA::solve_charge_transport_sor(void)
{
	cuReal2 normalized_max_error = cuReal2();

	pSTrans->iters_to_conv = 0;

	do {

		Zero_Errors();
		normalized_max_error = cuReal2();

		//1. solve V in each mesh separately (1 iteration each) - but do not set CMBND cells yet
		for (int idx = 0; idx < (int)pTransport.size(); idx++) {

			pTransport[idx]->IterateChargeSolver_SOR(SOR_damping_V, max_error, max_value);
		}

		//normalize error to maximum change in cpu memory
		normalized_max_error = cuReal2(max_error.to_cpu(), max_value.to_cpu());
		normalized_max_error.first = (normalized_max_error.second > 0 ? normalized_max_error.first / normalized_max_error.second : normalized_max_error.first);

		//2. now set CMBND cells
		set_cmbnd_charge_transport();

		pSTrans->iters_to_conv++;

	} while (normalized_max_error.first > pSTrans->errorMaxLaplace && pSTrans->iters_to_conv < pSTrans->maxLaplaceIterations);

	//continue next iteration if iterations timeout reached - with this timeout built in the program doesn't block if errorMaxLaplace cannot be reached. 
	if (pSTrans->iters_to_conv == pSTrans->maxLaplaceIterations) pSTrans->recalculate_transport = true;

	//2. update E in all meshes
	for (int idx = 0; idx < (int)pTransport.size(); idx++) {

		pTransport[idx]->CalculateElectricField();
	}

	//store the current max error in the energy term so it can be read if requested
	energy.from_cpu(normalized_max_error.first);
}

#endif

#endif