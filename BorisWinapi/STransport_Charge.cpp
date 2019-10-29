#include "stdafx.h"
#include "STransport.h"

#ifdef MODULE_TRANSPORT

#include "SuperMesh.h"

void STransport::solve_charge_transport_sor(void)
{
	DBL2 max_error = DBL2();

	iters_to_conv = 0;

	do {

		//get max error : the max change in V from one iteration to the next
		max_error = DBL2();

		//1. solve V in each mesh separately (1 iteration each) - but do not set CMBND cells yet
		
		for (int idx = 0; idx < (int)pTransport.size(); idx++) {

			DBL2 error;

			error = pTransport[idx]->IterateChargeSolver_SOR(SOR_damping.i);

			if (error.first > max_error.first) max_error.first = error.first;
			if (error.second > max_error.second) max_error.second = error.second;
		}

		//normalize error to maximum change
		max_error.first = (max_error.second > 0 ? max_error.first / max_error.second : max_error.first);

		//2. now set CMBND cells
		set_cmbnd_charge_transport();

		iters_to_conv++;

	} while (max_error.first > errorMaxLaplace && iters_to_conv < maxLaplaceIterations);

	//continue next iteration if iterations timeout reached - with this timeout built in the program doesn't block if errorMaxLaplace cannot be reached. 
	if (iters_to_conv == maxLaplaceIterations) recalculate_transport = true;

	//2. update E in all meshes
	for (int idx = 0; idx < (int)pTransport.size(); idx++) {

		pTransport[idx]->CalculateElectricField();
	}

	//store the current max error in the energy term so it can be read if requested
	energy = max_error.first;
}

//-------------------CMBND computation methods

void STransport::set_cmbnd_charge_transport(void)
{
	for (int idx1 = 0; idx1 < (int)CMBNDcontacts.size(); idx1++) {

		for (int idx2 = 0; idx2 < (int)CMBNDcontacts[idx1].size(); idx2++) {

			int idx_sec = CMBNDcontacts[idx1][idx2].mesh_idx.i;
			int idx_pri = CMBNDcontacts[idx1][idx2].mesh_idx.j;

			pV[idx_pri]->set_cmbnd_continuous<Transport>(
				*pV[idx_sec], CMBNDcontacts[idx1][idx2],
				&Transport::afunc_V_sec, &Transport::afunc_V_pri,
				&Transport::bfunc_V_sec, &Transport::bfunc_V_pri,
				&Transport::diff2_V_sec, &Transport::diff2_V_pri,
				*pTransport[idx_sec], *pTransport[idx_pri]);
		}
	}
}

#endif