#include "stdafx.h"
#include "STransport.h"

#ifdef MODULE_COMPILATION_TRANSPORT

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
	double total_current_out = 0.0;

	for (int idx = 0; idx < (int)pTransport.size(); idx++) {

		pTransport[idx]->CalculateElectricField();

		//in open potential mode count total thermoelectric current generated
		if (IsOpenPotential()) total_current_out += pTransport[idx]->mesh_thermoelectric_net_current;
	}

	//in open potential mode count total thermoelectric current generated and set potential depending on open potential resistance set
	if (IsOpenPotential()) SetPotential(total_current_out * open_potential_resistance);

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

			pV[idx_pri]->set_cmbnd_continuous<TransportBase>(
				*pV[idx_sec], CMBNDcontacts[idx1][idx2],
				&TransportBase::afunc_V_sec, &TransportBase::afunc_V_pri,
				&TransportBase::bfunc_V_sec, &TransportBase::bfunc_V_pri,
				&TransportBase::diff2_V_sec, &TransportBase::diff2_V_pri,
				*pTransport[idx_sec], *pTransport[idx_pri]);
		}
	}
}

#endif
