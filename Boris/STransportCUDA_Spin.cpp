#include "stdafx.h"
#include "STransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "STransport.h"
#include "SuperMesh.h"

void STransportCUDA::solve_spin_transport_sor(void)
{
	cuReal2 normalized_max_error = cuReal2();

	pSTrans->iters_to_conv = 0;

	bool start_iters = true;

	//Prime the spin solver for the charge part
	for (int idx = 0; idx < (int)pTransport.size(); idx++) {

		if (pTransport[idx]->Get_STSolveType() != STSOLVE_NONE) pTransport[idx]->PrimeSpinSolver_Charge();
	}

	//1. Solve V everywhere for current S until convergence criteria hit

	do {

		Zero_Errors();
		normalized_max_error = cuReal2();

		//1. solve V in each mesh separately (1 iteration each) - but do not set CMBND cells yet
		for (int idx = 0; idx < (int)pTransport.size(); idx++) {

			//use non-homogeneous Neumann boundary conditions for V? Only use them if iSHE is enabled and not a magnetic mesh
			bool use_NNeu = pTransport[idx]->pMeshBase->iSHA_nonzero() && !pTransport[idx]->pMeshBase->Magnetism_Enabled();

			if (pTransport[idx]->Get_STSolveType() != STSOLVE_NONE) pTransport[idx]->IterateSpinSolver_Charge_SOR(SOR_damping_V, max_error, max_value, use_NNeu);
			else pTransport[idx]->IterateChargeSolver_SOR(SOR_damping_V, max_error, max_value);
		}

		//normalize error to maximum change in cpu memory
		normalized_max_error = cuReal2(max_error.to_cpu(), max_value.to_cpu());
		normalized_max_error.first = (normalized_max_error.second > 0 ? normalized_max_error.first / normalized_max_error.second : normalized_max_error.first);

		start_iters = false;

		//now set CMBND cells for V
		set_cmbnd_spin_transport_V();

		pSTrans->iters_to_conv++;

	} while (normalized_max_error.first > pSTrans->errorMaxLaplace && pSTrans->iters_to_conv < pSTrans->maxLaplaceIterations);

	//2. update Jc in all meshes if a significant change occured
	for (int idx = 0; idx < (int)pTransport.size(); idx++) {

		pTransport[idx]->CalculateElectricField();
	}

	//--------------

	//3. Solve S everywhere based on obtained Jc until convergence criteria hit

	pSTrans->s_iters_to_conv = 0;

	start_iters = true;

	//Prime the spin solver for the spin part
	for (int idx = 0; idx < (int)pTransport.size(); idx++) {

		if (pTransport[idx]->Get_STSolveType() != STSOLVE_NONE) pTransport[idx]->PrimeSpinSolver_Spin();
	}

	do {

		Zero_Errors();
		normalized_max_error = cuReal2();

		//solve S in each mesh separately (1 iteration each) - but do not set CMBND cells yet
		for (int idx = 0; idx < (int)pTransport.size(); idx++) {

			//use non-homogeneous Neumann boundary conditions for S? Only use them if SHE is enabled and not a magnetic mesh
			bool use_NNeu = pTransport[idx]->pMeshBase->SHA_nonzero() && !pTransport[idx]->pMeshBase->Magnetism_Enabled();

			if (pTransport[idx]->Get_STSolveType() != STSOLVE_NONE) pTransport[idx]->IterateSpinSolver_Spin_SOR(SOR_damping_S, max_error, max_value, use_NNeu);
		}

		//normalize error to maximum change in cpu memory
		normalized_max_error = cuReal2(max_error.to_cpu(), max_value.to_cpu());
		normalized_max_error.first = (normalized_max_error.second > 0 ? normalized_max_error.first / normalized_max_error.second : normalized_max_error.first);

		start_iters = false;

		//now set CMBND cells for S
		set_cmbnd_spin_transport_S();

		pSTrans->s_iters_to_conv++;

	} while (normalized_max_error.first > pSTrans->s_errorMax && pSTrans->s_iters_to_conv < pSTrans->s_maxIterations);

	//always recalculate transport when using spin current solver - magnetization is bound to change, which will require updating s at least.
	pSTrans->recalculate_transport = true;

	//store the current max error in the energy term so it can be read if requested
	energy.from_cpu(normalized_max_error.first);
}

//------------

//Calculate interface spin accumulation torque (in magnetic meshes for NF interfaces with G interface conductance set)
void STransportCUDA::CalculateSAInterfaceField(void)
{
	//calculate and add Ts values in Ts_interf VECs
	for (int idx1 = 0; idx1 < (int)CMBNDcontacts.size(); idx1++) {

		for (int idx2 = 0; idx2 < (int)CMBNDcontacts[idx1].size(); idx2++) {

			int idx_sec = CMBNDcontacts[idx1][idx2].mesh_idx.i;
			int idx_pri = CMBNDcontacts[idx1][idx2].mesh_idx.j;

			//SA Interface Field currently only from micromagnetic to micromagnetic meshes
			if (!pTransport[idx_pri]->pMeshBase->is_atomistic() && !pTransport[idx_sec]->pMeshBase->is_atomistic()) {

				dynamic_cast<TransportCUDA*>(pTransport[idx_pri])->CalculateSAInterfaceField(dynamic_cast<TransportCUDA*>(pTransport[idx_sec]), CMBNDcontactsCUDA[idx1][idx2], CMBNDcontacts[idx1][idx2].IsPrimaryTop());
			}
		}
	}
}

#endif

#endif