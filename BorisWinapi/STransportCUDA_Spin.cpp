#include "stdafx.h"
#include "STransportCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include "STransport.h"
#include "SuperMesh.h"

void STransportCUDA::solve_spin_transport_sor(void)
{
	cuReal2 normalized_max_error = cuReal2();

	pSTrans->iters_to_conv = 0;

	bool start_iters = true;

	//1. Solve V everywhere for current S until convergence criteria hit

	do {

		Zero_Errors();
		normalized_max_error = cuReal2();

		//1. solve V in each mesh separately (1 iteration each) - but do not set CMBND cells yet
		for (int idx = 0; idx < (int)pTransport.size(); idx++) {

			//use non-homogeneous Neumann boundary conditions for V? Only use them if iSHE is enabled and not a magnetic mesh
			bool use_NNeu = IsNZ((double)pTransport[idx]->pMesh->iSHA) && !pTransport[idx]->pMesh->M.linear_size();

			if (pSTrans->fixed_SOR_damping) pTransport[idx]->IterateSpinSolver_Charge_SOR(SOR_damping_V, max_error, max_value, use_NNeu);
			else pTransport[idx]->IterateSpinSolver_Charge_aSOR(start_iters, pSTrans->s_errorMax, max_error, max_value, use_NNeu);
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

		pTransport[idx]->CalculateCurrentDensity();
	}

	//--------------

	//3. Solve S everywhere based on obtained Jc until convergence criteria hit

	pSTrans->s_iters_to_conv = 0;

	start_iters = true;

	do {

		pSTrans->aSOR_damping = DBL2();

		Zero_Errors();
		normalized_max_error = cuReal2();

		//solve S in each mesh separately (1 iteration each) - but do not set CMBND cells yet
		for (int idx = 0; idx < (int)pTransport.size(); idx++) {

			//use non-homogeneous Neumann boundary conditions for S? Only use them if SHE is enabled and not a magnetic mesh
			bool use_NNeu = IsNZ((double)pTransport[idx]->pMesh->SHA) && !pTransport[idx]->pMesh->M.linear_size();

			if (pSTrans->fixed_SOR_damping) pTransport[idx]->IterateSpinSolver_Spin_SOR(SOR_damping_S, max_error, max_value, use_NNeu);
			else pTransport[idx]->IterateSpinSolver_Spin_aSOR(start_iters, pSTrans->s_errorMax, max_error, max_value, use_NNeu);

			//store minimum and maximum damping values
			double damping = (*pS[idx])()->aSOR_get_damping_cpu();
			if (pSTrans->aSOR_damping == DBL2()) pSTrans->aSOR_damping = DBL2(damping, damping);
			else pSTrans->aSOR_damping = DBL2(min(pSTrans->aSOR_damping.i, damping), max(pSTrans->aSOR_damping.j, damping));
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

			pTransport[idx_pri]->CalculateSAInterfaceField(pTransport[idx_sec], CMBNDcontactsCUDA[idx1][idx2], CMBNDcontacts[idx1][idx2].IsPrimaryTop());
		}
	}
}

//-------------------Getters

//return interfacial spin torque in given mesh with matching transport module
cu_obj<cuVEC<cuReal3>>& STransportCUDA::GetInterfacialSpinTorque(TransportCUDA* pMeshTrans)
{
	if (!pMeshTrans->PrepareDisplayVEC(pMeshTrans->pMesh->h))
		return pMeshTrans->displayVEC;

	//calculate and add Ts values in Ts_interf VECs
	for (int idx1 = 0; idx1 < (int)CMBNDcontacts.size(); idx1++) {

		for (int idx2 = 0; idx2 < (int)CMBNDcontacts[idx1].size(); idx2++) {

			int idx_sec = CMBNDcontacts[idx1][idx2].mesh_idx.i;
			int idx_pri = CMBNDcontacts[idx1][idx2].mesh_idx.j;

			if (pTransport[idx_pri] == pMeshTrans)
				pTransport[idx_pri]->CalculateDisplaySAInterfaceTorque(pTransport[idx_sec], CMBNDcontactsCUDA[idx1][idx2], CMBNDcontacts[idx1][idx2].IsPrimaryTop());
		}
	}

	return pMeshTrans->displayVEC;
}

#endif

#endif