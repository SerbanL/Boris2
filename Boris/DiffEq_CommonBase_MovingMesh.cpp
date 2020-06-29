#include "stdafx.h"
#include "DiffEq_CommonBase.h"

#include "DiffEq_Common.h"
#include "Atom_DiffEq_Common.h"

#include "SuperMesh.h"

//----------------------------------- Moving Mesh Methods

void ODECommon_Base::SetMoveMeshTrigger(bool status, int meshId)
{
	//first turn it off
	for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

		podeSolver->pODE[idx]->pMesh->SetMoveMeshTrigger(false);
	}

	for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

		patom_odeSolver->pODE[idx]->paMesh->SetMoveMeshTrigger(false);
	}

	moving_mesh = false;

	//set to trigger on a specified mesh
	if (meshId >= 0 && status) {

		for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

			if (podeSolver->pODE[idx]->pMesh->get_id() == meshId) {

				podeSolver->pODE[idx]->pMesh->SetMoveMeshTrigger(true);
				moving_mesh = true;
				break;
			}
		}

		if (!moving_mesh) {

			for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

				if (patom_odeSolver->pODE[idx]->paMesh->get_id() == meshId) {

					patom_odeSolver->pODE[idx]->paMesh->SetMoveMeshTrigger(true);
					moving_mesh = true;
					break;
				}
			}
		}
	}
	//set to trigger on first ferromagnetic mesh
	else if (status) {

		if (podeSolver->pODE.size()) {

			podeSolver->pODE[0]->pMesh->SetMoveMeshTrigger(true);
			moving_mesh = true;
		}
		else if (patom_odeSolver->pODE.size()) {

			patom_odeSolver->pODE[0]->paMesh->SetMoveMeshTrigger(true);
			moving_mesh = true;
		}
	}
}

void ODECommon_Base::MovingMeshAlgorithm(SuperMesh* pSMesh)
{
	if (!moving_mesh) return;

	bool shift_done = false;

	for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

		double mesh_shift = podeSolver->pODE[idx]->pMesh->CheckMoveMesh();

		if (IsNZ(mesh_shift)) {

			moving_mesh_dwshift -= mesh_shift;

			//do the actual mesh shifts
			for (int idx_mesh = 0; idx_mesh < pSMesh->size(); idx_mesh++) {

				(*pSMesh)[idx_mesh]->MoveMesh(mesh_shift);
			}

			shift_done = true;

			//only one trigger used
			break;
		}
	}

	if (!shift_done) {

		for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

			double mesh_shift = patom_odeSolver->pODE[idx]->paMesh->CheckMoveMesh();

			if (IsNZ(mesh_shift)) {

				moving_mesh_dwshift -= mesh_shift;

				//do the actual mesh shifts
				for (int idx_mesh = 0; idx_mesh < pSMesh->size(); idx_mesh++) {

					(*pSMesh)[idx_mesh]->MoveMesh(mesh_shift);
				}

				//only one trigger used
				break;
			}
		}
	}
}

int ODECommon_Base::GetId_of_MoveMeshTrigger(void)
{
	if (!moving_mesh) return -1;

	for (int idx = 0; idx < podeSolver->pODE.size(); idx++) {

		if (podeSolver->pODE[idx]->pMesh->GetMoveMeshTrigger()) {

			return podeSolver->pODE[idx]->pMesh->get_id();
		}
	}

	for (int idx = 0; idx < patom_odeSolver->pODE.size(); idx++) {

		if (patom_odeSolver->pODE[idx]->paMesh->GetMoveMeshTrigger()) {

			return patom_odeSolver->pODE[idx]->paMesh->get_id();
		}
	}

	return -1;
}
