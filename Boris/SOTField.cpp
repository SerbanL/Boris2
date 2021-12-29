#include "stdafx.h"
#include "SOTField.h"

#ifdef MODULE_COMPILATION_SOTFIELD

#include "Mesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "SOTFieldCUDA.h"
#endif

SOTField::SOTField(Mesh *pMesh_) :
	Modules(),
	ProgramStateNames(this, {}, {})
{
	pMesh = pMesh_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

SOTField::~SOTField()
{
}

BError SOTField::Initialize(void)
{
	BError error(CLASS_STR(SOTField));

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect, 
		(MOD_)pMesh->Get_Module_Heff_Display() == MOD_SOTFIELD, 
		(MOD_)pMesh->Get_Module_Energy_Display() == MOD_SOTFIELD, 
		pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error)	initialized = true;

	return error;
}

BError SOTField::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(SOTField));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError SOTField::MakeCUDAModule(void)
{
	BError error(CLASS_STR(SOTField));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new SOTFieldCUDA(pMesh->pMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double SOTField::UpdateField(void)
{
	if (!pMesh->E.linear_size()) return 0.0;

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	
	if (pMesh->GetMeshType() == MESH_FERROMAGNETIC) {
		
#pragma omp parallel for
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			double grel = pMesh->grel;
			double Ms = pMesh->Ms;
			double SHA = pMesh->SHA;
			double flSOT = pMesh->flSOT;
			DBL3 STp = pMesh->STp;
			pMesh->update_parameters_mcoarse(idx, pMesh->grel, grel, pMesh->Ms, Ms, pMesh->SHA, SHA, pMesh->flSOT, flSOT, pMesh->STp, STp);

			if (IsNZ(grel)) {

				double a_const = -(SHA * MUB_E / (GAMMA * grel)) / (Ms * Ms * pMesh->GetMeshDimensions().z);

				int idx_E = pMesh->E.position_to_cellidx(pMesh->M.cellidx_to_position(idx));
				DBL3 p_vec = (STp ^ pMesh->E[idx_E]) * pMesh->elC[idx_E];

				DBL3 SOTField = a_const * ((pMesh->M[idx] ^ p_vec) + flSOT * Ms * p_vec);
				pMesh->Heff[idx] += SOTField;

				if (Module_Heff.linear_size()) Module_Heff[idx] = SOTField;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			DBL2 grel_AFM = pMesh->grel_AFM;
			DBL2 Ms_AFM = pMesh->Ms_AFM;
			double SHA = pMesh->SHA;
			double flSOT = pMesh->flSOT;
			DBL3 STp = pMesh->STp;
			pMesh->update_parameters_mcoarse(idx, pMesh->grel_AFM, grel_AFM, pMesh->Ms_AFM, Ms_AFM, pMesh->SHA, SHA, pMesh->flSOT, flSOT, pMesh->STp, STp);

			if (IsNZ(grel_AFM.i + grel_AFM.j)) {

				double a_const_A = -(SHA * MUB_E / (GAMMA * grel_AFM.i)) / (Ms_AFM.i * Ms_AFM.i * pMesh->GetMeshDimensions().z);
				double a_const_B = -(SHA * MUB_E / (GAMMA * grel_AFM.j)) / (Ms_AFM.j * Ms_AFM.j * pMesh->GetMeshDimensions().z);

				int idx_E = pMesh->E.position_to_cellidx(pMesh->M.cellidx_to_position(idx));
				DBL3 p_vec = (STp ^ pMesh->E[idx_E]) * pMesh->elC[idx_E];

				DBL3 SOTField_A = a_const_A * ((pMesh->M[idx] ^ p_vec) + flSOT * Ms_AFM.i * p_vec);
				DBL3 SOTField_B = a_const_B * ((pMesh->M2[idx] ^ p_vec) + flSOT * Ms_AFM.j * p_vec);

				pMesh->Heff[idx] += SOTField_A;
				pMesh->Heff2[idx] += SOTField_B;

				if (Module_Heff.linear_size()) Module_Heff[idx] = SOTField_A;
				if (Module_Heff2.linear_size()) Module_Heff2[idx] = SOTField_B;
			}
		}
	}

	//don't count this as a contribution to the total energy density
	return 0.0;
}

#endif