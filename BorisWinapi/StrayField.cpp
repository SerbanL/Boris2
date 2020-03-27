#include "stdafx.h"
#include "StrayField.h"

#ifdef MODULE_STRAYFIELD

#include "Mesh_Ferromagnetic.h"
#include "Mesh_Dipole.h"
#include "SuperMesh.h"
#include "DipoleTFunc.h"

StrayField::StrayField(SuperMesh *pSMesh_) :
	Modules(),
	ProgramStateNames(this, {}, {})
{
	pSMesh = pSMesh_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pSMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

BError StrayField::Initialize(void)
{
	BError error(CLASS_STR(StrayField));

	if (!initialized) {

		dipoles.clear();

		vector< VEC<DBL3>* > meshes_out;

		//identify all existing dipole meshes and output Heff meshes
		for (int idx = 0; idx < (int)pSMesh->size(); idx++) {

			if ((*pSMesh)[idx]->MComputation_Enabled()) {

				meshes_out.push_back(&((*pSMesh)[idx]->Heff));

				//for antiferromagnetic meshes we also need to add the stray field to the second sub-lattice
				if ((*pSMesh)[idx]->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

					meshes_out.push_back(&((*pSMesh)[idx]->Heff2));
				}
			}

			if ((*pSMesh)[idx]->GetMeshType() == MESH_DIPOLE) {

				DipoleMesh *pDipole = reinterpret_cast<DipoleMesh*>((*pSMesh)[idx]);

				//save pointer to dipole mesh
				dipoles.push_back(pDipole);
			}
		}

		//Initialize the mesh transfer object.
		if (!strayField.Initialize_MeshTransfer(vector< VEC<DBL3>* >(), meshes_out, MESHTRANSFERTYPE_WEIGHTED)) return error(BERROR_OUTOFMEMORY_CRIT);

		strayField.set(DBL3(0));

		//here need to calculate stray field from currently set dipole "meshes" (if any - there should be at least one otherwise there's no point using this module)
		CalculateStrayField();

		initialized = true;
	}

	return error;
}

BError StrayField::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(StrayField));

	Uninitialize();

	dipoles.clear();

	if (!strayField.resize(pSMesh->h_fm, pSMesh->sMeshRect_fm)) return error(BERROR_OUTOFMEMORY_CRIT);

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if(pModuleCUDA) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
#endif

	return error;
}

BError StrayField::MakeCUDAModule(void)
{
	BError error(CLASS_STR(StrayField));

#if COMPILECUDA == 1

	pModuleCUDA = new StrayFieldCUDA(pSMesh, this);
	error = pModuleCUDA->Error_On_Create();

#endif

	return error;
}

void StrayField::CalculateStrayField(void)
{
	//set field to zero so we can add strayfield contributions into it
	strayField.set(DBL3(0));

	//origin of super-mesh ferromagnetic rectangle
	DBL3 sMeshOrigin = pSMesh->sMeshRect_fm.get_s();

	for (int idx = 0; idx < (int)dipoles.size(); idx++) {

		//rectangle of dipole mesh
		Rect dipoleRect = dipoles[idx]->GetMeshRect();

		//dipole rectangular prism dimensions
		DBL3 abc = dipoleRect.size();

		//centre position of dipole mesh
		DBL3 dipole_centre = (dipoleRect.get_s() + dipoleRect.get_e()) / 2;

		//calculate stray field contribution from current dipole for all super-mesh cells
		for (int k = 0; k < strayField.size().k; k++) {
#pragma omp parallel for
			for (int j = 0; j < strayField.size().j; j++) {
				for (int i = 0; i < strayField.size().i; i++) {

					int cell_idx = i + j * int(strayField.size().i) + k * int(strayField.size().i * strayField.size().j);

					//distance from dipole_centre to current supermesh-cell
					DBL3 xyz = dipole_centre - DBL3((i + 0.5)*pSMesh->h_fm.x + sMeshOrigin.x, (j + 0.5)*pSMesh->h_fm.y + sMeshOrigin.y, (k + 0.5)*pSMesh->h_fm.z + sMeshOrigin.z);

					double p11, p22, p33, p12, p13, p23;

					p11 = DipoleTFunc::Nxx(xyz, abc);
					p22 = DipoleTFunc::Nyy(xyz, abc);
					p33 = DipoleTFunc::Nzz(xyz, abc);
					p12 = DipoleTFunc::Nxy(xyz, abc);
					p13 = DipoleTFunc::Nxz(xyz, abc);
					p23 = DipoleTFunc::Nyz(xyz, abc);

					//the Mdipole for this mesh (already correctly set for the mesh temperature)
					DBL3 Mdipole = dipoles[idx]->M[0];

					strayField[cell_idx].x -= p11 * Mdipole.x + p12 * Mdipole.y + p13 * Mdipole.z;
					strayField[cell_idx].y -= p12 * Mdipole.x + p22 * Mdipole.y + p23 * Mdipole.z;
					strayField[cell_idx].z -= p13 * Mdipole.x + p23 * Mdipole.y + p33 * Mdipole.z;
				}
			}
		}
	}
}

bool StrayField::CheckRecalculateStrayField(void)
{
	bool recalculateStrayField = false;

	for (int idx = 0; idx < (int)dipoles.size(); idx++) {

		//get recalculateStrayField flag for this dipole mesh then reset it
		recalculateStrayField |= dipoles[idx]->Check_recalculateStrayField();
		dipoles[idx]->Reset_recalculateStrayField();
	}

	return recalculateStrayField;
}

double StrayField::UpdateField(void)
{
	//recalculate stray field if needed (required when a dipole mesh changes, as indicated by its status flag)
	if(CheckRecalculateStrayField()) CalculateStrayField();

	//transfer stray field values to effective fields in the ferromagnetic meshes
	strayField.transfer_out();

	//not counting this to the total energy density for now
	return 0.0;
}

#endif