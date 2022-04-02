#include "stdafx.h"
#include "StrayField_BaseCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "SuperMesh.h"
#include "Mesh_DipoleCUDA.h"

#include "Mesh.h"

StrayField_BaseCUDA::StrayField_BaseCUDA(SuperMesh* pSMesh_)
{
	pSMesh = pSMesh_;

}

StrayField_BaseCUDA::~StrayField_BaseCUDA()
{
}

void StrayField_BaseCUDA::InitializeStrayField(void)
{
	dipoles.clear();
	managed_dipoles.clear();

	//identify all existing dipole meshes and output Heff meshes
	for (int idx = 0; idx < (int)pSMesh->size(); idx++) {

		//collect cuda pointers
		if ((*pSMesh)[idx]->GetMeshType() == MESH_DIPOLE) {

			DipoleMeshCUDA *pDipoleCUDA = dynamic_cast<DipoleMeshCUDA*>(dynamic_cast<Mesh*>((*pSMesh)[idx])->pMeshCUDA);
			dipoles.push_back(pDipoleCUDA);
			managed_dipoles.push_back(pDipoleCUDA->cuMesh.get_managed_object());
		}
	}

	strayField()->set(cuReal3(0.0));

	//here need to calculate stray field from currently set dipole "meshes" (if any - there should be at least one otherwise there's no point using this module)
	CalculateStrayField();
}

//check if stray field needs to be recalculated
bool StrayField_BaseCUDA::CheckRecalculateStrayField(void)
{
	for (int idx = 0; idx < (int)dipoles.size(); idx++) {

		if (dipoles[idx]->CheckRecalculateStrayField()) return true;
	}

	return false;
}

#endif

#endif