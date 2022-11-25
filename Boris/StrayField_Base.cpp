#include "stdafx.h"
#include "StrayField_Base.h"

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "Mesh_Dipole.h"
#include "SuperMesh.h"
#include "DipoleTFunc.h"

StrayField_Base::StrayField_Base(SuperMesh *pSMesh_)
{
	pSMesh = pSMesh_;
}

//call from Initialize before starting computations, so all dipole meshes can be collected
void StrayField_Base::InitializeStrayField(void)
{
	dipoles.clear();

	//identify all existing dipole meshes
	for (int idx = 0; idx < pSMesh->size(); idx++) {

		if ((*pSMesh)[idx]->GetMeshType() == MESH_DIPOLE) {

			DipoleMesh *pDipole = dynamic_cast<DipoleMesh*>((*pSMesh)[idx]);

			//save pointer to dipole mesh
			dipoles.push_back(pDipole);
		}
	}

	strayField.set(DBL3());

	//here need to calculate stray field from currently set dipole "meshes" (if any - there should be at least one otherwise there's no point using this module)
	CalculateStrayField();
}

//calculate stray field from all dipoles
void StrayField_Base::CalculateStrayField(void)
{
	if (!dipoles.size()) strayField.set(DBL3(0));
	else {

		//origin of this mesh
		DBL3 meshOrigin = strayField.rect.get_s();
		//cellsize of this mesh
		DBL3 h = strayField.h;

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
						DBL3 xyz = dipole_centre - DBL3((i + 0.5)*h.x + meshOrigin.x, (j + 0.5)*h.y + meshOrigin.y, (k + 0.5)*h.z + meshOrigin.z);

						double p11, p22, p33, p12, p13, p23;

						p11 = DipoleTFunc::Nxx(xyz, abc);
						p22 = DipoleTFunc::Nyy(xyz, abc);
						p33 = DipoleTFunc::Nzz(xyz, abc);
						p12 = DipoleTFunc::Nxy(xyz, abc);
						p13 = DipoleTFunc::Nxz(xyz, abc);
						p23 = DipoleTFunc::Nyz(xyz, abc);

						//the Mdipole for this mesh (already correctly set for the mesh temperature)
						DBL3 Mdipole = dipoles[idx]->M[0];

						if (idx == 0) {

							strayField[cell_idx].x = -(p11 * Mdipole.x + p12 * Mdipole.y + p13 * Mdipole.z);
							strayField[cell_idx].y = -(p12 * Mdipole.x + p22 * Mdipole.y - p23 * Mdipole.z);
							strayField[cell_idx].z = -(p13 * Mdipole.x - p23 * Mdipole.y + p33 * Mdipole.z);
						}
						else {

							strayField[cell_idx].x += -(p11 * Mdipole.x + p12 * Mdipole.y + p13 * Mdipole.z);
							strayField[cell_idx].y += -(p12 * Mdipole.x + p22 * Mdipole.y - p23 * Mdipole.z);
							strayField[cell_idx].z += -(p13 * Mdipole.x - p23 * Mdipole.y + p33 * Mdipole.z);
						}
					}
				}
			}
		}
	}
}

//check if stray field needs to be recalculated
bool StrayField_Base::CheckRecalculateStrayField(void)
{
	for (int idx = 0; idx < (int)dipoles.size(); idx++) {

		if (dipoles[idx]->CheckRecalculateStrayField()) return true;
	}

	return false;
}

#endif