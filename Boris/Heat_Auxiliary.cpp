#include "stdafx.h"
#include "Heat.h"

#ifdef MODULE_COMPILATION_HEAT

#include "Mesh.h"
#include "MeshParamsControl.h"

//set Temp uniformly to base temperature, unless a spatial variation is also specified through the cT mesh parameter
void Heat::SetBaseTemperature(double Temperature)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!pMesh->cT.is_sdep()) {

			dynamic_cast<HeatCUDA*>(pModuleCUDA)->SetBaseTemperature(Temperature);
		}
		else {

			dynamic_cast<HeatCUDA*>(pModuleCUDA)->SetBaseTemperature_Nonuniform(Temperature);
		}

		return;
	}
#endif

	if (!pMesh->cT.is_sdep()) {

		//uniform temperature
		pMesh->Temp.setnonempty(Temperature);

		if (pMesh->Temp_l.linear_size()) pMesh->Temp_l.setnonempty(Temperature);
	}
	else {

		//non-uniform temperature setting
#pragma omp parallel for
		for (int idx = 0; idx < pMesh->Temp.linear_size(); idx++) {

			if (pMesh->Temp.is_not_empty(idx)) {

				double cT = pMesh->cT;
				pMesh->update_parameters_tcoarse(idx, pMesh->cT, cT);

				pMesh->Temp[idx] = cT * Temperature;

				if (pMesh->Temp_l.linear_size()) {

					pMesh->Temp_l[idx] = cT * Temperature;
				}
			}
		}
	}
}

//-------------------Others

//called by MoveMesh method in this mesh - move relevant transport quantities
void Heat::MoveMesh_Heat(double x_shift)
{
	double mesh_end_size = pMesh->meshRect.size().x * MOVEMESH_ENDRATIO;

	Rect shift_rect = Rect(pMesh->meshRect.s + DBL3(mesh_end_size, 0, 0), pMesh->meshRect.e - DBL3(mesh_end_size, 0, 0));

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pMesh->pMeshCUDA->Temp()->shift_x(pMesh->Temp.linear_size(), x_shift, shift_rect);

		if (pMesh->Temp_l.linear_size()) pMesh->pMeshCUDA->Temp_l()->shift_x(pMesh->Temp_l.linear_size(), x_shift, shift_rect);

		return;
	}
#endif

	pMesh->Temp.shift_x(x_shift, shift_rect);

	if (pMesh->Temp_l.linear_size()) pMesh->Temp_l.shift_x(x_shift, shift_rect);
}

#endif