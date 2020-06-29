#include "stdafx.h"
#include "Atom_Heat.h"

#if defined(MODULE_COMPILATION_HEAT) && ATOMISTIC == 1

#include "Atom_Mesh.h"
#include "Atom_MeshParamsControl.h"

//set Temp uniformly to base temperature, unless a spatial variation is also specified through the cT mesh parameter
void Atom_Heat::SetBaseTemperature(double Temperature)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!paMesh->cT.is_sdep()) {

			dynamic_cast<Atom_HeatCUDA*>(pModuleCUDA)->SetBaseTemperature(Temperature);
		}
		else {

			dynamic_cast<Atom_HeatCUDA*>(pModuleCUDA)->SetBaseTemperature_Nonuniform(Temperature);
		}

		return;
	}
#endif

	if (!paMesh->cT.is_sdep()) {

		//uniform temperature
		paMesh->Temp.setnonempty(Temperature);

		if (paMesh->Temp_l.linear_size()) paMesh->Temp_l.setnonempty(Temperature);
	}
	else {

		//non-uniform temperature setting
#pragma omp parallel for
		for (int idx = 0; idx < paMesh->Temp.linear_size(); idx++) {

			if (paMesh->Temp.is_not_empty(idx)) {

				double cT = paMesh->cT;
				paMesh->update_parameters_tcoarse(idx, paMesh->cT, cT);

				paMesh->Temp[idx] = cT * Temperature;

				if (paMesh->Temp_l.linear_size()) {

					paMesh->Temp_l[idx] = cT * Temperature;
				}
			}
		}
	}
}

//-------------------Others

//called by MoveMesh method in this mesh - move relevant transport quantities
void Atom_Heat::MoveMesh_Heat(double x_shift)
{
	double mesh_end_size = paMesh->meshRect.size().x * MOVEMESH_ENDRATIO;

	Rect shift_rect = Rect(paMesh->meshRect.s + DBL3(mesh_end_size, 0, 0), paMesh->meshRect.e - DBL3(mesh_end_size, 0, 0));

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		paMesh->paMeshCUDA->Temp()->shift_x(paMesh->Temp.linear_size(), x_shift, shift_rect);

		if (paMesh->Temp_l.linear_size()) paMesh->paMeshCUDA->Temp_l()->shift_x(paMesh->Temp_l.linear_size(), x_shift, shift_rect);

		return;
	}
#endif

	paMesh->Temp.shift_x(x_shift, shift_rect);

	if (paMesh->Temp_l.linear_size()) paMesh->Temp_l.shift_x(x_shift, shift_rect);
}

#endif