#include "stdafx.h"
#include "Atom_Mesh.h"

//----------------------------------- VALUE GETTERS

DBL3 Atom_Mesh::GetAverageChargeCurrentDensity(Rect rectangle)
{
	//TO DO
	return DBL3();
	//return CallModuleMethod(&Transport::GetAverageChargeCurrent, rectangle);
}

DBL3 Atom_Mesh::GetAverageSpinCurrentX(Rect rectangle)
{
	//TO DO
	return DBL3();
	//return CallModuleMethod(&Transport::GetAverageSpinCurrent, 0, rectangle);
}

DBL3 Atom_Mesh::GetAverageSpinCurrentY(Rect rectangle)
{
	//TO DO
	return DBL3();
	//return CallModuleMethod(&Transport::GetAverageSpinCurrent, 1, rectangle);
}

DBL3 Atom_Mesh::GetAverageSpinCurrentZ(Rect rectangle)
{
	//TO DO
	return DBL3();
	//return CallModuleMethod(&Transport::GetAverageSpinCurrent, 2, rectangle);
}

double Atom_Mesh::GetAverageElectricalPotential(Rect rectangle)
{
#if COMPILECUDA == 1
	if (paMeshCUDA) {

		if (V.linear_size()) return paMeshCUDA->V()->average_nonempty(n_e.dim(), rectangle);
		else return 0.0;
	}
#endif

	if (V.linear_size()) return V.average_nonempty_omp(rectangle);
	else return 0.0;
}

DBL3 Atom_Mesh::GetAverageSpinAccumulation(Rect rectangle)
{
#if COMPILECUDA == 1
	if (paMeshCUDA) {

		if (S.linear_size()) return paMeshCUDA->S()->average_nonempty(n_e.dim(), rectangle);
		else return DBL3(0.0);
	}
#endif

	if (S.linear_size()) return S.average_nonempty_omp(rectangle);
	else return DBL3(0.0);
}

double Atom_Mesh::GetAverageElectricalConductivity(Rect rectangle)
{
#if COMPILECUDA == 1
	if (paMeshCUDA) {

		if (elC.linear_size()) return paMeshCUDA->elC()->average_nonempty(n_e.dim(), rectangle);
		else return 0.0;
	}
#endif

	if (elC.linear_size()) return elC.average_nonempty_omp(rectangle);
	else return 0.0;
}

double Atom_Mesh::GetAverageTemperature(Rect rectangle)
{
#if COMPILECUDA == 1
	if (paMeshCUDA) {

		if (Temp.linear_size()) return paMeshCUDA->Temp()->average_nonempty(n_t.dim(), rectangle);
		else return base_temperature;
	}
#endif

	if (Temp.linear_size()) return Temp.average_nonempty_omp(rectangle);
	else return base_temperature;
}

double Atom_Mesh::GetAverageLatticeTemperature(Rect rectangle)
{
#if COMPILECUDA == 1
	if (paMeshCUDA) {

		if (Temp_l.linear_size()) return paMeshCUDA->Temp_l()->average_nonempty(n_t.dim(), rectangle);
		else return base_temperature;
	}
#endif

	if (Temp_l.linear_size()) return Temp_l.average_nonempty_omp(rectangle);
	else return base_temperature;
}