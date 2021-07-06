#include "stdafx.h"
#include "Mesh.h"

#include "Transport.h"

//----------------------------------- VALUE GETTERS

//get average magnetization in given rectangle (entire mesh if none specified)
DBL3 Mesh::GetAverageMagnetization(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (M.linear_size()) return pMeshCUDA->M()->average_nonempty(n.dim(), rectangle);
		else return DBL3(0.0);
	}
#endif

	if (M.linear_size()) return M.average_nonempty_omp(rectangle);
	else return DBL3(0.0);
}

//get average magnetization in given rectangle (entire mesh if none specified); sub-lattice B
DBL3 Mesh::GetAveragemagnetization2(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (M2.linear_size()) return pMeshCUDA->M2()->average_nonempty(n.dim(), rectangle);
		else return DBL3(0.0);
	}
#endif

	if (M2.linear_size()) return M2.average_nonempty_omp(rectangle);
	else return DBL3(0.0);
}

//get magnetization magnitude min-max in given rectangle (entire mesh if none specified)
DBL2 Mesh::GetmagnetizationMinMax(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (M.linear_size()) return pMeshCUDA->M()->get_minmax(n.dim(), rectangle);
		else return DBL2(0.0);
	}
#endif

	if (M.linear_size()) return M.get_minmax(rectangle);
	else return DBL2(0.0);
}

//get magnetization component min-max in given rectangle (entire mesh if none specified)
DBL2 Mesh::GetmagnetizationXMinMax(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (M.linear_size()) return pMeshCUDA->M()->get_minmax_component_x(n.dim(), rectangle);
		else return DBL2(0.0);
	}
#endif

	if (M.linear_size()) return M.get_minmax_component_x(rectangle);
	else return DBL2(0.0);
}

//get magnetization component min-max in given rectangle (entire mesh if none specified)
DBL2 Mesh::GetmagnetizationYMinMax(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (M.linear_size()) return pMeshCUDA->M()->get_minmax_component_y(n.dim(), rectangle);
		else return DBL2(0.0);
	}
#endif

	if (M.linear_size()) return M.get_minmax_component_y(rectangle);
	else return DBL2(0.0);
}

//get magnetization component min-max in given rectangle (entire mesh if none specified)
DBL2 Mesh::GetmagnetizationZMinMax(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (M.linear_size()) return pMeshCUDA->M()->get_minmax_component_z(n.dim(), rectangle);
		else return DBL2(0.0);
	}
#endif

	if (M.linear_size()) return M.get_minmax_component_z(rectangle);
	else return DBL2(0.0);
}

//get topological charge using formula Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
double Mesh::GetTopologicalCharge(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) return pMeshCUDA->GetTopologicalCharge(rectangle);
#endif

	if (rectangle.IsNull()) rectangle = meshRect;

	if (M.linear_size()) {

		double Q = 0.0;

#pragma omp parallel for reduction(+:Q)
		for (int idx = 0; idx < M.linear_size(); idx++) {

			if (M.is_not_empty(idx)) {

				DBL3 pos = M.cellidx_to_position(idx);

				if (!rectangle.contains(pos)) continue;

				double M_mag = M[idx].norm();

				DBL33 M_grad = M.grad_neu(idx);

				DBL3 dm_dx = M_grad.x / M_mag;
				DBL3 dm_dy = M_grad.y / M_mag;

				Q += (M[idx] / M_mag) * (dm_dx ^ dm_dy) * M.h.x * M.h.y;
			}
		}

		return Q / (4 * PI * M.n.z);
	}
	else return 0.0;
}

double Mesh::GetAverageElectricalPotential(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (V.linear_size()) return pMeshCUDA->V()->average_nonempty(n_e.dim(), rectangle);
		else return 0.0;
	}
#endif

	if (V.linear_size()) return V.average_nonempty_omp(rectangle);
	else return 0.0;
}

DBL3 Mesh::GetAverageSpinAccumulation(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (S.linear_size()) return pMeshCUDA->S()->average_nonempty(n_e.dim(), rectangle);
		else return DBL3(0.0);
	}
#endif

	if (S.linear_size()) return S.average_nonempty_omp(rectangle);
	else return DBL3(0.0);
}

double Mesh::GetAverageElectricalConductivity(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (elC.linear_size()) return pMeshCUDA->elC()->average_nonempty(n_e.dim(), rectangle);
		else return 0.0;
	}
#endif

	if (elC.linear_size()) return elC.average_nonempty_omp(rectangle);
	else return 0.0;
}

double Mesh::GetAverageTemperature(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (Temp.linear_size()) return pMeshCUDA->Temp()->average_nonempty(n_t.dim(), rectangle);
		else return base_temperature;
	}
#endif

	if (Temp.linear_size()) return Temp.average_nonempty_omp(rectangle);
	else return base_temperature;
}

double Mesh::GetAverageLatticeTemperature(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		if (Temp_l.linear_size()) return pMeshCUDA->Temp_l()->average_nonempty(n_t.dim(), rectangle);
		else return base_temperature;
	}
#endif

	if (Temp_l.linear_size()) return Temp_l.average_nonempty_omp(rectangle);
	else return base_temperature;
}