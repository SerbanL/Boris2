#include "stdafx.h"
#include "Mesh_Ferromagnetic.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC

//----------------------------------- ODE CONTROL METHODS IN (ANTI)FERROMAGNETIC MESH : Mesh_..._ODEControl.cpp

//get rate of change of magnetization (overloaded by Ferromagnetic meshes)
DBL3 FMesh::dMdt(int idx)
{
	return meshODE.dMdt(idx);
}

//Save current magnetization in sM VECs (e.g. useful to reset dM / dt calculation)
void FMesh::SaveMagnetization(void)
{
	meshODE.SaveMagnetization();
}

//return average dm/dt in the given avRect (relative rect). Here m is the direction vector.
DBL3 FMesh::Average_dmdt(Rect avRect)
{
	if (avRect.IsNull()) avRect = meshRect - meshRect.s;
	Box box = M.box_from_rect_max(avRect + meshRect.s);

#if COMPILECUDA == 1
	if (pMeshCUDA) return reinterpret_cast<FMeshCUDA*>(pMeshCUDA)->Average_dmdt(box);
#endif

	OmpReduction<DBL3> reduction;
	reduction.new_average_reduction();

	for (int k = box.s.k; k < box.e.k; k++) {
#pragma omp parallel for
		for (int j = box.s.j; j < box.e.j; j++) {
			for (int i = box.s.i; i < box.e.i; i++) {

				int idx = i + j * n.x + k * n.x*n.y;

				if (M.is_not_empty(idx)) {

					reduction.reduce_average(meshODE.dMdt(idx) / M[idx].norm());
				}
			}
		}
	}

	return reduction.average();
}

//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
DBL3 FMesh::Average_mxdmdt(Rect avRect)
{
	if (avRect.IsNull()) avRect = meshRect - meshRect.s;
	Box box = M.box_from_rect_max(avRect + meshRect.s);

#if COMPILECUDA == 1
	if (pMeshCUDA) return reinterpret_cast<FMeshCUDA*>(pMeshCUDA)->Average_mxdmdt(box);
#endif

	OmpReduction<DBL3> reduction;
	reduction.new_average_reduction();

	for (int k = box.s.k; k < box.e.k; k++) {
#pragma omp parallel for
		for (int j = box.s.j; j < box.e.j; j++) {
			for (int i = box.s.i; i < box.e.i; i++) {

				int idx = i + j * n.x + k * n.x*n.y;

				if (M.is_not_empty(idx)) {

					double norm = M[idx].norm();
					reduction.reduce_average((M[idx] / norm) ^ (meshODE.dMdt(idx) / norm));
				}
			}
		}
	}

	return reduction.average();
}

#endif