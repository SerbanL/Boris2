#include "stdafx.h"
#include "Mesh_AntiFerromagnetic.h"

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

//----------------------------------- ODE CONTROL METHODS IN (ANTI)FERROMAGNETIC MESH : Mesh_..._ODEControl.cpp

//get rate of change of magnetization (overloaded by Ferromagnetic meshes)
DBL3 AFMesh::dMdt(int idx)
{
	return meshODE.dMdt(idx);
}

//return average dm/dt in the given avRect (relative rect). Here m is the direction vector.
DBL3 AFMesh::Average_dmdt(Rect avRect)
{
	if (avRect.IsNull()) avRect = meshRect;
	Box box = M.box_from_rect_max(avRect);

#if COMPILECUDA == 1
	if (pMeshCUDA) return reinterpret_cast<AFMeshCUDA*>(pMeshCUDA)->Average_dmdt(box);
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

DBL3 AFMesh::Average_dmdt2(Rect avRect)
{
	if (avRect.IsNull()) avRect = meshRect;
	Box box = M.box_from_rect_max(avRect);

#if COMPILECUDA == 1
	if (pMeshCUDA) return reinterpret_cast<AFMeshCUDA*>(pMeshCUDA)->Average_dmdt2(box);
#endif

	OmpReduction<DBL3> reduction;
	reduction.new_average_reduction();

	for (int k = box.s.k; k < box.e.k; k++) {
#pragma omp parallel for
		for (int j = box.s.j; j < box.e.j; j++) {
			for (int i = box.s.i; i < box.e.i; i++) {

				int idx = i + j * n.x + k * n.x*n.y;

				if (M2.is_not_empty(idx)) {

					reduction.reduce_average(meshODE.dMdt2(idx) / M2[idx].norm());
				}
			}
		}
	}

	return reduction.average();
}

//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
DBL3 AFMesh::Average_mxdmdt(Rect avRect)
{
	if (avRect.IsNull()) avRect = meshRect;
	Box box = M.box_from_rect_max(avRect);

#if COMPILECUDA == 1
	if (pMeshCUDA) return reinterpret_cast<AFMeshCUDA*>(pMeshCUDA)->Average_mxdmdt(box);
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

DBL3 AFMesh::Average_mxdmdt2(Rect avRect)
{
	if (avRect.IsNull()) avRect = meshRect;
	Box box = M.box_from_rect_max(avRect);

#if COMPILECUDA == 1
	if (pMeshCUDA) return reinterpret_cast<AFMeshCUDA*>(pMeshCUDA)->Average_mxdmdt2(box);
#endif

	OmpReduction<DBL3> reduction;
	reduction.new_average_reduction();

	for (int k = box.s.k; k < box.e.k; k++) {
#pragma omp parallel for
		for (int j = box.s.j; j < box.e.j; j++) {
			for (int i = box.s.i; i < box.e.i; i++) {

				int idx = i + j * n.x + k * n.x*n.y;

				if (M2.is_not_empty(idx)) {

					double norm = M2[idx].norm();
					reduction.reduce_average((M2[idx] / norm) ^ (meshODE.dMdt2(idx) / norm));
				}
			}
		}
	}

	return reduction.average();
}

#endif