#include "stdafx.h"
#include "Atom_Mesh_Cubic.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC

//----------------------------------- ODE METHODS : Atom_Mesh_Cubic_ODEControl.cpp

//get rate of change of magnetic moment (overloaded by Ferromagnetic meshes)
DBL3 Atom_Mesh_Cubic::dMdt(int idx)
{
	return meshODE.dMdt(idx);
}

//return average dm/dt in the given avRect (relative rect). Here m is the direction vector.
DBL3 Atom_Mesh_Cubic::Average_dmdt(Rect avRect)
{
	if (avRect.IsNull()) avRect = meshRect - meshRect.s;
	Box box = M1.box_from_rect_max(avRect + meshRect.s);

#if COMPILECUDA == 1
	if (paMeshCUDA) return reinterpret_cast<Atom_Mesh_CubicCUDA*>(paMeshCUDA)->Average_dmdt(box);
#endif

	OmpReduction<DBL3> reduction;
	reduction.new_average_reduction();

	for (int k = box.s.k; k < box.e.k; k++) {
#pragma omp parallel for
		for (int j = box.s.j; j < box.e.j; j++) {
			for (int i = box.s.i; i < box.e.i; i++) {

				int idx = i + j * n.x + k * n.x*n.y;

				if (M1.is_not_empty(idx)) {

					reduction.reduce_average(meshODE.dMdt(idx) / M1[idx].norm());
				}
			}
		}
	}

	return reduction.average();
}

//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
DBL3 Atom_Mesh_Cubic::Average_mxdmdt(Rect avRect)
{
	if (avRect.IsNull()) avRect = meshRect - meshRect.s;
	Box box = M1.box_from_rect_max(avRect + meshRect.s);

#if COMPILECUDA == 1
	if (paMeshCUDA) return reinterpret_cast<Atom_Mesh_CubicCUDA*>(paMeshCUDA)->Average_mxdmdt(box);
#endif

	OmpReduction<DBL3> reduction;
	reduction.new_average_reduction();

	for (int k = box.s.k; k < box.e.k; k++) {
#pragma omp parallel for
		for (int j = box.s.j; j < box.e.j; j++) {
			for (int i = box.s.i; i < box.e.i; i++) {

				int idx = i + j * n.x + k * n.x*n.y;

				if (M1.is_not_empty(idx)) {

					double norm = M1[idx].norm();
					reduction.reduce_average((M1[idx] / norm) ^ (meshODE.dMdt(idx) / norm));
				}
			}
		}
	}

	return reduction.average();
}

#endif