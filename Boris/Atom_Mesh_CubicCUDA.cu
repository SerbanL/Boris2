#include "Atom_Mesh_CubicCUDA.h"

#if COMPILECUDA == 1

#ifdef MESH_COMPILATION_ATOM_CUBIC

#include "BorisCUDALib.cuh"
#include "ManagedAtom_DiffEqCubicCUDA.h"

__global__ void Average_dmdt_Atom_Cubic_kernel(cuBox box, ManagedAtom_MeshCUDA& cuaMesh, ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, cuReal3& average, size_t& points_count)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;

	int idxbox = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal3 average_ = cuReal3();
	bool include_in_average = false;

	if (idxbox < box.size().dim()) {

		//indexes of this threads in box
		int ibox = idxbox % box.size().i;
		int jbox = (idxbox / box.size().i) % box.size().j;
		int kbox = idxbox / (box.size().i * box.size().j);

		//indexes of box start in mesh
		int i = box.s.i % M1.n.i;
		int j = (box.s.j / M1.n.i) % M1.n.j;
		int k = box.s.k / (M1.n.i * M1.n.j);

		//total index in mesh
		int idx = i + ibox + (j + jbox) * M1.n.x + (k + kbox) * M1.n.x*M1.n.y;

		if (M1.is_not_empty(idx)) {

			average_ = cuaDiffEq.dMdt(idx) / M1[idx].norm();
			include_in_average = true;
		}
	}

	reduction_avg(0, 1, &average_, average, points_count, include_in_average);
}

__global__ void Average_mxdmdt_Atom_Cubic_kernel(cuBox box, ManagedAtom_MeshCUDA& cuaMesh, ManagedAtom_DiffEqCubicCUDA& cuaDiffEq, cuReal3& average, size_t& points_count)
{
	cuVEC_VC<cuReal3>& M1 = *cuaMesh.pM1;

	int idxbox = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal3 average_ = cuReal3();
	bool include_in_average = false;

	if (idxbox < box.size().dim()) {

		//indexes of this threads in box
		int ibox = idxbox % box.size().i;
		int jbox = (idxbox / box.size().i) % box.size().j;
		int kbox = idxbox / (box.size().i * box.size().j);

		//indexes of box start in mesh
		int i = box.s.i % M1.n.i;
		int j = (box.s.j / M1.n.i) % M1.n.j;
		int k = box.s.k / (M1.n.i * M1.n.j);

		//total index in mesh
		int idx = i + ibox + (j + jbox) * M1.n.x + (k + kbox) * M1.n.x*M1.n.y;

		if (M1.is_not_empty(idx)) {

			cuBReal norm = M1[idx].norm();
			average_ = (M1[idx] / norm) ^ (cuaDiffEq.dMdt(idx) / norm);
			include_in_average = true;
		}
	}

	reduction_avg(0, 1, &average_, average, points_count, include_in_average);
}

//----------------------------------- ODE METHODS IN (ANTI)FERROMAGNETIC MESH : Mesh_FerromagneticCUDA.cu

//return average dm/dt in the given avRect (relative rect). Here m is the direction vector.
DBL3 Atom_Mesh_CubicCUDA::Average_dmdt(cuBox avBox)
{
	Zero_aux_values();

	Average_dmdt_Atom_Cubic_kernel <<< (avBox.size().dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (avBox, cuaMesh, Get_ManagedAtom_DiffEqCUDA(), aux_real3, aux_int);

	int num_points = aux_int.to_cpu();

	if (num_points) return aux_real3.to_cpu() / num_points;
	else return DBL3();
}

//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
DBL3 Atom_Mesh_CubicCUDA::Average_mxdmdt(cuBox avBox)
{
	Zero_aux_values();

	Average_mxdmdt_Atom_Cubic_kernel <<< (avBox.size().dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (avBox, cuaMesh, Get_ManagedAtom_DiffEqCUDA(), aux_real3, aux_int);

	int num_points = aux_int.to_cpu();

	if (num_points) return aux_real3.to_cpu() / num_points;
	else return DBL3();
}

#endif

#endif
