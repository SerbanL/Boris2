#include "Mesh_AntiFerromagneticCUDA.h"

#if COMPILECUDA == 1

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

#include "BorisCUDALib.cuh"
#include "ManagedDiffEqAFMCUDA.h"

__global__ void Average_dmdt_AFM_kernel(cuBox box, ManagedMeshCUDA& cuMesh, ManagedDiffEqAFMCUDA& cuDiffEq, cuReal3& average, size_t& points_count)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;

	int idxbox = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal3 average_ = cuReal3();
	bool include_in_average = false;

	if (idxbox < box.size().dim()) {

		//indexes of this threads in box
		int ibox = idxbox % box.size().i;
		int jbox = (idxbox / box.size().i) % box.size().j;
		int kbox = idxbox / (box.size().i * box.size().j);

		//indexes of box start in mesh
		int i = box.s.i % M.n.i;
		int j = (box.s.j / M.n.i) % M.n.j;
		int k = box.s.k / (M.n.i * M.n.j);

		//total index in mesh
		int idx = i + ibox + (j + jbox) * M.n.x + (k + kbox) * M.n.x*M.n.y;

		if (M.is_not_empty(idx)) {

			average_ = cuDiffEq.dMdt(idx) / M[idx].norm();
			include_in_average = true;
		}
	}

	reduction_avg(0, 1, &average_, average, points_count, include_in_average);
}

__global__ void Average_dmdt2_AFM_kernel(cuBox box, ManagedMeshCUDA& cuMesh, ManagedDiffEqAFMCUDA& cuDiffEq, cuReal3& average, size_t& points_count)
{
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;

	int idxbox = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal3 average_ = cuReal3();
	bool include_in_average = false;

	if (idxbox < box.size().dim()) {

		//indexes of this threads in box
		int ibox = idxbox % box.size().i;
		int jbox = (idxbox / box.size().i) % box.size().j;
		int kbox = idxbox / (box.size().i * box.size().j);

		//indexes of box start in mesh
		int i = box.s.i % M2.n.i;
		int j = (box.s.j / M2.n.i) % M2.n.j;
		int k = box.s.k / (M2.n.i * M2.n.j);

		//total index in mesh
		int idx = i + ibox + (j + jbox) * M2.n.x + (k + kbox) * M2.n.x*M2.n.y;

		if (M2.is_not_empty(idx)) {

			average_ = cuDiffEq.dMdt2(idx) / M2[idx].norm();
			include_in_average = true;
		}
	}

	reduction_avg(0, 1, &average_, average, points_count, include_in_average);
}

__global__ void Average_mxdmdt_AFM_kernel(cuBox box, ManagedMeshCUDA& cuMesh, ManagedDiffEqAFMCUDA& cuDiffEq, cuReal3& average, size_t& points_count)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;

	int idxbox = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal3 average_ = cuReal3();
	bool include_in_average = false;

	if (idxbox < box.size().dim()) {

		//indexes of this threads in box
		int ibox = idxbox % box.size().i;
		int jbox = (idxbox / box.size().i) % box.size().j;
		int kbox = idxbox / (box.size().i * box.size().j);

		//indexes of box start in mesh
		int i = box.s.i % M.n.i;
		int j = (box.s.j / M.n.i) % M.n.j;
		int k = box.s.k / (M.n.i * M.n.j);

		//total index in mesh
		int idx = i + ibox + (j + jbox) * M.n.x + (k + kbox) * M.n.x*M.n.y;

		if (M.is_not_empty(idx)) {

			cuBReal norm = M[idx].norm();
			average_ = (M[idx] / norm) ^ (cuDiffEq.dMdt(idx) / norm);
			include_in_average = true;
		}
	}

	reduction_avg(0, 1, &average_, average, points_count, include_in_average);
}

__global__ void Average_mxdmdt2_AFM_kernel(cuBox box, ManagedMeshCUDA& cuMesh, ManagedDiffEqAFMCUDA& cuDiffEq, cuReal3& average, size_t& points_count)
{
	cuVEC_VC<cuReal3>& M2 = *cuMesh.pM2;

	int idxbox = blockIdx.x * blockDim.x + threadIdx.x;

	cuReal3 average_ = cuReal3();
	bool include_in_average = false;

	if (idxbox < box.size().dim()) {

		//indexes of this threads in box
		int ibox = idxbox % box.size().i;
		int jbox = (idxbox / box.size().i) % box.size().j;
		int kbox = idxbox / (box.size().i * box.size().j);

		//indexes of box start in mesh
		int i = box.s.i % M2.n.i;
		int j = (box.s.j / M2.n.i) % M2.n.j;
		int k = box.s.k / (M2.n.i * M2.n.j);

		//total index in mesh
		int idx = i + ibox + (j + jbox) * M2.n.x + (k + kbox) * M2.n.x*M2.n.y;

		if (M2.is_not_empty(idx)) {

			cuBReal norm = M2[idx].norm();
			average_ = (M2[idx] / norm) ^ (cuDiffEq.dMdt2(idx) / norm);
			include_in_average = true;
		}
	}

	reduction_avg(0, 1, &average_, average, points_count, include_in_average);
}

//----------------------------------- ODE METHODS IN (ANTI)FERROMAGNETIC MESH : Mesh_AntiFerromagneticCUDA.cu

//return average dm/dt in the given avRect (relative rect). Here m is the direction vector.
DBL3 AFMeshCUDA::Average_dmdt(cuBox avBox)
{
	Zero_aux_values();

	Average_dmdt_AFM_kernel <<< (avBox.size().dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (avBox, cuMesh, Get_ManagedDiffEqCUDA(), aux_real3, aux_int);

	int num_points = aux_int.to_cpu();

	if (num_points) return aux_real3.to_cpu() / num_points;
	else return DBL3();
}

DBL3 AFMeshCUDA::Average_dmdt2(cuBox avBox)
{
	Zero_aux_values();

	Average_dmdt2_AFM_kernel <<< (avBox.size().dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (avBox, cuMesh, Get_ManagedDiffEqCUDA(), aux_real3, aux_int);

	int num_points = aux_int.to_cpu();

	if (num_points) return aux_real3.to_cpu() / num_points;
	else return DBL3();
}

//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
DBL3 AFMeshCUDA::Average_mxdmdt(cuBox avBox)
{
	Zero_aux_values();

	Average_mxdmdt_AFM_kernel <<< (avBox.size().dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (avBox, cuMesh, Get_ManagedDiffEqCUDA(), aux_real3, aux_int);

	int num_points = aux_int.to_cpu();

	if (num_points) return aux_real3.to_cpu() / num_points;
	else return DBL3();
}

DBL3 AFMeshCUDA::Average_mxdmdt2(cuBox avBox)
{
	Zero_aux_values();

	Average_mxdmdt2_AFM_kernel <<< (avBox.size().dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (avBox, cuMesh, Get_ManagedDiffEqCUDA(), aux_real3, aux_int);

	int num_points = aux_int.to_cpu();

	if (num_points) return aux_real3.to_cpu() / num_points;
	else return DBL3();
}

#endif

#endif
