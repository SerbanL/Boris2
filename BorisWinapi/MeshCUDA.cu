#include "MeshCUDA.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void Zero_Aux_MeshCUDA(cuBReal& aux_real)
{
	if (threadIdx.x == 0) aux_real = 0.0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void GetTopologicalCharge_Kernel(ManagedMeshCUDA& cuMesh, cuRect rectangle, cuBReal& Q)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal Q_ = 0.0;

	cuReal3 pos;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			pos = M.cellidx_to_position(idx);

			cuBReal M_mag = M[idx].norm();

			cuReal33 M_grad = M.grad_neu(idx);

			cuReal3 dm_dx = M_grad.x / M_mag;
			cuReal3 dm_dy = M_grad.y / M_mag;

			Q_ = (M[idx] / M_mag) * (dm_dx ^ dm_dy) * M.h.x * M.h.y / (4 * (cuBReal)PI);
		}
	}

	reduction_sum(0, 1, &Q_, Q, rectangle.contains(pos));
}

//get topological charge using formula Q = Integral(m.(dm/dx x dm/dy) dxdy) / 4PI
cuBReal MeshCUDA::GetTopologicalCharge(cuRect rectangle)
{
	if (rectangle.IsNull()) rectangle = meshRect;

	Zero_Aux_MeshCUDA <<< 1, CUDATHREADS >>> (aux_real);

	GetTopologicalCharge_Kernel <<< (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuMesh, rectangle, aux_real);

	return aux_real.to_cpu();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

__global__ void Compute_TopoChargeDensity_Kernel(ManagedMeshCUDA& cuMesh, cuVEC<cuBReal>& aux_vec_sca)
{
	cuVEC_VC<cuReal3>& M = *cuMesh.pM;

	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < M.linear_size()) {

		if (M.is_not_empty(idx)) {

			cuBReal M_mag = M[idx].norm();

			cuReal33 M_grad = M.grad_neu(idx);

			cuReal3 dm_dx = M_grad.x / M_mag;
			cuReal3 dm_dy = M_grad.y / M_mag;

			aux_vec_sca[idx] = (M[idx] / M_mag) * (dm_dx ^ dm_dy) * M.h.x * M.h.y / (4 * (cuBReal)PI);
		}
		else aux_vec_sca[idx] = 0.0;
	}
}

//compute topological charge density spatial dependence and have it available in aux_vec_sca
//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
void MeshCUDA::Compute_TopoChargeDensity(void)
{
	aux_vec_sca()->resize(h, meshRect);

	Compute_TopoChargeDensity_Kernel <<< (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cuMesh, aux_vec_sca);
}

#endif