#pragma once

#include "cuVEC_MeshTransfer.h"
#include "launchers.h"

//------------------------------------------------------------------- MESH TRANSFER IN

template <typename VType>
__global__ void transfer_in_kernel(size_t transfer_info_size, VType*& sMesh_quantity, cuVEC<VType>*& mesh_in, cuPair<cuINT3, cuReal>*& transfer_info)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < transfer_info_size) {

		cuINT3 full_index = transfer_info[idx].first;
		
		//weight to apply to external mesh values
		cuReal weight = transfer_info[idx].second;
		
		//obtain weighted value from external mesh
		//first get the external mesh
		cuVEC<VType>& cuvec_in = mesh_in[full_index.i];

		//now get its weighted value
		VType weighted_value = cuvec_in[full_index.j] * weight;
		
		//stored contribution in supermesh
		//use atomic operation since there may be multiple contributions to this cell
		atomicAdd(&sMesh_quantity[full_index.k], weighted_value);
	}
}

template void cuTransfer<float>::transfer_in(size_t size_transfer, float*& sMesh_quantity);
template void cuTransfer<double>::transfer_in(size_t size_transfer, double*& sMesh_quantity);

template void cuTransfer<cuFLT3>::transfer_in(size_t size_transfer, cuFLT3*& sMesh_quantity);
template void cuTransfer<cuDBL3>::transfer_in(size_t size_transfer, cuDBL3*& sMesh_quantity);

template void cuTransfer<cuFLT33>::transfer_in(size_t size_transfer, cuFLT33*& sMesh_quantity);
template void cuTransfer<cuDBL33>::transfer_in(size_t size_transfer, cuDBL33*& sMesh_quantity);

//transfer values from mesh_in meshes using transfer_info into sMesh_quantity which has given size
template <typename VType>
__host__ void cuTransfer<VType>::transfer_in(size_t size_transfer, VType*& sMesh_quantity)
{
	transfer_in_kernel <<< (size_transfer + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (size_transfer, sMesh_quantity, mesh_in, transfer_in_info);
}

//------------------------------------------------------------------- MESH TRANSFER OUT

template <typename VType>
__global__ void zero_mesh_out_kernel(int mesh_out_num, cuVEC<VType>*& mesh_out)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	for (int idxMesh = 0; idxMesh < mesh_out_num; idxMesh++) {

		cuVEC<VType>& cuvec_out = mesh_out[idxMesh];

		if (idx < cuvec_out.linear_size()) {

			cuvec_out[idx] = VType();
		}
	}
}

template <typename VType>
__global__ void transfer_out_kernel(size_t transfer_info_size, VType*& sMesh_quantity, cuVEC<VType>*& mesh_out, cuPair<cuINT3, cuReal>*& transfer_info)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < transfer_info_size) {

		cuINT3 full_index = transfer_info[idx].first;

		//weight to apply to supermesh values
		cuReal weight = transfer_info[idx].second;

		//obtained weighted value from supermesh
		VType weighted_value = sMesh_quantity[full_index.k] * weight;
		
		//stored contribution in external mesh

		//first get the external mesh
		cuVEC<VType>& cuvec_out = mesh_out[full_index.i];

		//use atomic operation since there may be multiple contributions to this cell
		atomicAdd(&cuvec_out[full_index.j], weighted_value);
	}
}

template void cuTransfer<float>::transfer_out(size_t size_transfer, float*& sMesh_quantity, int mesh_out_num);
template void cuTransfer<double>::transfer_out(size_t size_transfer, double*& sMesh_quantity, int mesh_out_num);

template void cuTransfer<cuFLT3>::transfer_out(size_t size_transfer, cuFLT3*& sMesh_quantity, int mesh_out_num);
template void cuTransfer<cuDBL3>::transfer_out(size_t size_transfer, cuDBL3*& sMesh_quantity, int mesh_out_num);

template void cuTransfer<cuFLT33>::transfer_out(size_t size_transfer, cuFLT33*& sMesh_quantity, int mesh_out_num);
template void cuTransfer<cuDBL33>::transfer_out(size_t size_transfer, cuDBL33*& sMesh_quantity, int mesh_out_num);

template <typename VType>
__host__ void cuTransfer<VType>::transfer_out(size_t size_transfer, VType*& sMesh_quantity, int mesh_out_num)
{
	if (mesh_out_num) {

		//clear the output meshes before transferring supermesh values in.
		//size_transfer is greater or equal to the total number of cells in all output meshes, so we can launch a kernel with size_transfer
		//not the most elegant way of zeroing all meshes in mesh_out but speed difference between zeroing al meshes separately with exactly sized kernels is negligible and this way we don't have to store mesh_out sizes
		//also this is not used either so far as all output contributions are to Heff where we need to add into it!
		zero_mesh_out_kernel <<< (size_transfer + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (mesh_out_num, mesh_out);
	}

	transfer_out_kernel <<< (size_transfer + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (size_transfer, sMesh_quantity, mesh_out, transfer_out_info);
}