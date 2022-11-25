#pragma once

#include "mGPU_transfer.h"

template <typename VType>
__global__ void setup_device_memory_handle_managedvalue_kernel(VType& cu_source_value, VType** pcu_source_value)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx == 0) pcu_source_value[0] = &cu_source_value;
}

template cudaError_t mGPU_Transfer<int>::setup_device_memory_handle_managedvalue(int device_idx, int& cu_source_value);
template cudaError_t mGPU_Transfer<float>::setup_device_memory_handle_managedvalue(int device_idx, float& cu_source_value);
template cudaError_t mGPU_Transfer<double>::setup_device_memory_handle_managedvalue(int device_idx, double& cu_source_value);

template cudaError_t mGPU_Transfer<cuINT3>::setup_device_memory_handle_managedvalue(int device_idx, cuINT3& cu_source_value);
template cudaError_t mGPU_Transfer<cuFLT3>::setup_device_memory_handle_managedvalue(int device_idx, cuFLT3& cu_source_value);
template cudaError_t mGPU_Transfer<cuDBL3>::setup_device_memory_handle_managedvalue(int device_idx, cuDBL3& cu_source_value);

template cudaError_t mGPU_Transfer<cuINT4>::setup_device_memory_handle_managedvalue(int device_idx, cuINT4& cu_source_value);
template cudaError_t mGPU_Transfer<cuFLT4>::setup_device_memory_handle_managedvalue(int device_idx, cuFLT4& cu_source_value);
template cudaError_t mGPU_Transfer<cuDBL4>::setup_device_memory_handle_managedvalue(int device_idx, cuDBL4& cu_source_value);

//configure device memory locations to use
//Here cu_source_value is stored in gpu memory (i.e. it is a managed memory location, and will normally be passed through reference from a managed cu object).
//In this case create a handle to it here : i.e. a location in cpu memory which contains the gpu memory address value.
template <typename VType>
cudaError_t mGPU_Transfer<VType>::setup_device_memory_handle_managedvalue(int device_idx, VType& cu_source_value)
{
	//the goal here is to obtain gpu address of cu_source_value (which is a value in gpu memory, so cannot use reference operator in host code), and store it in gpu_memory_handles[device_idx]
	//thus must launch kernel to use reference operator, and then we can finally use cudaMemcpy

	mGPU.select_device(device_idx);
	
	cu_arr<VType*> pcu_source_value(1);
	setup_device_memory_handle_managedvalue_kernel<VType> <<< (1 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (cu_source_value, pcu_source_value);

	return cudaMemcpy(&gpu_memory_handles[device_idx], pcu_source_value.get_array(), sizeof(VType*), cudaMemcpyDeviceToHost);
}

//---------------------------------------------------------------------------------------------------------------------------

template cudaError_t mGPU_Transfer<int>::setup_extra_device_memory_handle_managedvalue(int device_idx, int handle_idx, int& cu_source_value);
template cudaError_t mGPU_Transfer<float>::setup_extra_device_memory_handle_managedvalue(int device_idx, int handle_idx, float& cu_source_value);
template cudaError_t mGPU_Transfer<double>::setup_extra_device_memory_handle_managedvalue(int device_idx, int handle_idx, double& cu_source_value);

template cudaError_t mGPU_Transfer<cuINT3>::setup_extra_device_memory_handle_managedvalue(int device_idx, int handle_idx, cuINT3& cu_source_value);
template cudaError_t mGPU_Transfer<cuFLT3>::setup_extra_device_memory_handle_managedvalue(int device_idx, int handle_idx, cuFLT3& cu_source_value);
template cudaError_t mGPU_Transfer<cuDBL3>::setup_extra_device_memory_handle_managedvalue(int device_idx, int handle_idx, cuDBL3& cu_source_value);

template cudaError_t mGPU_Transfer<cuINT4>::setup_extra_device_memory_handle_managedvalue(int device_idx, int handle_idx, cuINT4& cu_source_value);
template cudaError_t mGPU_Transfer<cuFLT4>::setup_extra_device_memory_handle_managedvalue(int device_idx, int handle_idx, cuFLT4& cu_source_value);
template cudaError_t mGPU_Transfer<cuDBL4>::setup_extra_device_memory_handle_managedvalue(int device_idx, int handle_idx, cuDBL4& cu_source_value);

//set value of existing extra handle for given device
//NOTE handle_idx ranges from 1 for extra handles, so subtract 1 to index gpu_extra_memory_handles
template <typename VType>
cudaError_t mGPU_Transfer<VType>::setup_extra_device_memory_handle_managedvalue(int device_idx, int handle_idx, VType& cu_source_value)
{
	//the goal here is to obtain gpu address of cu_source_value (which is a value in gpu memory, so cannot use reference operator in host code), and store it in gpu_memory_handles[device_idx]
	//thus must launch kernel to use reference operator, and then we can finally use cudaMemcpy

	if (handle_idx == 0) return setup_device_memory_handle_managedvalue(device_idx, cu_source_value);
	else {

		mGPU.select_device(device_idx);

		cu_arr<VType*> pcu_source_value(1);
		setup_device_memory_handle_managedvalue_kernel<VType> << < (1 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (cu_source_value, pcu_source_value);

		return cudaMemcpy(&gpu_extra_memory_handles[device_idx][handle_idx - 1], pcu_source_value.get_array(), sizeof(VType*), cudaMemcpyDeviceToHost);
	}
}

//---------------------------------------------------------------------------------------------------------------------------

template <typename VType>
__global__ void reduce_basegpu_memory_kernel(VType* basegpu_aux_memory, size_t num_gpus, VType& basegpu_value)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx == 0) {

		for (int didx = 0; didx < num_gpus; didx++) {

			basegpu_value += basegpu_aux_memory[didx];
		}
	}
}

template int mGPU_Transfer<int>::reduce(void);
template float mGPU_Transfer<float>::reduce(void);
template double mGPU_Transfer<double>::reduce(void);

template cuINT3 mGPU_Transfer<cuINT3>::reduce(void);
template cuFLT3 mGPU_Transfer<cuFLT3>::reduce(void);
template cuDBL3 mGPU_Transfer<cuDBL3>::reduce(void);

template cuINT4 mGPU_Transfer<cuINT4>::reduce(void);
template cuFLT4 mGPU_Transfer<cuFLT4>::reduce(void);
template cuDBL4 mGPU_Transfer<cuDBL4>::reduce(void);

//reduce values from all gpus to a single output value
//typically each gpu will do its own reduction, resulting in values on each gpu to which we have handles gpu_memory_handles
//here we just need a final reduction of these values
template <typename VType>
VType mGPU_Transfer<VType>::reduce(void)
{
	for (int idx = 0; idx < mGPU.get_num_devices(); idx++) {

		if (mGPU.is_p2p(0, idx)) {

			//P2P transfer
			mGPU.select_device(0); //not strictly necessary, but if not in the right context already, use of UVA will result in additional significant delay (according to benchmarks!)
			cudaMemcpy(basegpu_aux_memory + idx, gpu_memory_handles[idx], sizeof(VType), cudaMemcpyDeviceToDevice);
		}
		else {

			//Indirect transfer through host (much slower than P2P)
			mGPU.select_device(0);
			cudaMemcpy(cpu_memory, gpu_memory_handles[idx], sizeof(VType), cudaMemcpyDeviceToHost);
			mGPU.select_device(0);
			cudaMemcpy(basegpu_aux_memory + idx, cpu_memory, sizeof(VType), cudaMemcpyHostToDevice);
		}
	}

	mGPU.select_device(0);
	pbasegpu_aux_value->from_cpu(VType());
	reduce_basegpu_memory_kernel<VType> <<< (1 + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (basegpu_aux_memory, mGPU.get_num_devices(), *pbasegpu_aux_value);
	return pbasegpu_aux_value->to_cpu();
}