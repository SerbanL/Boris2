#pragma once

#include "mGPU.h"

//Auxiliary class for performing memory transfers between devices, for configured linear memory spaces
//Uses mGPUConfig for configured memory access type.
template <typename VType>
class mGPU_Transfer
{

private:

	//use a pre-configured mGPUConfig object
	mGPUConfig& mGPU;

	//---- GPU Memory handles

	//BASE

	//gpu array memory handles used for transfers : gpu_memory_handles[idx] is a pointer stored in cpu memory, containing address value to gpu memory
	VType** gpu_memory_handles = nullptr;

	//size of transfer : same size for all devices, and user must ensure this does not exceed allocated memory on any of the devices
	size_t size_transfer = 0;

	//EXTRA

	//extra memory handles on each device : outer index is device (size is number of devices); inner index is number of handles allocated (initially empty)
	std::vector<std::vector<VType*>> gpu_extra_memory_handles;

	//---- Reduction

	//pointer in cpu memory to array on base gpu (base gpu is the first physical device)
	//these are used to implement reduction across multiple gpus
	VType* basegpu_aux_memory = nullptr;
	//base gpu single value for reduction of values in basegpu_aux_memory
	cu_obj<VType>* pbasegpu_aux_value = nullptr;

	//---- Indirect memory transfer

	//if indirect memory transfer is used, then transfer to this space first
	VType* cpu_memory = nullptr;

public:

	//------------------------------------------- CONSTRUCTOR

	mGPU_Transfer(mGPUConfig& mGPU_, size_t size_transfer_ = 0) :
		mGPU(mGPU_)
	{
		//BASE memory handles
		gpu_memory_handles = new VType*[mGPU.get_num_devices()];
		for (int idx = 0; idx < mGPU.get_num_devices(); idx++) {

			gpu_memory_handles[idx] = nullptr;
		}

		//EXTRA memory handles
		gpu_extra_memory_handles.resize(mGPU.get_num_devices());

		//transfer size and CPU buffer helper if P2P not available
		size_transfer = size_transfer_;
		if (!mGPU.is_all_p2p() && size_transfer) cpu_memory = new VType[size_transfer];

		//Auxiliary
		mGPU.select_device(0);
		gpu_alloc(basegpu_aux_memory, mGPU.get_num_devices());
		pbasegpu_aux_value = new cu_obj<VType>();
	}

	~mGPU_Transfer()
	{
		if (gpu_memory_handles) delete[] gpu_memory_handles;
		gpu_memory_handles = nullptr;

		if (cpu_memory) delete[] cpu_memory;
		cpu_memory = nullptr;

		mGPU.select_device(0);
		gpu_free(basegpu_aux_memory);

		delete pbasegpu_aux_value;
		pbasegpu_aux_value = nullptr;
	}

	//------------------------------------------- TRANSFER SIZE

	void set_transfer_size(size_t size_transfer_)
	{
		if (!mGPU.is_all_p2p() && size_transfer_ != size_transfer) {

			if (cpu_memory) delete[] cpu_memory;
			cpu_memory = nullptr;
			cpu_memory = new VType[size_transfer];
		}

		size_transfer = size_transfer_;
	}

	size_t get_transfer_size(void) const { return size_transfer; }

	//------------------------------------------- SETUP TRANSFERS

	//BASE

	//configure device memory locations to use
	//Here cu_source_pointer is stored in cpu memory, but contains an address to gpu memory : i.e. it is a handle already, so just copy it here.
	void setup_device_memory_handle(int device_idx, VType*& cu_source_pointer)
	{
		gpu_memory_handles[device_idx] = cu_source_pointer;
	}

	//configure device memory locations to use
	//Here cu_source_pointer is stored in gpu memory (i.e. it is a managed memory location, and will normally be passed through reference from a managed cu object).
	//In this case create a handle to it here : i.e. a location in cpu memory which contains the gpu memory address value.
	cudaError_t setup_device_memory_handle_managed(int device_idx, VType*& cu_source_pointer)
	{
		cudaError_t error = cudaMemcpy(&gpu_memory_handles[device_idx], &cu_source_pointer, sizeof(VType*), cudaMemcpyDeviceToHost);
		return error;
	}

	//configure device memory locations to use
	//Here cu_source_value is stored in gpu memory (i.e. it is a managed memory location, and will normally be passed through reference from a managed cu object).
	//In this case create a handle to it here : i.e. a location in cpu memory which contains the gpu memory address value.
	cudaError_t setup_device_memory_handle_managedvalue(int device_idx, VType& cu_source_value);

	//EXTRA

	//---Add new

	//set value of existing extra handle for given device
	void add_extra_device_memory_handle(int device_idx, VType*& cu_source_pointer)
	{
		gpu_extra_memory_handles[device_idx].push_back(nullptr);
		return setup_extra_device_memory_handle(device_idx, gpu_extra_memory_handles[device_idx].size(), cu_source_pointer);
	}

	//set value of existing extra handle for given device
	cudaError_t add_extra_device_memory_handle_managed(int device_idx, VType*& cu_source_pointer)
	{
		gpu_extra_memory_handles[device_idx].push_back(nullptr);
		return setup_extra_device_memory_handle_managed(device_idx, gpu_extra_memory_handles[device_idx].size(), cu_source_pointer);
	}

	//set value of existing extra handle for given device
	cudaError_t add_extra_device_memory_handle_managedvalue(int device_idx, VType& cu_source_value)
	{
		gpu_extra_memory_handles[device_idx].push_back(nullptr);
		return setup_extra_device_memory_handle_managedvalue(device_idx, gpu_extra_memory_handles[device_idx].size(), cu_source_value);
	}

	//---Set existing

	//set value of existing extra handle for given device. Can also change base handle (index with 0).
	//NOTE handle_idx ranges from 1 for extra handles, so subtract 1 to index gpu_extra_memory_handles
	void setup_extra_device_memory_handle(int device_idx, int handle_idx, VType*& cu_source_pointer)
	{
		if (handle_idx == 0) setup_device_memory_handle(device_idx, cu_source_pointer);
		else gpu_extra_memory_handles[device_idx][handle_idx - 1] = cu_source_pointer;
	}

	//set value of existing extra handle for given device. Can also change base handle (index with 0).
	//NOTE handle_idx ranges from 1 for extra handles, so subtract 1 to index gpu_extra_memory_handles
	cudaError_t setup_extra_device_memory_handle_managed(int device_idx, int handle_idx, VType*& cu_source_pointer)
	{
		if (handle_idx == 0) return setup_device_memory_handle_managed(device_idx, cu_source_pointer);
		else return cudaMemcpy(&gpu_extra_memory_handles[device_idx][handle_idx - 1], &cu_source_pointer, sizeof(VType*), cudaMemcpyDeviceToHost);
	}

	//set value of existing extra handle for given device. Can also change base handle (index with 0).
	//NOTE handle_idx ranges from 1 for extra handles, so subtract 1 to index gpu_extra_memory_handles
	cudaError_t setup_extra_device_memory_handle_managedvalue(int device_idx, int handle_idx, VType& cu_source_value);

	//------------------------------------------- DO TRANSFER

	//BASE

	//do a transfer between 2 devices of size_transfer, unless num_elements specified. Base handles only.
	void transfer(int device_to, int device_from, size_t num_elements = 0)
	{
		size_t size = (num_elements == 0 ? size_transfer : num_elements);

		if (mGPU.is_p2p(device_to, device_from)) {

			//P2P transfer
			mGPU.select_device(device_to);	//not strictly necessary, but if not in the right context already, use of UVA will result in additional significant delay (according to benchmarks!)
			cudaError_t error = cudaMemcpy(gpu_memory_handles[device_to], gpu_memory_handles[device_from], sizeof(VType)*size, cudaMemcpyDeviceToDevice);
		}
		else {

			//Indirect transfer through host (much slower than P2P)
			mGPU.select_device(device_from);
			cudaMemcpy(cpu_memory, gpu_memory_handles[device_from], sizeof(VType)*size, cudaMemcpyDeviceToHost);
			mGPU.select_device(device_to);
			cudaMemcpy(gpu_memory_handles[device_to], cpu_memory, sizeof(VType)*size, cudaMemcpyHostToDevice);
		}
	}

	//EXTRA

	//general transfer method where handle_to and handle_from are indexes, with 0 BASE and 1, 2, 3, ... are extra handle indexes
	//NOTE handle_idx ranges from 1 for extra handles, so subtract 1 to index gpu_extra_memory_handles
	void transfer(int device_to, int handle_to, int device_from, int handle_from, size_t num_elements = 0)
	{
		size_t size = (num_elements == 0 ? size_transfer : num_elements);

		if (mGPU.is_p2p(device_to, device_from)) {

			//P2P transfer

			if (handle_to == 0 && handle_from == 0) {
				
				mGPU.select_device(device_to);	//not strictly necessary, but if not in the right context already, use of UVA will result in additional significant delay (according to benchmarks!)
				cudaMemcpy(gpu_memory_handles[device_to], gpu_memory_handles[device_from], sizeof(VType)*size, cudaMemcpyDeviceToDevice);
			}
			else if (handle_to == 0) {

				mGPU.select_device(device_to);	//not strictly necessary, but if not in the right context already, use of UVA will result in additional significant delay (according to benchmarks!)
				cudaMemcpy(gpu_memory_handles[device_to], gpu_extra_memory_handles[device_from][handle_from - 1], sizeof(VType)*size, cudaMemcpyDeviceToDevice);
			}
			else if (handle_from == 0) {

				mGPU.select_device(device_to);	//not strictly necessary, but if not in the right context already, use of UVA will result in additional significant delay (according to benchmarks!)
				cudaMemcpy(gpu_extra_memory_handles[device_to][handle_to - 1], gpu_memory_handles[device_from], sizeof(VType)*size, cudaMemcpyDeviceToDevice);
			}
			else {

				mGPU.select_device(device_to);	//not strictly necessary, but if not in the right context already, use of UVA will result in additional significant delay (according to benchmarks!)
				cudaMemcpy(gpu_extra_memory_handles[device_to][handle_to - 1], gpu_extra_memory_handles[device_from][handle_from - 1], sizeof(VType)*size, cudaMemcpyDeviceToDevice);
			}
		}
		else {

			//Indirect transfer through host (much slower than P2P)
			
			mGPU.select_device(device_from);
			if (handle_from == 0) cudaMemcpy(cpu_memory, gpu_memory_handles[device_from], sizeof(VType)*size, cudaMemcpyDeviceToHost);
			else cudaMemcpy(cpu_memory, gpu_extra_memory_handles[device_from][handle_from - 1], sizeof(VType)*size, cudaMemcpyDeviceToHost);

			mGPU.select_device(device_to);
			if (handle_to == 0) cudaMemcpy(gpu_memory_handles[device_to], cpu_memory, sizeof(VType)*size, cudaMemcpyHostToDevice);
			else cudaMemcpy(gpu_extra_memory_handles[device_to][handle_to - 1], cpu_memory, sizeof(VType)*size, cudaMemcpyHostToDevice);
		}
	}

	//------------------------------------------- DO REDUCTION

	//reduce values from all gpus to a single output value
	//typically each gpu will do its own reduction, resulting in values on each gpu to which we have handles gpu_memory_handles
	//here we just need a final reduction of these values
	VType reduce(void);
};