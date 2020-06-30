#pragma once

#include <cuda_runtime.h>

#include "Funcs_Vectors.h"

//get total and free gpu memory in BYTEs
inline size_t cudaMemGetFree(void)
{
	size_t free_BYTEs, total_BYTEs;
	cudaMemGetInfo(&free_BYTEs, &total_BYTEs);
	return free_BYTEs;
}

inline size_t cudaMemGetTotal(void)
{
	size_t free_BYTEs, total_BYTEs;
	cudaMemGetInfo(&free_BYTEs, &total_BYTEs);
	return total_BYTEs;
}

//////////////////////////////////////////////// POINTERS TO GPU MEMORY STORED IN CPU MEMORY

//------------------------------------------------- GPU MEMORY MANAGEMENT using POINTERS STORED IN CPU MEMORY

//free gpu memory pointed to by cu_pointer, where cu_pointer itself is stored in cpu memory
template <typename Type>
cudaError_t gpu_free(Type*& cu_pointer)
{
	cudaError_t error = cudaSuccess;

	//free if memory allocated, then nullptr
	if (cu_pointer) error = cudaFree(cu_pointer);
	cu_pointer = nullptr;

	return error;
}

//allocate gpu memory (size elements) at location pointed to by cu_pointer, where cu_pointer itself is stored in cpu memory
template <typename Type>
cudaError_t gpu_alloc(Type*& cu_pointer, size_t size = 1)
{
	cudaError_t error = cudaSuccess;

	//cannot cudaMalloc with zero size
	if (!size) size = 1;

	//free if memory allocated (not nullptr)
	error = gpu_free(cu_pointer);
	if (error != cudaSuccess) return error;

	//allocate new memory - cudaMalloc always needs a cpu memory address for the first parameter
	error = cudaMalloc((void**)&cu_pointer, sizeof(Type) * size);
	if (error != cudaSuccess) cu_pointer = nullptr;

	return error;
}

//swap pointers so they swap memory pointed to, where both pointers are stored in cpu memory
template <typename Type>
void gpu_swap(Type*& cu_pointer1, Type*& cu_pointer2)
{
	//swap pointers
	Type* cu_pointer_swap = cu_pointer1;
	cu_pointer1 = cu_pointer2;
	cu_pointer2 = cu_pointer_swap;
}

//------------------------------------------------- MEMORY COPY TO/FROM GPU

//1. GPU to CPU

//copy size elements from gpu to cpu, where SType and DType are of same size and cu_source_pointer is stored in cpu memory
template <typename DType, typename SType, std::enable_if_t<sizeof(DType) == sizeof(SType)>* = nullptr>
cudaError_t gpu_to_cpu(DType* destination_pointer, const SType* cu_source_pointer, size_t size = 1)
{
	//nothing to copy? then do nothing (it's possible this can be called with zero size)
	if (size == 0) return cudaSuccess;

	return cudaMemcpy(destination_pointer, cu_source_pointer, sizeof(DType) * size, cudaMemcpyDeviceToHost);
}

//copy size elements from gpu to cpu, where SType and DType are not of same size, but convertible, and cu_source_pointer is stored in cpu memory
template <typename DType, typename SType, std::enable_if_t<sizeof(DType) != sizeof(SType)>* = nullptr>
cudaError_t gpu_to_cpu(DType* destination_pointer, const SType* cu_source_pointer, size_t size = 1)
{
	cudaError_t error = cudaSuccess;

	//nothing to copy? then do nothing (it's possible this can be called with zero size)
	if (size == 0) return cudaSuccess;

	//first copy from gpu to cpu for SType, then convert to DType
	std::vector<SType> convert_this;
	if (!malloc_vector(convert_this, size)) return cudaErrorMemoryAllocation;

	//transfer from gpu to cpu for SType
	error = gpu_to_cpu(convert_this.data(), cu_source_pointer, size);
	if (error != cudaSuccess) return error;

	//convert from SType to DType in cpu memory : conversion must be possible for SType to DType
#pragma omp parallel for
	for (int idx = 0; idx < size; idx++) {

		destination_pointer[idx] = (DType)convert_this[idx];
	}

	return error;
}

//return value in cpu memory from variable stored in gpu memory
template <typename Type>
Type get_gpu_value(const Type& cu_source_value)
{
	Type destination_value;

	if (cudaMemcpy(&destination_value, &cu_source_value, sizeof(Type), cudaMemcpyDeviceToHost) == cudaSuccess) return destination_value;
	else return Type();
}

//2. CPU to GPU

//copy size elements from cpu to gpu, where SType and DType are of same size and cu_destination_pointer is stored in cpu memory
template <typename DType, typename SType, std::enable_if_t<sizeof(DType) == sizeof(SType)>* = nullptr>
cudaError_t cpu_to_gpu(DType* cu_destination_pointer, const SType* source_pointer, size_t size = 1)
{
	//nothing to copy? then do nothing (it's possible this can be called with zero size)
	if (size == 0) return cudaSuccess;

	return cudaMemcpy(cu_destination_pointer, source_pointer, sizeof(DType) * size, cudaMemcpyHostToDevice);
}

//copy size elements from cpu to gpu, where SType and DType are not of same size, but convertible, and cu_destination_pointer is stored in cpu memory
template <typename DType, typename SType, std::enable_if_t<sizeof(DType) != sizeof(SType)>* = nullptr>
cudaError_t cpu_to_gpu(DType* cu_destination_pointer, const SType* source_pointer, size_t size = 1)
{
	//nothing to copy? then do nothing (it's possible this can be called with zero size)
	if (size == 0) return cudaSuccess;

	//first convert from SType to DType in cpu memory, then trasfer to GPU memory
	std::vector<DType> converted;
	if (!malloc_vector(converted, size)) return cudaErrorMemoryAllocation;

	//convert from SType to DType in cpu memory : conversion must be possible for SType to DType
#pragma omp parallel for
	for (int idx = 0; idx < size; idx++) {

		converted[idx] = (DType)source_pointer[idx];
	}

	//transfer from cpu to gpu for DSType
	return cpu_to_gpu(cu_destination_pointer, converted.data(), size);
}

//store cpu value in variable stored stored in gpu memory
template <typename Type>
cudaError_t set_gpu_value(Type& cu_destination_value, const Type& source_value)
{
	Type local_value = source_value;

	return cudaMemcpy(&cu_destination_value, &local_value, sizeof(Type), cudaMemcpyHostToDevice);
}

//3. GPU to GPU

//copy size elements from gpu to gpu, where SType and DType are of same size, and both cu_destination_pointer and cu_source_pointer are stored in cpu memory
template <typename DType, typename SType, std::enable_if_t<sizeof(DType) == sizeof(SType)>* = nullptr>
cudaError_t gpu_to_gpu(DType* cu_destination_pointer, const SType* cu_source_pointer, size_t size = 1)
{
	//nothing to copy? then do nothing (it's possible this can be called with zero size)
	if (size == 0) return cudaSuccess;

	return cudaMemcpy(cu_destination_pointer, cu_source_pointer, sizeof(DType) * size, cudaMemcpyDeviceToDevice);
}

//copy size elements from gpu to gpu, where SType and DType are not of same size, but convertible, and both cu_destination_pointer and cu_source_pointer are stored in cpu memory
template <typename DType, typename SType, std::enable_if_t<sizeof(DType) != sizeof(SType)>* = nullptr>
cudaError_t gpu_to_gpu(DType* cu_destination_pointer, const SType* cu_source_pointer, size_t size = 1)
{
	cudaError_t error = cudaSuccess;

	//nothing to copy? then do nothing (it's possible this can be called with zero size)
	if (size == 0) return cudaSuccess;

	//first copy from gpu to cpu for SType, then convert to DType
	std::vector<SType> convert_this;
	if (!malloc_vector(convert_this, size)) return cudaErrorMemoryAllocation;

	//transfer from gpu to cpu for SType
	error = gpu_to_cpu(convert_this.data(), cu_source_pointer, size);
	if (error != cudaSuccess) return error;

	//covert convert_this to DType before copying to destination
	std::vector<DType> converted;
	if (!malloc_vector(converted, size)) return cudaErrorMemoryAllocation;

	//convert from SType to DType in cpu memory : conversion must be possible for SType to DType
#pragma omp parallel for
	for (int idx = 0; idx < size; idx++) {

		converted[idx] = (DType)convert_this[idx];
	}

	return cpu_to_gpu(cu_destination_pointer, converted.data(), size);
}

//copy value from variable stored in gpu memory to variable stored in gpu memory
template <typename Type, std::enable_if_t<!std::is_pointer<Type>::value>* = nullptr>
cudaError_t gpu_to_gpu(Type& cu_destination_value, const Type& cu_source_value)
{
	return cudaMemcpy(&cu_destination_value, &cu_source_value, sizeof(Type), cudaMemcpyDeviceToDevice);
}

//------------------------------------------------- MEMORY INITIALIZATION

//this launches a kernel which initializes the cuda array to given value - specializations and kernel in "alloc_cpy.cuh"
//set value
template <typename DType>
void set_cuda_array_launcher(size_t size, DType* cu_dest_ptr, DType value);

//multiply by value
template <typename DType>
void mul_cuda_array_launcher(size_t size, DType* cu_dest_ptr, DType value);

//& with value
template <typename DType>
void and_cuda_array_launcher(size_t size, DType* cu_dest_ptr, DType value);

//| with value
template <typename DType>
void or_cuda_array_launcher(size_t size, DType* cu_dest_ptr, DType value);

//set value in array
template <typename DType>
void gpu_set(DType* cu_destination, DType value, size_t size)
{
	set_cuda_array_launcher(size, cu_destination, value);
}

//multiply value in array
template <typename DType>
void gpu_mul(DType* cu_destination, DType value, size_t size)
{
	mul_cuda_array_launcher(size, cu_destination, value);
}

//and value in array
template <typename DType>
void gpu_and(DType* cu_destination, DType value, size_t size)
{
	and_cuda_array_launcher(size, cu_destination, value);
}


//or value in array
template <typename DType>
void gpu_or(DType* cu_destination, DType value, size_t size)
{
	or_cuda_array_launcher(size, cu_destination, value);
}

//////////////////////////////////////////////// POINTERS TO GPU MEMORY STORED IN GPU MEMORY (referred to as managed since these pointers will be held in an object stored in gpu memory typically)

//------------------------------------------------- GPU NULLPTR

//null a managed gpu pointer, where cu_pointer itself is stored in gpu memory (is managed)
template <typename Type>
cudaError_t nullgpuptr(Type*& cu_pointer)
{
	Type* cu_pointer_handle = nullptr;
	return cudaMemcpy(&cu_pointer, &cu_pointer_handle, sizeof(Type*), cudaMemcpyHostToDevice);
}

//is the managed gpu pointer a nullptr? cu_pointer itself is stored in gpu memory (is managed)
template <typename Type>
bool isnullgpuptr(Type*& cu_pointer)
{
	Type* cu_pointer_handle = nullptr;

	cudaMemcpy(&cu_pointer_handle, &cu_pointer, sizeof(Type*), cudaMemcpyDeviceToHost);

	if (cu_pointer_handle) return false;
	return true;
}

//------------------------------------------------- GPU MEMORY MANAGEMENT for COMPLEX OBJECTS

//free gpu memory pointed to by cu_pointer, where cu_pointer itself is stored in gpu memory (is managed)
template <typename Type>
cudaError_t gpu_free_managed(Type*& cu_pointer)
{
	cudaError_t error = cudaSuccess;

	//only attempt to free if memory allocated
	if (isnullgpuptr(cu_pointer)) return error;

	//1. get handle to managed object
	Type* cu_pointer_handle = nullptr;

	error = cudaMemcpy(&cu_pointer_handle, &cu_pointer, sizeof(Type*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return error;

	//2. free memory using handle
	error = gpu_free(cu_pointer_handle);
	if (error != cudaSuccess) return error;

	//3. marked managed pointer as nullptr
	return nullgpuptr(cu_pointer);
}

//allocate gpu memory (size elements) at location pointed to by cu_pointer, where cu_pointer itself is stored in gpu memory (is managed)
template <typename Type>
cudaError_t gpu_alloc_managed(Type*& cu_pointer, size_t size = 1)
{
	cudaError_t error = cudaSuccess;

	//if memory allocated then free it first
	if (!isnullgpuptr(cu_pointer)) gpu_free_managed(cu_pointer);

	//1. allocate memory using a handle
	Type* cu_pointer_handle = nullptr;

	error = gpu_alloc(cu_pointer_handle, size);
	if (error != cudaSuccess) return error;

	//2. Finally bind allocated memory to managed object

	//cu_pointer_handle now contains the starting address of allocated gpu memory, but &cu_pointer_handle is the address of cu_pointer_handle stored in cpu memory
	//the following copies the starting address of allocated gpu memory to cu_pointer, so that cu_pointer then points to the allocated memory
	error = cudaMemcpy(&cu_pointer, &cu_pointer_handle, sizeof(Type*), cudaMemcpyHostToDevice);
	if (error != cudaSuccess) nullgpuptr(cu_pointer);

	return error;
}

//swap pointers where first parameter is for a managed array and the second is not. After this the first pointer will point to location of second array and vice-versa.
template <typename Type>
cudaError_t gpu_swap_managed(Type*& cu_pointer1_managed, Type*& cu_pointer2)
{
	cudaError_t error = cudaSuccess;

	//1. get handle to first array (the managed one)
	Type* cu_pointer1_handle = nullptr;

	error = cudaMemcpy(&cu_pointer1_handle, &cu_pointer1_managed, sizeof(Type*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return error;

	//2. swap pointers
	gpu_swap(cu_pointer1_handle, cu_pointer2);

	//3. bind handle back to managed array so it will point to memory location of the second array
	return cudaMemcpy(&cu_pointer1_managed, &cu_pointer1_handle, sizeof(Type*), cudaMemcpyHostToDevice);
}

//swap pointers where both parameter are for managed arrays. After this the first pointer will point to location of second array and vice-versa.
template <typename Type>
cudaError_t gpu_swap_managedx2(Type*& cu_pointer1_managed, Type*& cu_pointer2_managed)
{
	cudaError_t error = cudaSuccess;

	//1. get handle to first array
	Type* cu_pointer1_handle = nullptr;

	error = cudaMemcpy(&cu_pointer1_handle, &cu_pointer1_managed, sizeof(Type*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return error;

	//2. get handle to second array
	Type* cu_pointer2_handle = nullptr;

	error = cudaMemcpy(&cu_pointer2_handle, &cu_pointer2_managed, sizeof(Type*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return error;

	//3. swap pointers
	gpu_swap(cu_pointer1_handle, cu_pointer2_handle);

	//4. bind handle back to first array so it will point to memory location of the second array
	error = cudaMemcpy(&cu_pointer1_managed, &cu_pointer1_handle, sizeof(Type*), cudaMemcpyHostToDevice);
	if (error != cudaSuccess) return error;

	//5. bind handle back to second array so it will point to memory location of the first array
	return cudaMemcpy(&cu_pointer2_managed, &cu_pointer2_handle, sizeof(Type*), cudaMemcpyHostToDevice);
}

//------------------------------------------------- MEMORY COPY TO/FROM GPU for COMPLEX OBJECTS

//1. GPU to CPU

//copy size elements from gpu to cpu, where SType and DType are of same size and cu_source_pointer is stored in gpu memory (is managed)
template <typename DType, typename SType, std::enable_if_t<sizeof(DType) == sizeof(SType)>* = nullptr>
cudaError_t gpu_to_cpu_managed(DType* destination_pointer, SType*& cu_source_pointer, size_t size = 1, size_t offset = 0)
{
	cudaError_t error = cudaSuccess;

	//1. get handle to source object
	SType* cu_source_pointer_handle = nullptr;

	error = cudaMemcpy(&cu_source_pointer_handle, &cu_source_pointer, sizeof(SType*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return error;

	//2. copy from gpu using the handle we just obtained to cpu
	return cudaMemcpy(destination_pointer, cu_source_pointer_handle + offset, sizeof(DType) * size, cudaMemcpyDeviceToHost);
}

//copy size elements from gpu to cpu, where SType and DType are not of same size, but convertible, and cu_source_pointer is stored in gpu memory (is managed)
template <typename DType, typename SType, std::enable_if_t<sizeof(DType) != sizeof(SType)>* = nullptr>
cudaError_t gpu_to_cpu_managed(DType* destination_pointer, SType*& cu_source_pointer, size_t size = 1, size_t offset = 0)
{
	cudaError_t error = cudaSuccess;

	//1. get handle to source object
	SType* cu_source_pointer_handle = nullptr;

	error = cudaMemcpy(&cu_source_pointer_handle, &cu_source_pointer, sizeof(SType*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return error;

	//2. copy to cpu using the handle we just obtained (conversion done in gpu_to_cpu function)
	error = gpu_to_cpu(destination_pointer, cu_source_pointer_handle + offset, size);

	return error;
}

//2. CPU to GPU

//copy size elements from cpu to gpu, where SType and DType are of same size and cu_destination_pointer is stored in gpu memory (is managed)
template <typename DType, typename SType, std::enable_if_t<sizeof(DType) == sizeof(SType)>* = nullptr>
cudaError_t cpu_to_gpu_managed(DType*& cu_destination_pointer, SType* source_pointer, size_t size = 1, size_t offset = 0)
{
	cudaError_t error = cudaSuccess;

	//copy from cpu to gpu : get handle then copy to gpu memory using the handle

	//1. get handle to destination object
	DType* cu_destination_pointer_handle = nullptr;

	error = cudaMemcpy(&cu_destination_pointer_handle, &cu_destination_pointer, sizeof(DType*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return error;

	//2. copy from cpu to the handle we just obtained
	return cudaMemcpy(cu_destination_pointer_handle + offset, source_pointer, sizeof(DType) * size, cudaMemcpyHostToDevice);
}

//copy size elements from cpu to gpu, where SType and DType are not of same size, but convertible, and cu_destination_pointer is stored in gpu memory (is managed)
template <typename DType, typename SType, std::enable_if_t<sizeof(DType) != sizeof(SType)>* = nullptr>
cudaError_t cpu_to_gpu_managed(DType*& cu_destination_pointer, SType* source_pointer, size_t size = 1, size_t offset = 0)
{
	cudaError_t error = cudaSuccess;

	//1. get handle to destination object
	DType* cu_destination_pointer_handle = nullptr;

	error = cudaMemcpy(&cu_destination_pointer_handle, &cu_destination_pointer, sizeof(DType*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return error;

	//2. copy to gpu using the handle we just obtained (conversion done in cpu_to_gpu function)
	error = cpu_to_gpu(cu_destination_pointer_handle + offset, source_pointer, size);

	return error;
}

//3. GPU to GPU

//copy from gpu-addressable memory to gpu-addressable memory, where both source and destination are managed.
template <typename DType, typename SType, std::enable_if_t<sizeof(DType) == sizeof(SType)>* = nullptr>
cudaError_t gpu_to_gpu_managed(DType*& cu_destination_pointer, SType*& cu_source_pointer, size_t size = 1, size_t offset_dst = 0, size_t offset_src = 0)
{
	cudaError_t error = cudaSuccess;

	//copy from gpu to gpu : get handles then copy from gpu to gpu memory using the handles

	//1. get handle to destination object
	DType* cu_destination_pointer_handle = nullptr;

	error = cudaMemcpy(&cu_destination_pointer_handle, &cu_destination_pointer, sizeof(DType*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return error;

	//2. get handle to source object
	SType* cu_source_pointer_handle = nullptr;

	error = cudaMemcpy(&cu_source_pointer_handle, &cu_source_pointer, sizeof(SType*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return error;

	//3. copy from gpu to gpu using the handles we just obtained
	return cudaMemcpy(cu_destination_pointer_handle + offset_dst, cu_source_pointer_handle + offset_src, sizeof(DType) * size, cudaMemcpyDeviceToDevice);
}

//copy from gpu-addressable memory to gpu-addressable memory, where only source is managed.
template <typename DType, typename SType, std::enable_if_t<sizeof(DType) == sizeof(SType)>* = nullptr>
cudaError_t gpu_to_gpu_managed2nd(DType* cu_destination_pointer, SType*& cu_source_pointer, size_t size = 1, size_t offset = 0)
{
	cudaError_t error = cudaSuccess;

	//copy from gpu to gpu : get handles then copy from gpu to gpu memory using the handles

	//1. get handle to source object
	SType* cu_source_pointer_handle = nullptr;

	error = cudaMemcpy(&cu_source_pointer_handle, &cu_source_pointer, sizeof(SType*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return error;

	//2. copy from gpu to gpu using the handle we just obtained
	return cudaMemcpy(cu_destination_pointer, cu_source_pointer_handle + offset, sizeof(DType) * size, cudaMemcpyDeviceToDevice);
}

//copy from gpu-addressable memory to gpu-addressable memory, where only destination is managed.
template <typename DType, typename SType, std::enable_if_t<sizeof(DType) == sizeof(SType)>* = nullptr>
cudaError_t gpu_to_gpu_managed1st(DType*&  cu_destination_pointer, SType* cu_source_pointer, size_t size = 1, size_t offset = 0)
{
	cudaError_t error = cudaSuccess;

	//copy from gpu to gpu : get handles then copy from gpu to gpu memory using the handles

	//1. get handle to destination object
	DType* cu_destination_pointer_handle = nullptr;

	error = cudaMemcpy(&cu_destination_pointer_handle, &cu_destination_pointer, sizeof(DType*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return error;

	//2. copy from gpu to gpu using the handle we just obtained
	return cudaMemcpy(cu_destination_pointer_handle + offset, cu_source_pointer, sizeof(DType) * size, cudaMemcpyDeviceToDevice);
}

//------------------------------------------------- MEMORY INITIALIZATION

//set value in array : managed version
template <typename DType>
void gpu_set_managed(DType*& cu_destination_pointer, DType value, size_t size, size_t offset = 0)
{
	cudaError_t error = cudaSuccess;

	//1. get handle to destination object
	DType* cu_destination_pointer_handle = nullptr;

	error = cudaMemcpy(&cu_destination_pointer_handle, &cu_destination_pointer, sizeof(DType*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return;

	//2. set values using handle
	gpu_set(cu_destination_pointer_handle + offset, value, size);
}

//multiply value in array : managed version
template <typename DType>
void gpu_mul_managed(DType*& cu_destination_pointer, DType value, size_t size, size_t offset = 0)
{
	cudaError_t error = cudaSuccess;

	//1. get handle to destination object
	DType* cu_destination_pointer_handle = nullptr;

	error = cudaMemcpy(&cu_destination_pointer_handle, &cu_destination_pointer, sizeof(DType*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return;

	//2. multiply values using handle
	gpu_mul(cu_destination_pointer_handle + offset, value, size);
}

//and value in array : managed version
template <typename DType>
void gpu_and_managed(DType*& cu_destination_pointer, DType value, size_t size, size_t offset = 0)
{
	cudaError_t error = cudaSuccess;

	//1. get handle to destination object
	DType* cu_destination_pointer_handle = nullptr;

	error = cudaMemcpy(&cu_destination_pointer_handle, &cu_destination_pointer, sizeof(DType*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return;

	//2. set values using handle
	gpu_and(cu_destination_pointer_handle + offset, value, size);
}

//or value in array : managed version
template <typename DType>
void gpu_or_managed(DType*& cu_destination_pointer, DType value, size_t size, size_t offset = 0)
{
	cudaError_t error = cudaSuccess;

	//1. get handle to destination object
	DType* cu_destination_pointer_handle = nullptr;

	error = cudaMemcpy(&cu_destination_pointer_handle, &cu_destination_pointer, sizeof(DType*), cudaMemcpyDeviceToHost);
	if (error != cudaSuccess) return;

	//2. set values using handle
	gpu_or(cu_destination_pointer_handle + offset, value, size);
}