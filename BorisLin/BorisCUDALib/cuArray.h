#pragma once

#include "alloc_cpy.h"

//cu_arr is used to manage an array in gpu memory. cu_arr can reside on the cpu, or it can be cu_obj managed, in which case it resides on the gpu and can be passed into a cuda kernel where it can be used almost like a std::vector

//NOTE : currently only using cu_arr NOT cu_obj managed (this latter part of the code not done as I never had to use it -> could be useful/make some code neater but in practice you can get by without)

//EXAMPLE USAGE 1 - cu_arr is not managed (WORKS).
//
//cu_arr<float> arr_flt(1e6);
//the above declares an array of 1e6 floats in gpu memory.

//set all values to 1.0:
//arr_flt.set(1.0);

//if we have:
//__global__ cuda_kernel(size_t size, float* array);
//
//launch it as:
//
//cuda_kernel<<<...>>>(arr_flt.size(), arr_flt);
//
//arr_flt will be converted to a float*. NOTE : do not declare float*& array in the kernel parameters, as the pointer in the cu_arr object is in cpu memory (but it points to allocated gpu memory)
//
//EXAMPLE USAGE 2 - cu_arr is cu_obj managed (NOT AVAILABLE YET).
//
//cu_obj<cu_arr<float>> arr_flt(1e6);
//arr_flt()->set(1.0);
//
//if we have:
//__global__ cuda_kernel(cu_arr<float>& array);
//
//launch it as:
//
//cuda_kernel<<<...>>>(arr_flt);
//
//arr_flt can then be used almost like a std::vector : index it, get its size, etc.
//
//EXAMPLE USAGE 3 - cu_arr is cu_obj managed and stores pointers to other cu_obj managed objects (NOT AVAILABLE YET).
//
//Suppose we want to pass into a cuda kernel an array of cuVECs which are cu_obj managed.
//
//e.g. we have 2 cuVECs:

//cu_obj<cuVEC<float>> M(cuSZ3(100));
//cu_obj<cuVEC<float>> M2(cuSZ3(200));
//
//Then use:
//
//cu_obj<cu_arr<cuVEC<float>>> vec_arr;
//
//Store pointers to the cuVECs as:
//vec_arr()->push_back(M());
//vec_arr()->push_back(M2());
//
//Note, this differs from std::vector usage. the vec_arr actually stores pointers to the cuVECs, not copies of them.
//
//if we have:
//__global__ cuda_kernel(cu_arr<cuVEC<float>>& vec_arr);
//
//launch it as:
//cuda_kernel<<<...>>>(vec_arr);
//
//Inside the kernel access the cuVEC simply by indexing, e.g. vec_arr[0] is M, vec_arr[1] is M2.

//EXAMPLE USAGE 4 : cu_arr is not managed and stores pointers to other cu_obj managed objects (WORKS).
//
//Combine usage as in example 1 and 3:
//
//cu_obj<cuVEC<float>> M(cuSZ3(100));
//cu_obj<cuVEC<float>> M2(cuSZ3(200));
//
//cu_arr<cuVEC<float>> vec_arr;
//
//vec_arr.push_back(M.get_managed_object();
//vec_arr.push_back(M2.get_managed_object();
//

//EXAMPLE USAGE 5 : cu_arr stores pointers to arrays in gpu memory (WORKS).
//
//cu_arr<float> arr1(1000), arr2(2000);
//
//cu_arr<float*> arr_col;
//
//arr_col.push_back(arr1.get_managed_array());
//arr_col.push_back(arr2.get_managed_array());
//
//if we have __global__ cuda_kernel(float** arr_col);
//
//launch it as:
//cuda_kernel<<<...>>>(arr_col);
//
//inside the kernel you can now access arr1 elements as arr_col[0][idx], and arr2 elements as arr_col[1][idx]
//you might also need to pass info to the kernel about the number of stored arrays and dimensions of each array, etc.

template <typename VType>
class cu_arr
{

private:

	//array in gpu memory : cu_array itself is stored in cpu memory, but contains an address to gpu memory, i.e. you can only dereference it in device code
	VType* cu_array;

	//size value stored in cpu memory - only used if stand-alone object (as opposed to cu_obj managed)
	size_t arr_size;

	//pcu_array is stored in cpu memory, but contains an address to gpu memory, such that *pcu_array = cu_array, i.e. contains the value of cu_array but in gpu, not cpu, memory.
	VType** pcu_array;

private:

	//------------------------------------------- GPU MEMORY MANAGEMENT HELPERS : cuArray_mng_gpu.h


	//------------------------------------------- RESIZING cu_obj managed : cuArray_sizing.h

public:

	//--------------------------------------------CONSTRUCTORS : cu_obj "managed constructors" only : cuArray_mng_gpu.h


	//------------------------------------------- CONSTRUCTOR

	//void constructor
	__host__ cu_arr(void);

	//size constructor
	__host__ cu_arr(size_t size_);

	//------------------------------------------- DESTRUCTOR

	//destructor
	__host__ ~cu_arr();

	//------------------------------------------- GET ARRAY

	//use this when not cu_obj managed, e.g.:
	//cu_arr<float> arr;
	//if we have __global__ void cuda_kernel(size_t size, float* arr); then launch kernel as:
	//cuda_kernel<<<...>>>(arr.size(), arr); 
	__host__ operator VType*&() { return cu_array; }

	//get gpu address of stored array in cpu memory, i.e. get_managed_array() returns a pointer stored in cpu memory, which contains an address to gpu memory, such that *get_managed_array() is the stored array in gpu memory
	//This is useful to build a cu_arr of arrays.
	//Example if we have cu_arr<float> arr1, cu_arr<float*> arr_col, then you can use arr_col.push_back(arr1.get_managed_array()); arr_col can now be passed to a __global__, where arr_col[0][0] accesses first element in arr1, etc.
	VType**& get_managed_array(void)
	{
		return pcu_array;
	}

	//get gpu address of stored array in cpu memory directly; this can be used with thrust device pointers, e.g. thrust::device_ptr<VType> dev_ptr(arr.get_array());
	VType*& get_array(void)
	{
		return cu_array;
	}

	//------------------------------------------- INDEXING

	__device__ VType operator[](int idx) { return cu_array[idx]; }

	//------------------------------------------- RESIZING : cuArray_sizing.h

	__host__ bool resize(size_t size_);

	__host__ void clear(void);

	//------------------------------------------- STORE ENTRIES : cuArray_sizing.h

	//new_entry is a pointer in cpu memory to an object in gpu memory
	__host__ void push_back(VType*& new_entry);

	//new_entry_cpu resides in cpu memory
	//__host__ void push_back_from_cpu(VType& new_entry_cpu);

	//------------------------------------------- GPU <-> CPU TRANSFER : cuArray_transfer.h

	//copy values from a std::vector into gpu memory. Doesn't set size, but copies up to currently allocated size.
	template <typename Type>
	__host__ void copy_from_vector(std::vector<Type>& cpuvec);

	//copy values to a std::vector into cpu memory. Doesn't set size, but copies up to currently allocated size, starting at given offset in cpuvec
	template <typename Type>
	__host__ void copy_to_vector(std::vector<Type>& cpuvec, size_t offset = 0);

	//------------------------------------------- SET VALUE : cuArray_aux.h

	//set all entries to given value
	__host__ void set(VType value);

	//set single value from cpu memory at given index
	__host__ void setvalue(int index, VType value);

	//------------------------------------------- GET SIZE : cuArray_aux.h

	__host__ size_t size(void);

	//--------------------
};