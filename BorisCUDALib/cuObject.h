#pragma once

#include "alloc_cpy.h"

#include <type_traits>

//cu_obj resides in cpu memory and manages objects in gpu memory.

//EXAMPLES:

//cu_obj<cuVEC<float>> vec(cuSZ3(100));

//The above makes a cuVEC<float> in gpu memory and is initialized by calling the cuVEC constructor with a cuSZ3(100), i.e. with 100*100*100 cells.

//vec()->set(1.0);

//The above calls the cuVEC host method set, which sets all cell values to 1.0.

//To pass the cuVEC in a cuda kernel launch then:

//__global__ cuda_kernel(cuVEC<float>& vec);			//this is the cuda kernel declaration
//cuda_kernel<<<...>>>(vec);

//In the above call the conversion operator passes by reference the dereferenced managed type (cuVEC<float>).
//In the cuda kernel use vec just as you might use it in cpu code, i.e. call any of its __device__ methods directly and use any of its accessible data members directly without any further complications (they are all stored in gpu memory)

template <typename MType>
class cu_obj
{

private:

	MType * managed_cuda_object = nullptr;

public:

	//------------------------------------------- CONSTRUCTOR

	//void constructor
	cu_obj(void)
	{
		//allocate memory for managed object
		gpu_alloc(managed_cuda_object);

		//call "managed" constructor - managed object has to implement this, where it will allocate gpu memory for its stored types
		call_constructor(std::bool_constant< std::is_fundamental<MType>::value >());
	}

	//parameters to pass to managed constructor
	template <typename ... PType>
	cu_obj(PType ... params)
	{
		//allocate memory for managed object
		gpu_alloc(managed_cuda_object);

		//call "managed" constructor - managed object has to implement this, where it will allocate gpu memory for its stored types
		call_constructor(std::bool_constant< std::is_fundamental<MType>::value >(), params...);
	}

	//------------------------------------------- DESTRUCTOR

	//destructor
	~cu_obj()
	{
		//free memory allocated for managed object
		if (managed_cuda_object) {

			//call "managed" destructor - managed object has to implement this, where it will free gpu memory for its stored types
			call_destructor(std::bool_constant< std::is_fundamental<MType>::value >());

			//free object itself before destructing the cu_obj
			gpu_free(managed_cuda_object);
		}
	}

	//------------------------------------------- COPY CONSTRUCTOR and ASSIGNMENT

	//implement MType::construct_cu_obj(const MType& coppyThis); in managed object
	cu_obj(const cu_obj& copyThis)
	{
		//allocate new memory for managed object
		gpu_alloc(managed_cuda_object);

		//call "managed" copy constructor - managed object has to implement this
		call_constructor_ref(std::bool_constant< std::is_fundamental<MType>::value >(), *(copyThis.managed_cuda_object));
	}

	//this makes an identical but independent managed cuda object : (implement void MType::assign_cu_obj(const MType& copyThis); in managed object)
	cu_obj& operator=(const cu_obj& copyThis)
	{
		//free memory allocated for managed object
		if (managed_cuda_object) {

			//call "managed" destructor - managed object has to implement this, where it will free gpu memory for its stored types
			call_destructor(std::bool_constant< std::is_fundamental<MType>::value >());

			//free object itself before destructing the cu_obj
			gpu_free(managed_cuda_object);
		}

		//allocate new memory for managed object
		gpu_alloc(managed_cuda_object);

		//call assignment operator on managed object - managed object has to implement this (implement void MType::assign_cu_obj(const MType& copyThis); in managed object)
		call_assignment(std::bool_constant< std::is_fundamental<MType>::value >(), *(copyThis.managed_cuda_object));

		return *this;
	}

	//------------------------------------------- CONVERSION

	//convertor : used to pass managed object by reference when launching a cuda kernel
	operator MType&()
	{
		return *managed_cuda_object;
	}

	//as above but with function call
	MType& get_dereferenced()
	{
		return *managed_cuda_object;
	}

	//get reference to managed object pointer : remember this is still in gpu memory; useful if we used with a gpu_to_... method from alloc_cpy.h to store in a cu_arr
	MType*& get_managed_object(void)
	{
		return managed_cuda_object;
	}

	//you can also use this instead of get_managed_object
	MType*& operator()(void)
	{
		return managed_cuda_object;
	}

	//------------------------------------------- SIMPLE READ / WRITE

	//for simple managed types (e.g a managed cuFLT3) can just read off the value to cpu directly using this
	MType to_cpu(void) const
	{
		MType cpu_value;

		gpu_to_cpu(&cpu_value, managed_cuda_object);
		
		return cpu_value;
	}

	//for simple managed types (e.g a managed cuFLT3) can just write the value to gpu directly using this
	void from_cpu(const MType& value)
	{
		MType local_value = value;
		cpu_to_gpu(managed_cuda_object, &local_value);
	}

	//------------------------------------------- CONSTRUCTOR / DESTRUCTOR CALLS

	//constructor for non-fundamental types
	template <typename ... PType>
	void call_constructor(std::false_type, PType ... params)
	{
		managed_cuda_object->construct_cu_obj(params...);
	}

	//constructor for fundamental types : nothing to construct
	void call_constructor(std::true_type)
	{
	}

	//constructor with parameter const references for non-fundamental types
	template <typename ... PType>
	void call_constructor_ref(std::false_type, const PType& ... params)
	{
		managed_cuda_object->construct_cu_obj(params...);
	}

	//constructor with references for fundamental types : nothing to construct
	void call_constructor_ref(std::true_type)
	{
	}

	//destructor (non-fundamental types)
	void call_destructor(std::false_type)
	{
		managed_cuda_object->destruct_cu_obj();
	}

	//destructor for fundamental types : nothing to destruct
	void call_destructor(std::true_type)
	{
	}

	//assignment with const reference : non-fundamental types
	void call_assignment(std::false_type, const MType& copyThis)
	{
		managed_cuda_object->assign_cu_obj(copyThis);
	}

	//assignment with const reference : fundamental types
	void call_assignment(std::true_type, const MType& copyThis)
	{
		gpu_to_gpu(managed_cuda_object, copyThis.managed_cuda_object);
	}
};
