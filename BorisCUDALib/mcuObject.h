#pragma once

#include "alloc_cpy.h"

#include <type_traits>
#include "mGPU.h"

//mcu_obj resides in cpu memory and manages objects in multiple gpu memories (i.e. multiple devices)
//Requires a Policy class which 'overloads' methods of managed object (call these directly) so as to implement multiple device management.

//EXAMPLES:

//1. First need to configure a multi-GPU config object, e.g. for 3 devices:
//mGPUConfig mGPU({0, 1, 2});

//2. Next define the mcu_obj managed object, e.g. a cuVEC_VC across multiple devices as configured through mGPU (mcuVEC is the policy class):
//mcu_obj<cuVEC_VC<cuReal3>, mcuVEC<cuReal3, cuVEC_VC<cuReal3>>> M(mGPU);

//3. Now can access methods at the policy class level on M, which should mirror __host__ methods defined in cuVEC_VC (the policy class manages behaviour of these methods across multiple devices), e.g.:
//M.assign(h, rect, cuReal3(1, 2, 3));

//e.g. the above method sets rectangle rect and cellsize h, but the policy class actually handles 3 separate cuVEC_VC objects for the respective 3 devices, each ~1/3 of the required size.
//the user doesn't see all these details, to the user it simply appears as a single large cuVEC_VC object.

//4. M can be passed to a kernel (which accepts cuVEC_VC<cuReal3>&) for computations on a particular device, e.g.:
//
//for (auto& mGPU = M.device_begin(); mGPU != mGPU.device_end(); mGPU++) {
//
//	cuda_func_kernel <<< (M.device_size(mGPU) + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>> (M.get_deviceobject(mGPU));
//}
//
//NOTE: It's not necessary to use a "#pragma omp parallel for num_threads(mGPU.get_num_devices())" directive, since cuda kernel launch is asynchronous
//Thus the above for loop iterates over all devices, so kernels are launched simultaneously for all configured devices and work is done in parallel (in theory 3 times faster than for a single device configured).
//
//5. It may be necessary to exchange information between devices, e.g. for differential operator computations. In this case halos are used behind the scene. All the user has to do is exchange halo data before launching kernels where this would be needed (differential operators used):
//M.exchange_halos();
//
//6. Further advanced memory transfers between devices can be managed using a mGPU_Transfer object (see mGPU_Transfer.h).

//PRINCIPLES:

//1.
//All code works for any mGPU configuration and any number of devices (if devices are physically present), e.g. if mGPU{0}) is configured, then code produces same result, but 3 times slower, etc.
//Can even configure mGPU({0, 0, 0}) etc., but this is not very useful except for testing "multi-GPU" code when only a single device is present.

//2.
//The managed object should not be aware of any multi-GPU details, so any existing code can be converted with relatively little effort to multi-GPU code.
//i.e. all the managed object code executes on a single device. multi-device management is only done at the Policy class level.

//USAGE TEMPLATE :

/*
//class for object on device, managed by mcu_obj (i.e. the mcu_obj-managed object)
class PoissonRHS
{
public:

	//required
	__host__ void construct_cu_obj(void) {}
	//required
	__host__ void destruct_cu_obj(void) {}

	//optional
	__host__ void assign_cu_obj(const PoissonRHS& copyThis);

	//specific for PoissonRHS as example, but not required for this template
	__device__ cuBReal Poisson_RHS(int idx) { return 0.0; }
};

//policy class for access to mcu_obj-managed object host public methods, which implements multi-device coordination for managed object
//as defined, the policy class above is "empty", i.e. minimal definition required.
//to implement policy, it should redefine public __host__ methods from PoissonRHS, meant to be called directly, instead of those in PoissonRHS.
//this allows opportunity for multi-device management if needed. Note, PoissonRHS could also be cu_obj managed (i.e. single device code), where its methods would be called through cu_obj instead.
class Policy_PoissonRHS
{
private:

	//////////////////////////////////////////////////////////////////////////////////
	//
	// SPECIAL OBJECTS FROM CONSTRUCTOR (all Policy classes)

	//reference to mcu_obj manager for which this is a policy class
	mcu_obj<PoissonRHS, Policy_PoissonRHS>& mng;

	//multi-GPU configuration (list of physical devices configured with memory transfer type configuration set)
	mGPUConfig& mGPU;

private:

	//////////////////////////////////////////////////////////////////////////////////
	//
	// AUXILIARY (all Policy classes)

	//clear all allocated memory
	void clear_memory_aux(void) {}

public:

	//////////////////////////////////////////////////////////////////////////////////
	//
	// CONSTRUCTORS (all Policy classes)

	//--------------------------------------------CONSTRUCTORS : mcuVEC_mng.h

	//constructor for this policy class : void
	Policy_PoissonRHS(mcu_obj<PoissonRHS, Policy_PoissonRHS>& mng_, mGPUConfig& mGPU_) :
		mng(mng_),
		mGPU(mGPU_)
	{}

	void construct_policy(void) {}

	//assignment operator
	Policy_PoissonRHS& operator=(const Policy_PoissonRHS& copyThis) {}

	//destructor
	virtual ~Policy_PoissonRHS() { clear_memory_aux(); }
};

//Declare object as:
//mcu_obj<PoissonRHS, Policy_PoissonRHS> mcuPoissonRHS(mGPU);
*/

//MType : managed type
//Policy : policy class for managed type
template <typename MType, typename Policy>
class mcu_obj :
	public Policy
{
	friend Policy;

private:

	//used to iterate over configured devices, so user doesn't have to worry about selecting them, or which device numbers have been set.
	mGPUConfig& mGPU;

	//array of pointers to gpu memory locations of managed cu objects
	//array size is set to number of available gpus
	//each managed object is for the corresponding gpu, e.g. at index 0 we have device 0, etc.
	MType ** managed_cuda_objects = nullptr;

private:

	//clear all allocated memory
	void clear(void)
	{
		//free memory allocated for managed object
		if (managed_cuda_objects) {

			for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {
			
				//call "managed" destructor - managed object has to implement this, where it will free gpu memory for its stored types
				call_destructor(mGPU, std::integral_constant<bool, std::is_fundamental<MType>::value >());

				//free object itself before destructing the cu_obj
				gpu_free(managed_cuda_objects[mGPU]);
				managed_cuda_objects[mGPU] = nullptr;
			}

			//free array made in constructor
			delete[] managed_cuda_objects;
			managed_cuda_objects = nullptr;
		}
	}

public:

	//------------------------------------------- CONSTRUCTOR

	//void constructor
	mcu_obj(mGPUConfig& mGPU_) :
		mGPU(mGPU_),
		Policy(*this, mGPU_)
	{
		//array of pointers set to number of required gpus
		managed_cuda_objects = new MType*[mGPU.get_num_devices()];

		//construct each object on the respective device
		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			//allocate memory for managed object
			managed_cuda_objects[mGPU] = nullptr;
			gpu_alloc(managed_cuda_objects[mGPU]);

			//call "managed" constructor - managed object has to implement this, where it will allocate gpu memory for its stored types
			call_constructor(mGPU, std::integral_constant<bool, std::is_fundamental<MType>::value >());
		}

		Policy::construct_policy();
	}

	//parameters to pass to managed constructor
	//this needs to be handled carefuly. Cannot pass parameters to Policy class constructor, since this is executed before mcu_obj constructor, and it requires mcu_obj to be finished.
	//also cannot call managed object constructor with parameters from here, the policy class must handle this.
	//instead, call Policy class void constructor first, then create managed objects with void constructors.
	//after this, call Policy class constructor with parameters (construct_policy), which will coordinate the managed objects as needed.
	template <typename ... PType>
	mcu_obj(mGPUConfig& mGPU_, PType ... params) :
		mGPU(mGPU_),
		Policy(*this, mGPU_)
	{
		//array of pointers set to number of required gpus
		managed_cuda_objects = new MType*[mGPU.get_num_devices()];

		//construct each object on the respective device
		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			//allocate memory for managed object
			managed_cuda_objects[mGPU] = nullptr;
			gpu_alloc(managed_cuda_objects[mGPU]);

			//call "managed" constructor (without parameters) - managed object has to implement this, where it will allocate gpu memory for its stored types
			call_constructor(mGPU, std::integral_constant<bool, std::is_fundamental<MType>::value >());
		}

		//now that everything has been constructed (memory allocated for objects), call Policy class "constructor" with parameters.
		Policy::construct_policy(params...);
	}

	//------------------------------------------- DESTRUCTOR

	//destructor
	~mcu_obj()
	{
		clear();
	}

	//------------------------------------------- COPY CONSTRUCTOR and ASSIGNMENT

	//implement MType::construct_cu_obj(const MType& copyThis); in managed object
	mcu_obj(const mcu_obj& copyThis) :
		mGPU(copyThis.mGPU),
		Policy(*this, copyThis.mGPU)
	{
		//array of pointers set to number of required gpus
		managed_cuda_objects = new MType*[mGPU.get_num_devices()];

		//construct each object on the respective device
		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			//allocate memory for managed object
			managed_cuda_objects[mGPU] = nullptr;
			gpu_alloc(managed_cuda_objects[mGPU]);

			//call "managed" copy constructor - managed object has to implement this
			call_constructor_ref(mGPU, std::integral_constant<bool, std::is_fundamental<MType>::value >(), *(copyThis.managed_cuda_objects[mGPU]));
		}

		//invoke assignment operator on Policy
		*static_cast<Policy*>(this) = *static_cast<const Policy*>(&copyThis);
	}
	
	//this makes an identical but independent managed cuda object : (implement void MType::assign_cu_obj(const MType& copyThis); in managed object)
	//copyThis must have same number of gpu devices set
	mcu_obj& operator=(const mcu_obj& copyThis)
	{
		//free memory allocated for managed object first
		clear();

		//copy from copyThis
		managed_cuda_objects = new MType*[mGPU.get_num_devices()];

		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			//allocate new memory for managed object
			managed_cuda_objects[mGPU] = nullptr;
			gpu_alloc(managed_cuda_objects[mGPU]);

			//call assignment operator on managed object - managed object has to implement this (implement void MType::assign_cu_obj(const MType& copyThis); in managed object)
			call_assignment(mGPU, std::integral_constant<bool, std::is_fundamental<MType>::value >(), *(copyThis.managed_cuda_objects[mGPU]));
		}

		//invoke assignment operator on Policy
		*static_cast<Policy*>(this) = *static_cast<const Policy*>(&copyThis);

		return *this;
	}

	//------------------------------------------- CONVERSION

	//as above but with function call
	//useful for passing managed device to __global__ kernel when iterating over devices, e.g. M.get_deviceobject(mGPU)
	MType& get_deviceobject(int idx)
	{
		return *managed_cuda_objects[idx];
	}

	//get reference to managed object pointer : remember this is still in gpu memory; useful if we used with a gpu_to_... method from alloc_cpy.h to store in a cu_arr
	MType*& get_managed_object(int idx)
	{
		return managed_cuda_objects[idx];
	}

	//you can also use this instead of get_managed_object
	//Useful for accessing __host__ methods in managed device, e.g. M(mGPU)->...
	MType*& operator()(int idx)
	{
		return managed_cuda_objects[idx];
	}

	//------------------------------------------- SIMPLE READ / WRITE

	//mcu_obj not meant to manage simple types. instead use mcu_val for simplicity 
	//(mcu_obj could be made to manage simple types by defining a separate policy class, but it seems overly complicated in this case; instead, all we need is a simple container of cu_obj, i.e. mcu_val).

	//------------------------------------------- CONSTRUCTOR / DESTRUCTOR CALLS

	//constructor for non-fundamental types
	template <typename ... PType>
	void call_constructor(int idx, std::false_type, PType ... params)
	{
		managed_cuda_objects[idx]->construct_cu_obj(params...);
	}

	//constructor for fundamental types : nothing to construct
	void call_constructor(int idx, std::true_type)
	{
	}

	//constructor with parameter const references for non-fundamental types
	template <typename ... PType>
	void call_constructor_ref(int idx, std::false_type, const PType& ... params)
	{
		managed_cuda_objects[idx]->construct_cu_obj(params...);
	}

	//constructor with references for fundamental types : nothing to construct
	void call_constructor_ref(int idx, std::true_type)
	{
	}

	//destructor (non-fundamental types)
	void call_destructor(int idx, std::false_type)
	{
		managed_cuda_objects[idx]->destruct_cu_obj();
	}

	//destructor for fundamental types : nothing to destruct
	void call_destructor(int idx, std::true_type)
	{
	}

	//assignment with const reference : non-fundamental types
	void call_assignment(int idx, std::false_type, const MType& copyThis)
	{
		managed_cuda_objects[idx]->assign_cu_obj(copyThis);
	}

	//assignment with const reference : fundamental types
	void call_assignment(int idx, std::true_type, const MType& copyThis)
	{
		gpu_to_gpu(managed_cuda_objects[idx], copyThis.managed_cuda_objects[idx]);
	}

	//------------------------------------------- DEVICE ITERATOR

	mGPUConfig& device_begin(void) { return mGPU.device_begin(); }
};
