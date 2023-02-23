#pragma once

#include "mcuObject.h"

#include "cu_prng.h"

//Policy class for managing cu_prng across multiple devices

class mcuBorisRand
{
private:

	//////////////////////////////////////////////////////////////////////////////////
	//
	// SPECIAL OBJECTS FROM CONSTRUCTOR (all Policy classes)

	//reference to mcu_obj manager for which this is a policy class
	mcu_obj<cuBorisRand<>, mcuBorisRand>& mng;

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
	mcuBorisRand(mcu_obj<cuBorisRand<>, mcuBorisRand>& mng_, mGPUConfig& mGPU_) :
		mng(mng_),
		mGPU(mGPU_)
	{}

	void construct_policy(void) {}

	void construct_policy(unsigned seed, size_t mem_size) { initialize(seed, mem_size); }

	//assignment operator
	mcuBorisRand& operator=(const mcuBorisRand& copyThis) { return *this; }

	//destructor
	virtual ~mcuBorisRand() { clear_memory_aux(); }

	//----------------------------------------- resize / reseed

	//make a prng with starting seed and memory size : in general recommended mem_size is the kernel size the prng will be used in divided by 128.
	//Larger memory size might not result in speed increase, lower memory size starts to slow the prng due to atomic operations used.
	cudaError_t initialize(unsigned seed, size_t mem_size)
	{
		cudaError_t error = cudaSuccess;
		
		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			//initialize each device, making sure they have different seeds (don't want same sequence generate in each!)
			mng(mGPU)->initialize(seed + mGPU, mem_size);
		}

		return error;
	}

	//resize array to just 1 element with a default seed of 0
	void clear(void)
	{
		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			mng(mGPU)->clear();
		}
	}
};

//Macro to simplify declaration
#define mcu_BorisRand mcu_obj<cuBorisRand, mcuBorisRand>
