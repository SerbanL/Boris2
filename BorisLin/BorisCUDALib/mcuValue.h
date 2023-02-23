#pragma once

#include "cuObject.h"
#include "mGPU.h"

//mcu_val is used to simple types, and is just a container of cu_obj, one for each configured device

template <typename VType>
class mcu_val
{

private:

	//used to iterate over configured devices, so user doesn't have to worry about selecting them, or which device numbers have been set.
	mGPUConfig& mGPU;

	//contained cu_obj objects, intended for simple types, one for each gpu configured.
	cu_obj<VType> ** pcuval = nullptr;

public:

	//------------------------------------------- CONSTRUCTOR

	//void constructor
	mcu_val(mGPUConfig& mGPU_) :
		mGPU(mGPU_)
	{
		//array of pointers set to number of required gpus
		pcuval = new cu_obj<VType>*[mGPU.get_num_devices()];

		//construct each cu_obj on the respective device
		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			pcuval[mGPU] = nullptr;
			pcuval[mGPU] = new cu_obj<VType>();
		}
	}

	//value constructor (each cu_obj gets same value)
	mcu_val(mGPUConfig& mGPU_, VType value) :
		mGPU(mGPU_)
	{
		//array of pointers set to number of required gpus
		pcuval = new cu_obj<VType>*[mGPU.get_num_devices()];

		//construct each cu_arr on the respective device
		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			pcuval[mGPU] = nullptr;
			pcuval[mGPU] = new cu_obj<VType>();
			pcuval[mGPU]->from_cpu(value);
		}
	}

	//------------------------------------------- DESTRUCTOR

	//destructor
	~mcu_val()
	{
		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			delete pcuval[mGPU];
			pcuval[mGPU] = nullptr;
		}

		delete[] pcuval;
		pcuval = nullptr;
	}

	//------------------------------------------- CONVERSION

	//used to pass managed object by reference when launching a cuda kernel
	VType& operator()(int idx)
	{
		return *pcuval[idx];
	}

	//------------------------------------------- SIMPLE READ / WRITE

	//read value in cpu memory for given device
	VType to_cpu(int idx) const
	{
		mGPU.select_device(idx);
		return pcuval[idx]->to_cpu();
	}

	//read value in cpu memory, as average from all devices
	VType to_cpu(void) const
	{
		VType value = VType();

		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			value += pcuval[mGPU]->to_cpu();
		}

		return value / mGPU.get_num_devices();
	}

	//read value in cpu memory, as sum from all devices
	VType to_cpu_sum(void) const
	{
		VType value = VType();

		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			value += pcuval[mGPU]->to_cpu();
		}

		return value;
	}

	//read value in cpu memory, as maximum from all devices
	VType to_cpu_max(void) const
	{
		VType value = VType(-1e38);

		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			VType value_device = pcuval[mGPU]->to_cpu();
			if (value_device > value) value = value_device;
		}

		return value;
	}

	//read value in cpu memory, as minimum from all devices
	VType to_cpu_min(void) const
	{
		VType value = VType(+1e38);

		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			VType value_device = pcuval[mGPU]->to_cpu();
			if (value_device < value) value = value_device;
		}

		return value;
	}

	//write the value to gpu directly using this (set value for indexed device)
	void from_cpu(int idx, const VType& value)
	{
		mGPU.select_device(idx);
		pcuval[idx]->from_cpu(value);
	}

	//write the value to gpu directly using this (set same value for all)
	void from_cpu(const VType& value)
	{
		for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

			pcuval[mGPU]->from_cpu(value);
		}
	}
};