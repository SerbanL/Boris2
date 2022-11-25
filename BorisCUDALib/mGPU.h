#pragma once

#include <cuda_runtime.h>

//EXAMPLE :

//To iterate over devices:

//for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {
//
// ... //mGPU is convertible to a linear device index
//}

//multi-GPU memory transfer configuration (p2p or via host)
//also used to iterate over devices
class mGPUConfig
{
private:

	//----------------- Memory transfer data

	//set to true if all memory transfers are p2p
	bool all_p2p = false;

	//special mode for testing : disable p2p access between all devices
	bool disable_p2p = false;

	//matrix to indicate type of transfer between any 2 devices. 0 : indirect (through host), 1 : p2p
	//NOTE : if all transfers are p2p, then clear this vector and set all_p2p = true
	//NOTE : first index is for device to transfer TO, second index is for device to transfer FROM
	int** transfer_type = nullptr;

	//----------------- Device iterator data

	int iterator_idx = 0;

	//configured device numbers this object manages - must ensure devices are actually available before configuring this object
	//index this to get the device number to set with cudaSetDevice, and nothing else: DO NOT USE pdevices[idx] AS AN INDEX!
	int* pdevices = nullptr;
	int num_devices = 0;

public:

	/////////////CONSTRUCTORS

	mGPUConfig(std::initializer_list<int> devices, bool disable_p2p_ = false)
	{
		disable_p2p = disable_p2p_;
		num_devices = devices.size();

		//store devices configured
		pdevices = new int[num_devices];
		std::copy(devices.begin(), devices.end(), pdevices);

		if (!disable_p2p) {

			//transfer_type size : number of devices squared
			transfer_type = new int*[num_devices];
			for (int idx = 0; idx < num_devices; idx++) transfer_type[idx] = new int[num_devices];

			all_p2p = true;

			for (int idx = 0; idx < num_devices; idx++) {

				cudaSetDevice(pdevices[idx]);
				cudaDeviceReset();
			}

			for (int idx_device_to = 0; idx_device_to < num_devices; idx_device_to++) {

				//the current device
				cudaSetDevice(pdevices[idx_device_to]);

				//check all other devices, which this current device can access in p2p mode
				for (int idx_device_from = 0; idx_device_from < num_devices; idx_device_from++) {

					transfer_type[idx_device_to][idx_device_from] = 1;
					if (pdevices[idx_device_to] == pdevices[idx_device_from]) continue;

					int can_access_peer;
					cudaDeviceCanAccessPeer(&can_access_peer, pdevices[idx_device_to], pdevices[idx_device_from]);

					if (can_access_peer == 0) {

						all_p2p = false;
						transfer_type[idx_device_to][idx_device_from] = 0;
					}
					//enable p2p from current device to device_from
					else cudaDeviceEnablePeerAccess(pdevices[idx_device_from], 0);
				}
			}

			//if all p2p (ideal case) then transfer_type matrix not needed
			if (all_p2p) {

				for (int idx = 0; idx < num_devices; idx++) {

					if (transfer_type[idx]) delete[] transfer_type[idx];
					transfer_type[idx] = nullptr;
				}

				delete[] transfer_type;
				transfer_type = nullptr;
			}
		}
	}

	~mGPUConfig()
	{
		if (pdevices) {

			delete[] pdevices;
			pdevices = nullptr;
		}

		if (transfer_type) {

			for (int idx = 0; idx < num_devices; idx++) {

				if (transfer_type[idx]) delete[] transfer_type[idx];
				transfer_type[idx] = nullptr;
			}

			delete[] transfer_type;
			transfer_type = nullptr;
		}
	}

	/////////////ITERATOR

	//start iterator
	mGPUConfig& device_begin(void)
	{
		iterator_idx = 0;
		cudaSetDevice(pdevices[iterator_idx]);
		return *this;
	}

	//check for last device
	int device_end(void) const { return num_devices; }

	//check iterator index
	bool operator!=(int value) { return iterator_idx != value; }

	//increment iterator
	mGPUConfig& operator++(int)
	{
		iterator_idx++;
		if (iterator_idx < num_devices) cudaSetDevice(pdevices[iterator_idx]);
		return *this;
	}

	/////////////INFO ITERATOR

	//number of devices available
	int get_num_devices(void) const { return num_devices; }

	//conversion operator to get linear device index (iterator value)
	operator int() const { return iterator_idx; }

	/////////////INFO MEMORY TRANSFER

	//check for a particular device combination if p2p is enabled
	bool is_p2p(int device_to, int device_from)
	{
		if (all_p2p) return true;
		else if (disable_p2p) return false;
		else return transfer_type[device_to][device_from];
	}

	bool is_all_p2p(void) { return all_p2p; }

	/////////////MANUAL DEVICE SELECTION

	//select a particular device with an index from 0 to num_devices
	void select_device(int idx) const { cudaSetDevice(pdevices[idx]); }

	bool select_previous_device(void) const
	{
		if (iterator_idx > 0) cudaSetDevice(pdevices[iterator_idx - 1]);
		else return false;
		return true;
	}

	bool select_next_device(void) const
	{
		if (iterator_idx < num_devices - 1) cudaSetDevice(pdevices[iterator_idx + 1]);
		else return false;
		return true;
	}

	void select_current_device(void) const { cudaSetDevice(pdevices[iterator_idx]); }
};