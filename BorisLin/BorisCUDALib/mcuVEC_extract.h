#pragma once

#include "mcuVEC.h"

//--------------------------------------------EXTRACT A LINE PROFILE : mcuVEC_extract.h

//cordinates are relative

//1. Profile values only, without stencil operation

//extract profile to a mcu_arr : extract size points starting at start in the direction step for the given number of points (size); use weighted average to extract profile with h stencil only
//profile_gpu resized as needed for each device
template <typename VType, typename MType>
void mcuVEC<VType, MType>::extract_profilevalues(size_t size, mcu_arr<VType>& profile_gpu, cuReal3 start, cuReal3 step)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		mng(0)->extract_profilevalues(size, profile_gpu.get_cu_arr(0), start, step);
		return;
	}

	cuReal3 end = (start + (size - 1) * step);
	cuReal3 intersection_s;
	cuReal3 intersection_e;
	
	//linear step
	cuBReal lstep = step.norm();
	if (!lstep) return;

	//normalized vector step
	cuReal3 nstep = step.normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//for given start and end points, find intersection values with each sub-rect (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / lstep) * lstep) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / lstep) * lstep) * nstep;

			//extract profile on this device with correct size
			size_t sub_size = round((intersection_e - intersection_s).norm() / lstep) + 1;
			//coordinates now must be relative to subrect
			mng(mGPU)->extract_profilevalues(sub_size, profile_gpu.get_cu_arr(mGPU), intersection_s + rect.s - prect_d[mGPU].s, step);
		}
		else {

			//cannot extract profile in this sub-rect, so make sure to clear respective cu_arr
			if (profile_gpu.size(mGPU)) profile_gpu.clear(mGPU);
		}
	}
}

//these specifically apply for VType == cuReal3, allowing extraction of the x, y, z components separately
template <typename VType, typename MType>
void mcuVEC<VType, MType>::extract_profilevalues_component_x(size_t size, mcu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		mng(0)->extract_profilevalues_component_x(size, profile_gpu.get_cu_arr(0), start, step);
		return;
	}

	cuReal3 end = (start + (size - 1) * step);
	cuReal3 intersection_s;
	cuReal3 intersection_e;

	//linear step
	cuBReal lstep = step.norm();
	if (!lstep) return;

	//normalized vector step
	cuReal3 nstep = step.normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//for given start and end points, find intersection values with each sub-rect  (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / lstep) * lstep) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / lstep) * lstep) * nstep;

			//extract profile on this device with correct size
			size_t sub_size = round((intersection_e - intersection_s).norm() / lstep) + 1;
			//coordinates now must be relative to subrect
			mng(mGPU)->extract_profilevalues_component_x(sub_size, profile_gpu.get_cu_arr(mGPU), intersection_s + rect.s - prect_d[mGPU].s, step);
		}
		else {

			//cannot extract profile in this sub-rect, so make sure to clear respective cu_arr
			if (profile_gpu.size(mGPU)) profile_gpu.clear(mGPU);
		}
	}
}

template <typename VType, typename MType>
void mcuVEC<VType, MType>::extract_profilevalues_component_y(size_t size, mcu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		mng(0)->extract_profilevalues_component_y(size, profile_gpu.get_cu_arr(0), start, step);
		return;
	}

	cuReal3 end = (start + (size - 1) * step);
	cuReal3 intersection_s;
	cuReal3 intersection_e;

	//linear step
	cuBReal lstep = step.norm();
	if (!lstep) return;

	//normalized vector step
	cuReal3 nstep = step.normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//for given start and end points, find intersection values with each sub-rect  (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / lstep) * lstep) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / lstep) * lstep) * nstep;

			//extract profile on this device with correct size
			size_t sub_size = round((intersection_e - intersection_s).norm() / lstep) + 1;
			//coordinates now must be relative to subrect
			mng(mGPU)->extract_profilevalues_component_y(sub_size, profile_gpu.get_cu_arr(mGPU), intersection_s + rect.s - prect_d[mGPU].s, step);
		}
		else {

			//cannot extract profile in this sub-rect, so make sure to clear respective cu_arr
			if (profile_gpu.size(mGPU)) profile_gpu.clear(mGPU);
		}
	}
}

template <typename VType, typename MType>
void mcuVEC<VType, MType>::extract_profilevalues_component_z(size_t size, mcu_arr<cuBReal>& profile_gpu, cuReal3 start, cuReal3 step)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		mng(0)->extract_profilevalues_component_z(size, profile_gpu.get_cu_arr(0), start, step);
		return;
	}

	cuReal3 end = (start + (size - 1) * step);
	cuReal3 intersection_s;
	cuReal3 intersection_e;

	//linear step
	cuBReal lstep = step.norm();
	if (!lstep) return;

	//normalized vector step
	cuReal3 nstep = step.normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//for given start and end points, find intersection values with each sub-rect  (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / lstep) * lstep) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / lstep) * lstep) * nstep;

			//extract profile on this device with correct size
			size_t sub_size = round((intersection_e - intersection_s).norm() / lstep) + 1;
			//coordinates now must be relative to subrect
			mng(mGPU)->extract_profilevalues_component_z(sub_size, profile_gpu.get_cu_arr(mGPU), intersection_s + rect.s - prect_d[mGPU].s, step);
		}
		else {

			//cannot extract profile in this sub-rect, so make sure to clear respective cu_arr
			if (profile_gpu.size(mGPU)) profile_gpu.clear(mGPU);
		}
	}
}

//--------------------------------------------

//2. Profile values only, with stencil operation around profile point

//extract profile components: extract starting at start in the direction end - step, with given step; use weighted average to extract profile with given stencil
//all coordinates are relative positions. Return profile values in profile_gpu (which is resized as needed for each device)
template <typename VType, typename MType>
bool mcuVEC<VType, MType>::extract_profile(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, mcu_arr<VType>& profile_gpu)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		return mng(0)->extract_profile(start, end, step, stencil, profile_gpu.get_cu_arr(0));
	}

	cuReal3 intersection_s;
	cuReal3 intersection_e;

	//normalized vector step (profile direction)
	cuReal3 nstep = (end - start).normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//for given start and end points, find intersection values with each sub-rect  (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / step) * step) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / step) * step) * nstep;

			//extract profile on this device
			//coordinates now must be relative to subrect
			if (!mng(mGPU)->extract_profile(intersection_s + rect.s - prect_d[mGPU].s, intersection_e + rect.s - prect_d[mGPU].s, step, stencil, profile_gpu.get_cu_arr(mGPU))) return false;
		}
		else {

			//cannot extract profile in this sub-rect, so make sure to clear respective cu_arr
			if (profile_gpu.size(mGPU)) profile_gpu.clear(mGPU);
		}
	}

	return true;
}

//as above, but only store profile in internal memory (line_profile) so we can read it out later as needed
template <typename VType, typename MType>
bool mcuVEC<VType, MType>::extract_profile(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		return mng(0)->extract_profile(start, end, step, stencil);
	}

	cuReal3 intersection_s;
	cuReal3 intersection_e;

	//normalized vector step (profile direction)
	cuReal3 nstep = (end - start).normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//for given start and end points, find intersection values with each sub-rect  (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / step) * step) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / step) * step) * nstep;

			//extract profile on this device
			//coordinates now must be relative to subrect
			if (!mng(mGPU)->extract_profile(intersection_s + rect.s - prect_d[mGPU].s, intersection_e + rect.s - prect_d[mGPU].s, step, stencil)) return false;
		}
	}

	return true;
}

//extract profile components: extract starting at start in the direction end - step, with given step; use weighted average to extract profile with given stencil
//all coordinates are relative positions. Return profile values in profile_cpu.
template <typename VType, typename MType>
template <typename SType>
bool mcuVEC<VType, MType>::extract_profile(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<SType>& profile_cpu)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		return mng(0)->extract_profile(start, end, step, stencil, profile_cpu);
	}

	cuReal3 intersection_s;
	cuReal3 intersection_e;

	//normalized vector step (profile direction)
	cuReal3 nstep = (end - start).normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//for given start and end points, find intersection values with each sub-rect  (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / step) * step) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / step) * step) * nstep;

			//extract profile on this device to auxiliary storage first
			//coordinates now must be relative to subrect
			if (!mng(mGPU)->extract_profile(intersection_s + rect.s - prect_d[mGPU].s, intersection_e + rect.s - prect_d[mGPU].s, step, stencil, profile_aux.get_cu_arr(mGPU))) return false;
		}
		else {

			//cannot extract profile in this sub-rect, so make sure to clear respective auxiliar storage cu_arr
			if (profile_aux.size(mGPU)) profile_aux.clear(mGPU);
		}
	}

	//now copy profile to profile_cpu, making sure it has the correct size
	size_t size = round((end - start).norm() / step) + 1;
	if (profile_cpu.size() != size) if (!malloc_vector(profile_cpu, size)) return false;
	profile_aux.copy_to_vector(profile_cpu);

	return true;
}

//--------------------------------------------

//3. Profile values for individual components together with profile position returned in cuReal2/DBL2, with stencil operation around profile point (capped to mesh size)

//as above but only component x and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_gpu (which is resized as needed for each device)
template <typename VType, typename MType>
bool mcuVEC<VType, MType>::extract_profile_component_x(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, mcu_arr<cuReal2>& profile_gpu)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		return mng(0)->extract_profile_component_x(start, end, step, stencil, profile_gpu.get_cu_arr(0));
	}

	cuReal3 intersection_s;
	cuReal3 intersection_e;

	//normalized vector step (profile direction)
	cuReal3 nstep = (end - start).normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {
		 
		//for given start and end points, find intersection values with each sub-rect  (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / step) * step) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / step) * step) * nstep;

			//extract profile on this device
			//coordinates now must be relative to subrect
			if (!mng(mGPU)->extract_profile_component_x(intersection_s + rect.s - prect_d[mGPU].s, intersection_e + rect.s - prect_d[mGPU].s, step, stencil, profile_gpu.get_cu_arr(mGPU))) return false;
		}
		else {

			//cannot extract profile in this sub-rect, so make sure to clear respective cu_arr
			if (profile_gpu.size(mGPU)) profile_gpu.clear(mGPU);
		}
	}

	return true;
}

//as above but only component x and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_cpu. profile_cpu resized as needed.
template <typename VType, typename MType>
bool mcuVEC<VType, MType>::extract_profile_component_x(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		return mng(0)->extract_profile_component_x(start, end, step, stencil, profile_cpu);
	}

	cuReal3 intersection_s;
	cuReal3 intersection_e;

	//normalized vector step (profile direction)
	cuReal3 nstep = (end - start).normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//for given start and end points, find intersection values with each sub-rect (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / step) * step) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / step) * step) * nstep;

			//extract profile on this device to auxiliary storage first
			//coordinates now must be relative to subrect
			if (!mng(mGPU)->extract_profile_component_x(intersection_s + rect.s - prect_d[mGPU].s, intersection_e + rect.s - prect_d[mGPU].s, step, stencil, profile_component_aux.get_cu_arr(mGPU))) return false;
		}
		else {

			//cannot extract profile in this sub-rect, so make sure to clear respective auxiliar storage cu_arr
			if (profile_component_aux.size(mGPU)) profile_component_aux.clear(mGPU);
		}
	}

	//now copy profile to profile_cpu, making sure it has the correct size
	size_t size = round((end - start).norm() / step) + 1;
	if (profile_cpu.size() != size) if (!malloc_vector(profile_cpu, size)) return false;
	profile_component_aux.copy_to_vector(profile_cpu);

	return true;
}


//as above but only component y and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_gpu (which is resized as needed for each device)
template <typename VType, typename MType>
bool mcuVEC<VType, MType>::extract_profile_component_y(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, mcu_arr<cuReal2>& profile_gpu)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		return mng(0)->extract_profile_component_y(start, end, step, stencil, profile_gpu.get_cu_arr(0));
	}

	cuReal3 intersection_s;
	cuReal3 intersection_e;

	//normalized vector step (profile direction)
	cuReal3 nstep = (end - start).normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//for given start and end points, find intersection values with each sub-rect  (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / step) * step) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / step) * step) * nstep;

			//extract profile on this device
			//coordinates now must be relative to subrect
			if (!mng(mGPU)->extract_profile_component_y(intersection_s + rect.s - prect_d[mGPU].s, intersection_e + rect.s - prect_d[mGPU].s, step, stencil, profile_gpu.get_cu_arr(mGPU))) return false;
		}
		else {

			//cannot extract profile in this sub-rect, so make sure to clear respective cu_arr
			if (profile_gpu.size(mGPU)) profile_gpu.clear(mGPU);
		}
	}

	return true;
}

//as above but only component y and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_cpu. profile_cpu resized as needed.
template <typename VType, typename MType>
bool mcuVEC<VType, MType>::extract_profile_component_y(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		return mng(0)->extract_profile_component_y(start, end, step, stencil, profile_cpu);
	}

	cuReal3 intersection_s;
	cuReal3 intersection_e;

	//normalized vector step (profile direction)
	cuReal3 nstep = (end - start).normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//for given start and end points, find intersection values with each sub-rect  (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / step) * step) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / step) * step) * nstep;

			//extract profile on this device to auxiliary storage first
			//coordinates now must be relative to subrect
			if (!mng(mGPU)->extract_profile_component_y(intersection_s + rect.s - prect_d[mGPU].s, intersection_e + rect.s - prect_d[mGPU].s, step, stencil, profile_component_aux.get_cu_arr(mGPU))) return false;
		}
		else {

			//cannot extract profile in this sub-rect, so make sure to clear respective auxiliar storage cu_arr
			if (profile_component_aux.size(mGPU)) profile_component_aux.clear(mGPU);
		}
	}

	//now copy profile to profile_cpu, making sure it has the correct size
	size_t size = round((end - start).norm() / step) + 1;
	if (profile_cpu.size() != size) if (!malloc_vector(profile_cpu, size)) return false;
	profile_component_aux.copy_to_vector(profile_cpu);

	return true;
}

//as above but only component z and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_gpu (which is resized as needed for each device)
template <typename VType, typename MType>
bool mcuVEC<VType, MType>::extract_profile_component_z(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, mcu_arr<cuReal2>& profile_gpu)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		return mng(0)->extract_profile_component_z(start, end, step, stencil, profile_gpu.get_cu_arr(0));
	}

	cuReal3 intersection_s;
	cuReal3 intersection_e;

	//normalized vector step (profile direction)
	cuReal3 nstep = (end - start).normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//for given start and end points, find intersection values with each sub-rect  (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / step) * step) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / step) * step) * nstep;

			//extract profile on this device
			//coordinates now must be relative to subrect
			if (!mng(mGPU)->extract_profile_component_z(intersection_s + rect.s - prect_d[mGPU].s, intersection_e + rect.s - prect_d[mGPU].s, step, stencil, profile_gpu.get_cu_arr(mGPU))) return false;
		}
		else {

			//cannot extract profile in this sub-rect, so make sure to clear respective cu_arr
			if (profile_gpu.size(mGPU)) profile_gpu.clear(mGPU);
		}
	}

	return true;
}

//as above but only component z and pack in profile position too (for VAL3 floating types only). Return data as extracted component in profile_cpu. profile_cpu resized as needed.
template <typename VType, typename MType>
bool mcuVEC<VType, MType>::extract_profile_component_z(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		return mng(0)->extract_profile_component_z(start, end, step, stencil, profile_cpu);
	}

	cuReal3 intersection_s;
	cuReal3 intersection_e;

	//normalized vector step (profile direction)
	cuReal3 nstep = (end - start).normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//for given start and end points, find intersection values with each sub-rect  (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / step) * step) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / step) * step) * nstep;

			//extract profile on this device to auxiliary storage first
			//coordinates now must be relative to subrect
			if (!mng(mGPU)->extract_profile_component_z(intersection_s + rect.s - prect_d[mGPU].s, intersection_e + rect.s - prect_d[mGPU].s, step, stencil, profile_component_aux.get_cu_arr(mGPU))) return false;
		}
		else {

			//cannot extract profile in this sub-rect, so make sure to clear respective auxiliar storage cu_arr
			if (profile_component_aux.size(mGPU)) profile_component_aux.clear(mGPU);
		}
	}

	//now copy profile to profile_cpu, making sure it has the correct size
	size_t size = round((end - start).norm() / step) + 1;
	if (profile_cpu.size() != size) if (!malloc_vector(profile_cpu, size)) return false;
	profile_component_aux.copy_to_vector(profile_cpu);

	return true;
}

//as above but only component which has largest value for the first point and pack in profile position too (after stencil averaging) (for VAL3 floating types only). Return data in profile_gpu (which is resized as needed for each device)
template <typename VType, typename MType>
bool mcuVEC<VType, MType>::extract_profile_component_max(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, mcu_arr<cuReal2>& profile_gpu)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		return mng(0)->extract_profile_component_max(start, end, step, stencil, profile_gpu.get_cu_arr(0));
	}

	cuReal3 intersection_s;
	cuReal3 intersection_e;

	//normalized vector step (profile direction)
	cuReal3 nstep = (end - start).normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//for given start and end points, find intersection values with each sub-rect  (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / step) * step) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / step) * step) * nstep;

			//extract profile on this device
			//coordinates now must be relative to subrect
			if (!mng(mGPU)->extract_profile_component_max(intersection_s + rect.s - prect_d[mGPU].s, intersection_e + rect.s - prect_d[mGPU].s, step, stencil, profile_gpu.get_cu_arr(mGPU))) return false;
		}
		else {

			//cannot extract profile in this sub-rect, so make sure to clear respective cu_arr
			if (profile_gpu.size(mGPU)) profile_gpu.clear(mGPU);
		}
	}

	return true;
}

//as above but only component which has largest value for the first point and pack in profile position too (after stencil averaging) (for VAL3 floating types only). Return data in profile_cpu. profile_cpu resized as needed.
template <typename VType, typename MType>
bool mcuVEC<VType, MType>::extract_profile_component_max(cuReal3 start, cuReal3 end, cuBReal step, cuReal3 stencil, std::vector<DBL2>& profile_cpu)
{
	if (mGPU.get_num_devices() == 1) {

		//if using a single device then go directy to profile extraction
		return mng(0)->extract_profile_component_max(start, end, step, stencil, profile_cpu);
	}

	cuReal3 intersection_s;
	cuReal3 intersection_e;

	//normalized vector step (profile direction)
	cuReal3 nstep = (end - start).normalized();

	for (mGPU.device_begin(); mGPU != mGPU.device_end(); mGPU++) {

		//for given start and end points, find intersection values with each sub-rect  (absolute coordinates)
		if (prect_d[mGPU].intersection_test(start + rect.s, end + rect.s, &intersection_s, &intersection_e)) {

			//make relative to rect
			intersection_s -= rect.s;
			intersection_e -= rect.s;

			//snap coordinates to step, inside each sub-rect
			intersection_s = start + (ceil_epsilon((intersection_s - start).norm() / step) * step) * nstep;
			intersection_e = start + (floor_epsilon((intersection_e - start).norm() / step) * step) * nstep;

			//extract profile on this device to auxiliary storage first
			//coordinates now must be relative to subrect
			if (!mng(mGPU)->extract_profile_component_max(intersection_s + rect.s - prect_d[mGPU].s, intersection_e + rect.s - prect_d[mGPU].s, step, stencil, profile_component_aux.get_cu_arr(mGPU))) return false;
		}
		else {

			//cannot extract profile in this sub-rect, so make sure to clear respective auxiliar storage cu_arr
			if (profile_component_aux.size(mGPU)) profile_component_aux.clear(mGPU);
		}
	}

	//now copy profile to profile_cpu, making sure it has the correct size
	size_t size = round((end - start).norm() / step) + 1;
	if (profile_cpu.size() != size) if (!malloc_vector(profile_cpu, size)) return false;
	profile_component_aux.copy_to_vector(profile_cpu);

	return true;
}