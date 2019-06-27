#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_SDEMAG

#include "BorisCUDALib.h"

using namespace std;

//collect various kernel types in this struct used in demag kernel collections
//cuKerType is cu_obj managed, thus resides entirely in gpu memory
class cuKerType {

public:

	//2D off-diagonal for self-demag (only need real version)
	cuVEC<cuReal> K2D_odiag;

	//3D Complex
	cuVEC<cuReIm3> Kdiag_cmpl, Kodiag_cmpl;

	//3D Real
	cuVEC<cuReal3> Kdiag_real, Kodiag_real;

	//shift, source and destination cellsizes for which this kernel applies (not normalized)
	cuReal3 shift, h_src, h_dst;

	//The dimensions of the FFT space (In and Out in the kernel multiplication methods) - determines the size of kernels depending on their type. Set through AllocateKernels method.
	cuSZ3 N;

	//internal_demag == true -> use real kernels with reduced memory usage, 2D or 3D.
	//internal_demag == false -> use complex kernels 3D only.
	bool internal_demag;

	//if z shifted only then kernel multiplications and storage can be made more efficient for 2D. For 3D we can make storage more efficient.
	bool zshifted;

	//has this kernel been calculated and ready to use?
	bool kernel_calculated;

private:

	//------------------------------ AUXILIARY
	
	__host__ void FreeKernels(void);

public:
	//------------------------------ CONSTRUCTORS / DESTRUCTOR

	//void constructor
	__host__ void construct_cu_obj(void);

	//destructor
	__host__ void destruct_cu_obj(void);

	//------------------------------ MEMORY MANAGEMENT

	//allocate required gpu memory for kernels
	__host__ bool AllocateKernels(cuRect from_rect, cuRect this_rect, cuSZ3 N_cpu);

	//------------------------------ SETTERS

	__host__ void Set_Shift_and_Cellsizes(cuReal3 shift_, cuReal3 h_src_, cuReal3 h_dst_);

	//Flags
	__host__ void SetFlag_InternalDemag(bool status) { set_gpu_value(internal_demag, status); }

	__host__ void SetFlag_zShifted(bool status) { set_gpu_value(zshifted, status); }

	__host__ void SetFlag_Calculated(bool status) { set_gpu_value(kernel_calculated, status); }

	//Kernels
	__host__ void Set_K2D_odiag(std::vector<double>& K2D_odiag_cpu);

	template <typename cpuVEC>
	__host__ void Set_Kdiag_cmpl(cpuVEC& Kdiag_cmpl_cpu);

	template <typename cpuVEC>
	__host__ void Set_Kodiag_cmpl(cpuVEC& Kodiag_cmpl_cpu);

	template <typename cpuVEC>
	__host__ void Set_Kdiag_real(cpuVEC& Kdiag_real_cpu);

	template <typename cpuVEC>
	__host__ void Set_Kodiag_real(cpuVEC& Kodiag_real_cpu);

	//------------------------------ GETTERS

	__host__ cuReal3 Get_shift(void) { return get_gpu_value(shift); }

	__host__ cuReal3 Get_h_src(void) { return get_gpu_value(h_src); }

	__host__ cuReal3 Get_h_dst(void) { return get_gpu_value(h_dst); }

	__host__ cuSZ3 Get_N(void) { return get_gpu_value(N); }

	__host__ bool GetFlag_InternalDemag(void) { return get_gpu_value(internal_demag); }

	__host__ bool GetFlag_zShifted(void) { return get_gpu_value(zshifted); }

	__host__ bool GetFlag_Calculated(void) { return get_gpu_value(kernel_calculated); }

	//------------------------------ RUN-TIME KERNEL MULTIPLICATION

	//Kernel multiplication device methods
	//They take an input, output (dimensions N when not transposed)
	//Also take an idx : this is the idx in In and Out

	///////////////////////////////////////// 2D

	__device__ void cu_KernelMultiplication_2D_Self_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);
	__device__ void cu_KernelMultiplication_2D_Self_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);

	__device__ void cu_KernelMultiplication_2D_Self_transpose_xy_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);
	__device__ void cu_KernelMultiplication_2D_Self_transpose_xy_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);

	__device__ void cu_KernelMultiplication_2D_zShifted_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);
	__device__ void cu_KernelMultiplication_2D_zShifted_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);

	__device__ void cu_KernelMultiplication_2D_zShifted_transpose_xy_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);
	__device__ void cu_KernelMultiplication_2D_zShifted_transpose_xy_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);

	__device__ void cu_KernelMultiplication_2D_inversezShifted_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);
	__device__ void cu_KernelMultiplication_2D_inversezShifted_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);

	__device__ void cu_KernelMultiplication_2D_inversezShifted_transpose_xy_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);
	__device__ void cu_KernelMultiplication_2D_inversezShifted_transpose_xy_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);

	__device__ void cu_KernelMultiplication_2D_Regular_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);
	__device__ void cu_KernelMultiplication_2D_Regular_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);

	///////////////////////////////////////// 3D

	__device__ void cu_KernelMultiplication_3D_Self_transpose_xy_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);
	__device__ void cu_KernelMultiplication_3D_Self_transpose_xy_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);

	__device__ void cu_KernelMultiplication_3D_zShifted_transpose_xy_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);
	__device__ void cu_KernelMultiplication_3D_zShifted_transpose_xy_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);

	__device__ void cu_KernelMultiplication_3D_Regular_Set(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);
	__device__ void cu_KernelMultiplication_3D_Regular_Add(int idx, cuComplex* Inx, cuComplex* Iny, cuComplex* Inz, cuComplex* Outx, cuComplex* Outy, cuComplex* Outz);
};

//------------------------------ AUXILIARY

__host__ inline void cuKerType::FreeKernels(void)
{
	K2D_odiag.clear();

	Kdiag_cmpl.clear();
	Kodiag_cmpl.clear();

	Kdiag_real.clear();
	Kodiag_real.clear();
}

//------------------------------ CONSTRUCTORS / DESTRUCTOR

//void constructor
__host__ inline void cuKerType::construct_cu_obj(void)
{
	K2D_odiag.construct_cu_obj();

	Kdiag_cmpl.construct_cu_obj();
	Kodiag_cmpl.construct_cu_obj();

	Kdiag_real.construct_cu_obj();
	Kodiag_real.construct_cu_obj();

	set_gpu_value(internal_demag, false);
	set_gpu_value(zshifted, false);
	set_gpu_value(kernel_calculated, false);

	set_gpu_value(shift, cuReal3());
	set_gpu_value(h_src, cuReal3());
	set_gpu_value(h_dst, cuReal3());
}

//destructor
__host__ inline void cuKerType::destruct_cu_obj(void)
{
	K2D_odiag.destruct_cu_obj();

	Kdiag_cmpl.destruct_cu_obj();
	Kodiag_cmpl.destruct_cu_obj();

	Kdiag_real.destruct_cu_obj();
	Kodiag_real.destruct_cu_obj();
}

//------------------------------ MEMORY MANAGEMENT

//allocate required gpu memory for kernels
__host__ inline bool cuKerType::AllocateKernels(cuRect from_rect, cuRect this_rect, cuSZ3 N_cpu)
{
	FreeKernels();

	//Set N in gpu memory for later use with kernel multiplication methods
	set_gpu_value(N, N_cpu);

	//check if internal demag
	if (from_rect == this_rect) {

		set_gpu_value(internal_demag, true);
		set_gpu_value(zshifted, false);

		if (N_cpu.z == 1) {

			if (!Kdiag_real.resize(cuSZ3(N_cpu.x / 2 + 1, N_cpu.y / 2 + 1, 1))) return false;
			if (!K2D_odiag.resize(cuSZ3(N_cpu.x / 2 + 1, N_cpu.y / 2 + 1, 1))) return false;
		}
		else {

			if (!Kdiag_real.resize(cuSZ3(N_cpu.x / 2 + 1, N_cpu.y / 2 + 1, N_cpu.z / 2 + 1))) return false;
			if (!Kodiag_real.resize(cuSZ3(N_cpu.x / 2 + 1, N_cpu.y / 2 + 1, N_cpu.z / 2 + 1))) return false;
		}
	}
	//not internal demag
	else {

		set_gpu_value(internal_demag, false);

		//see if we can use z-shifted kernels instead of full kernels
		cuReal3 shift = from_rect.s - this_rect.s;

		if (cuIsZ(shift.x) && cuIsZ(shift.y)) {

			set_gpu_value(zshifted, true);

			if (N_cpu.z == 1) {

				//z shifted (for 2D) : can use real kernels of reduced dimensions
				//Kxx, Kyy, Kzz : real
				//Kxy : real
				//Kxz, Kyz : imaginary (so must adjust for i multiplication when multiplying kernels)

				if (!Kdiag_real.resize(cuSZ3(N_cpu.x / 2 + 1, N_cpu.y / 2 + 1, N_cpu.z))) return false;
				if (!Kodiag_real.resize(cuSZ3(N_cpu.x / 2 + 1, N_cpu.y / 2 + 1, N_cpu.z))) return false;
			}
			else {

				//z shifted for 3D : can use kernels of reduced dimensions but must be complex
				//
				//Kxx : y - symmetrical (+), z - Re part symmetrical (+), Im part asymmetrical (-)
				//Kyy : y - symmetrical (+), z - Re part symmetrical (+), Im part asymmetrical (-)
				//Kzz : y - symmetrical (+), z - Re part symmetrical (+), Im part asymmetrical (-)
				//
				//Kxy : y - asymmetrical (-), z - Re part symmetrical  (+), Im part asymmetrical (-)
				//Kxz : y - symmetrical  (+), z - Re part asymmetrical (-), Im part symmetrical  (+)
				//Kyz : y - asymmetrical (-), z - Re part asymmetrical (-), Im part symmetrical  (+)

				if (!Kdiag_cmpl.resize(cuSZ3(N_cpu.x / 2 + 1, N_cpu.y / 2 + 1, N_cpu.z / 2 + 1))) return false;
				if (!Kodiag_cmpl.resize(cuSZ3(N_cpu.x / 2 + 1, N_cpu.y / 2 + 1, N_cpu.z / 2 + 1))) return false;
			}
		}

		else {

			//full complex kernels

			set_gpu_value(zshifted, false);

			if (!Kdiag_cmpl.resize(cuSZ3(N_cpu.x / 2 + 1, N_cpu.y, N_cpu.z))) return false;
			if (!Kodiag_cmpl.resize(cuSZ3(N_cpu.x / 2 + 1, N_cpu.y, N_cpu.z))) return false;
		}
	}

	return true;
}

//------------------------------ SETTERS

__host__ inline void cuKerType::Set_Shift_and_Cellsizes(cuReal3 shift_, cuReal3 h_src_, cuReal3 h_dst_)
{
	set_gpu_value(shift, shift_);
	set_gpu_value(h_src, h_src_);
	set_gpu_value(h_dst, h_dst_);
}

//Kernels
__host__ inline void cuKerType::Set_K2D_odiag(std::vector<double>& K2D_odiag_cpu)
{
	K2D_odiag.copy_from_vector(K2D_odiag_cpu);
}

template <typename cpuVEC>
__host__ void cuKerType::Set_Kdiag_cmpl(cpuVEC& Kdiag_cmpl_cpu)
{
	Kdiag_cmpl.copy_from_cpuvec(Kdiag_cmpl_cpu);
}

template <typename cpuVEC>
__host__ void cuKerType::Set_Kodiag_cmpl(cpuVEC& Kodiag_cmpl_cpu)
{
	Kodiag_cmpl.copy_from_cpuvec(Kodiag_cmpl_cpu);
}

template <typename cpuVEC>
__host__ void cuKerType::Set_Kdiag_real(cpuVEC& Kdiag_real_cpu)
{
	Kdiag_real.copy_from_cpuvec(Kdiag_real_cpu);
}

template <typename cpuVEC>
__host__ void cuKerType::Set_Kodiag_real(cpuVEC& Kodiag_real_cpu)
{
	Kodiag_real.copy_from_cpuvec(Kodiag_real_cpu);
}

#endif

#endif