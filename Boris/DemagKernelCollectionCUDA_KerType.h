#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SDEMAG

#include "BorisCUDALib.h"



//collect various kernel types in this struct used in demag kernel collections
//cuKerType is cu_obj managed, thus resides entirely in gpu memory
class cuKerType {

public:

	//2D off-diagonal for self-demag (only need real version)
	cuVEC<cuBReal> K2D_odiag;

	//3D Complex
	cuVEC<cuReIm3> Kdiag_cmpl, Kodiag_cmpl;

	//3D Real
	cuVEC<cuReal3> Kdiag_real, Kodiag_real;

	//shift, source and destination cellsizes for which this kernel applies (not normalized)
	cuReal3 shift, h_src, h_dst;

	//The dimensions of the FFT space (In and Out in the kernel multiplication methods) - determines the size of kernels depending on their type. Set through AllocateKernels method.
	cuSZ3 N;

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
	__host__ void SetFlag_zShifted(bool status) { set_gpu_value(zshifted, status); }

	__host__ void SetFlag_Calculated(bool status) { set_gpu_value(kernel_calculated, status); }

	//Kernels
	__host__ void Set_K2D_odiag(std::vector<double>& K2D_odiag_cpu);
	//all-GPU version
	void Set_K2D_odiag(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int transpose_xy);

	template <typename cpuVEC>
	__host__ void Set_Kdiag_cmpl(cpuVEC& Kdiag_cmpl_cpu);
	//all-GPU version (full size : (N.x/2 + 1) * N.y * N.z
	void Set_Kdiag_cmpl(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int component, int transpose_xy);
	//all-GPU version (reduced size : (N.x/2 + 1) * (N.y/2 + 1) * (N.z/2 + 1)
	void Set_Kdiag_cmpl_reduced(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int component, int transpose_xy);

	template <typename cpuVEC>
	__host__ void Set_Kodiag_cmpl(cpuVEC& Kodiag_cmpl_cpu);
	//all-GPU version (full size : (N.x/2 + 1) * N.y * N.z
	void Set_Kodiag_cmpl(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int component, int transpose_xy);
	//all-GPU version (reduced size : (N.x/2 + 1) * (N.y/2 + 1) * (N.z/2 + 1)
	void Set_Kodiag_cmpl_reduced(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int component, int transpose_xy);

	template <typename cpuVEC>
	__host__ void Set_Kdiag_real(cpuVEC& Kdiag_real_cpu);
	//all-GPU version
	void Set_Kdiag_real(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int component, int transpose_xy);

	template <typename cpuVEC>
	__host__ void Set_Kodiag_real(cpuVEC& Kodiag_real_cpu);
	//all-GPU version  (2D : -Im, Re, Im parts)
	void Set_Kodiag_real_2D(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int component, int transpose_xy);
	//all-GPU version  (3D : -Re, -Im, -Im parts)
	void Set_Kodiag_real_3D(size_t size, cu_arr<cufftDoubleComplex>& cuOut, int component, int transpose_xy);

	//------------------------------ GETTERS

	__host__ cuReal3 Get_shift(void) { return get_gpu_value(shift); }

	__host__ cuReal3 Get_h_src(void) { return get_gpu_value(h_src); }

	__host__ cuReal3 Get_h_dst(void) { return get_gpu_value(h_dst); }

	__host__ cuSZ3 Get_N(void) { return get_gpu_value(N); }

	__host__ bool GetFlag_zShifted(void) { return get_gpu_value(zshifted); }

	__host__ bool GetFlag_Calculated(void) { return get_gpu_value(kernel_calculated); }
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
				//Kxx : y - symmetrical (+), z - Re part symmetrical (+), Im part inv. symmetric (-)
				//Kyy : y - symmetrical (+), z - Re part symmetrical (+), Im part inv. symmetric (-)
				//Kzz : y - symmetrical (+), z - Re part symmetrical (+), Im part inv. symmetric (-)
				//
				//Kxy : y - inv. symmetric (-), z - Re part symmetrical  (+), Im part inv. symmetric (-)
				//Kxz : y - symmetrical  (+), z - Re part inv. symmetric (-), Im part symmetrical  (+)
				//Kyz : y - inv. symmetric (-), z - Re part inv. symmetric (-), Im part symmetrical  (+)

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
