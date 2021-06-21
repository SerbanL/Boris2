#include "stdafx.h"
#include "DemagKernelCollectionCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_SDEMAG

#include "DemagTFuncCUDA.h"

#include "DemagTFunc.h"

//-------------------------- KERNEL CALCULATION

//2D kernels (Kdiag_real, and K2D_odiag, with full use of kernel symmetries)
BError DemagKernelCollectionCUDA::Calculate_Demag_Kernels_2D_Self_onGPU(int index)
{
	BError error(__FUNCTION__);

	//-------------------------------------------

	cu_arr<cufftDoubleReal> cuInx, cuIny, cuInz;
	if (!cuInx.resize(N.x*N.y)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuIny.resize(N.x*N.y)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuInz.resize(N.x*N.y)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	cu_arr<cufftDoubleComplex> cuOut;
	if (!cuOut.resize((N.x / 2 + 1)*N.y)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	//-------------- CALCULATE DEMAG TENSOR

	//Demag tensor components
	//
	// D11 D12 0
	// D12 D22 0
	// 0   0   D33

	DemagTFuncCUDA dtf_gpu;

	if (pbc_images.IsNull()) {

		if (!dtf_gpu.CalcDiagTens2D(cuInx, cuIny, cuInz, n, N, h / h_max)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}
	else {

		if (!dtf_gpu.CalcDiagTens2D_PBC(
			cuInx, cuIny, cuInz, N, h / h_max,
			true, ASYMPTOTIC_DISTANCE, pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}

	//-------------- SETUP FFT

	cufftHandle planTens2D_x, planTens2D_y;

	int embed[1] = { 0 };
	int ndims_x[1] = { (int)N.x };
	int ndims_y[1] = { (int)N.y };

	//Forward fft along x direction (out-of-place):
	cufftResult cuffterr = cufftPlanMany(&planTens2D_x, 1, ndims_x,
		embed, 1, N.x,
		embed, 1, (N.x / 2 + 1),
		CUFFT_D2Z, (int)N.y);

	//Forward fft along y direction (out-of-place):
	cuffterr = cufftPlanMany(&planTens2D_y, 1, ndims_y,
		embed, (N.x / 2 + 1), 1,
		embed, (N.x / 2 + 1), 1,
		CUFFT_D2Z, (int)(N.x / 2 + 1));

	//-------------- FFT REAL TENSOR INTO REAL KERNELS

	auto fftcomponent = [&](cu_arr<cufftDoubleReal>& cuIn, int component) -> void
	{
		//Now FFT along x from input to output
		cufftExecD2Z(planTens2D_x, cuIn, cuOut);

		//extract real parts from cuOut and place in cuIn
		cuOut_to_cuIn_Re((N.x / 2 + 1) * N.y, cuIn, cuOut);

		//Now FFT along y from input to output
		cufftExecD2Z(planTens2D_y, cuIn, cuOut);

		//extract real parts from cuOut and place in kernel
		(*kernels[index])()->Set_Kdiag_real((N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1), cuOut, component, transpose_xy);
	};

	//diagonal components
	fftcomponent(cuInx, 1);
	fftcomponent(cuIny, 2);
	fftcomponent(cuInz, 3);

	//off-diagonal component
	if (pbc_images.IsNull()) {

		if (!dtf_gpu.CalcOffDiagTens2D(cuInx, n, N, h / h_max)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}
	else {

		if (!dtf_gpu.CalcOffDiagTens2D_PBC(
			cuInx, N, h / h_max,
			true, ASYMPTOTIC_DISTANCE, pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}

	//Now FFT along x from input to output
	cufftExecD2Z(planTens2D_x, cuInx, cuOut);

	//extract imaginary parts from cuOut and place in cuIn
	cuOut_to_cuIn_Im((N.x / 2 + 1) * N.y, cuInx, cuOut);

	//Now FFT along y from input to output
	cufftExecD2Z(planTens2D_y, cuInx, cuOut);

	//extract imaginary part from cuOut and place in kernel, also changing sign
	(*kernels[index])()->Set_K2D_odiag((N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1), cuOut, transpose_xy);

	//-------------- CLEANUP

	cufftDestroy(planTens2D_x);
	cufftDestroy(planTens2D_y);

	//Done
	return error;
}

//2D layers, z shift only : Kernels can be stored as real with use of kernel symmetries. Kxx, Kyy, Kzz, Kxy real, Kxz, Kyz imaginary
BError DemagKernelCollectionCUDA::Calculate_Demag_Kernels_2D_zShifted_onGPU(int index)
{
	BError error(__FUNCTION__);

	//-------------------------------------------

	cu_arr<cufftDoubleReal> cuInx, cuIny, cuInz;
	if (!cuInx.resize(N.x*N.y)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuIny.resize(N.x*N.y)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuInz.resize(N.x*N.y)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	cu_arr<cufftDoubleComplex> cuOut;
	if (!cuOut.resize((N.x / 2 + 1)*N.y)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	//-------------- CALCULATE DEMAG TENSOR

	//Demag tensor components
	//
	// D11 D12 0
	// D12 D22 0
	// 0   0   D33

	DemagTFuncCUDA dtf_gpu;

	if (pbc_images.IsNull()) {

		if (!dtf_gpu.CalcDiagTens2D_Shifted_Irregular(
			cuInx, cuIny, cuInz, 
			n, N, (*kernels[index])()->Get_h_src() / h_max, (*kernels[index])()->Get_h_dst() / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}
	else {
		
		if (!dtf_gpu.CalcDiagTens2D_Shifted_Irregular_PBC(
			cuInx, cuIny, cuInz, 
			N, (*kernels[index])()->Get_h_src() / h_max, (*kernels[index])()->Get_h_dst() / h_max, (*kernels[index])()->Get_shift() / h_max,
			true, ASYMPTOTIC_DISTANCE, pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}

	//-------------- SETUP FFT

	cufftHandle planTens2D_x, planTens2D_y;

	int embed[1] = { 0 };
	int ndims_x[1] = { (int)N.x };
	int ndims_y[1] = { (int)N.y };

	//Forward fft along x direction (out-of-place):
	cufftResult cuffterr = cufftPlanMany(&planTens2D_x, 1, ndims_x,
		embed, 1, N.x,
		embed, 1, (N.x / 2 + 1),
		CUFFT_D2Z, (int)N.y);

	//Forward fft along y direction (out-of-place):
	cuffterr = cufftPlanMany(&planTens2D_y, 1, ndims_y,
		embed, (N.x / 2 + 1), 1,
		embed, (N.x / 2 + 1), 1,
		CUFFT_D2Z, (int)(N.x / 2 + 1));

	//-------------- FFT REAL TENSOR INTO REAL KERNELS

	auto fftdiagcomponent = [&](cu_arr<cufftDoubleReal>& cuIn, int component) -> void
	{
		//Now FFT along x from input to output
		cufftExecD2Z(planTens2D_x, cuIn, cuOut);

		//extract real parts from cuOut and place in cuIn
		cuOut_to_cuIn_Re((N.x / 2 + 1) * N.y, cuIn, cuOut);

		//Now FFT along y from input to output
		cufftExecD2Z(planTens2D_y, cuIn, cuOut);

		//extract real parts from cuOut and place in kernel
		(*kernels[index])()->Set_Kdiag_real((N.x / 2 + 1) * (N.y / 2 + 1), cuOut, component, transpose_xy);
	};

	//diagonal components
	fftdiagcomponent(cuInx, 1);
	fftdiagcomponent(cuIny, 2);
	fftdiagcomponent(cuInz, 3);

	auto fftodiagcomponent = [&](cu_arr<cufftDoubleReal>& cuIn, int component) -> void
	{
		//Now FFT along x from input to output
		cufftExecD2Z(planTens2D_x, cuIn, cuOut);

		//extract Im, Im, Re parts from cuOut and place in cuIn
		if (component == 1) cuOut_to_cuIn_Im((N.x / 2 + 1) * N.y, cuIn, cuOut);
		else if (component == 2) cuOut_to_cuIn_Im((N.x / 2 + 1) * N.y, cuIn, cuOut);
		else if (component == 3) cuOut_to_cuIn_Re((N.x / 2 + 1) * N.y, cuIn, cuOut);

		//Now FFT along y from input to output
		cufftExecD2Z(planTens2D_y, cuIn, cuOut);

		//extract -Im, Re, Im parts from cuOut and place in kernel
		(*kernels[index])()->Set_Kodiag_real_2D((N.x / 2 + 1) * (N.y / 2 + 1), cuOut, component, transpose_xy);
	};

	//off-diagonal component
	if (pbc_images.IsNull()) {

		if (!dtf_gpu.CalcOffDiagTens2D_Shifted_Irregular(
			cuInx, cuIny, cuInz, 
			n, N, (*kernels[index])()->Get_h_src() / h_max, (*kernels[index])()->Get_h_dst() / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}
	else {

		if (!dtf_gpu.CalcOffDiagTens2D_Shifted_Irregular_PBC(
			cuInx, cuIny, cuInz, 
			N, (*kernels[index])()->Get_h_src() / h_max, (*kernels[index])()->Get_h_dst() / h_max, (*kernels[index])()->Get_shift() / h_max,
			true, ASYMPTOTIC_DISTANCE, pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}

	//off-diagonal components
	fftodiagcomponent(cuInx, 1);
	fftodiagcomponent(cuIny, 2);
	fftodiagcomponent(cuInz, 3);

	//-------------- CLEANUP

	cufftDestroy(planTens2D_x);
	cufftDestroy(planTens2D_y);

	//Done
	return error;
}

//2D layers, complex kernels most general case (Kdiag_cmpl, and Kodiag_cmpl, without any kernel symmetries)
BError DemagKernelCollectionCUDA::Calculate_Demag_Kernels_2D_Complex_Full_onGPU(int index)
{
	BError error(__FUNCTION__);

	//-------------------------------------------

	cu_arr<cufftDoubleReal> cuInx, cuIny, cuInz;
	if (!cuInx.resize(N.x*N.y)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuIny.resize(N.x*N.y)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuInz.resize(N.x*N.y)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	cu_arr<cufftDoubleComplex> cuOut;
	if (!cuOut.resize((N.x / 2 + 1)*N.y)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	//-------------- CALCULATE DEMAG TENSOR

	//Demag tensor components
	//
	// D11 D12 D13
	// D12 D22 D23
	// D13 D23 D33

	DemagTFuncCUDA dtf_gpu;

	if (pbc_images.IsNull()) {

		if (!dtf_gpu.CalcDiagTens2D_Shifted_Irregular(
			cuInx, cuIny, cuInz,
			n, N, (*kernels[index])()->Get_h_src() / h_max, (*kernels[index])()->Get_h_dst() / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}
	else {

		if (!dtf_gpu.CalcDiagTens2D_Shifted_Irregular_PBC(
			cuInx, cuIny, cuInz,
			N, (*kernels[index])()->Get_h_src() / h_max, (*kernels[index])()->Get_h_dst() / h_max, (*kernels[index])()->Get_shift() / h_max,
			true, ASYMPTOTIC_DISTANCE, pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}

	//-------------- SETUP FFT

	cufftHandle planTens2D_x, planTens2D_y;

	int embed[1] = { 0 };
	int ndims_x[1] = { (int)N.x };
	int ndims_y[1] = { (int)N.y };

	//Forward fft along x direction (out-of-place):
	cufftResult cuffterr = cufftPlanMany(&planTens2D_x, 1, ndims_x,
		embed, 1, N.x,
		embed, 1, (N.x / 2 + 1),
		CUFFT_D2Z, (int)N.y);

	//Forward fft along y direction (out-of-place):
	cuffterr = cufftPlanMany(&planTens2D_y, 1, ndims_y,
		embed, (N.x / 2 + 1), 1,
		embed, (N.x / 2 + 1), 1,
		CUFFT_Z2Z, (int)(N.x / 2 + 1));

	//-------------- FFT REAL TENSOR INTO REAL KERNELS

	auto fftdiagcomponent = [&](cu_arr<cufftDoubleReal>& cuIn, int component, bool odiag) -> void
	{
		//Now FFT along x from input to output
		cufftExecD2Z(planTens2D_x, cuIn, cuOut);

		//extract real parts from cuOut and place in cuIn
		//cuOut_to_cuIn_Re((N.x / 2 + 1) * N.y, cuIn, cuOut);

		//Now FFT along y from input to output
		cufftExecZ2Z(planTens2D_y, cuOut, cuOut, CUFFT_FORWARD);

		//extract real parts from cuOut and place in kernel
		if (!odiag) (*kernels[index])()->Set_Kdiag_cmpl((N.x / 2 + 1) * N.y, cuOut, component, transpose_xy);
		else (*kernels[index])()->Set_Kodiag_cmpl((N.x / 2 + 1) * N.y, cuOut, component, transpose_xy);
	};

	//diagonal components
	fftdiagcomponent(cuInx, 1, false);
	fftdiagcomponent(cuIny, 2, false);
	fftdiagcomponent(cuInz, 3, false);

	//off-diagonal component
	if (pbc_images.IsNull()) {

		if (!dtf_gpu.CalcOffDiagTens2D_Shifted_Irregular(
			cuInx, cuIny, cuInz,
			n, N, (*kernels[index])()->Get_h_src() / h_max, (*kernels[index])()->Get_h_dst() / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}
	else {

		if (!dtf_gpu.CalcOffDiagTens2D_Shifted_Irregular_PBC(
			cuInx, cuIny, cuInz,
			N, (*kernels[index])()->Get_h_src() / h_max, (*kernels[index])()->Get_h_dst() / h_max, (*kernels[index])()->Get_shift() / h_max,
			true, ASYMPTOTIC_DISTANCE, pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}

	//off-diagonal components
	fftdiagcomponent(cuInx, 1, true);
	fftdiagcomponent(cuIny, 2, true);
	fftdiagcomponent(cuInz, 3, true);

	//-------------- CLEANUP

	cufftDestroy(planTens2D_x);
	cufftDestroy(planTens2D_y);

	//Done
	return error;
}

//3D real kernels (Kdiag_real, and Kodiag_real, with full use of kernel symmetries)
BError DemagKernelCollectionCUDA::Calculate_Demag_Kernels_3D_Self_onGPU(int index)
{
	BError error(__FUNCTION__);

	//-------------------------------------------

	cu_arr<cufftDoubleReal> cuInx, cuIny, cuInz;
	if (!cuInx.resize(N.x*N.y*N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuIny.resize(N.x*N.y*N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuInz.resize(N.x*N.y*N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	cu_arr<cufftDoubleComplex> cuOut;
	if (!cuOut.resize((N.x / 2 + 1)*N.y*N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	//-------------- DEMAG TENSOR

	//Demag tensor components
	//
	// D11 D12 D13
	// D12 D22 D23
	// D13 D23 D33

	DemagTFuncCUDA dtf_gpu;

	if (pbc_images.IsNull()) {

		if (!dtf_gpu.CalcDiagTens3D(cuInx, cuIny, cuInz, n, N, h / h_max)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}
	else {

		if (!dtf_gpu.CalcDiagTens3D_PBC(
			cuInx, cuIny, cuInz, N, h / h_max,
			true, ASYMPTOTIC_DISTANCE, pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}

	//-------------- SETUP FFT

	cufftHandle planTens3D_x, planTens3D_y, planTens3D_z;

	int embed[1] = { 0 };
	int ndims_x[1] = { (int)N.x };
	int ndims_y[1] = { (int)N.y };
	int ndims_z[1] = { (int)N.z };

	//Forward fft along x direction (out-of-place):
	cufftResult cuffterr = cufftPlanMany(&planTens3D_x, 1, ndims_x,
		embed, 1, N.x,
		embed, 1, (N.x / 2 + 1),
		CUFFT_D2Z, (int)N.y*N.z);

	//Forward fft along y direction (out-of-place):
	cuffterr = cufftPlanMany(&planTens3D_y, 1, ndims_y,
		embed, (N.x / 2 + 1), 1,
		embed, (N.x / 2 + 1), 1,
		CUFFT_D2Z, (int)(N.x / 2 + 1));

	//Forward fft along z direction (out-of-place):
	cuffterr = cufftPlanMany(&planTens3D_z, 1, ndims_z,
		embed, (N.x / 2 + 1) * (N.y / 2 + 1), 1,
		embed, (N.x / 2 + 1) * (N.y / 2 + 1), 1,
		CUFFT_D2Z, (int)(N.x / 2 + 1)*(N.y / 2 + 1));

	//-------------- FFT REAL TENSOR INTO REAL KERNELS

	auto fftdiagcomponent = [&](cu_arr<cufftDoubleReal>& cuIn, int component) -> void
	{
		//Now FFT along x from input to output
		cufftExecD2Z(planTens3D_x, cuIn, cuOut);

		//extract real parts from cuOut and place in cuIn
		cuOut_to_cuIn_Re((N.x / 2 + 1) * N.y * N.z, cuIn, cuOut);

		//Now FFT along y from input to output
		for (int k = 0; k < N.z; k++) {

			cufftExecD2Z(planTens3D_y, cuIn + k * (N.x / 2 + 1) * N.y, cuOut + k * (N.x / 2 + 1) * (N.y / 2 + 1));
		}

		//extract real parts from cuOut and place in cuIn
		cuOut_to_cuIn_Re((N.x / 2 + 1) * (N.y / 2 + 1) * N.z, cuIn, cuOut);

		//Now FFT along z from input to output
		cufftExecD2Z(planTens3D_z, cuIn, cuOut);

		//extract real parts from cuOut and place in kernel
		(*kernels[index])()->Set_Kdiag_real((N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1), cuOut, component, transpose_xy);
	};

	auto fftodiagcomponent = [&](cu_arr<cufftDoubleReal>& cuIn, int component) -> void
	{
		//Now FFT along x from input to output
		cufftExecD2Z(planTens3D_x, cuIn, cuOut);

		//extract real or imaginary parts from cuOut and place in cuIn, depending on the computed component
		if (component == 1) cuOut_to_cuIn_Im((N.x / 2 + 1) * N.y * N.z, cuIn, cuOut);
		else if (component == 2) cuOut_to_cuIn_Im((N.x / 2 + 1) * N.y * N.z, cuIn, cuOut);
		else if (component == 3) cuOut_to_cuIn_Re((N.x / 2 + 1) * N.y * N.z, cuIn, cuOut);

		//Now FFT along y from input to output
		for (int k = 0; k < N.z; k++) {

			cufftExecD2Z(planTens3D_y, cuIn + k * (N.x / 2 + 1) * N.y, cuOut + k * (N.x / 2 + 1) * (N.y / 2 + 1));
		}

		//extract real or imaginary parts from cuOut and place in cuIn, depending on the computed component
		if (component == 1) cuOut_to_cuIn_Im((N.x / 2 + 1) * (N.y / 2 + 1) * N.z, cuIn, cuOut);
		else if (component == 2) cuOut_to_cuIn_Re((N.x / 2 + 1) * (N.y / 2 + 1) * N.z, cuIn, cuOut);
		else if (component == 3) cuOut_to_cuIn_Im((N.x / 2 + 1) * (N.y / 2 + 1) * N.z, cuIn, cuOut);

		//Now FFT along z from input to output
		cufftExecD2Z(planTens3D_z, cuIn, cuOut);

		//extract -Re, -Im, -Im parts from cuOut and place in kernel
		(*kernels[index])()->Set_Kodiag_real_3D((N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1), cuOut, component, transpose_xy);
	};

	//-------------- CALCULATE DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	//diagonal components

	fftdiagcomponent(cuInx, 1);
	fftdiagcomponent(cuIny, 2);
	fftdiagcomponent(cuInz, 3);

	//-------------- CALCULATE OFF-DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	if (pbc_images.IsNull()) {

		//no pbcs in any dimension -> use these
		if (!dtf_gpu.CalcOffDiagTens3D(cuInx, cuIny, cuInz, n, N, h / h_max)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}
	else {

		//pbcs used in at least one dimension
		if (!dtf_gpu.CalcOffDiagTens3D_PBC(
			cuInx, cuIny, cuInz, N, h / h_max,
			true, ASYMPTOTIC_DISTANCE, pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}

	//off-diagonal components
	fftodiagcomponent(cuInx, 1);
	fftodiagcomponent(cuIny, 2);
	fftodiagcomponent(cuInz, 3);

	//-------------- CLEANUP

	cufftDestroy(planTens3D_x);
	cufftDestroy(planTens3D_y);
	cufftDestroy(planTens3D_z);

	return error;
}

//3D layers, z shift only : Kernels can be stored with use of kernel symmetries (but still complex).
BError DemagKernelCollectionCUDA::Calculate_Demag_Kernels_3D_zShifted_onGPU(int index)
{
	BError error(__FUNCTION__);

	//-------------------------------------------

	cu_arr<cufftDoubleReal> cuInx, cuIny, cuInz;
	if (!cuInx.resize(N.x*N.y*N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuIny.resize(N.x*N.y*N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuInz.resize(N.x*N.y*N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	cu_arr<cufftDoubleComplex> cuOut;
	if (!cuOut.resize((N.x / 2 + 1)*N.y*N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	//-------------- DEMAG TENSOR

	//Demag tensor components
	//
	// D11 D12 D13
	// D12 D22 D23
	// D13 D23 D33

	DemagTFuncCUDA dtf_gpu;

	if (pbc_images.IsNull()) {

		if (!dtf_gpu.CalcDiagTens3D_Shifted(cuInx, cuIny, cuInz, n, N, h / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}
	else {

		if (!dtf_gpu.CalcDiagTens3D_Shifted_PBC(
			cuInx, cuIny, cuInz, N, h / h_max, (*kernels[index])()->Get_shift() / h_max,
			true, ASYMPTOTIC_DISTANCE, pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}

	//-------------- SETUP FFT

	cufftHandle planTens3D_x, planTens3D_y, planTens3D_z;

	int embed[1] = { 0 };
	int ndims_x[1] = { (int)N.x };
	int ndims_y[1] = { (int)N.y };
	int ndims_z[1] = { (int)N.z };

	//Forward fft along x direction (out-of-place):
	cufftResult cuffterr = cufftPlanMany(&planTens3D_x, 1, ndims_x,
		embed, 1, N.x,
		embed, 1, (N.x / 2 + 1),
		CUFFT_D2Z, (int)N.y*N.z);

	//Forward fft along y direction (out-of-place):
	cuffterr = cufftPlanMany(&planTens3D_y, 1, ndims_y,
		embed, (N.x / 2 + 1), 1,
		embed, (N.x / 2 + 1), 1,
		CUFFT_Z2Z, (int)(N.x / 2 + 1));

	//Forward fft along z direction (out-of-place):
	cuffterr = cufftPlanMany(&planTens3D_z, 1, ndims_z,
		embed, (N.x / 2 + 1) * (N.y / 2 + 1), 1,
		embed, (N.x / 2 + 1) * (N.y / 2 + 1), 1,
		CUFFT_Z2Z, (int)(N.x / 2 + 1)*(N.y / 2 + 1));

	//-------------- FFT REAL TENSOR INTO REAL KERNELS

	auto fftcomponent = [&](cu_arr<cufftDoubleReal>& cuIn, int component, bool odiag) -> void
	{
		//Now FFT along x from input to output
		cufftExecD2Z(planTens3D_x, cuIn, cuOut);

		//Now FFT along y from input to output
		for (int k = 0; k < N.z; k++) {

			cufftExecZ2Z(planTens3D_y, cuOut + k * (N.x / 2 + 1) * N.y, cuOut + k * (N.x / 2 + 1) * (N.y / 2 + 1), CUFFT_FORWARD);
		}

		//Now FFT along z from input to output
		cufftExecZ2Z(planTens3D_z, cuOut, cuOut, CUFFT_FORWARD);

		//extract and place in kernel
		if (!odiag) (*kernels[index])()->Set_Kdiag_cmpl_reduced((N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1), cuOut, component, transpose_xy);
		else (*kernels[index])()->Set_Kodiag_cmpl_reduced((N.x / 2 + 1) * (N.y / 2 + 1) * (N.z / 2 + 1), cuOut, component, transpose_xy);
	};

	//-------------- CALCULATE DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	//diagonal components

	fftcomponent(cuInx, 1, false);
	fftcomponent(cuIny, 2, false);
	fftcomponent(cuInz, 3, false);

	//-------------- CALCULATE OFF-DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	if (pbc_images.IsNull()) {

		//no pbcs in any dimension -> use these
		if (!dtf_gpu.CalcOffDiagTens3D_Shifted(cuInx, cuIny, cuInz, n, N, h / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}
	else {

		//pbcs used in at least one dimension
		if (!dtf_gpu.CalcOffDiagTens3D_Shifted_PBC(
			cuInx, cuIny, cuInz, N, h / h_max, (*kernels[index])()->Get_shift() / h_max,
			true, ASYMPTOTIC_DISTANCE, pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}

	//off-diagonal components
	fftcomponent(cuInx, 1, true);
	fftcomponent(cuIny, 2, true);
	fftcomponent(cuInz, 3, true);

	//-------------- CLEANUP

	cufftDestroy(planTens3D_x);
	cufftDestroy(planTens3D_y);
	cufftDestroy(planTens3D_z);

	return error;
}

//3D complex kernels (Kdiag_cmpl, and Kodiag_cmpl, without any kernel symmetries)
BError DemagKernelCollectionCUDA::Calculate_Demag_Kernels_3D_Complex_Full_onGPU(int index)
{
	BError error(__FUNCTION__);

	//-------------------------------------------

	cu_arr<cufftDoubleReal> cuInx, cuIny, cuInz;
	if (!cuInx.resize(N.x*N.y*N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuIny.resize(N.x*N.y*N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);
	if (!cuInz.resize(N.x*N.y*N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	cu_arr<cufftDoubleComplex> cuOut;
	if (!cuOut.resize((N.x / 2 + 1)*N.y*N.z)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

	//-------------- DEMAG TENSOR

	//Demag tensor components
	//
	// D11 D12 D13
	// D12 D22 D23
	// D13 D23 D33

	DemagTFuncCUDA dtf_gpu;

	if (pbc_images.IsNull()) {

		if (!dtf_gpu.CalcDiagTens3D_Shifted(cuInx, cuIny, cuInz, n, N, h / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}
	else {

		if (!dtf_gpu.CalcDiagTens3D_Shifted_PBC(
			cuInx, cuIny, cuInz, N, h / h_max, (*kernels[index])()->Get_shift() / h_max,
			true, ASYMPTOTIC_DISTANCE, pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}

	//-------------- SETUP FFT

	cufftHandle planTens3D_x, planTens3D_y, planTens3D_z;

	int embed[1] = { 0 };
	int ndims_x[1] = { (int)N.x };
	int ndims_y[1] = { (int)N.y };
	int ndims_z[1] = { (int)N.z };

	//Forward fft along x direction (out-of-place):
	cufftResult cuffterr = cufftPlanMany(&planTens3D_x, 1, ndims_x,
		embed, 1, N.x,
		embed, 1, (N.x / 2 + 1),
		CUFFT_D2Z, (int)N.y*N.z);

	//Forward fft along y direction (out-of-place):
	cuffterr = cufftPlanMany(&planTens3D_y, 1, ndims_y,
		embed, (N.x / 2 + 1), 1,
		embed, (N.x / 2 + 1), 1,
		CUFFT_Z2Z, (int)(N.x / 2 + 1));

	//Forward fft along z direction (out-of-place):
	cuffterr = cufftPlanMany(&planTens3D_z, 1, ndims_z,
		embed, (N.x / 2 + 1) * N.y, 1,
		embed, (N.x / 2 + 1) * N.y, 1,
		CUFFT_Z2Z, (int)(N.x / 2 + 1)*N.y);

	//-------------- FFT REAL TENSOR INTO REAL KERNELS

	auto fftcomponent = [&](cu_arr<cufftDoubleReal>& cuIn, int component, bool odiag) -> void
	{
		//Now FFT along x from input to output
		cufftExecD2Z(planTens3D_x, cuIn, cuOut);

		//Now FFT along y from input to output
		for (int k = 0; k < N.z; k++) {

			cufftExecZ2Z(planTens3D_y, cuOut + k * (N.x / 2 + 1) * N.y, cuOut + k * (N.x / 2 + 1) * N.y, CUFFT_FORWARD);
		}

		//Now FFT along z from input to output
		cufftExecZ2Z(planTens3D_z, cuOut, cuOut, CUFFT_FORWARD);

		//extract and place in kernel
		if (!odiag) (*kernels[index])()->Set_Kdiag_cmpl((N.x / 2 + 1) * N.y * N.z, cuOut, component, transpose_xy);
		else (*kernels[index])()->Set_Kodiag_cmpl((N.x / 2 + 1) * N.y * N.z, cuOut, component, transpose_xy);
	};

	//-------------- CALCULATE DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	//diagonal components

	fftcomponent(cuInx, 1, false);
	fftcomponent(cuIny, 2, false);
	fftcomponent(cuInz, 3, false);

	//-------------- CALCULATE OFF-DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	if (pbc_images.IsNull()) {

		//no pbcs in any dimension -> use these
		if (!dtf_gpu.CalcOffDiagTens3D_Shifted(cuInx, cuIny, cuInz, n, N, h / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}
	else {

		//pbcs used in at least one dimension
		if (!dtf_gpu.CalcOffDiagTens3D_Shifted_PBC(
			cuInx, cuIny, cuInz, N, h / h_max, (*kernels[index])()->Get_shift() / h_max,
			true, ASYMPTOTIC_DISTANCE, pbc_images.x, pbc_images.y, pbc_images.z)) return error(BERROR_OUTOFGPUMEMORY_NCRIT);
	}

	//off-diagonal components
	fftcomponent(cuInx, 1, true);
	fftcomponent(cuIny, 2, true);
	fftcomponent(cuInz, 3, true);

	//-------------- CLEANUP

	cufftDestroy(planTens3D_x);
	cufftDestroy(planTens3D_y);
	cufftDestroy(planTens3D_z);

	return error;
}

#endif

#endif
