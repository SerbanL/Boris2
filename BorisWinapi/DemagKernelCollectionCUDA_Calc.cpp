#include "stdafx.h"
#include "DemagKernelCollectionCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_SDEMAG

#include "DemagTFunc.h"

//-------------------------- KERNEL CALCULATION

//this initializes all the convolution kernels for the given mesh dimensions. 2D is for n.z == 1.
BError DemagKernelCollectionCUDA::Calculate_Demag_Kernels(vector<DemagKernelCollectionCUDA*>& kernelCollection)
{
	BError error(__FUNCTION__);

	for (int index = 0; index < Rect_collection.size(); index++) {

		if (!error) {

			//before allocating and computing kernel, check to see if any other DemagKernelCollection module has not already calculated one like it

			//Rect_collection is in one-to-one correspondence with kernelCollection
			//For a demag kernel to be identical to the one we need, it must have the same shift, same h source and same h destination

			//shift for source as Rect_collection[index] and destination as this_rect
			cuReal3 shift = (this_rect.s - Rect_collection[index].s);
			if (use_multiple_outputs) shift *= -1.0;

			cuReal3 h_src = kernelCollection[index]->h;
			cuReal3 h_dst = h;
			
			if (use_multiple_outputs) {

				h_src = h;
				h_dst = kernelCollection[index]->h;
			}

			for (int idx = 0; idx < kernelCollection.size(); idx++) {

				shared_ptr<cu_obj<cuKerType>> existing_kernel = kernelCollection[idx]->KernelAlreadyComputed(shift, h_src, h_dst);

				if (existing_kernel != nullptr) {

					//found one : just increase ref count
					kernels[index] = existing_kernel;

					//is it z shifted? keep copy of flag in cpu memory
					zshifted[index] = (*kernels[index])()->GetFlag_zShifted();

					//is it inverse z-shifted?
					if (cuIsZ(shift.x) && cuIsZ(shift.y)) {

						if (cuIsNZ(shift.z) && cuIsZ(shift.z + (*kernels[index])()->Get_shift().z)) {

							//yes it is. mark it here so we can adjust the kernel multiplications
							inverse_shifted[index] = true;
						}
					}

					break;
				}
			}
			
			if (kernels[index] == nullptr) {

				//no -> allocate then compute it
				kernels[index] = shared_ptr<cu_obj<cuKerType>>(new cu_obj<cuKerType>());
				if (!(*kernels[index])()->AllocateKernels(Rect_collection[index], this_rect, N)) return error(BERROR_OUTOFGPUMEMORY_CRIT);

				//is it z shifted? keep copy of flag in cpu memory
				zshifted[index] = (*kernels[index])()->GetFlag_zShifted();

				//now compute it
				if ((*kernels[index])()->GetFlag_InternalDemag()) {

					(*kernels[index])()->Set_Shift_and_Cellsizes(cuReal3(), h, h);

					//use self versions
					if (n.z == 1) error = Calculate_Demag_Kernels_2D_Self(index);
					else error = Calculate_Demag_Kernels_3D_Self(index);
				}
				else {

					//shift for source as Rect_collection[index] and destination as this_rect
					if (use_multiple_outputs) {

						(*kernels[index])()->Set_Shift_and_Cellsizes(Rect_collection[index].s - this_rect.s, h, kernelCollection[index]->h);
					}
					else {

						(*kernels[index])()->Set_Shift_and_Cellsizes(this_rect.s - Rect_collection[index].s, kernelCollection[index]->h, h);
					}

					if (n.z == 1) {
						
						if (cuIsZ((*kernels[index])()->Get_shift().x) && cuIsZ((*kernels[index])()->Get_shift().y)) {

							//z-shifted kernels for 2D
							error = Calculate_Demag_Kernels_2D_zShifted(index);
						}
						
						else {

							//general 2D kernels (not z-shifted)
							error = Calculate_Demag_Kernels_2D_Complex_Full(index);
						}
					}
					else {
						
						if (cuIsZ((*kernels[index])()->Get_shift().x) && cuIsZ((*kernels[index])()->Get_shift().y)) {

							//z-shifted kernels for 3D
							error = Calculate_Demag_Kernels_3D_zShifted(index);
						}
						
						else {

							//general 3D kernels (not z-shifted)
							error = Calculate_Demag_Kernels_3D_Complex_Full(index);
						}
					}
				}

				//set flag to say it's been computed so it could be reused if needed
				if (!error) (*kernels[index])()->SetFlag_Calculated(true);
			}
		}
	}

	//now build kernels_gpu from kernels
	//kernels_gpu is an array in gpu memory
	//each element in kernels_gpu is a pointer to a cuKerType kernel - 1-2-1 relation with kernels
	//also build inverse_shifted_gpu
	for (int idx = 0; idx < kernels.size(); idx++) {

		kernels_gpu.push_back((cuKerType*&)(*kernels[idx]).get_managed_object());

		//transfer bool value from cpu to gpu memory in a cu_obj, so we can then place it in the cu_arr
		//the push_back method in cu_arr is meant to take a pointer in gpu memory
		//TO DO : make a push_back_from_cpu method that places in cu_arr a value directly from cpu memory
		cu_obj<bool> gpu_value;
		gpu_value.from_cpu(inverse_shifted[idx]);
		inverse_shifted_gpu.push_back((bool*&)gpu_value.get_managed_object());
	}

	num_kernels.from_cpu((int)kernels.size());

	return error;
}

//search to find a matching kernel that has already been computed and return pointer to it -> kernel can be identified from shift, source and destination discretisation
shared_ptr<cu_obj<cuKerType>> DemagKernelCollectionCUDA::KernelAlreadyComputed(cuReal3 shift, cuReal3 h_src, cuReal3 h_dst)
{
	//kernels[index] must not be nullptr, must have kernel_calculated = true and shift, h_src, h_dst must match the corresponding values in kernels[index] 
	
	for (int idx = 0; idx < kernels.size(); idx++) {

		if (kernels[idx] && (*kernels[idx])()->GetFlag_Calculated()) {

			//match in source and destination cellsizes?
			if ((*kernels[idx])()->Get_h_src() == h_src && (*kernels[idx])()->Get_h_dst() == h_dst) {

				//do the shifts match?
				if ((*kernels[idx])()->Get_shift() == shift) {

					return kernels[idx];
				}

				//are the shifts z shifts that differ only in sign? Only applies in 2D mode.
				if (N.z == 1 && cuIsZ(shift.x) && cuIsZ(shift.y) && cuIsZ((*kernels[idx])()->Get_shift().x) && cuIsZ((*kernels[idx])()->Get_shift().y) && cuIsZ(shift.z + (*kernels[idx])()->Get_shift().z)) {

					return kernels[idx];
				}
			}
		}
	}
	
	return nullptr;
}

//2D kernels (Kdiag_real, and K2D_odiag, with full use of kernel symmetries)
BError DemagKernelCollectionCUDA::Calculate_Demag_Kernels_2D_Self(int index)
{
	BError error(__FUNCTION__);

	//-------------- CALCULATE DEMAG TENSOR

	//Demag tensor components
	//
	// D11 D12 0
	// D12 D22 0
	// 0   0   D33

	//Ddiag : D11, D22, D33 are the diagonal tensor elements
	VEC<DBL3> Ddiag;
	if (!Ddiag.resize(N)) return error(BERROR_OUTOFMEMORY_NCRIT);

	//off-diagonal tensor elements
	vector<double> Dodiag;
	if (!malloc_vector(Dodiag, N.x*N.y)) return error(BERROR_OUTOFMEMORY_NCRIT);

	//use ratios instead of cellsizes directly - same result but better in terms of floating point errors
	DemagTFunc dtf;

	if (!dtf.CalcDiagTens2D(Ddiag, n, N, h / h_max)) return error(BERROR_OUTOFMEMORY_NCRIT);
	if (!dtf.CalcOffDiagTens2D(Dodiag, n, N, h / h_max)) return error(BERROR_OUTOFMEMORY_NCRIT);

	//-------------- DEMAG KERNELS ON CPU-ADDRESSABLE MEMORY

	VEC<DBL3> Kdiag_cpu;
	vector<double> K2D_odiag_cpu;

	if (!transpose_xy) {

		if (!Kdiag_cpu.resize(SZ3(N.x / 2 + 1, N.y / 2 + 1, 1))) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		if (!Kdiag_cpu.resize(SZ3(N.y / 2 + 1, N.x / 2 + 1, 1))) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	if (!malloc_vector(K2D_odiag_cpu, (N.x / 2 + 1) * (N.y / 2 + 1))) return error(BERROR_OUTOFMEMORY_CRIT);

	//-------------- SETUP FFT

	//setup fft object with fft computation lines
	FFTMethods_Cpp<double> fft;

	size_t maxN = maximum(N.x / 2 + 1, N.y / 2 + 1);

	vector<ReIm3> fft_line(maxN);
	vector<ReIm3> fft_line2(maxN);
	vector<ReIm> fft_line2d(maxN);
	vector<ReIm> fft_line2d2(maxN);

	//-------------- FFT REAL TENSOR INTO REAL KERNELS

	//NOTE : don't use parallel for loops as it will mess up the packing in the D tensor
	//If you want parallel loops you'll need to allocate additional temporary spaces, so not worth it for initialization
	//rather have slightly slower initialization than fail due to running out of memory for large problem sizes

	//FFT into Kernel forms ready for convolution multiplication
	for (int j = 0; j < N.y; j++) {

		fft.CopyRealShuffle(Ddiag.data() + j * N.x, fft_line.data(), N.x / 2);
		fft.FFT_Radix4_DIT(fft_line.data(), log2(N.x) - 1, N.x / 2);
		fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.x / 2);

		fft.CopyRealShuffle(Dodiag.data() + j * N.x, fft_line2d.data(), N.x / 2);
		fft.FFT_Radix4_DIT(fft_line2d.data(), log2(N.x) - 1, N.x / 2);
		fft.RealfromComplexFFT(fft_line2d.data(), fft_line2d2.data(), N.x / 2);

		//pack into Ddiag and Dodiag for next step
		for (int i = 0; i < N.x / 2 + 1; i++) {

			//even w.r.t. to x so output is purely real
			Ddiag[i + j * (N.x / 2 + 1)] = DBL3(fft_line2[i].x.Re, fft_line2[i].y.Re, fft_line2[i].z.Re);

			//odd w.r.t. to x so output is purely imaginary
			Dodiag[i + j * (N.x / 2 + 1)] = fft_line2d2[i].Im;
		}
	}

	for (int i = 0; i < N.x / 2 + 1; i++) {

		fft.CopyRealShuffle(Ddiag.data() + i, fft_line.data(), N.x / 2 + 1, N.y / 2);
		fft.FFT_Radix4_DIT(fft_line.data(), log2(N.y) - 1, N.y / 2);
		fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.y / 2);

		//the input sequence should actually be purely imaginary not purely real (but see below)
		fft.CopyRealShuffle(Dodiag.data() + i, fft_line2d.data(), N.x / 2 + 1, N.y / 2);
		fft.FFT_Radix4_DIT(fft_line2d.data(), log2(N.y) - 1, N.y / 2);
		fft.RealfromComplexFFT(fft_line2d.data(), fft_line2d2.data(), N.y / 2);

		//pack into output real kernels with reduced strides
		for (int j = 0; j < N.y / 2 + 1; j++) {

			if (!transpose_xy) {

				//even w.r.t. y so output is purely real
				Kdiag_cpu[i + j * (N.x / 2 + 1)] = DBL3(fft_line2[j].x.Re, fft_line2[j].y.Re, fft_line2[j].z.Re);

				//odd w.r.t. y so the purely imaginary input becomes purely real
				//however since we used CopyRealShuffle and RealfromComplexFFT, i.e. treating the input as purely real rather than purely imaginary we need to account for the i * i = -1 term, hence the - sign below
				K2D_odiag_cpu[i + j * (N.x / 2 + 1)] = -fft_line2d2[j].Im;
			}
			else {

				//even w.r.t. y so output is purely real
				Kdiag_cpu[j + i * (N.y / 2 + 1)] = DBL3(fft_line2[j].x.Re, fft_line2[j].y.Re, fft_line2[j].z.Re);

				//odd w.r.t. y so the purely imaginary input becomes purely real
				//however since we used CopyRealShuffle and RealfromComplexFFT, i.e. treating the input as purely real rather than purely imaginary we need to account for the i * i = -1 term, hence the - sign below
				K2D_odiag_cpu[j + i * (N.y / 2 + 1)] = -fft_line2d2[j].Im;
			}
		}
	}

	//-------------- TRANSFER TO GPU KERNELS

	(*kernels[index])()->Set_Kdiag_real(Kdiag_cpu);
	(*kernels[index])()->Set_K2D_odiag(K2D_odiag_cpu);

	//Done
	return error;
}

//2D layers, z shift only : Kernels can be stored as real with use of kernel symmetries. Kxx, Kyy, Kzz, Kxy real, Kxz, Kyz imaginary
BError DemagKernelCollectionCUDA::Calculate_Demag_Kernels_2D_zShifted(int index)
{
	BError error(__FUNCTION__);

	//-------------- DEMAG TENSOR

	//Demag tensor components
	//
	// D11 D12 D13
	// D12 D22 D23
	// D13 D23 D33

	//D11, D22, D33 are the diagonal tensor elements
	//D12, D13, D23 are the off-diagonal tensor elements
	VEC<DBL3> D;

	if (!D.resize(N)) return error(BERROR_OUTOFMEMORY_NCRIT);

	//object used to compute tensor elements
	DemagTFunc dtf;

	//-------------- DEMAG KERNEL ON CPU-ADDRESSABLE MEMORY

	VEC<DBL3> K_cpu;

	if (!transpose_xy) {

		if (!K_cpu.resize(SZ3(N.x / 2 + 1, N.y / 2 + 1, N.z))) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		if (!K_cpu.resize(SZ3(N.y / 2 + 1, N.x / 2 + 1, N.z))) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	//-------------- FFT SETUP

	//setup fft object with fft computation lines
	FFTMethods_Cpp<double> fft;

	size_t maxN = maximum(N.x / 2 + 1, N.y / 2 + 1);

	vector<ReIm3> fft_line(maxN);
	vector<ReIm3> fft_line2(maxN);

	//lambda used to transform an input tensor into an output kernel
	auto tensor_to_kernel = [&](VEC<DBL3>& tensor, VEC<DBL3>& kernel, bool off_diagonal) -> void {

		//-------------- FFT REAL TENSOR

		//NOTE : don't use parallel for loops as it will mess up the packing in the D tensor
		//If you want parallel loops you'll need to allocate additional temporary spaces, so not worth it for initialization
		//rather have slightly slower initialization than fail due to running out of memory for large problem sizes

		//FFT into Kernel forms ready for convolution multiplication - diagonal components
		for (int j = 0; j < N.y; j++) {

			fft.CopyRealShuffle(tensor.data() + j * N.x, fft_line.data(), N.x / 2);
			fft.FFT_Radix4_DIT(fft_line.data(), log2(N.x) - 1, N.x / 2);
			fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.x / 2);

			if (!off_diagonal) {

				//diagonal elements even in x

				//pack into tensor
				for (int i = 0; i < N.x / 2 + 1; i++) {

					tensor[i + j * (N.x / 2 + 1)] = DBL3(fft_line2[i].x.Re, fft_line2[i].y.Re, fft_line2[i].z.Re);
				}
			}
			else {

				//Nxy, Nxz odd in x, Nyz even in x

				//pack into tensor
				for (int i = 0; i < N.x / 2 + 1; i++) {

					tensor[i + j * (N.x / 2 + 1)] = DBL3(fft_line2[i].x.Im, fft_line2[i].y.Im, fft_line2[i].z.Re);
				}
			}
		}

		for (int i = 0; i < N.x / 2 + 1; i++) {

			fft.CopyRealShuffle(tensor.data() + i, fft_line.data(), N.x / 2 + 1, N.y / 2);
			fft.FFT_Radix4_DIT(fft_line.data(), log2(N.y) - 1, N.y / 2);
			fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.y / 2);

			if (!off_diagonal) {

				if (!transpose_xy) {

					//pack into output real kernels
					for (int j = 0; j < N.y / 2 + 1; j++) {

						kernel[i + j * (N.x / 2 + 1)] = DBL3(fft_line2[j].x.Re, fft_line2[j].y.Re, fft_line2[j].z.Re);
					}
				}
				else {

					//pack into output real kernels
					for (int j = 0; j < N.y / 2 + 1; j++) {

						kernel[j + i * (N.y / 2 + 1)] = DBL3(fft_line2[j].x.Re, fft_line2[j].y.Re, fft_line2[j].z.Re);
					}
				}
			}
			else {

				//Nxy odd in y, Nxz even in y, Nyz odd in y

				if (!transpose_xy) {

					//pack into output real kernels
					for (int j = 0; j < N.y / 2 + 1; j++) {

						//adjust for i * i = -1 in Nxy element
						kernel[i + j * (N.x / 2 + 1)] = DBL3(-fft_line2[j].x.Im, fft_line2[j].y.Re, fft_line2[j].z.Im);
					}
				}
				else {

					//pack into output real kernels
					for (int j = 0; j < N.y / 2 + 1; j++) {

						//adjust for i * i = -1 in Nxy element
						kernel[j + i * (N.y / 2 + 1)] = DBL3(-fft_line2[j].x.Im, fft_line2[j].y.Re, fft_line2[j].z.Im);
					}
				}
			}
		}
	};

	//-------------- CALCULATE DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	//no need to pass the actual cellsize values, just normalized values will do

	if (!dtf.CalcDiagTens2D_Shifted_Irregular(D, n, N, (*kernels[index])()->Get_h_src() / h_max, (*kernels[index])()->Get_h_dst() / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFMEMORY_NCRIT);

	tensor_to_kernel(D, K_cpu, false);

	//transfer to GPU
	(*kernels[index])()->Set_Kdiag_real(K_cpu);

	//-------------- CALCULATE OFF-DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	if (!dtf.CalcOffDiagTens2D_Shifted_Irregular(D, n, N, (*kernels[index])()->Get_h_src() / h_max, (*kernels[index])()->Get_h_dst() / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFMEMORY_NCRIT);

	tensor_to_kernel(D, K_cpu, true);

	//transfer to GPU
	(*kernels[index])()->Set_Kodiag_real(K_cpu);

	//Done
	return error;
}

//2D layers, complex kernels most general case (Kdiag_cmpl, and Kodiag_cmpl, without any kernel symmetries)
BError DemagKernelCollectionCUDA::Calculate_Demag_Kernels_2D_Complex_Full(int index)
{
	BError error(__FUNCTION__);

	//-------------- DEMAG TENSOR

	//Demag tensor components
	//
	// D11 D12 D13
	// D12 D22 D23
	// D13 D23 D33

	//D11, D22, D33 are the diagonal tensor elements
	//D12, D13, D23 are the off-diagonal tensor elements
	VEC<DBL3> D;

	if (!D.resize(N)) return error(BERROR_OUTOFMEMORY_NCRIT);

	//object used to compute tensor elements
	DemagTFunc dtf;

	//-------------- DEMAG KERNEL ON CPU-ADDRESSABLE MEMORY

	VEC<ReIm3> K_cpu, Scratch;

	if (!Scratch.resize(SZ3(N.x / 2 + 1, N.y, N.z))) return error(BERROR_OUTOFMEMORY_CRIT);

	if (!transpose_xy) {

		if (!K_cpu.resize(SZ3(N.x / 2 + 1, N.y, N.z))) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		if (!K_cpu.resize(SZ3(N.y, N.x / 2 + 1, N.z))) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	//-------------- FFT SETUP

	//setup fft object with fft computation lines
	FFTMethods_Cpp<double> fft;

	size_t maxN = maximum(N.x / 2 + 1, N.y);

	vector<ReIm3> fft_line(maxN);
	vector<ReIm3> fft_line2(maxN);

	//lambda used to transform an input tensor into an output kernel
	auto tensor_to_kernel = [&](VEC<DBL3>& tensor, VEC<ReIm3>& kernel, VEC<ReIm3>& Scratch) -> void {

		//-------------- FFT REAL TENSOR

		//NOTE : don't use parallel for loops as it will mess up the packing in the D tensor
		//If you want parallel loops you'll need to allocate additional temporary spaces, so not worth it for initialization
		//rather have slightly slower initialization than fail due to running out of memory for large problem sizes

		//FFT into Kernel forms ready for convolution multiplication - diagonal components
		for (int j = 0; j < N.y; j++) {

			fft.CopyRealShuffle(tensor.data() + j * N.x, fft_line.data(), N.x / 2);
			fft.FFT_Radix4_DIT(fft_line.data(), log2(N.x) - 1, N.x / 2);
			fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.x / 2);

			//pack into scratch space
			for (int i = 0; i < N.x / 2 + 1; i++) {

				Scratch[i + j * (N.x / 2 + 1)] = fft_line2[i];
			}
		}

		for (int i = 0; i < N.x / 2 + 1; i++) {

			fft.CopyShuffle(Scratch.data() + i, fft_line.data(), N.x / 2 + 1, N.y);
			fft.FFT_Radix4_DIT(fft_line.data(), log2(N.y), N.y);

			if (!transpose_xy) {

				//pack into output kernel
				for (int j = 0; j < N.y; j++) {

					kernel[i + j * (N.x / 2 + 1)] = fft_line[j];
				}
			}
			else {

				//pack into output kernel
				for (int j = 0; j < N.y; j++) {

					kernel[j + i * N.y] = fft_line[j];
				}
			}
		}
	};

	//-------------- CALCULATE DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	//no need to pass the actual cellsize values, just normalized values will do

	if (!dtf.CalcDiagTens2D_Shifted_Irregular(D, n, N, (*kernels[index])()->Get_h_src() / h_max, (*kernels[index])()->Get_h_dst() / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFMEMORY_NCRIT);

	tensor_to_kernel(D, K_cpu, Scratch);

	//transfer to GPU
	(*kernels[index])()->Set_Kdiag_cmpl(K_cpu);

	//-------------- CALCULATE OFF-DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	if (!dtf.CalcOffDiagTens2D_Shifted_Irregular(D, n, N, (*kernels[index])()->Get_h_src() / h_max, (*kernels[index])()->Get_h_dst() / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFMEMORY_NCRIT);

	tensor_to_kernel(D, K_cpu, Scratch);

	//transfer to GPU
	(*kernels[index])()->Set_Kodiag_cmpl(K_cpu);

	//Done
	return error;
}

//3D real kernels (Kdiag_real, and Kodiag_real, with full use of kernel symmetries)
BError DemagKernelCollectionCUDA::Calculate_Demag_Kernels_3D_Self(int index)
{
	BError error(__FUNCTION__);

	//-------------- DEMAG TENSOR

	//Demag tensor components
	//
	// D11 D12 D13
	// D12 D22 D23
	// D13 D23 D33

	//D11, D22, D33 are the diagonal tensor elements
	//D12, D13, D23 are the off-diagonal tensor elements
	VEC<DBL3> D;

	if (!D.resize(N)) return error(BERROR_OUTOFMEMORY_NCRIT);

	//object used to compute tensor elements
	DemagTFunc dtf;

	//-------------- DEMAG KERNEL ON CPU-ADDRESSABLE MEMORY

	VEC<DBL3> K_cpu;

	if (!transpose_xy) {

		if (!K_cpu.resize(SZ3(N.x / 2 + 1, N.y / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		if (!K_cpu.resize(SZ3(N.y / 2 + 1, N.x / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	//-------------- FFT SETUP

	//setup fft object with fft computation lines
	FFTMethods_Cpp<double> fft;

	size_t maxN = maximum(N.x / 2 + 1, N.y / 2 + 1, N.z / 2 + 1);

	vector<ReIm3> fft_line(maxN);
	vector<ReIm3> fft_line2(maxN);

	//lambda used to transform an input tensor into an output kernel
	auto tensor_to_kernel = [&](VEC<DBL3>& tensor, VEC<DBL3>& kernel, bool off_diagonal) -> void {

		//-------------- FFT REAL TENSOR

		//NOTE : don't use parallel for loops as it will mess up the packing in the D tensor
		//If you want parallel loops you'll need to allocate additional temporary spaces, so not worth it for initialization
		//rather have slightly slower initialization than fail due to running out of memory for large problem sizes

		//FFT into Kernel forms ready for convolution multiplication - diagonal components
		for (int k = 0; k < N.z; k++) {
			for (int j = 0; j < N.y; j++) {

				fft.CopyRealShuffle(tensor.data() + j * N.x + k * N.x*N.y, fft_line.data(), N.x / 2);
				fft.FFT_Radix4_DIT(fft_line.data(), log2(N.x) - 1, N.x / 2);
				fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.x / 2);

				//pack into lower half of tensor row for next step (keep same row and plane strides)
				for (int i = 0; i < N.x / 2 + 1; i++) {

					if (!off_diagonal) {

						//even w.r.t. to x so output is purely real
						tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = DBL3(fft_line2[i].x.Re, fft_line2[i].y.Re, fft_line2[i].z.Re);
					}
					else {

						//Dxy : odd x, Dxz : odd x, Dyz : even x
						tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = DBL3(fft_line2[i].x.Im, fft_line2[i].y.Im, fft_line2[i].z.Re);
					}
				}
			}
		}

		for (int k = 0; k < N.z; k++) {
			for (int i = 0; i < N.x / 2 + 1; i++) {

				fft.CopyRealShuffle(tensor.data() + i + k * (N.x / 2 + 1)*N.y, fft_line.data(), N.x / 2 + 1, N.y / 2);
				fft.FFT_Radix4_DIT(fft_line.data(), log2(N.y) - 1, N.y / 2);
				fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.y / 2);

				//pack into lower half of tensor column for next step (keep same row and plane strides)
				for (int j = 0; j < N.y / 2 + 1; j++) {

					if (!off_diagonal) {

						//even w.r.t. to y so output is purely real
						tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = DBL3(fft_line2[j].x.Re, fft_line2[j].y.Re, fft_line2[j].z.Re);
					}
					else {

						//Dxy : odd y, Dxz : even y, Dyz : odd y
						tensor[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = DBL3(fft_line2[j].x.Im, fft_line2[j].y.Re, fft_line2[j].z.Im);
					}
				}
			}
		}

		for (int j = 0; j < N.y / 2 + 1; j++) {
			for (int i = 0; i < N.x / 2 + 1; i++) {

				fft.CopyRealShuffle(tensor.data() + i + j * (N.x / 2 + 1), fft_line.data(), (N.x / 2 + 1)*N.y, N.z / 2);
				fft.FFT_Radix4_DIT(fft_line.data(), log2(N.z) - 1, N.z / 2);
				fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.z / 2);

				//pack into output kernels with reduced strides
				for (int k = 0; k < N.z / 2 + 1; k++) {

					if (!off_diagonal) {

						//even w.r.t. to z so output is purely real

						if (!transpose_xy) {

							kernel[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(fft_line2[k].x.Re, fft_line2[k].y.Re, fft_line2[k].z.Re);
						}
						else {

							kernel[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(fft_line2[k].x.Re, fft_line2[k].y.Re, fft_line2[k].z.Re);
						}
					}
					else {

						//Dxy : even z, Dxz : odd z, Dyz : odd z
						//Also multiply by -1 since all off-diagonal tensor elements have been odd twice
						//The final output is thus purely real but we always treated the input as purely real even when it should have been purely imaginary
						//This means we need to account for i * i = -1 at the end

						if (!transpose_xy) {

							kernel[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(-fft_line2[k].x.Re, -fft_line2[k].y.Im, -fft_line2[k].z.Im);
						}
						else {

							kernel[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = DBL3(-fft_line2[k].x.Re, -fft_line2[k].y.Im, -fft_line2[k].z.Im);
						}
					}
				}
			}
		}
	};

	//-------------- CALCULATE DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	//no need to pass the actual cellsize values, just normalized values will do
	if (!dtf.CalcDiagTens3D(D, n, N, h / h_max)) return error(BERROR_OUTOFMEMORY_NCRIT);

	tensor_to_kernel(D, K_cpu, false);

	//transfer to GPU
	(*kernels[index])()->Set_Kdiag_real(K_cpu);

	//-------------- CALCULATE OFF-DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	if (!dtf.CalcOffDiagTens3D(D, n, N, h / h_max)) return error(BERROR_OUTOFMEMORY_NCRIT);

	tensor_to_kernel(D, K_cpu, true);

	//transfer to GPU
	(*kernels[index])()->Set_Kodiag_real(K_cpu);

	//Done
	return error;
}

//3D layers, z shift only : Kernels can be stored with use of kernel symmetries (but still complex).
BError DemagKernelCollectionCUDA::Calculate_Demag_Kernels_3D_zShifted(int index)
{
	BError error(__FUNCTION__);

	//-------------- DEMAG TENSOR

	//Demag tensor components
	//
	// D11 D12 D13
	// D12 D22 D23
	// D13 D23 D33

	//D11, D22, D33 are the diagonal tensor elements
	//D12, D13, D23 are the off-diagonal tensor elements
	VEC<DBL3> D;

	if (!D.resize(N)) return error(BERROR_OUTOFMEMORY_NCRIT);

	//object used to compute tensor elements
	DemagTFunc dtf;

	//-------------- DEMAG KERNEL ON CPU-ADDRESSABLE MEMORY

	VEC<ReIm3> K_cpu, Scratch;

	if (!Scratch.resize(SZ3(N.x / 2 + 1, N.y, N.z))) return error(BERROR_OUTOFMEMORY_CRIT);

	if (!transpose_xy) {

		if (!K_cpu.resize(SZ3(N.x / 2 + 1, N.y / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		if (!K_cpu.resize(SZ3(N.y / 2 + 1, N.x / 2 + 1, N.z / 2 + 1))) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	//-------------- FFT SETUP

	//setup fft object with fft computation lines
	FFTMethods_Cpp<double> fft;

	size_t maxN = maximum(N.x / 2 + 1, N.y, N.z);

	vector<ReIm3> fft_line(maxN);
	vector<ReIm3> fft_line2(maxN);

	//lambda used to transform an input tensor into an output kernel
	auto tensor_to_kernel = [&](VEC<DBL3>& tensor, VEC<ReIm3>& kernel, VEC<ReIm3>& Scratch) -> void {

		//-------------- FFT REAL TENSOR

		//NOTE : don't use parallel for loops as it will mess up the packing in the D tensor
		//If you want parallel loops you'll need to allocate additional temporary spaces, so not worth it for initialization
		//rather have slightly slower initialization than fail due to running out of memory for large problem sizes

		//FFT into Kernel forms ready for convolution multiplication - diagonal components
		for (int k = 0; k < N.z; k++) {
			for (int j = 0; j < N.y; j++) {

				fft.CopyRealShuffle(tensor.data() + j * N.x + k * N.x*N.y, fft_line.data(), N.x / 2);
				fft.FFT_Radix4_DIT(fft_line.data(), log2(N.x) - 1, N.x / 2);
				fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.x / 2);

				//pack into scratch space
				for (int i = 0; i < N.x / 2 + 1; i++) {

					Scratch[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = fft_line2[i];
				}
			}
		}

		for (int k = 0; k < N.z; k++) {
			for (int i = 0; i < N.x / 2 + 1; i++) {

				fft.CopyShuffle(Scratch.data() + i + k * (N.x / 2 + 1)*N.y, fft_line.data(), N.x / 2 + 1, N.y);
				fft.FFT_Radix4_DIT(fft_line.data(), log2(N.y), N.y);

				//pack into scratch space
				for (int j = 0; j < N.y; j++) {

					Scratch[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = fft_line[j];
				}
			}
		}

		for (int j = 0; j < N.y / 2 + 1; j++) {
			for (int i = 0; i < N.x / 2 + 1; i++) {

				fft.CopyShuffle(Scratch.data() + i + j * (N.x / 2 + 1), fft_line.data(), (N.x / 2 + 1)*N.y, N.z);
				fft.FFT_Radix4_DIT(fft_line.data(), log2(N.z), N.z);

				if (!transpose_xy) {

					//pack into output kernel with reduced strides
					for (int k = 0; k < N.z / 2 + 1; k++) {

						kernel[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = fft_line[k];
					}
				}
				else {

					//pack into output kernel with reduced strides
					for (int k = 0; k < N.z / 2 + 1; k++) {

						kernel[j + i * (N.y / 2 + 1) + k * (N.x / 2 + 1) * (N.y / 2 + 1)] = fft_line[k];
					}
				}
			}
		}
	};

	//-------------- CALCULATE DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	//no need to pass the actual cellsize values, just normalized values will do
	if (!dtf.CalcDiagTens3D_Shifted(D, n, N, h / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFMEMORY_NCRIT);

	tensor_to_kernel(D, K_cpu, Scratch);

	//transfer to GPU
	(*kernels[index])()->Set_Kdiag_cmpl(K_cpu);

	//-------------- CALCULATE OFF-DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	if (!dtf.CalcOffDiagTens3D_Shifted(D, n, N, h / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFMEMORY_NCRIT);

	tensor_to_kernel(D, K_cpu, Scratch);

	//transfer to GPU
	(*kernels[index])()->Set_Kodiag_cmpl(K_cpu);

	//Done
	return error;
}

//3D complex kernels (Kdiag_cmpl, and Kodiag_cmpl, without any kernel symmetries)
BError DemagKernelCollectionCUDA::Calculate_Demag_Kernels_3D_Complex_Full(int index)
{
	BError error(__FUNCTION__);

	//-------------- DEMAG TENSOR

	//Demag tensor components
	//
	// D11 D12 D13
	// D12 D22 D23
	// D13 D23 D33

	//D11, D22, D33 are the diagonal tensor elements
	//D12, D13, D23 are the off-diagonal tensor elements
	VEC<DBL3> D;

	if (!D.resize(N)) return error(BERROR_OUTOFMEMORY_NCRIT);

	//object used to compute tensor elements
	DemagTFunc dtf;

	//-------------- DEMAG KERNEL ON CPU-ADDRESSABLE MEMORY

	VEC<ReIm3> K_cpu, Scratch;

	if (!Scratch.resize(SZ3(N.x / 2 + 1, N.y, N.z))) return error(BERROR_OUTOFMEMORY_CRIT);

	if (!transpose_xy) {

		if (!K_cpu.resize(SZ3(N.x / 2 + 1, N.y, N.z))) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		if (!K_cpu.resize(SZ3(N.y, N.x / 2 + 1, N.z))) return error(BERROR_OUTOFMEMORY_CRIT);
	}

	//-------------- FFT SETUP

	//setup fft object with fft computation lines
	FFTMethods_Cpp<double> fft;

	size_t maxN = maximum(N.x / 2 + 1, N.y, N.z);

	vector<ReIm3> fft_line(maxN);
	vector<ReIm3> fft_line2(maxN);

	//lambda used to transform an input tensor into an output kernel
	auto tensor_to_kernel = [&](VEC<DBL3>& tensor, VEC<ReIm3>& kernel, VEC<ReIm3>& Scratch) -> void {

		//-------------- FFT REAL TENSOR

		//NOTE : don't use parallel for loops as it will mess up the packing in the D tensor
		//If you want parallel loops you'll need to allocate additional temporary spaces, so not worth it for initialization
		//rather have slightly slower initialization than fail due to running out of memory for large problem sizes

		//FFT into Kernel forms ready for convolution multiplication - diagonal components
		for (int k = 0; k < N.z; k++) {
			for (int j = 0; j < N.y; j++) {

				fft.CopyRealShuffle(tensor.data() + j * N.x + k * N.x*N.y, fft_line.data(), N.x / 2);
				fft.FFT_Radix4_DIT(fft_line.data(), log2(N.x) - 1, N.x / 2);
				fft.RealfromComplexFFT(fft_line.data(), fft_line2.data(), N.x / 2);

				//pack into Scratch space
				for (int i = 0; i < N.x / 2 + 1; i++) {

					Scratch[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = fft_line2[i];
				}
			}
		}

		for (int k = 0; k < N.z; k++) {
			for (int i = 0; i < N.x / 2 + 1; i++) {

				fft.CopyShuffle(Scratch.data() + i + k * (N.x / 2 + 1)*N.y, fft_line.data(), N.x / 2 + 1, N.y);
				fft.FFT_Radix4_DIT(fft_line.data(), log2(N.y), N.y);

				//pack into Scratch space
				for (int j = 0; j < N.y; j++) {

					Scratch[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = fft_line[j];
				}
			}
		}

		for (int j = 0; j < N.y; j++) {
			for (int i = 0; i < N.x / 2 + 1; i++) {

				fft.CopyShuffle(Scratch.data() + i + j * (N.x / 2 + 1), fft_line.data(), (N.x / 2 + 1)*N.y, N.z);
				fft.FFT_Radix4_DIT(fft_line.data(), log2(N.z), N.z);

				if (!transpose_xy) {

					//pack into output kernel
					for (int k = 0; k < N.z; k++) {

						kernel[i + j * (N.x / 2 + 1) + k * (N.x / 2 + 1) * N.y] = fft_line[k];
					}
				}
				else {

					//pack into output kernel
					for (int k = 0; k < N.z; k++) {

						kernel[j + i * N.y + k * (N.x / 2 + 1) * N.y] = fft_line[k];
					}
				}
			}
		}
	};

	//-------------- CALCULATE DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	

	//no need to pass the actual cellsize values, just normalized values will do
	if (!dtf.CalcDiagTens3D_Shifted(D, n, N, h / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFMEMORY_NCRIT);

	tensor_to_kernel(D, K_cpu, Scratch);

	//transfer to GPU
	(*kernels[index])()->Set_Kdiag_cmpl(K_cpu);

	//-------------- CALCULATE OFF-DIAGONAL TENSOR ELEMENTS THEN TRANSFORM INTO KERNEL

	if (!dtf.CalcOffDiagTens3D_Shifted(D, n, N, h / h_max, (*kernels[index])()->Get_shift() / h_max)) return error(BERROR_OUTOFMEMORY_NCRIT);

	tensor_to_kernel(D, K_cpu, Scratch);

	//transfer to GPU
	(*kernels[index])()->Set_Kodiag_cmpl(K_cpu);

	//Done
	return error;
}

#endif

#endif
