#include "ConvolutionDataCUDA.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"

//---------------------------------------------------Copy input arrays from cuVEC or cuVEC_VC ( all <cuReal3> )

template <typename cuVECIn>
__global__ void In_to_cuFFTArrays_forInPlace(cuVECIn& In, cuReal* cuSx, cuReal* cuSy, cuReal* cuSz, cuSZ3& N)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = In.n;
	cuINT3 ijk = cuINT3(idx % (N.x + 2), (idx / (N.x + 2)) % N.y, idx / ((N.x + 2)*N.y));

	if (idx < (N.x + 2)*N.y*N.z) {

		if (ijk.i < n.x && ijk.j < n.y && ijk.k < n.z) {

			int idx_cell = ijk.i + ijk.j * n.x + ijk.k * n.x*n.y;

			cuSx[idx] = In[idx_cell].x;
			cuSy[idx] = In[idx_cell].y;
			cuSz[idx] = In[idx_cell].z;
		}
		else {

			cuSx[idx] = 0.0;
			cuSy[idx] = 0.0;
			cuSz[idx] = 0.0;
		}
	}
}

template <typename cuVECIn>
__global__ void In_to_cuFFTArrays_forOutOfPlace(cuVECIn& In, cuReal* cuSx, cuReal* cuSy, cuReal* cuSz, cuSZ3& N)
{	
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = In.n;

	if (idx < n.dim()) {

		cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));
		int idx_out = ijk.i + ijk.j * N.x + ijk.k * N.x*n.y;

		cuSx[idx_out] = In[idx].x;
		cuSy[idx_out] = In[idx].y;
		cuSz[idx_out] = In[idx].z;
	}
}

//---------------------------------------------------Set/Add to output cuVEC_VC to cuVEC or cuVEC_VC ( all <cuReal3> )

template <typename cuVECOut>
__global__ void cuFFTArrays_to_Out_Set_forInPlace(cuVEC_VC<cuReal3>& In, cuVECOut& Out, cuReal* cuSx, cuReal* cuSy, cuReal* cuSz, cuSZ3& N, cuReal& energy, bool do_reduction)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = Out.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuReal energy_ = 0.0;

	if (idx < Out.linear_size()) {

		cuReal3 Heff_value;

		Heff_value.x = cuSx[ijk.i + ijk.j * (N.x + 2) + ijk.k * (N.x + 2) * N.y] / N.dim();
		Heff_value.y = cuSy[ijk.i + ijk.j * (N.x + 2) + ijk.k * (N.x + 2) * N.y] / N.dim();
		Heff_value.z = cuSz[ijk.i + ijk.j * (N.x + 2) + ijk.k * (N.x + 2) * N.y] / N.dim();

		if (do_reduction) {

			int non_empty_cells = In.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuReal)MU0 * In[idx] * Heff_value / (2 * non_empty_cells);
		}

		Out[idx] = Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

template <typename cuVECOut>
__global__ void cuFFTArrays_to_Out_Add_forInPlace(cuVEC_VC<cuReal3>& In, cuVECOut& Out, cuReal* cuSx, cuReal* cuSy, cuReal* cuSz, cuSZ3& N, cuReal& energy, bool do_reduction)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = Out.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuReal energy_ = 0.0;

	if (idx < Out.linear_size()) {

		cuReal3 Heff_value;

		Heff_value.x = cuSx[ijk.i + ijk.j * (N.x + 2) + ijk.k * (N.x + 2) * N.y] / N.dim();
		Heff_value.y = cuSy[ijk.i + ijk.j * (N.x + 2) + ijk.k * (N.x + 2) * N.y] / N.dim();
		Heff_value.z = cuSz[ijk.i + ijk.j * (N.x + 2) + ijk.k * (N.x + 2) * N.y] / N.dim();

		if (do_reduction) {

			int non_empty_cells = In.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuReal)MU0 * In[idx] * Heff_value / (2 * non_empty_cells);
		}

		Out[idx] += Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

template <typename cuVECOut>
__global__ void cuFFTArrays_to_Out_Set_forOutOfPlace(cuVEC_VC<cuReal3>& In, cuVECOut& Out, cuReal* cuSx, cuReal* cuSy, cuReal* cuSz, cuSZ3& N, cuReal& energy, bool do_reduction)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = Out.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuReal energy_ = 0.0;

	if (idx < Out.linear_size()) {

		cuReal3 Heff_value;

		Heff_value.x = cuSx[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();
		Heff_value.y = cuSy[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();
		Heff_value.z = cuSz[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();

		if (do_reduction) {

			int non_empty_cells = In.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuReal)MU0 * In[idx] * Heff_value / (2 * non_empty_cells);
		}

		Out[idx] = Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

template <typename cuVECOut>
__global__ void cuFFTArrays_to_Out_Set_forOutOfPlace_batched_xy(cuVEC_VC<cuReal3>& In, cuVECOut& Out, cuReal* cuSx, cuReal* cuSy, cuReal* cuSz, cuSZ3& N, cuReal& energy, bool do_reduction)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = Out.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuReal energy_ = 0.0;

	if (idx < Out.linear_size()) {

		cuReal3 Heff_value;

		Heff_value.x = cuSx[ijk.i + ijk.j * N.x + ijk.k * N.x * N.y] / N.dim();
		Heff_value.y = cuSy[ijk.i + ijk.j * N.x + ijk.k * N.x * N.y] / N.dim();
		Heff_value.z = cuSz[ijk.i + ijk.j * N.x + ijk.k * N.x * N.y] / N.dim();

		if (do_reduction) {

			int non_empty_cells = In.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuReal)MU0 * In[idx] * Heff_value / (2 * non_empty_cells);
		}

		Out[idx] = Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

template <typename cuVECOut>
__global__ void cuFFTArrays_to_Out_Add_forOutOfPlace(cuVEC_VC<cuReal3>& In, cuVECOut& Out, cuReal* cuSx, cuReal* cuSy, cuReal* cuSz, cuSZ3& N, cuReal& energy, bool do_reduction)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = Out.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuReal energy_ = 0.0;

	if (idx < Out.linear_size()) {

		cuReal3 Heff_value;

		Heff_value.x = cuSx[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();
		Heff_value.y = cuSy[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();
		Heff_value.z = cuSz[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();

		if (do_reduction) {

			int non_empty_cells = In.get_nonempty_cells();
			if (non_empty_cells) energy_ = -(cuReal)MU0 * In[idx] * Heff_value / (2 * non_empty_cells);
		}

		Out[idx] += Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

//---------------------------------------------------Set/Add to output cuVEC to cuVEC or cuVEC_VC ( all <cuReal3> )

template <typename cuVECOut>
__global__ void cuFFTArrays_to_Out_Set_forInPlace(cuVEC<cuReal3>& In, cuVECOut& Out, cuReal* cuSx, cuReal* cuSy, cuReal* cuSz, cuSZ3& N, cuReal& energy, bool do_reduction)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = Out.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuReal energy_ = 0.0;

	if (idx < Out.linear_size()) {

		cuReal3 Heff_value;

		Heff_value.x = cuSx[ijk.i + ijk.j * (N.x + 2) + ijk.k * (N.x + 2) * N.y] / N.dim();
		Heff_value.y = cuSy[ijk.i + ijk.j * (N.x + 2) + ijk.k * (N.x + 2) * N.y] / N.dim();
		Heff_value.z = cuSz[ijk.i + ijk.j * (N.x + 2) + ijk.k * (N.x + 2) * N.y] / N.dim();

		if (do_reduction) {

			size_t non_empty_cells = In.get_aux_integer();
			if (non_empty_cells) energy_ = -(cuReal)MU0 * In[idx] * Heff_value / (2 * non_empty_cells);
		}

		Out[idx] = Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

template <typename cuVECOut>
__global__ void cuFFTArrays_to_Out_Add_forInPlace(cuVEC<cuReal3>& In, cuVECOut& Out, cuReal* cuSx, cuReal* cuSy, cuReal* cuSz, cuSZ3& N, cuReal& energy, bool do_reduction)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = Out.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuReal energy_ = 0.0;

	if (idx < Out.linear_size()) {

		cuReal3 Heff_value;

		Heff_value.x = cuSx[ijk.i + ijk.j * (N.x + 2) + ijk.k * (N.x + 2) * N.y] / N.dim();
		Heff_value.y = cuSy[ijk.i + ijk.j * (N.x + 2) + ijk.k * (N.x + 2) * N.y] / N.dim();
		Heff_value.z = cuSz[ijk.i + ijk.j * (N.x + 2) + ijk.k * (N.x + 2) * N.y] / N.dim();

		if (do_reduction) {

			size_t non_empty_cells = In.get_aux_integer();
			if (non_empty_cells) energy_ = -(cuReal)MU0 * In[idx] * Heff_value / (2 * non_empty_cells);
		}

		Out[idx] += Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

template <typename cuVECOut>
__global__ void cuFFTArrays_to_Out_Set_forOutOfPlace(cuVEC<cuReal3>& In, cuVECOut& Out, cuReal* cuSx, cuReal* cuSy, cuReal* cuSz, cuSZ3& N, cuReal& energy, bool do_reduction)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = Out.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuReal energy_ = 0.0;

	if (idx < Out.linear_size()) {

		cuReal3 Heff_value;

		Heff_value.x = cuSx[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();
		Heff_value.y = cuSy[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();
		Heff_value.z = cuSz[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();

		if (do_reduction) {

			size_t non_empty_cells = In.get_aux_integer();
			if (non_empty_cells) energy_ = -(cuReal)MU0 * In[idx] * Heff_value / (2 * non_empty_cells);
		}

		Out[idx] = Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

template <typename cuVECOut>
__global__ void cuFFTArrays_to_Out_Add_forOutOfPlace(cuVEC<cuReal3>& In, cuVECOut& Out, cuReal* cuSx, cuReal* cuSy, cuReal* cuSz, cuSZ3& N, cuReal& energy, bool do_reduction)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = Out.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuReal energy_ = 0.0;

	if (idx < Out.linear_size()) {

		cuReal3 Heff_value;

		Heff_value.x = cuSx[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();
		Heff_value.y = cuSy[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();
		Heff_value.z = cuSz[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();

		if (do_reduction) {

			size_t non_empty_cells = In.get_aux_integer();
			if (non_empty_cells) energy_ = -(cuReal)MU0 * In[idx] * Heff_value / (2 * non_empty_cells);
		}

		Out[idx] += Heff_value;
	}

	if (do_reduction) reduction_sum(0, 1, &energy_, energy);
}

template <typename cuVECOut>
__global__ void cuFFTArrays_to_Out_Set_weighted_forOutOfPlace(cuVEC<cuReal3>& In, cuVECOut& Out, cuReal* cuSx, cuReal* cuSy, cuReal* cuSz, cuSZ3& N, cuReal& energy, cuReal& energy_weight)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = Out.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuReal energy_ = 0.0;

	if (idx < Out.linear_size()) {

		cuReal3 Heff_value;

		Heff_value.x = cuSx[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();
		Heff_value.y = cuSy[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();
		Heff_value.z = cuSz[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();

		size_t non_empty_cells = In.get_aux_integer();
		if (non_empty_cells) energy_ = -(cuReal)MU0 * In[idx] * Heff_value * energy_weight / (2 * non_empty_cells);

		Out[idx] = Heff_value;
	}

	reduction_sum(0, 1, &energy_, energy);
}

template <typename cuVECOut>
__global__ void cuFFTArrays_to_Out_Add_weighted_forOutOfPlace(cuVEC<cuReal3>& In, cuVECOut& Out, cuReal* cuSx, cuReal* cuSy, cuReal* cuSz, cuSZ3& N, cuReal& energy, cuReal& energy_weight)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuSZ3 n = Out.n;
	cuINT3 ijk = cuINT3(idx % n.x, (idx / n.x) % n.y, idx / (n.x*n.y));

	cuReal energy_ = 0.0;

	if (idx < Out.linear_size()) {

		cuReal3 Heff_value;

		Heff_value.x = cuSx[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();
		Heff_value.y = cuSy[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();
		Heff_value.z = cuSz[ijk.i + ijk.j * N.x + ijk.k * N.x * n.y] / N.dim();

		size_t non_empty_cells = In.get_aux_integer();
		if (non_empty_cells) energy_ = -(cuReal)MU0 * In[idx] * Heff_value * energy_weight / (2 * non_empty_cells);

		Out[idx] += Heff_value;
	}

	reduction_sum(0, 1, &energy_, energy);
}

//-------------------------- RUN-TIME METHODS

template void ConvolutionDataCUDA::CopyInputData(cu_obj<cuVEC<cuReal3>>& In);
template void ConvolutionDataCUDA::CopyInputData(cu_obj<cuVEC_VC<cuReal3>>& In);

//Copy data to cuSx, cuSy, cuSz arrays at start of convolution iteration
template <typename cuVECIn>
void ConvolutionDataCUDA::CopyInputData(cu_obj<cuVECIn>& In)
{
	//fine to dereference In() here even though it points to gpu memory : it gets passed by reference to a cuda kernel so the dereferencing is actually done in gpu code
	In_to_cuFFTArrays_forOutOfPlace << < (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (*In(), cuIn_x, cuIn_y, cuIn_z, cuN);

	//we need to zero the parts above n.z and n.y before a new convolution
		
	if (!transpose_xy) {
		
		//zero pad upper y region (from n.y up to N.y) but only up to n.z
		//with pbc enabled n.y will be equal to N.y so no need to launch zero padding kernel
		if (N.y - n.y) {

			cu_zeropad((N.x / 2 + 1)*(N.y - n.y)*n.z, cuS_x, cuS_y, cuS_z, cuNc_xy, cuUpper_y_region);
		}
	}
	//else {

		//in the transpose_xy mode the upper y region zero padding is done after transposition in the forward_fft methods
	//}

	//for 3D problems we also need to zero pad the upper z region (from n.z up to N.z)
	if (n.z > 1 && !q2D_level) {

		//with pbc enabled n.z will be equal to N.z so no need to launch zero padding kernel
		if (N.z - n.z) {

			cu_zeropad((N.x / 2 + 1)*N.y*(N.z - n.z), cuS_x, cuS_y, cuS_z, cuNc_xy, cuUpper_z_region);
		}
	}
	else if (n.z > 1 && q2D_level && n.z != N.z / 2) {

		//if using q2D level, make sure to zero pad from n.z up to N.z / 2 as these values might not be the same
		cu_zeropad((N.x / 2 + 1)*N.y*(N.z / 2 - n.z), cuS_x, cuS_y, cuS_z, cuNc_xy_q2d, cuUpper_z_region_q2d);
	}
}

void ConvolutionDataCUDA::forward_fft_2D(void)
{
	//Forward 2D FFT
#if SINGLEPRECISION == 1

	if (!transpose_xy) {

		cufftExecR2C(plan2D_fwd_x, cuIn_x, cuS_x);
		cufftExecC2C(plan2D_y, cuS_x, cuS_x, CUFFT_FORWARD);

		cufftExecR2C(plan2D_fwd_x, cuIn_y, cuS_y);
		cufftExecC2C(plan2D_y, cuS_y, cuS_y, CUFFT_FORWARD);

		cufftExecR2C(plan2D_fwd_x, cuIn_z, cuS_z);
		cufftExecC2C(plan2D_y, cuS_z, cuS_z, CUFFT_FORWARD);
	}
	else {

		if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {

			//interleave transpose_xy operations

			cufftExecR2C(plan2D_fwd_x, cuIn_x, cuSquart_x);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuSquart_x, cuS_x, cuNcquart_xy, cuNc_xy, cuNc_xy);

			cufftExecR2C(plan2D_fwd_x, cuIn_y, cuSquart_y);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuSquart_y, cuS_y, cuNcquart_xy, cuNc_xy, cuNc_xy);

			cufftExecR2C(plan2D_fwd_x, cuIn_z, cuSquart_z);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuSquart_z, cuS_z, cuNcquart_xy, cuNc_xy, cuNc_xy);

			//with pbc enabled n.y will be equal to N.y so no need to launch zero padding kernel
			if (N.y - n.y) {

				cu_zeropad((N.x / 2 + 1)*(N.y - n.y), cuS_x, cuS_y, cuS_z, cuNc_yx, cuUpper_y_transposed_region);
			}

			cufftExecC2C(plan2D_y, cuS_x, cuS_x, CUFFT_FORWARD);
			cufftExecC2C(plan2D_y, cuS_y, cuS_y, CUFFT_FORWARD);
			cufftExecC2C(plan2D_y, cuS_z, cuS_z, CUFFT_FORWARD);
		}
		else {

			//transpose_xy operation in a single call

			cufftExecR2C(plan2D_fwd_x, cuIn_x, cuSquart_x);
			cufftExecR2C(plan2D_fwd_x, cuIn_y, cuSquart_y);
			cufftExecR2C(plan2D_fwd_x, cuIn_z, cuSquart_z);

			cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_x, cuSquart_y, cuSquart_z, cuS_x, cuS_y, cuS_z, cuNcquart_xy, cuNc_xy, cuNc_xy);
			
			//with pbc enabled n.y will be equal to N.y so no need to launch zero padding kernel
			if (N.y - n.y) {

				cu_zeropad((N.x / 2 + 1)*(N.y - n.y)*n.z, cuS_x, cuS_y, cuS_z, cuNc_yx, cuUpper_y_transposed_region);
			}

			cufftExecC2C(plan2D_y, cuS_x, cuS_x, CUFFT_FORWARD);
			cufftExecC2C(plan2D_y, cuS_y, cuS_y, CUFFT_FORWARD);
			cufftExecC2C(plan2D_y, cuS_z, cuS_z, CUFFT_FORWARD);
		}
	}

#else

	if (!transpose_xy) {

		cufftExecD2Z(plan2D_fwd_x, cuIn_x, cuS_x);
		cufftExecZ2Z(plan2D_y, cuS_x, cuS_x, CUFFT_FORWARD);

		cufftExecD2Z(plan2D_fwd_x, cuIn_y, cuS_y);
		cufftExecZ2Z(plan2D_y, cuS_y, cuS_y, CUFFT_FORWARD);

		cufftExecD2Z(plan2D_fwd_x, cuIn_z, cuS_z);
		cufftExecZ2Z(plan2D_y, cuS_z, cuS_z, CUFFT_FORWARD);
	}
	else {

		if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {

			//interleave transpose_xy operations

			cufftExecD2Z(plan2D_fwd_x, cuIn_x, cuSquart_x);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuSquart_x, cuS_x, cuNcquart_xy, cuNc_xy, cuNc_xy);

			cufftExecD2Z(plan2D_fwd_x, cuIn_y, cuSquart_y);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuSquart_y, cuS_y, cuNcquart_xy, cuNc_xy, cuNc_xy);

			cufftExecD2Z(plan2D_fwd_x, cuIn_z, cuSquart_z);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuSquart_z, cuS_z, cuNcquart_xy, cuNc_xy, cuNc_xy);

			//with pbc enabled n.y will be equal to N.y so no need to launch zero padding kernel
			if (N.y - n.y) {

				cu_zeropad((N.x / 2 + 1)*(N.y - n.y), cuS_x, cuS_y, cuS_z, cuNc_yx, cuUpper_y_transposed_region);
			}

			cufftExecZ2Z(plan2D_y, cuS_x, cuS_x, CUFFT_FORWARD);
			cufftExecZ2Z(plan2D_y, cuS_y, cuS_y, CUFFT_FORWARD);
			cufftExecZ2Z(plan2D_y, cuS_z, cuS_z, CUFFT_FORWARD);
}
		else {

			//transpose_xy operation in a single call

			cufftExecD2Z(plan2D_fwd_x, cuIn_x, cuSquart_x);
			cufftExecD2Z(plan2D_fwd_x, cuIn_y, cuSquart_y);
			cufftExecD2Z(plan2D_fwd_x, cuIn_z, cuSquart_z);

			cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_x, cuSquart_y, cuSquart_z, cuS_x, cuS_y, cuS_z, cuNcquart_xy, cuNc_xy, cuNc_xy);
			
			//with pbc enabled n.y will be equal to N.y so no need to launch zero padding kernel
			if (N.y - n.y) {

				cu_zeropad((N.x / 2 + 1)*(N.y - n.y)*n.z, cuS_x, cuS_y, cuS_z, cuNc_yx, cuUpper_y_transposed_region);
			}

			cufftExecZ2Z(plan2D_y, cuS_x, cuS_x, CUFFT_FORWARD);
			cufftExecZ2Z(plan2D_y, cuS_y, cuS_y, CUFFT_FORWARD);
			cufftExecZ2Z(plan2D_y, cuS_z, cuS_z, CUFFT_FORWARD);
		}
	}

#endif
}

void ConvolutionDataCUDA::inverse_fft_2D(void)
{
	//Inverse 2D FFT
#if SINGLEPRECISION == 1

	if (!transpose_xy) {

		cufftExecC2C(plan2D_y, cuS_x, cuS_x, CUFFT_INVERSE);
		cufftExecC2R(plan2D_inv_x, cuS_x, cuOut_x);

		cufftExecC2C(plan2D_y, cuS_y, cuS_y, CUFFT_INVERSE);
		cufftExecC2R(plan2D_inv_x, cuS_y, cuOut_y);

		cufftExecC2C(plan2D_y, cuS_z, cuS_z, CUFFT_INVERSE);
		cufftExecC2R(plan2D_inv_x, cuS_z, cuOut_z);
	}
	else {
		
		if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {

			//interleave transpose_xy operations

			cufftExecC2C(plan2D_y, cuS_x, cuS_x, CUFFT_INVERSE);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS_x, cuSquart_x, cuNcquart_yx, cuNc_yx, cuNc_yx);
			cufftExecC2R(plan2D_inv_x, cuSquart_x, cuOut_x);

			cufftExecC2C(plan2D_y, cuS_y, cuS_y, CUFFT_INVERSE);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS_y, cuSquart_y, cuNcquart_yx, cuNc_yx, cuNc_yx);
			cufftExecC2R(plan2D_inv_x, cuSquart_y, cuOut_y);

			cufftExecC2C(plan2D_y, cuS_z, cuS_z, CUFFT_INVERSE);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS_z, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNc_yx);
			cufftExecC2R(plan2D_inv_x, cuSquart_z, cuOut_z);
		}
		else {

			//transpose_xy operation in a single call

			cufftExecC2C(plan2D_y, cuS_x, cuS_x, CUFFT_INVERSE);
			cufftExecC2C(plan2D_y, cuS_y, cuS_y, CUFFT_INVERSE);
			cufftExecC2C(plan2D_y, cuS_z, cuS_z, CUFFT_INVERSE);

			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS_x, cuS_y, cuS_z, cuSquart_x, cuSquart_y, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNc_yx);

			cufftExecC2R(plan2D_inv_x, cuSquart_x, cuOut_x);
			cufftExecC2R(plan2D_inv_x, cuSquart_y, cuOut_y);
			cufftExecC2R(plan2D_inv_x, cuSquart_z, cuOut_z);
		}
	}

#else

	if (!transpose_xy) {

		cufftExecZ2Z(plan2D_y, cuS_x, cuS_x, CUFFT_INVERSE);
		cufftExecZ2D(plan2D_inv_x, cuS_x, cuOut_x);

		cufftExecZ2Z(plan2D_y, cuS_y, cuS_y, CUFFT_INVERSE);
		cufftExecZ2D(plan2D_inv_x, cuS_y, cuOut_y);

		cufftExecZ2Z(plan2D_y, cuS_z, cuS_z, CUFFT_INVERSE);
		cufftExecZ2D(plan2D_inv_x, cuS_z, cuOut_z);
		}
	else {

		if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {

			//interleave transpose_xy operations

			cufftExecZ2Z(plan2D_y, cuS_x, cuS_x, CUFFT_INVERSE);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS_x, cuSquart_x, cuNcquart_yx, cuNc_yx, cuNc_yx);
			cufftExecZ2D(plan2D_inv_x, cuSquart_x, cuOut_x);

			cufftExecZ2Z(plan2D_y, cuS_y, cuS_y, CUFFT_INVERSE);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS_y, cuSquart_y, cuNcquart_yx, cuNc_yx, cuNc_yx);
			cufftExecZ2D(plan2D_inv_x, cuSquart_y, cuOut_y);

			cufftExecZ2Z(plan2D_y, cuS_z, cuS_z, CUFFT_INVERSE);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS_z, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNc_yx);
			cufftExecZ2D(plan2D_inv_x, cuSquart_z, cuOut_z);
		}
		else {

			//transpose_xy operation in a single call

			cufftExecZ2Z(plan2D_y, cuS_x, cuS_x, CUFFT_INVERSE);
			cufftExecZ2Z(plan2D_y, cuS_y, cuS_y, CUFFT_INVERSE);
			cufftExecZ2Z(plan2D_y, cuS_z, cuS_z, CUFFT_INVERSE);

			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS_x, cuS_y, cuS_z, cuSquart_x, cuSquart_y, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNc_yx);

			cufftExecZ2D(plan2D_inv_x, cuSquart_x, cuOut_x);
			cufftExecZ2D(plan2D_inv_x, cuSquart_y, cuOut_y);
			cufftExecZ2D(plan2D_inv_x, cuSquart_z, cuOut_z);
		}
	}

#endif
}

void ConvolutionDataCUDA::inverse_fft_2D_2(void)
{
	//Inverse 2D FFT
#if SINGLEPRECISION == 1

	if (!transpose_xy) {

		cufftExecC2C(plan2D_y, cuS2_x, cuS2_x, CUFFT_INVERSE);
		cufftExecC2R(plan2D_inv_x, cuS2_x, cuOut_x);

		cufftExecC2C(plan2D_y, cuS2_y, cuS2_y, CUFFT_INVERSE);
		cufftExecC2R(plan2D_inv_x, cuS2_y, cuOut_y);

		cufftExecC2C(plan2D_y, cuS2_z, cuS2_z, CUFFT_INVERSE);
		cufftExecC2R(plan2D_inv_x, cuS2_z, cuOut_z);
	}
	else {

		if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {

			//interleave transpose_xy operations

			cufftExecC2C(plan2D_y, cuS2_x, cuS2_x, CUFFT_INVERSE);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS2_x, cuSquart_x, cuNcquart_yx, cuNc_yx, cuNc_yx);
			cufftExecC2R(plan2D_inv_x, cuSquart_x, cuOut_x);

			cufftExecC2C(plan2D_y, cuS2_y, cuS2_y, CUFFT_INVERSE);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS2_y, cuSquart_y, cuNcquart_yx, cuNc_yx, cuNc_yx);
			cufftExecC2R(plan2D_inv_x, cuSquart_y, cuOut_y);

			cufftExecC2C(plan2D_y, cuS2_z, cuS2_z, CUFFT_INVERSE);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS2_z, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNc_yx);
			cufftExecC2R(plan2D_inv_x, cuSquart_z, cuOut_z);
		}
		else {

			//transpose_xy operation in a single call

			cufftExecC2C(plan2D_y, cuS2_x, cuS2_x, CUFFT_INVERSE);
			cufftExecC2C(plan2D_y, cuS2_y, cuS2_y, CUFFT_INVERSE);
			cufftExecC2C(plan2D_y, cuS2_z, cuS2_z, CUFFT_INVERSE);

			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS2_x, cuS2_y, cuS2_z, cuSquart_x, cuSquart_y, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNc_yx);

			cufftExecC2R(plan2D_inv_x, cuSquart_x, cuOut_x);
			cufftExecC2R(plan2D_inv_x, cuSquart_y, cuOut_y);
			cufftExecC2R(plan2D_inv_x, cuSquart_z, cuOut_z);
		}
	}

#else

	if (!transpose_xy) {

		cufftExecZ2Z(plan2D_y, cuS2_x, cuS2_x, CUFFT_INVERSE);
		cufftExecZ2D(plan2D_inv_x, cuS2_x, cuOut_x);

		cufftExecZ2Z(plan2D_y, cuS2_y, cuS2_y, CUFFT_INVERSE);
		cufftExecZ2D(plan2D_inv_x, cuS2_y, cuOut_y);

		cufftExecZ2Z(plan2D_y, cuS2_z, cuS2_z, CUFFT_INVERSE);
		cufftExecZ2D(plan2D_inv_x, cuS2_z, cuOut_z);
	}
	else {

		if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {

			//interleave transpose_xy operations

			cufftExecZ2Z(plan2D_y, cuS2_x, cuS2_x, CUFFT_INVERSE);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS2_x, cuSquart_x, cuNcquart_yx, cuNc_yx, cuNc_yx);
			cufftExecZ2D(plan2D_inv_x, cuSquart_x, cuOut_x);

			cufftExecZ2Z(plan2D_y, cuS2_y, cuS2_y, CUFFT_INVERSE);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS2_y, cuSquart_y, cuNcquart_yx, cuNc_yx, cuNc_yx);
			cufftExecZ2D(plan2D_inv_x, cuSquart_y, cuOut_y);

			cufftExecZ2Z(plan2D_y, cuS2_z, cuS2_z, CUFFT_INVERSE);
			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS2_z, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNc_yx);
			cufftExecZ2D(plan2D_inv_x, cuSquart_z, cuOut_z);
		}
		else {

			//transpose_xy operation in a single call

			cufftExecZ2Z(plan2D_y, cuS2_x, cuS2_x, CUFFT_INVERSE);
			cufftExecZ2Z(plan2D_y, cuS2_y, cuS2_y, CUFFT_INVERSE);
			cufftExecZ2Z(plan2D_y, cuS2_z, cuS2_z, CUFFT_INVERSE);

			cu_transpose_xy((N.x / 2 + 1)*n.y, cuS2_x, cuS2_y, cuS2_z, cuSquart_x, cuSquart_y, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNc_yx);

			cufftExecZ2D(plan2D_inv_x, cuSquart_x, cuOut_x);
			cufftExecZ2D(plan2D_inv_x, cuSquart_y, cuOut_y);
			cufftExecZ2D(plan2D_inv_x, cuSquart_z, cuOut_z);
		}
	}

#endif
}

void ConvolutionDataCUDA::forward_fft_3D(void)
{
	//Forward 3D FFT
#if SINGLEPRECISION == 1
	
	if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {

		//interleave transpose_xy operations

		cufftExecR2C(plan3D_fwd_x, cuIn_x, cuSquart_x);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_x, cuS_x, cuNcquart_xy, cuNcquart_xy, cuNc_xy);

		cufftExecR2C(plan3D_fwd_x, cuIn_y, cuSquart_y);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_y, cuS_y, cuNcquart_xy, cuNcquart_xy, cuNc_xy);

		cufftExecR2C(plan3D_fwd_x, cuIn_z, cuSquart_z);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_z, cuS_z, cuNcquart_xy, cuNcquart_xy, cuNc_xy);

		//with pbc enabled n.y will be equal to N.y so no need to launch zero padding kernel
		if (N.y - n.y) {

			cu_zeropad((N.x / 2 + 1)*(N.y - n.y)*n.z, cuS_x, cuS_y, cuS_z, cuNc_yx, cuUpper_y_transposed_region);
		}

		cufftExecC2C(plan3D_y, cuS_x, cuS_x, CUFFT_FORWARD);
		cufftExecC2C(plan3D_z, cuS_x, cuS_x, CUFFT_FORWARD);

		cufftExecC2C(plan3D_y, cuS_y, cuS_y, CUFFT_FORWARD);
		cufftExecC2C(plan3D_z, cuS_y, cuS_y, CUFFT_FORWARD);

		cufftExecC2C(plan3D_y, cuS_z, cuS_z, CUFFT_FORWARD);
		cufftExecC2C(plan3D_z, cuS_z, cuS_z, CUFFT_FORWARD);
	}
	else {

		//transpose_xy operation in a single call

		cufftExecR2C(plan3D_fwd_x, cuIn_x, cuSquart_x);
		cufftExecR2C(plan3D_fwd_x, cuIn_y, cuSquart_y);
		cufftExecR2C(plan3D_fwd_x, cuIn_z, cuSquart_z);

		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_x, cuSquart_y, cuSquart_z, cuS_x, cuS_y, cuS_z, cuNcquart_xy, cuNcquart_xy, cuNc_xy);
		
		//with pbc enabled n.y will be equal to N.y so no need to launch zero padding kernel
		if (N.y - n.y) {

			cu_zeropad((N.x / 2 + 1)*(N.y - n.y)*n.z, cuS_x, cuS_y, cuS_z, cuNc_yx, cuUpper_y_transposed_region);
		}

		cufftExecC2C(plan3D_y, cuS_x, cuS_x, CUFFT_FORWARD);
		cufftExecC2C(plan3D_z, cuS_x, cuS_x, CUFFT_FORWARD);

		cufftExecC2C(plan3D_y, cuS_y, cuS_y, CUFFT_FORWARD);
		cufftExecC2C(plan3D_z, cuS_y, cuS_y, CUFFT_FORWARD);

		cufftExecC2C(plan3D_y, cuS_z, cuS_z, CUFFT_FORWARD);
		cufftExecC2C(plan3D_z, cuS_z, cuS_z, CUFFT_FORWARD);
	}

#else

	if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {

		//interleave transpose_xy operations

		cufftExecD2Z(plan3D_fwd_x, cuIn_x, cuSquart_x);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_x, cuS_x, cuNcquart_xy, cuNcquart_xy, cuNc_xy);

		cufftExecD2Z(plan3D_fwd_x, cuIn_y, cuSquart_y);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_y, cuS_y, cuNcquart_xy, cuNcquart_xy, cuNc_xy);

		cufftExecD2Z(plan3D_fwd_x, cuIn_z, cuSquart_z);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_z, cuS_z, cuNcquart_xy, cuNcquart_xy, cuNc_xy);

		//with pbc enabled n.y will be equal to N.y so no need to launch zero padding kernel
		if (N.y - n.y) {

			cu_zeropad((N.x / 2 + 1)*(N.y - n.y)*n.z, cuS_x, cuS_y, cuS_z, cuNc_yx, cuUpper_y_transposed_region);
		}

		cufftExecZ2Z(plan3D_y, cuS_x, cuS_x, CUFFT_FORWARD);
		cufftExecZ2Z(plan3D_z, cuS_x, cuS_x, CUFFT_FORWARD);

		cufftExecZ2Z(plan3D_y, cuS_y, cuS_y, CUFFT_FORWARD);
		cufftExecZ2Z(plan3D_z, cuS_y, cuS_y, CUFFT_FORWARD);

		cufftExecZ2Z(plan3D_y, cuS_z, cuS_z, CUFFT_FORWARD);
		cufftExecZ2Z(plan3D_z, cuS_z, cuS_z, CUFFT_FORWARD);
	}
	else {

		//transpose_xy operation in a single call

		cufftExecD2Z(plan3D_fwd_x, cuIn_x, cuSquart_x);
		cufftExecD2Z(plan3D_fwd_x, cuIn_y, cuSquart_y);
		cufftExecD2Z(plan3D_fwd_x, cuIn_z, cuSquart_z);

		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_x, cuSquart_y, cuSquart_z, cuS_x, cuS_y, cuS_z, cuNcquart_xy, cuNcquart_xy, cuNc_xy);
		
		//with pbc enabled n.y will be equal to N.y so no need to launch zero padding kernel
		if (N.y - n.y) {

			cu_zeropad((N.x / 2 + 1)*(N.y - n.y)*n.z, cuS_x, cuS_y, cuS_z, cuNc_yx, cuUpper_y_transposed_region);
		}

		cufftExecZ2Z(plan3D_y, cuS_x, cuS_x, CUFFT_FORWARD);
		cufftExecZ2Z(plan3D_z, cuS_x, cuS_x, CUFFT_FORWARD);

		cufftExecZ2Z(plan3D_y, cuS_y, cuS_y, CUFFT_FORWARD);
		cufftExecZ2Z(plan3D_z, cuS_y, cuS_y, CUFFT_FORWARD);

		cufftExecZ2Z(plan3D_y, cuS_z, cuS_z, CUFFT_FORWARD);
		cufftExecZ2Z(plan3D_z, cuS_z, cuS_z, CUFFT_FORWARD);
	}

#endif
}

void ConvolutionDataCUDA::forward_fft_q2D(void)
{
	//Forward 3D FFT
#if SINGLEPRECISION == 1
	
	if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {
		
		//interleave transpose_xy operations
		
		cufftExecR2C(plan3D_fwd_x, cuIn_x, cuSquart_x);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_x, cuS_x, cuNcquart_xy, cuNcquart_xy, cuNc_xy);

		cufftExecR2C(plan3D_fwd_x, cuIn_y, cuSquart_y);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_y, cuS_y, cuNcquart_xy, cuNcquart_xy, cuNc_xy);

		cufftExecR2C(plan3D_fwd_x, cuIn_z, cuSquart_z);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_z, cuS_z, cuNcquart_xy, cuNcquart_xy, cuNc_xy);

		//with pbc enabled n.y will be equal to N.y so no need to launch zero padding kernel
		if (N.y - n.y) {

			cu_zeropad((N.x / 2 + 1)*(N.y - n.y)*n.z, cuS_x, cuS_y, cuS_z, cuNc_yx, cuUpper_y_transposed_region);
		}

		cufftExecC2C(plan3D_y, cuS_x, cuS_x, CUFFT_FORWARD);
		cufftExecC2C(plan3D_y, cuS_y, cuS_y, CUFFT_FORWARD);
		cufftExecC2C(plan3D_y, cuS_z, cuS_z, CUFFT_FORWARD);
	}
	else {

		//transpose_xy operation in a single call

		cufftExecR2C(plan3D_fwd_x, cuIn_x, cuSquart_x);
		cufftExecR2C(plan3D_fwd_x, cuIn_y, cuSquart_y);
		cufftExecR2C(plan3D_fwd_x, cuIn_z, cuSquart_z);

		//it's fine to use cuNc_xy even though its z-dimension is N.z and cuS has n.z dimension for q2D mode
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_x, cuSquart_y, cuSquart_z, cuS_x, cuS_y, cuS_z, cuNcquart_xy, cuNcquart_xy, cuNc_xy);
		
		//with pbc enabled n.y will be equal to N.y so no need to launch zero padding kernel
		if (N.y - n.y) {

			cu_zeropad((N.x / 2 + 1)*(N.y - n.y)*n.z, cuS_x, cuS_y, cuS_z, cuNc_yx, cuUpper_y_transposed_region);
		}

		cufftExecC2C(plan3D_y, cuS_x, cuS_x, CUFFT_FORWARD);
		cufftExecC2C(plan3D_y, cuS_y, cuS_y, CUFFT_FORWARD);
		cufftExecC2C(plan3D_y, cuS_z, cuS_z, CUFFT_FORWARD);
	}

#else

	if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {

		//interleave transpose_xy operations

		cufftExecD2Z(plan3D_fwd_x, cuIn_x, cuSquart_x);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_x, cuS_x, cuNcquart_xy, cuNcquart_xy, cuNc_xy);

		cufftExecD2Z(plan3D_fwd_x, cuIn_y, cuSquart_y);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_y, cuS_y, cuNcquart_xy, cuNcquart_xy, cuNc_xy);

		cufftExecD2Z(plan3D_fwd_x, cuIn_z, cuSquart_z);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_z, cuS_z, cuNcquart_xy, cuNcquart_xy, cuNc_xy);

		//with pbc enabled n.y will be equal to N.y so no need to launch zero padding kernel
		if (N.y - n.y) {

			cu_zeropad((N.x / 2 + 1)*(N.y - n.y)*n.z, cuS_x, cuS_y, cuS_z, cuNc_yx, cuUpper_y_transposed_region);
		}

		cufftExecZ2Z(plan3D_y, cuS_x, cuS_x, CUFFT_FORWARD);
		cufftExecZ2Z(plan3D_y, cuS_y, cuS_y, CUFFT_FORWARD);
		cufftExecZ2Z(plan3D_y, cuS_z, cuS_z, CUFFT_FORWARD);
	}
	else {

		//transpose_xy operation in a single call

		cufftExecD2Z(plan3D_fwd_x, cuIn_x, cuSquart_x);
		cufftExecD2Z(plan3D_fwd_x, cuIn_y, cuSquart_y);
		cufftExecD2Z(plan3D_fwd_x, cuIn_z, cuSquart_z);

		//it's fine to use cuNc_xy even though its z-dimension is N.z and cuS has n.z dimension for q2D mode
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuSquart_x, cuSquart_y, cuSquart_z, cuS_x, cuS_y, cuS_z, cuNcquart_xy, cuNcquart_xy, cuNc_xy);
		
		//with pbc enabled n.y will be equal to N.y so no need to launch zero padding kernel
		if (N.y - n.y) {

			cu_zeropad((N.x / 2 + 1)*(N.y - n.y)*n.z, cuS_x, cuS_y, cuS_z, cuNc_yx, cuUpper_y_transposed_region);
		}

		cufftExecZ2Z(plan3D_y, cuS_x, cuS_x, CUFFT_FORWARD);
		cufftExecZ2Z(plan3D_y, cuS_y, cuS_y, CUFFT_FORWARD);
		cufftExecZ2Z(plan3D_y, cuS_z, cuS_z, CUFFT_FORWARD);
	}

#endif
}

void ConvolutionDataCUDA::inverse_fft_3D(void)
{
	//Inverse 3D FFT
#if SINGLEPRECISION == 1
	
	if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {

		//interleave transpose_xy operations

		cufftExecC2C(plan3D_z, cuS_x, cuS_x, CUFFT_INVERSE);
		cufftExecC2C(plan3D_y, cuS_x, cuS_x, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_x, cuSquart_x, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecC2R(plan3D_inv_x, cuSquart_x, cuOut_x);

		cufftExecC2C(plan3D_z, cuS_y, cuS_y, CUFFT_INVERSE);
		cufftExecC2C(plan3D_y, cuS_y, cuS_y, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_y, cuSquart_y, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecC2R(plan3D_inv_x, cuSquart_y, cuOut_y);

		cufftExecC2C(plan3D_z, cuS_z, cuS_z, CUFFT_INVERSE);
		cufftExecC2C(plan3D_y, cuS_z, cuS_z, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_z, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecC2R(plan3D_inv_x, cuSquart_z, cuOut_z);
	}
	else {

		//transpose_xy operation in a single call

		cufftExecC2C(plan3D_z, cuS_x, cuS_x, CUFFT_INVERSE);
		cufftExecC2C(plan3D_y, cuS_x, cuS_x, CUFFT_INVERSE);

		cufftExecC2C(plan3D_z, cuS_y, cuS_y, CUFFT_INVERSE);
		cufftExecC2C(plan3D_y, cuS_y, cuS_y, CUFFT_INVERSE);

		cufftExecC2C(plan3D_z, cuS_z, cuS_z, CUFFT_INVERSE);
		cufftExecC2C(plan3D_y, cuS_z, cuS_z, CUFFT_INVERSE);

		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_x, cuS_y, cuS_z, cuSquart_x, cuSquart_y, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNcquart_yx);

		cufftExecC2R(plan3D_inv_x, cuSquart_x, cuOut_x);
		cufftExecC2R(plan3D_inv_x, cuSquart_y, cuOut_y);
		cufftExecC2R(plan3D_inv_x, cuSquart_z, cuOut_z);
	}
	
#else

	if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {

		//interleave transpose_xy operations

		cufftExecZ2Z(plan3D_z, cuS_x, cuS_x, CUFFT_INVERSE);
		cufftExecZ2Z(plan3D_y, cuS_x, cuS_x, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_x, cuSquart_x, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecZ2D(plan3D_inv_x, cuSquart_x, cuOut_x);

		cufftExecZ2Z(plan3D_z, cuS_y, cuS_y, CUFFT_INVERSE);
		cufftExecZ2Z(plan3D_y, cuS_y, cuS_y, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_y, cuSquart_y, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecZ2D(plan3D_inv_x, cuSquart_y, cuOut_y);

		cufftExecZ2Z(plan3D_z, cuS_z, cuS_z, CUFFT_INVERSE);
		cufftExecZ2Z(plan3D_y, cuS_z, cuS_z, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_z, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecZ2D(plan3D_inv_x, cuSquart_z, cuOut_z);
	}
	else {

		//transpose_xy operation in a single call

		cufftExecZ2Z(plan3D_z, cuS_x, cuS_x, CUFFT_INVERSE);
		cufftExecZ2Z(plan3D_y, cuS_x, cuS_x, CUFFT_INVERSE);

		cufftExecZ2Z(plan3D_z, cuS_y, cuS_y, CUFFT_INVERSE);
		cufftExecZ2Z(plan3D_y, cuS_y, cuS_y, CUFFT_INVERSE);

		cufftExecZ2Z(plan3D_z, cuS_z, cuS_z, CUFFT_INVERSE);
		cufftExecZ2Z(plan3D_y, cuS_z, cuS_z, CUFFT_INVERSE);

		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_x, cuS_y, cuS_z, cuSquart_x, cuSquart_y, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNcquart_yx);

		cufftExecZ2D(plan3D_inv_x, cuSquart_x, cuOut_x);
		cufftExecZ2D(plan3D_inv_x, cuSquart_y, cuOut_y);
		cufftExecZ2D(plan3D_inv_x, cuSquart_z, cuOut_z);
	}

#endif
}

void ConvolutionDataCUDA::inverse_fft_3D_2(void)
{
	//Inverse 3D FFT
#if SINGLEPRECISION == 1

	if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {

		//interleave transpose_xy operations

		cufftExecC2C(plan3D_z, cuS2_x, cuS2_x, CUFFT_INVERSE);
		cufftExecC2C(plan3D_y, cuS2_x, cuS2_x, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS2_x, cuSquart_x, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecC2R(plan3D_inv_x, cuSquart_x, cuOut_x);

		cufftExecC2C(plan3D_z, cuS2_y, cuS2_y, CUFFT_INVERSE);
		cufftExecC2C(plan3D_y, cuS2_y, cuS2_y, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS2_y, cuSquart_y, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecC2R(plan3D_inv_x, cuSquart_y, cuOut_y);

		cufftExecC2C(plan3D_z, cuS2_z, cuS2_z, CUFFT_INVERSE);
		cufftExecC2C(plan3D_y, cuS2_z, cuS2_z, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS2_z, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecC2R(plan3D_inv_x, cuSquart_z, cuOut_z);
	}
	else {

		//transpose_xy operation in a single call

		cufftExecC2C(plan3D_z, cuS2_x, cuS2_x, CUFFT_INVERSE);
		cufftExecC2C(plan3D_y, cuS2_x, cuS2_x, CUFFT_INVERSE);

		cufftExecC2C(plan3D_z, cuS2_y, cuS2_y, CUFFT_INVERSE);
		cufftExecC2C(plan3D_y, cuS2_y, cuS2_y, CUFFT_INVERSE);

		cufftExecC2C(plan3D_z, cuS2_z, cuS2_z, CUFFT_INVERSE);
		cufftExecC2C(plan3D_y, cuS2_z, cuS2_z, CUFFT_INVERSE);

		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS2_x, cuS2_y, cuS2_z, cuSquart_x, cuSquart_y, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNcquart_yx);

		cufftExecC2R(plan3D_inv_x, cuSquart_x, cuOut_x);
		cufftExecC2R(plan3D_inv_x, cuSquart_y, cuOut_y);
		cufftExecC2R(plan3D_inv_x, cuSquart_z, cuOut_z);
	}

#else

	if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {

		//interleave transpose_xy operations

		cufftExecZ2Z(plan3D_z, cuS2_x, cuS2_x, CUFFT_INVERSE);
		cufftExecZ2Z(plan3D_y, cuS2_x, cuS2_x, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS2_x, cuSquart_x, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecZ2D(plan3D_inv_x, cuSquart_x, cuOut_x);

		cufftExecZ2Z(plan3D_z, cuS2_y, cuS2_y, CUFFT_INVERSE);
		cufftExecZ2Z(plan3D_y, cuS2_y, cuS2_y, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS2_y, cuSquart_y, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecZ2D(plan3D_inv_x, cuSquart_y, cuOut_y);

		cufftExecZ2Z(plan3D_z, cuS2_z, cuS2_z, CUFFT_INVERSE);
		cufftExecZ2Z(plan3D_y, cuS2_z, cuS2_z, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS2_z, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecZ2D(plan3D_inv_x, cuSquart_z, cuOut_z);
	}
	else {

		//transpose_xy operation in a single call

		cufftExecZ2Z(plan3D_z, cuS2_x, cuS2_x, CUFFT_INVERSE);
		cufftExecZ2Z(plan3D_y, cuS2_x, cuS2_x, CUFFT_INVERSE);

		cufftExecZ2Z(plan3D_z, cuS2_y, cuS2_y, CUFFT_INVERSE);
		cufftExecZ2Z(plan3D_y, cuS2_y, cuS2_y, CUFFT_INVERSE);

		cufftExecZ2Z(plan3D_z, cuS2_z, cuS2_z, CUFFT_INVERSE);
		cufftExecZ2Z(plan3D_y, cuS2_z, cuS2_z, CUFFT_INVERSE);

		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS2_x, cuS2_y, cuS2_z, cuSquart_x, cuSquart_y, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNcquart_yx);

		cufftExecZ2D(plan3D_inv_x, cuSquart_x, cuOut_x);
		cufftExecZ2D(plan3D_inv_x, cuSquart_y, cuOut_y);
		cufftExecZ2D(plan3D_inv_x, cuSquart_z, cuOut_z);
	}

#endif
}

void ConvolutionDataCUDA::inverse_fft_q2D(void)
{
	//Inverse 3D FFT
#if SINGLEPRECISION == 1

	if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {
		
		//interleave transpose_xy operations

		cufftExecC2C(plan3D_y, cuS_x, cuS_x, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_x, cuSquart_x, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecC2R(plan3D_inv_x, cuSquart_x, cuOut_x);

		cufftExecC2C(plan3D_y, cuS_y, cuS_y, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_y, cuSquart_y, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecC2R(plan3D_inv_x, cuSquart_y, cuOut_y);

		cufftExecC2C(plan3D_y, cuS_z, cuS_z, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_z, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecC2R(plan3D_inv_x, cuSquart_z, cuOut_z);
	}
	else {

		//transpose_xy operation in a single call

		cufftExecC2C(plan3D_y, cuS_x, cuS_x, CUFFT_INVERSE);
		cufftExecC2C(plan3D_y, cuS_y, cuS_y, CUFFT_INVERSE);
		cufftExecC2C(plan3D_y, cuS_z, cuS_z, CUFFT_INVERSE);

		//it's fine to use cuNc_yx even though its z-dimension is N.z and cuS has n.z dimension for q2D mode
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_x, cuS_y, cuS_z, cuSquart_x, cuSquart_y, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNcquart_yx);

		cufftExecC2R(plan3D_inv_x, cuSquart_x, cuOut_x);
		cufftExecC2R(plan3D_inv_x, cuSquart_y, cuOut_y);
		cufftExecC2R(plan3D_inv_x, cuSquart_z, cuOut_z);
	}

#else

	if (N.x / 2 >= TRANSPOSE_XY_YTHRESHOLD) {

		//interleave transpose_xy operations

		cufftExecZ2Z(plan3D_y, cuS_x, cuS_x, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_x, cuSquart_x, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecZ2D(plan3D_inv_x, cuSquart_x, cuOut_x);

		cufftExecZ2Z(plan3D_y, cuS_y, cuS_y, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_y, cuSquart_y, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecZ2D(plan3D_inv_x, cuSquart_y, cuOut_y);

		cufftExecZ2Z(plan3D_y, cuS_z, cuS_z, CUFFT_INVERSE);
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_z, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNcquart_yx);
		cufftExecZ2D(plan3D_inv_x, cuSquart_z, cuOut_z);
	}
	else {

		//transpose_xy operation in a single call

		cufftExecZ2Z(plan3D_y, cuS_x, cuS_x, CUFFT_INVERSE);
		cufftExecZ2Z(plan3D_y, cuS_y, cuS_y, CUFFT_INVERSE);
		cufftExecZ2Z(plan3D_y, cuS_z, cuS_z, CUFFT_INVERSE);

		//it's fine to use cuNc_yx even though its z-dimension is N.z and cuS has n.z dimension for q2D mode
		cu_transpose_xy((N.x / 2 + 1)*n.y*n.z, cuS_x, cuS_y, cuS_z, cuSquart_x, cuSquart_y, cuSquart_z, cuNcquart_yx, cuNc_yx, cuNcquart_yx);

		cufftExecZ2D(plan3D_inv_x, cuSquart_x, cuOut_x);
		cufftExecZ2D(plan3D_inv_x, cuSquart_y, cuOut_y);
		cufftExecZ2D(plan3D_inv_x, cuSquart_z, cuOut_z);
	}

#endif
}

template void ConvolutionDataCUDA::FinishConvolution_Set(cu_obj<cuVEC<cuReal3>>& In, cu_obj<cuVEC<cuReal3>>& Out, cu_obj<cuReal>& energy, bool get_energy);
template void ConvolutionDataCUDA::FinishConvolution_Set(cu_obj<cuVEC<cuReal3>>& In, cu_obj<cuVEC_VC<cuReal3>>& Out, cu_obj<cuReal>& energy, bool get_energy);
template void ConvolutionDataCUDA::FinishConvolution_Set(cu_obj<cuVEC_VC<cuReal3>>& In, cu_obj<cuVEC<cuReal3>>& Out, cu_obj<cuReal>& energy, bool get_energy);
template void ConvolutionDataCUDA::FinishConvolution_Set(cu_obj<cuVEC_VC<cuReal3>>& In, cu_obj<cuVEC_VC<cuReal3>>& Out, cu_obj<cuReal>& energy, bool get_energy);

//Copy convolution result (in cuS arrays) to output and obtain energy value : product of In with Out times -MU0 / (2 * non_empty_points), where non_empty_points = In.get_nonempty_points();
template <typename cuVECIn, typename cuVECOut>
void ConvolutionDataCUDA::FinishConvolution_Set(cu_obj<cuVECIn>& In, cu_obj<cuVECOut>& Out, cu_obj<cuReal>& energy, bool get_energy)
{
	if (get_energy) {

		//set aux_integer in In
		In()->count_nonempty_cells(n.dim());
	}

	cuFFTArrays_to_Out_Set_forOutOfPlace << < (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (*In(), *Out(), cuOut_x, cuOut_y, cuOut_z, cuN, energy, get_energy);
}

template void ConvolutionDataCUDA::FinishConvolution_Add(cu_obj<cuVEC<cuReal3>>& In, cu_obj<cuVEC<cuReal3>>& Out, cu_obj<cuReal>& energy, bool get_energy);
template void ConvolutionDataCUDA::FinishConvolution_Add(cu_obj<cuVEC<cuReal3>>& In, cu_obj<cuVEC_VC<cuReal3>>& Out, cu_obj<cuReal>& energy, bool get_energy);
template void ConvolutionDataCUDA::FinishConvolution_Add(cu_obj<cuVEC_VC<cuReal3>>& In, cu_obj<cuVEC<cuReal3>>& Out, cu_obj<cuReal>& energy, bool get_energy);
template void ConvolutionDataCUDA::FinishConvolution_Add(cu_obj<cuVEC_VC<cuReal3>>& In, cu_obj<cuVEC_VC<cuReal3>>& Out, cu_obj<cuReal>& energy, bool get_energy);

template <typename cuVECIn, typename cuVECOut>
void ConvolutionDataCUDA::FinishConvolution_Add(cu_obj<cuVECIn>& In, cu_obj<cuVECOut>& Out, cu_obj<cuReal>& energy, bool get_energy)
{
	if (get_energy) {

		//set aux_integer in In
		In()->count_nonempty_cells(n.dim());
	}

	cuFFTArrays_to_Out_Add_forOutOfPlace << < (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (*In(), *Out(), cuOut_x, cuOut_y, cuOut_z, cuN, energy, get_energy);
}

template void ConvolutionDataCUDA::FinishConvolution_Set(cu_obj<cuVEC<cuReal3>>& In, cu_obj<cuVEC<cuReal3>>& Out, cu_obj<cuReal>& energy, cu_obj<cuReal>& energy_weight);
template void ConvolutionDataCUDA::FinishConvolution_Set(cu_obj<cuVEC<cuReal3>>& In, cu_obj<cuVEC_VC<cuReal3>>& Out, cu_obj<cuReal>& energy, cu_obj<cuReal>& energy_weight);
template void ConvolutionDataCUDA::FinishConvolution_Set(cu_obj<cuVEC_VC<cuReal3>>& In, cu_obj<cuVEC<cuReal3>>& Out, cu_obj<cuReal>& energy, cu_obj<cuReal>& energy_weight);
template void ConvolutionDataCUDA::FinishConvolution_Set(cu_obj<cuVEC_VC<cuReal3>>& In, cu_obj<cuVEC_VC<cuReal3>>& Out, cu_obj<cuReal>& energy, cu_obj<cuReal>& energy_weight);

//Copy convolution result (in cuS arrays) to output and obtain energy value : weighted product of In with Out times -MU0 / (2 * non_empty_points), where non_empty_points = In.get_nonempty_points();
template <typename cuVECIn, typename cuVECOut>
void ConvolutionDataCUDA::FinishConvolution_Set(cu_obj<cuVECIn>& In, cu_obj<cuVECOut>& Out, cu_obj<cuReal>& energy, cu_obj<cuReal>& energy_weight)
{
	//set aux_integer in In
	In()->count_nonempty_cells(n.dim());

	cuFFTArrays_to_Out_Set_weighted_forOutOfPlace << < (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (*In(), *Out(), cuOut_x, cuOut_y, cuOut_z, cuN, energy, energy_weight);
}

template void ConvolutionDataCUDA::FinishConvolution_Add(cu_obj<cuVEC<cuReal3>>& In, cu_obj<cuVEC<cuReal3>>& Out, cu_obj<cuReal>& energy, cu_obj<cuReal>& energy_weight);
template void ConvolutionDataCUDA::FinishConvolution_Add(cu_obj<cuVEC<cuReal3>>& In, cu_obj<cuVEC_VC<cuReal3>>& Out, cu_obj<cuReal>& energy, cu_obj<cuReal>& energy_weight);
template void ConvolutionDataCUDA::FinishConvolution_Add(cu_obj<cuVEC_VC<cuReal3>>& In, cu_obj<cuVEC<cuReal3>>& Out, cu_obj<cuReal>& energy, cu_obj<cuReal>& energy_weight);
template void ConvolutionDataCUDA::FinishConvolution_Add(cu_obj<cuVEC_VC<cuReal3>>& In, cu_obj<cuVEC_VC<cuReal3>>& Out, cu_obj<cuReal>& energy, cu_obj<cuReal>& energy_weight);

//Add convolution result (in cuS arrays) to output and obtain energy value : weighted product of In with Out times -MU0 / (2 * non_empty_points), where non_empty_points = In.get_nonempty_points();
template <typename cuVECIn, typename cuVECOut>
void ConvolutionDataCUDA::FinishConvolution_Add(cu_obj<cuVECIn>& In, cu_obj<cuVECOut>& Out, cu_obj<cuReal>& energy, cu_obj<cuReal>& energy_weight)
{
	//set aux_integer in In
	In()->count_nonempty_cells(n.dim());

	cuFFTArrays_to_Out_Add_weighted_forOutOfPlace << < (n.dim() + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> > (*In(), *Out(), cuOut_x, cuOut_y, cuOut_z, cuN, energy, energy_weight);
}

#endif