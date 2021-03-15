#include "DWRunTimeFitCUDA.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.cuh"

///////////////////////////////////////////////////////////////
// AUXILIARY

__global__ void Zero_DWRunTimeFitCUDA_Values(cuBReal& As, cuBReal& Ae, cuBReal& x0, cuBReal& dw, size_t& av_points_x0, cuBReal& weight)
{
	if (threadIdx.x == 0) As = 0.0;
	else if (threadIdx.x == 1) Ae = 0.0;
	else if (threadIdx.x == 2) x0 = 0.0;
	else if (threadIdx.x == 3) dw = 0.0;
	else if (threadIdx.x == 4) av_points_x0 = 0;
	else if (threadIdx.x == 5) weight = 0.0;
}

///////////////////////////////////////////////////////////////
// PARALLEL MONTE-CARLO METROPOLIS - WITH REDUCTION

//reduce start and end values : finish average at the end by dividing by av_points_s and av_points_e
__global__ void DWRunTimeFitCUDA_endpoints_kernel(size_t size, cuReal2* pxy_data, cuBReal& As, cuBReal& Ae)
{
	//kernel launched with at least 2*num_end_points : lower for As, upper for Ae. Better than launching 2 separate kernels as these are usually very small.
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal value_s = 0.0, value_e = 0.0;

	int num_end_points = size * (cuBReal)DWPOS_ENDSTENCIL;

	if (idx < num_end_points) {

		value_s = pxy_data[idx].j / num_end_points;
	}
	else if (idx < 2 * num_end_points) {

		value_e = pxy_data[size - 2 * num_end_points + idx].j / num_end_points;
	}

	reduction_sum(0, 1, &value_s, As);
	reduction_sum(0, 1, &value_e, Ae);
}

__global__ void DWRunTimeFitCUDA_smoothing_kernel(cuReal2* pxy_data, size_t size_smoothed, int num_stencil_points, cuReal2* pxy_data_smoothed)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < size_smoothed) {

		idx += num_stencil_points;

		cuBReal value = 0.0;
		for (int sidx = idx - num_stencil_points; sidx < idx + num_stencil_points; sidx++) {

			value += pxy_data[sidx].j / (2 * num_stencil_points + 1);
		}

		pxy_data_smoothed[idx - num_stencil_points].i = pxy_data[idx].i;
		pxy_data_smoothed[idx - num_stencil_points].j = value;
	}
}

__global__ void DWRunTimeFitCUDA_x0_kernel(size_t size_smoothed, cuReal2* pxy_data_smoothed, cuBReal& As, cuBReal& Ae, cuBReal& x0, size_t& av_points_x0)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	cuBReal x0value = 0.0;
	bool include_in_average = false;

	if (idx < size_smoothed - 1) {

		cuBReal ycentre = (As + Ae) / 2;

		if ((pxy_data_smoothed[idx].j - ycentre) * (pxy_data_smoothed[idx + 1].j - ycentre) < 0) {

			x0value = cu_interpolate(
				cuReal2(pxy_data_smoothed[idx].j, pxy_data_smoothed[idx].i),
				cuReal2(pxy_data_smoothed[idx + 1].j, pxy_data_smoothed[idx + 1].i),
				ycentre);

			include_in_average = true;
		}
	}

	reduction_avg(0, 1, &x0value, x0, av_points_x0, include_in_average);
}

__global__ void DWRunTimeFitCUDA_dw_kernel(size_t size, cuReal2* pxy_data, cuBReal& As, cuBReal& Ae, cuBReal& x0, size_t& av_points_x0, cuBReal& dw, cuBReal& weight)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	//domain wall width value times weight
	cuBReal w_i_dw_i = 0.0;
	//weight only (need to reduce to find total weight)
	cuBReal w_i = 0.0;
	bool include_in_average = false;

	int num_end_points = size * (cuBReal)DWPOS_ENDSTENCIL;

	if (idx < size - num_end_points && idx >= num_end_points) {
		
		cuBReal c = (As + Ae) / 2;
		cuBReal m = (As - Ae) / 2;

		//function is f(x) = [ (As - Ae) * tanh(-PI * (x - x0) / dw) + (As + Ae) ] / 2 = m * tanh(-PI*(x - x0) / dw) + c
		//at each point find f(x), then solve for a dw value and an attached weight obtained from least squares equation. Obtain final domain wall width using weighted average.
		cuBReal nval = (pxy_data[idx].j - c) / m;
		if (abs(nval) < (cuBReal)DWPOS_YTHRESHOLD_MAX && abs(nval) > (cuBReal)DWPOS_YTHRESHOLD_MIN) {

			//domain wall width for this point
			cuBReal dw_i = fabs(PI * (pxy_data[idx].i - x0) / atanh(nval));

			if (dw_i) {

				//function evaluated for dw_i
				cuBReal f_i = m * tanh(-PI * (pxy_data[idx].i - x0) / dw_i) + c;
				//weight for dw_i (obtained from least squares equation)
				w_i = fabs((m*m - f_i * f_i) * (pxy_data[idx].i - x0));

				//total domain wall as weighted sum
				w_i_dw_i = w_i * dw_i;
				include_in_average = true;
			}
		}
	}

	//reduce numerator and denominator for weighted sum
	reduction_sum(0, 1, &w_i_dw_i, dw, include_in_average);
	reduction_sum(0, 1, &w_i, weight, include_in_average);
}

//Fit the extracted profile for position and width, assuming the magnetization component follows f(x) = [ (As - Ae) * tanh(-PI * (x - x0) / dw) + (As + Ae) ] / 2
//Here As, Ae are the start and end values - profile must be long enough to include at least DWPOS_ENDSTENCIL (length ratio) flat parts of the tanh profile
//x0 is the centre position relative to start of profile
//dw is the domain wall width
//xy_data contains profile as x coordinates and corresponding y values
//return x0, dw
cuReal2 DWPosWidthCUDA::FitDomainWallCUDA(double length)
{
	//1. Find end values

	size_t size = xy_data.size();
	if (size < DWPOS_MINPROFILEPOINTS) return cuReal2();

	//zero all reduction values
	Zero_DWRunTimeFitCUDA_Values <<< 1, CUDATHREADS >>> (As, Ae, x0, dw, av_points_x0, weight);

	int num_end_points = size * DWPOS_ENDSTENCIL;
	DWRunTimeFitCUDA_endpoints_kernel <<< (2 * num_end_points + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(size, xy_data, As, Ae);

	//2. Find centre position using calculated ycentre and amplitude values

	int num_stencil_points = size * DWPOS_STENCIL;

	//produce smoothed xy data using nearest neighbor average
	if (xy_data_smoothed_size != size - 2 * num_stencil_points) {

		if (!xy_data_smoothed.resize(size - 2 * num_stencil_points)) return cuReal2();
		else xy_data_smoothed_size = size - 2 * num_stencil_points;
	}

	//smooth
	DWRunTimeFitCUDA_smoothing_kernel <<< (xy_data_smoothed_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(xy_data, xy_data_smoothed_size, num_stencil_points, xy_data_smoothed);

	//find x0
	DWRunTimeFitCUDA_x0_kernel <<< (xy_data_smoothed_size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >> >
		(xy_data_smoothed_size, xy_data_smoothed, As, Ae, x0, av_points_x0);

	//If x0 is not within reasonable bounds then fail
	cuBReal x0_cpu = x0.to_cpu();
	if (x0_cpu < length * DWPOS_ENDSTENCIL || x0_cpu > length * (1.0 - DWPOS_ENDSTENCIL)) return cuReal2();

	//3. Find DW width using x0, As, Ae values

	//find dw
	DWRunTimeFitCUDA_dw_kernel <<< (size + CUDATHREADS) / CUDATHREADS, CUDATHREADS >>>
		(size, xy_data, As, Ae, x0, av_points_x0, dw, weight);

	cuBReal weight_cpu = weight.to_cpu();
	cuBReal dw_cpu = dw.to_cpu();

	if (weight_cpu) dw_cpu /= weight_cpu;
	else return cuReal2();

	//if dw width is not within reasonable bounds then fail
	if (dw_cpu > length * (1.0 - 2 * DWPOS_ENDSTENCIL) || dw_cpu < 0) return cuReal2();

	//remember x0 is relative to start of profile, so caller will have to adjust for this to make it relative to start of mesh
	return cuReal2(x0_cpu, dw_cpu);
}

#endif