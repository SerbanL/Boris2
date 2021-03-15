#include "stdafx.h"
#include "DWRunTimeFit.h"

//Fit the extracted profile for position and width, assuming the magnetization component follows f(x) = [ (As - Ae) * tanh(-PI * (x - x0) / dw) + (As + Ae) ] / 2
//Here As, Ae are the start and end values - profile must be long enough to include at least DWPOS_ENDSTENCIL (length ratio) flat parts of the tanh profile
//x0 is the centre position relative to start of profile
//dw is the domain wall width
//xy_data contains profile as x coordinates and corresponding y values
//return x0, dw
DBL2 DWPosWidth::FitDomainWall(std::vector<DBL2>& xy_data)
{
	if (xy_data.size() < DWPOS_MINPROFILEPOINTS) return DBL2();

	//1. Find end values

	double As = 0.0, Ae = 0.0;
	int num_end_points = xy_data.size() * DWPOS_ENDSTENCIL;

#pragma omp parallel for reduction(+:As)
	for (int idx = 0; idx < num_end_points; idx++) {

		As += xy_data[idx].j / num_end_points;
	}

#pragma omp parallel for reduction(+:Ae)
	for (int idx = xy_data.size() - num_end_points; idx < xy_data.size(); idx++) {

		Ae += xy_data[idx].j / num_end_points;
	}

	double ycentre = (As + Ae) / 2;
	double amplitude = abs(As - Ae);
	if (amplitude == 0.0) return DBL2();

	//2. Find centre position using calculated ycentre and amplitude values

	auto stencil_average = [](std::vector<DBL2>& xy_data, int pidx, int num_stencil_points) -> double
	{
		double value = 0.0;
		for (int idx = pidx - num_stencil_points; idx < pidx + num_stencil_points; idx++) {

			value += xy_data[idx].j / (2 * num_stencil_points + 1);
		}
		return value;
	};

	int num_stencil_points = xy_data.size() * DWPOS_STENCIL;

	//produce smoothed xy data using nearest neighbor average
	if (xy_data_smoothed.size() != xy_data.size() - 2 * num_stencil_points) if (!malloc_vector(xy_data_smoothed, xy_data.size() - 2 * num_stencil_points)) return DBL2();

#pragma omp parallel for
	for (int idx = 0; idx < xy_data_smoothed.size(); idx++) {

		double v1 = stencil_average(xy_data, idx + num_stencil_points, num_stencil_points);
		xy_data_smoothed[idx] = DBL2(xy_data[idx + num_stencil_points].i, v1);
	}

	//now find all crossing points from smoothed data and average them to find x0
	//Note, working on smoothed data to find x0 works well, because if the profile follows the expected tanh function, smoothing will not change the x0 point so it's safe to use a large stencil (in fact needed to extract good x0 values at high temperatures).
	//You do not want to use the smoothed data to extract dw however since it doesn't follow the correct tanh profile.
	double x0 = 0.0;
	int x0_points = 0;
#pragma omp parallel for reduction(+:x0, x0_points)
	for (int idx = 0; idx < xy_data_smoothed.size() - 1; idx++) {

		if ((xy_data_smoothed[idx].j - ycentre) * (xy_data_smoothed[idx + 1].j - ycentre) < 0) {

			x0 += interpolate(
				DBL2(xy_data_smoothed[idx].j, xy_data_smoothed[idx].i), 
				DBL2(xy_data_smoothed[idx + 1].j, xy_data_smoothed[idx + 1].i), 
				ycentre);

			x0_points++;
		}
	}

	if (x0_points) x0 /= x0_points;
	else return DBL2();

	double length = xy_data.back().i - xy_data.front().i;

	//If x0 is not within reasonable bounds then fail
	if (x0 < length * DWPOS_ENDSTENCIL || x0 > length * (1.0 - DWPOS_ENDSTENCIL)) return DBL2();

	//3. Find DW width using x0, As, Ae values
	/*
	//OLD METHOD WITH EQUAL WEIGHTS
	double dw = 0.0;
	double c = ycentre, a = (As - Ae) / 2;
	int num_dw_points = 0;

#pragma omp parallel for reduction(+:dw, num_dw_points)
	for (int idx = num_end_points; idx < xy_data.size() - num_end_points; idx++) {

		//function is f(x) = [ (As - Ae) * tanh(-PI * (x - x0) / dw) + (As + Ae) ] / 2
		//at each point find f(x), then solve for a dw value. at the end average all dw values. 
		//Works very well even at high temperatures, same answer as LMA fitting, but far more reliable and fast.
		//Note at high temperatures fitting can fail sometimes (detected through bounds checks), so just don't include the failed fits in the final average.
		double nval = (xy_data[idx].j - c) / a;
		if (abs(nval) < DWPOS_YTHRESHOLD_MAX && abs(nval) > DWPOS_YTHRESHOLD_MIN) {

			double value = -PI * (xy_data[idx].i - x0) / atanh(nval);
			if (value > 0.0) {

				dw += value;
				num_dw_points++;
			}
		}
	}
	
	if (num_dw_points) dw /= num_dw_points;
	else return DBL2();
	*/
	
	double dw = 0.0, w = 0.0;
	double c = ycentre; 
	double m = (As - Ae) / 2;

#pragma omp parallel for reduction(+:dw, w)
	for (int idx = num_end_points; idx < xy_data.size() - num_end_points; idx++) {

		//function is f(x) = [ (As - Ae) * tanh(-PI * (x - x0) / dw) + (As + Ae) ] / 2 = m * tanh(-PI*(x - x0) / dw) + c
		//at each point find f(x), then solve for a dw value and an attached weight obtained from least squares equation. Obtain final domain wall width using weighted average.
		double nval = (xy_data[idx].j - c) / m;
		if (abs(nval) < DWPOS_YTHRESHOLD_MAX && abs(nval) > DWPOS_YTHRESHOLD_MIN) {

			//domain wall width for this point
			double dw_i = abs(PI * (xy_data[idx].i - x0) / atanh(nval));

			if (dw_i) {

				//function evaluated for dw_i
				double f_i = m * tanh(-PI * (xy_data[idx].i - x0) / dw_i) + c;
				//weight for dw_i (obtained from least squares equation)
				double w_i = abs((m*m - f_i * f_i) * (xy_data[idx].i - x0));

				//total domain wall as weighted sum
				dw += w_i * dw_i;
				//total weight
				w += w_i;
			}
		}
	}

	if (w) dw /= w;
	else return DBL2();
	
	//if dw width is not within reasonable bounds then fail
	if (dw > length * (1.0 - 2 * DWPOS_ENDSTENCIL) || dw < 0) return DBL2();

	//remember x0 is relative to start of profile, so caller will have to adjust for this to make it relative to start of mesh
	return DBL2(x0, dw);
}