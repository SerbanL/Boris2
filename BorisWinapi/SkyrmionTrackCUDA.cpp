#include "stdafx.h"
#include "SkyrmionTrack.h"
#include "SimulationData.h"

//-------------------------------------------------- Get_skyshift : CUDA 

#if COMPILECUDA == 1

DBL2 SkyrmionTrack::Get_skyshiftCUDA(size_t size, DBL3 h, cu_obj<cuVEC_VC<cuReal3>>& M, Rect skyRect)
{
	//must have a set rectangle
	if (skyRect.IsNull()) return DBL2();

	int skyTrack_idx = Get_skyTrack_index(skyRect);

	//shift the skyrmion rectangle to current tracking position
	skyRect += DBL3(skyTrack_ShiftLast[skyTrack_idx].x, skyTrack_ShiftLast[skyTrack_idx].y, 0.0);

	//Find M averages in the 4 skyRect xy-plane quadrants
	DBL3 bottom_left = (DBL3)M()->average_nonempty(size, skyRect.get_quadrant_bl());
	DBL3 bottom_right = (DBL3)M()->average_nonempty(size, skyRect.get_quadrant_br());
	DBL3 top_left = (DBL3)M()->average_nonempty(size, skyRect.get_quadrant_tl());
	DBL3 top_right = (DBL3)M()->average_nonempty(size, skyRect.get_quadrant_tr());

	//the new shift value
	DBL2 skyTrack_newShift = DBL2();

	//if left half z-component value modulus has increased compared to the right half then it contains less of the skyrmion ring -> skyrmion must have shifted along +x so follow it
	if (mod(bottom_left.z + top_left.z) > mod(bottom_right.z + top_right.z)) {

		skyTrack_newShift.x += h.x;
	}
	else {

		skyTrack_newShift.x -= h.x;
	}

	//if bottom half z-component value modulus has increased compared to the top half then it contains less of the skyrmion ring -> skyrmion must have shifted along +y so follow it
	if (mod(bottom_left.z + bottom_right.z) > mod(top_left.z + top_right.z)) {

		skyTrack_newShift.y += h.y;
	}
	else {

		skyTrack_newShift.y -= h.y;
	}

	//set actual total shift as average of new and last total shift values - this eliminates tracking oscillations
	skyTrack_Shift[skyTrack_idx] = skyTrack_ShiftLast[skyTrack_idx] + skyTrack_newShift / 2;

	//save current total shift for next time
	skyTrack_ShiftLast[skyTrack_idx] += skyTrack_newShift;

	return skyTrack_Shift[skyTrack_idx];
}

DBL4 SkyrmionTrack::Get_skypos_diametersCUDA(size_t size, DBL3 h, Rect M_rect, cu_obj<cuVEC_VC<cuReal3>>& M, Rect skyRect)
{
	//must have a set rectangle
	if (skyRect.IsNull()) return DBL4();

	int skyTrack_idx = Get_skyTrack_index(skyRect);

	CurveFitting fit;

	vector<double> params;

	//skyRect was the original skyrmion rectangle, used to identify it. We now want to work with the updated skyrmion rectangle.
	skyRect = skyTrack_rect[skyTrack_idx];

	//get x axis data then fit to find skyrmion diameter and center position
	auto fit_x_axis = [&](void) -> DBL2 {

		int points_x = (skyRect.e.x - skyRect.s.x) / h.x;

		if (xy_data.size() < points_x) xy_data.resize(points_x);
		if (data_gpu.size() < points_x) data_gpu.resize(points_x);
		if (data_cpu.size() < points_x) data_cpu.resize(points_x);
		
		//Extract profile from M (gpu to gpu)
		M()->extract_profile_component_z(points_x, data_gpu, cuReal3(skyRect.s.x, (skyRect.e.y + skyRect.s.y) / 2, h.z / 2), cuReal3(h.x, 0, 0));

		//Transfer extracted profile from gpu to cpu
		data_gpu.copy_to_cpuvector(data_cpu);

#pragma omp parallel for
		for (int idx = 0; idx < points_x; idx++) {

			xy_data[idx].x = ((double)idx + 0.5) * h.x + skyRect.s.x;
			xy_data[idx].y = data_cpu[idx];
		}

		fit.FitSkyrmion_LMA(xy_data, params, points_x);

		return DBL2(params[0] * 2, params[1]);
	};

	//get y axis data then fit to find skyrmion diameter and center position
	auto fit_y_axis = [&](void) -> DBL2 {

		int points_y = (skyRect.e.y - skyRect.s.y) / h.y;

		if (xy_data.size() < points_y) xy_data.resize(points_y);
		if (data_gpu.size() < points_y) data_gpu.resize(points_y);
		if (data_cpu.size() < points_y) data_cpu.resize(points_y);

		//Extract profile from M (gpu to gpu)
		M()->extract_profile_component_z(points_y, data_gpu, cuReal3((skyRect.e.x + skyRect.s.x) / 2, skyRect.s.y, h.z / 2), cuReal3(0, h.y, 0));

		//Transfer extracted profile from gpu to cpu
		data_gpu.copy_to_cpuvector(data_cpu);

#pragma omp parallel for
		for (int idx = 0; idx < points_y; idx++) {

			xy_data[idx].x = ((double)idx + 0.5) * h.y + skyRect.s.y;
			xy_data[idx].y = data_cpu[idx];
		}

		fit.FitSkyrmion_LMA(xy_data, params, points_y);

		return DBL2(params[0] * 2, params[1]);
	};

	//1. First fitting along x direction - only need skyrmion center x position now so we can center the rect along x (the "radius" returned will not be the actual radius).

	DBL2 dia_pos = fit_x_axis();

	//center rectangle along x
	skyRect += DBL3(dia_pos.j - (skyRect.e.x + skyRect.s.x) / 2, 0.0, 0.0);

	//2. Fit along y direction - this gives us the correct y axis diameter, and also allows us to center the rectangle along y

	dia_pos = fit_y_axis();
	double diameter_y = dia_pos.i;

	//center rectangle along y
	skyRect += DBL3(0.0, dia_pos.j - (skyRect.e.y + skyRect.s.y) / 2, 0.0);

	//3. Finally fit along x direction again to obtain correct x diameter.

	dia_pos = fit_x_axis();
	double diameter_x = dia_pos.i;

	double position_x = (skyRect.s.x + skyRect.e.x) / 2;
	double position_y = (skyRect.s.y + skyRect.e.y) / 2;

	//Update the skyrmion rectangle for next time - center it on the skyrmion with dimensions 50% larger than the diameter.
	//Also make sure to cap the rectangle to the mesh dimensions so we don't attempt to read data outside of M
	double start_x = position_x - diameter_x * 0.75;
	if (start_x < 0.0) start_x = 0.0;

	double start_y = position_y - diameter_y * 0.75;
	if (start_y < 0.0) start_y = 0.0;

	double end_x = position_x + diameter_x * 0.75;
	if (end_x > M_rect.e.x) end_x = M_rect.e.x;

	double end_y = position_y + diameter_y * 0.75;
	if (end_y > M_rect.e.y) end_y = M_rect.e.y;

	//Update the skyrmion rectangle for next time - center it on the skyrmion with dimensions 50% larger than the diameter.
	skyTrack_rect[skyTrack_idx] = Rect(DBL3(start_x, start_y, 0.0), DBL3(end_x, end_y, h.z));

	return DBL4(position_x, position_y, diameter_x, diameter_y);
}

#endif