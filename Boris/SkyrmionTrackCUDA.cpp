#include "stdafx.h"
#include "SkyrmionTrack.h"
#include "SimulationData.h"

//-------------------------------------------------- Get_skyshift : CUDA 

#if COMPILECUDA == 1

DBL2 SkyrmionTrack::Get_skyshiftCUDA(size_t size, DBL3 h, Rect M_rect, cu_obj<cuVEC_VC<cuReal3>>& M, Rect skyRect)
{
	//must have a set rectangle
	if (skyRect.IsNull()) return DBL2();

	int skyTrack_idx = Get_skyTrack_index(skyRect, M_rect);
	if (skyTrack_idx < 0) return DBL2();

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

	int skyTrack_idx = Get_skyTrack_index(skyRect, M_rect);
	if (skyTrack_idx < 0) return DBL4();

	CurveFitting fit;

	std::vector<double> params;

	//skyRect was the original skyrmion rectangle, used to identify it. We now want to work with the updated skyrmion rectangle.
	skyRect = skyTrack_rect[skyTrack_idx];

	//in-plane mesh dimensions
	DBL2 meshDim = DBL2(M_rect.e.x - M_rect.s.x, M_rect.e.y - M_rect.s.y);

	/*

	//OLD METHOD : NOT RESILIENT ENOUGH

	//maximum number of points along any dimension in the skyrmion tracking window; multiply by 2 since the skyrmion window is adjusted to be 2 times the skyrmion diameter so in the extreme it could be two times the mesh size.
	int max_points = maximum((M_rect.e.x - M_rect.s.x) / h.x, (M_rect.e.y - M_rect.s.y) / h.y) * 2;
	if (max_points <= 0) return DBL4();

	//get x axis data then fit to find skyrmion diameter and center position
	auto fit_x_axis = [&](void) -> DBL2 {

		int points_x = (skyRect.e.x - skyRect.s.x) / h.x;
		
		//skyrmion too large for current mesh size so fail the fitting
		if (points_x > max_points) return DBL2(-1);

		//to improve fitting accuracy extended the number of points either side by filling in start / end values.
		int end_points = 100;

		//keep working arrays to maximum possible size (memory use is not significant and this avoids allocating memory often).
		if (xy_data.size() < max_points + 2 * end_points) xy_data.resize(max_points + 2 * end_points);
		if (data_gpu.size() < max_points) data_gpu.resize(max_points);
		if (data_cpu.size() != data_gpu.size()) data_cpu.resize(data_gpu.size());

		//Extract profile from M (gpu to gpu)
		if (data_gpu.size() >= points_x) M()->extract_profile_component_x(points_x, data_gpu, cuReal3(skyRect.s.x, (skyRect.e.y + skyRect.s.y) / 2, h.z / 2), cuReal3(h.x, 0, 0));
		//just in case
		else return DBL2(-1);

		//Transfer extracted profile from gpu to cpu
		data_gpu.copy_to_vector(data_cpu);

		double window_start = h.x / 2 + skyRect.s.x;
		double window_end = ((double)points_x - 0.5) * h.x + skyRect.s.x;

#pragma omp parallel for
		for (int idx = 0; idx < points_x; idx++) {

			double position = window_start + idx * h.x;
			xy_data[idx + end_points].x = position;
			xy_data[idx + end_points].y = data_cpu[idx];
		}

#pragma omp parallel for
		for (int idx = 0; idx < end_points; idx++) {

			//left end
			xy_data[idx].x = window_start - (end_points - idx) * h.x;
			xy_data[idx].y = 0.0;

			//right end
			xy_data[idx + points_x + end_points].x = window_end + (idx + 1) * h.x;
			xy_data[idx + points_x + end_points].y = 0.0;
		}

		fit.FitSkyrmion_Longitudinal_LMA(xy_data, params, points_x + 2 * end_points);

		return DBL2(params[0] * 2, params[1]);
	};

	//get y axis data then fit to find skyrmion diameter and center position
	auto fit_y_axis = [&](void) -> DBL2 {

		int points_y = (skyRect.e.y - skyRect.s.y) / h.y;
		if (points_y > max_points) return DBL2(-1);

		//to improve fitting accuracy extended the number of points either side by filling in start / end values.
		int end_points = 100;

		//keep working arrays to maximum possible size (memory use is not significant and this avoids allocating memory often).
		if (xy_data.size() < max_points + 2 * end_points) xy_data.resize(max_points + 2 * end_points);
		if (data_gpu.size() < max_points) data_gpu.resize(max_points);
		if (data_cpu.size() != data_gpu.size()) data_cpu.resize(data_gpu.size());

		//Extract profile from M (gpu to gpu)
		if (data_gpu.size() >= points_y) M()->extract_profile_component_y(points_y, data_gpu, cuReal3((skyRect.e.x + skyRect.s.x) / 2, skyRect.s.y, h.z / 2), cuReal3(0, h.y, 0));
		else return DBL2(-1);

		//Transfer extracted profile from gpu to cpu
		data_gpu.copy_to_vector(data_cpu);

		double window_start = h.y / 2 + skyRect.s.y;
		double window_end = ((double)points_y - 0.5) * h.y + skyRect.s.y;

#pragma omp parallel for
		for (int idx = 0; idx < points_y; idx++) {

			double position = window_start + idx * h.y;
			xy_data[idx + end_points].x = position;
			xy_data[idx + end_points].y = data_cpu[idx];
		}

#pragma omp parallel for
		for (int idx = 0; idx < end_points; idx++) {

			//left end
			xy_data[idx].x = window_start - (end_points - idx) * h.y;
			xy_data[idx].y = 0.0;

			//right end
			xy_data[idx + points_y + end_points].x = window_end + (idx + 1) * h.y;
			xy_data[idx + points_y + end_points].y = 0.0;
		}

		fit.FitSkyrmion_Longitudinal_LMA(xy_data, params, points_y + 2 * end_points);

		return DBL2(params[0] * 2, params[1]);
	};

	*/

	//maximum number of points along any dimension in the skyrmion tracking window; multiply by 2 since the skyrmion window is adjusted to be 2 times the skyrmion diameter so in the extreme it could be two times the mesh size.
	int max_points = maximum((M_rect.e.x - M_rect.s.x) / h.x, (M_rect.e.y - M_rect.s.y) / h.y) * 2;
	if (max_points <= 0) return DBL4();

	auto fit_x_axis_zerocrossing = [&](void) -> DBL2 {

		auto search_line = [&](double pos_y) -> DBL2 {

			int points_x = (skyRect.e.x - skyRect.s.x) / h.x;

			//skyrmion too large for current mesh size so fail the fitting
			if (points_x > max_points) return DBL2(-1);

			//keep working arrays to maximum possible size (memory use is not significant and this avoids allocating memory often).
			if (data_gpu.size() < max_points) data_gpu.resize(max_points);
			if (data_cpu.size() != data_gpu.size()) data_cpu.resize(data_gpu.size());

			//Extract profile from M (gpu to gpu)
			double position = skyRect.s.x;
			position -= floor_epsilon(position / meshDim.x) * meshDim.x;

			if (data_gpu.size() >= points_x) M()->extract_profile_component_z(points_x, data_gpu, cuReal3(position, pos_y, h.z / 2), cuReal3(h.x, 0, 0));
			//just in case
			else return DBL2(-1);

			//Transfer extracted profile from gpu to cpu
			data_gpu.copy_to_vector(data_cpu);

			bool plus_sign = data_cpu[0] > 0;
			double first_crossing = 0.0, second_crossing = 0.0;

			for (int idx = 0; idx < (skyRect.e.x - skyRect.s.x) / h.x; idx++) {

				double value = data_cpu[idx];

				if ((plus_sign && value < 0) || (!plus_sign && value > 0)) {

					plus_sign = !plus_sign;

					if (!first_crossing) {

						first_crossing = skyRect.s.x + idx * h.x;
					}
					else {

						second_crossing = skyRect.s.x + idx * h.x;
						break;
					}
				}
			}

			if (first_crossing && second_crossing) return DBL2(second_crossing - first_crossing, (first_crossing + second_crossing) / 2);
			else return DBL2(-1);
		};

		//initially search through the center of the tracker rectangle
		double pos_y = (skyRect.e.y + skyRect.s.y) / 2;
		pos_y -= floor_epsilon(pos_y / meshDim.y) * meshDim.y;

		DBL2 dia_pos = search_line(pos_y);

		if (dia_pos.i > 0) return dia_pos;
		else {

			DBL2 max_dia_pos = DBL2(-1);

			//bounds couldn't be found, so search line by line for the largest bounds distance
			for (int idx_y = 0; idx_y < (skyRect.e.y - skyRect.s.y) / h.y; idx_y++) {

				double pos_y = skyRect.s.y + idx_y * h.y;
				pos_y -= floor_epsilon(pos_y / meshDim.y) * meshDim.y;
				dia_pos = search_line(pos_y);
				if (dia_pos >= 0) {

					if (dia_pos.i > max_dia_pos.i) max_dia_pos = dia_pos;
				}
			}

			if (max_dia_pos >= 0) return max_dia_pos;

			//searched everything and still couldn't find bounds : no skyrmion present in current rectangle, or skyrmion too small for current cellsize
			//finally try to set the tracker rectangle to the entire mesh and search again

			//set tracker rectangle to entire mesh, remembering we need a relative rect
			skyRect = M_rect - M_rect.s;

			for (int idx_y = 0; idx_y < (skyRect.e.y - skyRect.s.y) / h.y; idx_y++) {

				dia_pos = search_line(skyRect.s.y + idx_y * h.y);
				if (dia_pos >= 0) {

					if (dia_pos.i > max_dia_pos.i) max_dia_pos = dia_pos;
				}
			}

			return max_dia_pos;
		}
	};

	auto fit_y_axis_zerocrossing = [&](void) -> DBL2 {

		auto search_line = [&](double pos_x) -> DBL2 {

			int points_y = (skyRect.e.y - skyRect.s.y) / h.y;

			//skyrmion too large for current mesh size so fail the fitting
			if (points_y > max_points) return DBL2(-1);

			//keep working arrays to maximum possible size (memory use is not significant and this avoids allocating memory often).
			if (data_gpu.size() < max_points) data_gpu.resize(max_points);
			if (data_cpu.size() != data_gpu.size()) data_cpu.resize(data_gpu.size());

			//Extract profile from M (gpu to gpu)
			double position = skyRect.s.y;
			position -= floor_epsilon(position / meshDim.y) * meshDim.y;

			if (data_gpu.size() >= points_y) M()->extract_profile_component_z(points_y, data_gpu, cuReal3(pos_x, position, h.z / 2), cuReal3(0, h.y, 0));
			else return DBL2(-1);

			//Transfer extracted profile from gpu to cpu
			data_gpu.copy_to_vector(data_cpu);

			bool plus_sign = data_cpu[0] > 0;
			double first_crossing = 0.0, second_crossing = 0.0;

			for (int idx = 0; idx < (skyRect.e.y - skyRect.s.y) / h.y; idx++) {

				double value = data_cpu[idx];

				if ((plus_sign && value < 0) || (!plus_sign && value > 0)) {

					plus_sign = !plus_sign;

					if (!first_crossing) {

						first_crossing = skyRect.s.y + idx * h.y;
					}
					else {

						second_crossing = skyRect.s.y + idx * h.y;
						break;
					}
				}
			}

			if (first_crossing && second_crossing) return DBL2(second_crossing - first_crossing, (first_crossing + second_crossing) / 2);
			else return DBL2(-1);
		};

		//initially search through the center of the tracker rectangle
		double pos_x = (skyRect.e.x + skyRect.s.x) / 2;
		//wrap around if needed
		pos_x -= floor_epsilon(pos_x / meshDim.x) * meshDim.x;

		DBL2 dia_pos = search_line(pos_x);

		if (dia_pos.i > 0) return dia_pos;
		else {

			DBL2 max_dia_pos = DBL2(-1);

			//bounds couldn't be found, so search line by line for the largest bounds distance
			for (int idx_x = 0; idx_x < (skyRect.e.x - skyRect.s.x) / h.x; idx_x++) {

				double pos_x = skyRect.s.x + idx_x * h.x;
				//wrap around if needed
				pos_x -= floor_epsilon(pos_x / meshDim.x) * meshDim.x;
				dia_pos = search_line(pos_x);
				if (dia_pos >= 0) {

					if (dia_pos.i > max_dia_pos.i) max_dia_pos = dia_pos;
				}
			}

			if (max_dia_pos >= 0) return max_dia_pos;

			//searched everything and still couldn't find bounds : no skyrmion present in current rectangle, or skyrmion too small for current cellsize
			//finally try to set the tracker rectangle to the entire mesh and search again

			//set tracker rectangle to entire mesh, remembering we need a relative rect
			skyRect = M_rect - M_rect.s;

			for (int idx_x = 0; idx_x < (skyRect.e.x - skyRect.s.x) / h.x; idx_x++) {

				dia_pos = search_line(skyRect.s.x + idx_x * h.x);
				if (dia_pos >= 0) {

					if (dia_pos.i > max_dia_pos.i) max_dia_pos = dia_pos;
				}
			}

			return max_dia_pos;
		}
	};

	//1. Fitting along x direction

	DBL2 dia_pos = fit_x_axis_zerocrossing();

	//need these checks just in case the fitting fails
	if (dia_pos.i < 0) return DBL4();

	double diameter_x = dia_pos.i;
	double position_x = dia_pos.j;

	//center rectangle along x
	skyRect += DBL3(position_x - (skyRect.e.x + skyRect.s.x) / 2, 0.0, 0.0);

	//2. Fit along y direction - this gives us the correct y axis diameter and y center position, and also allows us to center the rectangle along y

	dia_pos = fit_y_axis_zerocrossing();

	//need these checks just in case the fitting fails
	if (dia_pos.i < 0) return DBL4();

	double diameter_y = dia_pos.i;
	double position_y = dia_pos.j;

	//center rectangle along y
	skyRect += DBL3(0.0, position_y - (skyRect.e.y + skyRect.s.y) / 2, 0.0);
	
	//3. Fitting along x direction again

	dia_pos = fit_x_axis_zerocrossing();

	//need these checks just in case the fitting fails
	if (dia_pos.i < 0) return DBL4();

	diameter_x = dia_pos.i;
	position_x = dia_pos.j;

	//center rectangle along x
	skyRect += DBL3(position_x - (skyRect.e.x + skyRect.s.x) / 2, 0.0, 0.0);

	//Update the skyrmion rectangle for next time - center it on the skyrmion with dimensions dia_mul times larger than the diameter.
	double start_x = position_x - diameter_x * dia_mul / 2;
	double start_y = position_y - diameter_y * dia_mul / 2;
	double end_x = position_x + diameter_x * dia_mul / 2;
	double end_y = position_y + diameter_y * dia_mul / 2;

	//Update the skyrmion rectangle for next time - center it on the skyrmion with dimensions 2 times larger than the diameter.
	skyTrack_rect[skyTrack_idx] = Rect(DBL3(start_x, start_y, 0.0), DBL3(end_x, end_y, h.z));

	return DBL4(position_x, position_y, diameter_x, diameter_y);
}

#endif