#include "stdafx.h"
#include "SkyrmionTrack.h"
#include "SimulationData.h"

//-------------------------------------------------- Get_skyshift : non CUDA 

//Calculate skyrmion shift for a skyrmion in M initially in the given rectangle (relative coordinates)
//The bottom-left coordinates of the rectangle are used to uniquely identify the skyrmion (new entry made if this is first call for this Rect), and a shift is calculated
//so that the shifted rectangle is centered over the skyrmion - the total shift for this entry is returned. This method should be called often enough so tracking is not lost.
//This method only creates entries, does not clean any.
DBL2 SkyrmionTrack::Get_skyshift(VEC_VC<DBL3>& M, Rect skyRect)
{
	//must have a set rectangle
	if (skyRect.IsNull()) return DBL2();

	int skyTrack_idx = Get_skyTrack_index(skyRect);

	//shift the skyrmion rectangle to current tracking position
	skyRect += DBL3(skyTrack_ShiftLast[skyTrack_idx].x, skyTrack_ShiftLast[skyTrack_idx].y, 0.0);

	//Find M averages in the 4 skyRect xy-plane quadrants
	DBL3 bottom_left = M.average_nonempty_omp(skyRect.get_quadrant_bl());
	DBL3 bottom_right = M.average_nonempty_omp(skyRect.get_quadrant_br());
	DBL3 top_left = M.average_nonempty_omp(skyRect.get_quadrant_tl());
	DBL3 top_right = M.average_nonempty_omp(skyRect.get_quadrant_tr());

	//the new shift value
	DBL2 skyTrack_newShift = DBL2();

	//if left half z-component value modulus has increased compared to the right half then it contains less of the skyrmion ring -> skyrmion must have shifted along +x so follow it
	if (mod(bottom_left.z + top_left.z) > mod(bottom_right.z + top_right.z)) {

		skyTrack_newShift.x += M.h.x;
	}
	else {

		skyTrack_newShift.x -= M.h.x;
	}

	//if bottom half z-component value modulus has increased compared to the top half then it contains less of the skyrmion ring -> skyrmion must have shifted along +y so follow it
	if (mod(bottom_left.z + bottom_right.z) > mod(top_left.z + top_right.z)) {

		skyTrack_newShift.y += M.h.y;
	}
	else {

		skyTrack_newShift.y -= M.h.y;
	}

	//set actual total shift as average of new and last total shift values - this eliminates tracking oscillations
	skyTrack_Shift[skyTrack_idx] = skyTrack_ShiftLast[skyTrack_idx] + skyTrack_newShift / 2;
	
	//save current total shift for next time
	skyTrack_ShiftLast[skyTrack_idx] += skyTrack_newShift;

	return skyTrack_Shift[skyTrack_idx];
}

//additionally return the x and y diameters
DBL4 SkyrmionTrack::Get_skypos_diameters(VEC_VC<DBL3>& M, Rect skyRect)
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

		int points_x = (skyRect.e.x - skyRect.s.x) / M.h.x;
		
		//to improve fitting accuracy extended the number of points either side by filling in start / end values.
		int end_points = 100;

		if (xy_data.size() < points_x + 2 * end_points) xy_data.resize(points_x + 2 * end_points);

		double window_start = M.h.x / 2 + skyRect.s.x;
		double window_end = ((double)points_x - 0.5) * M.h.x + skyRect.s.x;

#pragma omp parallel for
		for (int idx = 0; idx < points_x; idx++) {

			double position = window_start + idx * M.h.x;
			xy_data[idx + end_points].x = position;

			DBL3 value = M.weighted_average(DBL3(position, (skyRect.e.y + skyRect.s.y) / 2, M.h.z / 2), M.h);
			xy_data[idx + end_points].y = value.x;
		}

#pragma omp parallel for
		for (int idx = 0; idx < end_points; idx++) {

			//left end
			xy_data[idx].x = window_start - (end_points - idx) * M.h.x;
			xy_data[idx].y = 0.0;

			//right end
			xy_data[idx + points_x + end_points].x = window_end + (idx + 1) * M.h.x;
			xy_data[idx + points_x + end_points].y = 0.0;
		}

		fit.FitSkyrmion_Longitudinal_LMA(xy_data, params, points_x + 2 * end_points);

		//return diameter and center x
		return DBL2(params[0] * 2, params[1]);
	};

	//get y axis data then fit to find skyrmion diameter and center position
	auto fit_y_axis = [&](void) -> DBL2 {

		int points_y = (skyRect.e.y - skyRect.s.y) / M.h.y;
		
		//to improve fitting accuracy extended the number of points either side by filling in start / end values.
		int end_points = 100;

		if (xy_data.size() < points_y + 2 * end_points) xy_data.resize(points_y + 2 * end_points);

		double window_start = M.h.y / 2 + skyRect.s.y;
		double window_end = ((double)points_y - 0.5) * M.h.y + skyRect.s.y;

#pragma omp parallel for
		for (int idx = 0; idx < points_y; idx++) {

			double position = window_start + idx * M.h.y;
			xy_data[idx + end_points].x = position;

			DBL3 value = M.weighted_average(DBL3((skyRect.e.x + skyRect.s.x) / 2, position, M.h.z / 2), M.h);
			xy_data[idx + end_points].y = value.y;
		}

#pragma omp parallel for
		for (int idx = 0; idx < end_points; idx++) {

			//left end
			xy_data[idx].x = window_start - (end_points - idx) * M.h.y;
			xy_data[idx].y = 0.0;

			//right end
			xy_data[idx + points_y + end_points].x = window_end + (idx + 1) * M.h.y;
			xy_data[idx + points_y + end_points].y = 0.0;
		}

		//fit.FitSkyrmion_LMA(xy_data, params, points_y + 2 * end_points);
		fit.FitSkyrmion_Longitudinal_LMA(xy_data, params, points_y + 2 * end_points);

		//return diameter and center y
		return DBL2(params[0] * 2, params[1]);
	};

	//1. Fitting along x direction

	DBL2 dia_pos = fit_x_axis();
	
	//need these checks just in case the fitting fails : otherwise we can crash the program
	if (isnan(dia_pos.i) || isnan(dia_pos.j) || dia_pos.i < 0 || dia_pos.j < 0) return DBL4();

	double diameter_x = dia_pos.i;
	double position_x = dia_pos.j;

	//center rectangle along x
	skyRect += DBL3(dia_pos.j - (skyRect.e.x + skyRect.s.x) / 2, 0.0, 0.0);

	//make sure the rectangle doesn't go out of bounds
	if (skyRect.s.x < 0.0) skyRect.s.x = 0.0;
	if (skyRect.e.x > M.rect.e.x) skyRect.e.x = M.rect.e.x;

	//2. Fit along y direction - this gives us the correct y axis diameter and y center position, and also allows us to center the rectangle along y

	dia_pos = fit_y_axis();

	//need these checks just in case the fitting fails : otherwise we can crash the program
	if (isnan(dia_pos.i) || isnan(dia_pos.j) || dia_pos.i < 0 || dia_pos.j < 0) return DBL4();

	double diameter_y = dia_pos.i;
	double position_y = dia_pos.j;

	//center rectangle along y
	skyRect += DBL3(0.0, dia_pos.j - (skyRect.e.y + skyRect.s.y) / 2, 0.0);

	//make sure the rectangle doesn't go out of bounds
	if (skyRect.s.y < 0.0) skyRect.s.y = 0.0;
	if (skyRect.e.y > M.rect.e.y) skyRect.e.y = M.rect.e.y;

	//Update the skyrmion rectangle for next time - center it on the skyrmion with dimensions 2 times larger than the diameter.
	//Also make sure to cap the rectangle to the mesh dimensions so we don't attempt to read data outside of M
	double start_x = position_x - diameter_x;
	if (start_x < 0.0) start_x = 0.0;

	double start_y = position_y - diameter_y;
	if (start_y < 0.0) start_y = 0.0;

	double end_x = position_x + diameter_x;
	if (end_x > M.rect.e.x) end_x = M.rect.e.x;

	double end_y = position_y + diameter_y;
	if (end_y > M.rect.e.y) end_y = M.rect.e.y;

	//Update the skyrmion rectangle for next time - center it on the skyrmion with dimensions 2 times larger than the diameter.
	skyTrack_rect[skyTrack_idx] = Rect(DBL3(start_x, start_y, 0.0), DBL3(end_x, end_y, M.h.z));

	return DBL4(position_x, position_y, diameter_x, diameter_y);
}

//-------------------------------------------------- AUXILIARY

//for given skyrmion identifying rectangle obtain an index in skyTrack_Shift - either an existing entry if skyRectOrigin found in skyTrack_Id, else make a new entry
int SkyrmionTrack::Get_skyTrack_index(Rect skyRect)
{
	DBL2 skyRectOrigin = DBL2(skyRect.s.x, skyRect.s.y);

	for (int idx = 0; idx < skyTrack_Id.size(); idx++) {

		if (skyRectOrigin == skyTrack_Id[idx]) return idx;
	}

	//nothing found so create new entry
	skyTrack_Id.push_back(skyRectOrigin);
	skyTrack_Shift.push_back(DBL2());
	skyTrack_ShiftLast.push_back(DBL2());
	skyTrack_rect.push_back(skyRect);

	return skyTrack_Id.size() - 1;
}

//Clean any skyrmion tracker entries by comparing the entries held here against external lists (Simulation::dataBoxList and Simulation::saveDataList) - if not found in external lists then delete them here
//This method is called by (A)FMesh::UpdateConfiguration
void SkyrmionTrack::UpdateConfiguration(vector_lut<DatumConfig>& saveDataList)
{
	//check saveDataList to see if there is a DATA_SKYSHIFT or DATA_SKYPOS entry with matching rectangle xy origin.
	auto FindEntry = [&](DBL2 skyRectOrigin) -> bool {

		for (int idx = 0; idx < saveDataList.size(); idx++) {

			if (saveDataList[idx].datumId == DATA_SKYSHIFT || saveDataList[idx].datumId == DATA_SKYPOS) {

				//found it
				if (DBL2(saveDataList[idx].rectangle.s.x, saveDataList[idx].rectangle.s.y) == skyRectOrigin) return true;
			}
		}

		//couldn't find any
		return false;
	};

	int entry_idx = 0;

	while (entry_idx < skyTrack_Id.size()) {

		if (!FindEntry(skyTrack_Id[entry_idx])) {

			//didn't find entry so delete it
			skyTrack_Id.erase(skyTrack_Id.begin() + entry_idx);
			skyTrack_Shift.erase(skyTrack_Shift.begin() + entry_idx);
			skyTrack_ShiftLast.erase(skyTrack_ShiftLast.begin() + entry_idx);
			skyTrack_rect.erase(skyTrack_rect.begin() + entry_idx);

			//check next entry but don't increment entry_idx as we just erased entry at entry_idx
		}
		//found entry, check next one
		else entry_idx++;
	}
}