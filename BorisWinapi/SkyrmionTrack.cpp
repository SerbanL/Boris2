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

	int skyTrack_idx = Get_skyTrack_index(DBL2(skyRect.s.x, skyRect.s.y));

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

//-------------------------------------------------- Get_skyshift : CUDA 

#if COMPILECUDA == 1
DBL2 SkyrmionTrack::Get_skyshiftCUDA(size_t size, DBL3 h, cu_obj<cuVEC_VC<cuReal3>>& M, Rect skyRect)
{
	//must have a set rectangle
	if (skyRect.IsNull()) return DBL2();

	int skyTrack_idx = Get_skyTrack_index(DBL2(skyRect.s.x, skyRect.s.y));

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
#endif

//-------------------------------------------------- AUXILIARY

//for given skyrmion identifying rectangle obtain an index in skyTrack_Shift - either an existing entry if skyRectOrigin found in skyTrack_Id, else make a new entry
int SkyrmionTrack::Get_skyTrack_index(DBL2 skyRectOrigin)
{
	for (int idx = 0; idx < skyTrack_Id.size(); idx++) {

		if (skyRectOrigin == skyTrack_Id[idx]) return idx;
	}

	//nothing found so create new entry
	skyTrack_Id.push_back(skyRectOrigin);
	skyTrack_Shift.push_back(DBL2());
	skyTrack_ShiftLast.push_back(DBL2());

	return skyTrack_Id.size() - 1;
}

//Clean any skyrmion tracker entries by comparing the entries held here against external lists (Simulation::dataBoxList and Simulation::saveDataList) - if not found in external lists then delete them here
//This method is called by FMesh::UpdateConfiguration
void SkyrmionTrack::UpdateConfiguration(vector_lut<DatumConfig>& saveDataList)
{
	//check saveDataList to see if there is a DATA_SKYSHIFT entry with matching rectangle xy origin.
	auto FindEntry = [&](DBL2 skyRectOrigin) -> bool {

		for (int idx = 0; idx < saveDataList.size(); idx++) {

			if (saveDataList[idx].datumId == DATA_SKYSHIFT) {

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

			//check next entry but don't increment entry_idx as we just erased entry at entry_idx
		}
		//found entry, check next one
		else entry_idx++;
	}
}