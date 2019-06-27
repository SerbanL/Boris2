#pragma once

#include "Boris_Enums_Defs.h"

#include "BorisLib.h"
#include "ErrorHandler.h"

#if COMPILECUDA == 1
#include "BorisCUDALib.h"
#endif

class DatumConfig;

//Class used to track one or more skyrmions in a given ferromagnetic mesh
class SkyrmionTrack :
	public ProgramState<SkyrmionTrack,
	tuple<std::vector<DBL2>, std::vector<DBL2>, std::vector<DBL2>>,
	tuple<> >
{

private:

	//identify skyrmion tracker by bottom-left Rect coordinates
	vector<DBL2> skyTrack_Id;

	//the skyrmion total shift - same vector size as skyTrack_Id
	vector<DBL2> skyTrack_Shift;

	//the skyrmion total shift from last pass; use it to eliminate shift oscillations - same vector size as skyTrack_Id
	vector<DBL2> skyTrack_ShiftLast;

private:

	//for given skyrmion identifying rectangle obtain an index in skyTrack_Shift - either an existing entry if skyRectOrigin found in skyTrack_Id, else make a new entry
	int Get_skyTrack_index(DBL2 skyRectOrigin);

public:

	SkyrmionTrack(void) :
		ProgramStateNames(this,
			{ VINFO(skyTrack_Id), VINFO(skyTrack_Shift), VINFO(skyTrack_ShiftLast) }, {})
	{}

	~SkyrmionTrack() {}

	//-----------------------------------

	//implement pure virtual method from ProgramState
	void RepairObjectState(void) 
	{
		if (skyTrack_Shift.size() != skyTrack_Id.size()) skyTrack_Shift.resize(skyTrack_Id.size());
		if (skyTrack_ShiftLast.size() != skyTrack_Id.size()) skyTrack_ShiftLast.resize(skyTrack_Id.size());
	}

	//----------------------------------- PUBLIC METHODS

	//Calculate skyrmion shift for a skyrmion in M initially in the given rectangle (relative coordinates)
	//The bottom-left coordinates of the rectangle are used to uniquely identify the skyrmion (new entry made if this is first call for this Rect), and a shift is calculated
	//so that the shifted rectangle is centered over the skyrmion - the total shift for this entry is returned. This method should be called often enough so tracking is not lost.
	//This method only creates entries, does not clean any.
	DBL2 Get_skyshift(VEC_VC<DBL3>& M, Rect skyRect);

#if COMPILECUDA == 1
	//as for the non-CUDA version but also pass in M.n.dim() and M.h so we don't have to get them from gpu memory
	DBL2 Get_skyshiftCUDA(size_t size, DBL3 h, cu_obj<cuVEC_VC<cuReal3>>& M, Rect skyRect);
#endif

	//Clean any skyrmion tracker entries by comparing the entries held here against external lists (Simulation::dataBoxList and Simulation::saveDataList) - if not found in external lists then delete them here
	//This method is called by FMesh::UpdateConfiguration
	void UpdateConfiguration(vector_lut<DatumConfig>& saveDataList);
};

