#pragma once

#include "Boris_Enums_Defs.h"

#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "DWRunTimeFitDefs.h"

//Class used to analyse domain wall and extract position and width - efficient and robust for runtime, especially atomistic simulations
class DWPosWidthCUDA
{

private:

	//auxiliary data
	cu_obj<cuBReal> As, Ae, x0, dw, weight;
	cu_obj<size_t> av_points_x0;

	//auxiliary vector containing smoothed xy data
	cu_arr<cuReal2> xy_data_smoothed;
	size_t xy_data_smoothed_size = 0;

	//profile data
	cu_arr<cuReal2> xy_data;
	size_t xy_data_size = 0;

public:

	DWPosWidthCUDA(void) {}
	~DWPosWidthCUDA() {}

	//----------------------------------- PUBLIC METHODS

	//Fit the extracted profile (xy_data) for position and width, assuming the magnetization component follows f(x) = [ (As - Ae) * tanh(-PI * (x - x0) / dw) + (As + Ae) ] / 2
	//Here As, Ae are the start and end values - profile must be long enough to include at least DWPOS_ENDSTENCIL (length ratio) flat parts of the tanh profile
	//x0 is the centre position relative to start of profile
	//dw is the domain wall width
	//xy_data contains profile as x coordinates and corresponding y values
	//return x0, dw
	cuReal2 FitDomainWallCUDA(double length);

	//prepare xy_data profile storage with correct dimensions and return through reference so it can be filled in
	cu_arr<cuReal2>& get_xy_data_ref(size_t size)
	{
		if (xy_data_size != size) {

			xy_data.resize(size);
			xy_data_size = size;
		}

		return xy_data;
	}
};

#endif

