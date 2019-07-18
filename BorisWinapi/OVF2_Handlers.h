#pragma once

#include <string>

#include "ErrorHandler.h"

#include "BorisLib.h"

class OVF2 {

	enum header_id { MESHTYPE, MESHUNIT, VALUEDIM, XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX, XNODES, YNODES, ZNODES, XSTEP, YSTEP, ZSTEP, BEGIN_DATA_UC, BEGIN_DATA_LC };

	vector_lut<string> headers;

public:

	OVF2(void);

	//read an OOMMF OVF2 file containing uniform vector data, and set the data VEC from it
	BError Read_OVF2_VEC(std::string fileName, VEC<DBL3>& data);
};