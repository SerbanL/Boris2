#pragma once

#include <string>

#include "ErrorHandler.h"

#include "BorisLib.h"

class OVF2 {

	enum header_id { OVF2_HEADER, MESHTYPE, MESHUNIT, VALUEDIM, XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX, XNODES, YNODES, ZNODES, XSTEP, YSTEP, ZSTEP, BEGIN_DATA, END_DATA, BEGIN_SEGMENT, END_SEGMENT, BEGIN_HEADER, END_HEADER	};

	vector_lut<std::string> headers;

	enum data_type { DATA_BINARY4, DATA_BINARY8, DATA_TEXT };

	vector_lut<std::string> data_headers;

public:

	OVF2(void);

	//read an OOMMF OVF2 file containing uniform vector data, and set the data VEC from it
	BError Read_OVF2_VEC(std::string fileName, VEC<DBL3>& data);

	//read an OOMMF OVF2 file containing uniform scalar data, and set the data VEC from it
	BError Read_OVF2_VEC(std::string fileName, VEC<double>& data);

	//write an OOMMF OVF2 file containing uniform vector data
	//you can write normalized data to norm (divide by it, no normalization by default)
	//you can also choose the type of data output : data_type = bin4 for single precision binary, data_type = bin8 for double precision binary, or data_type = text
	BError Write_OVF2_VEC(std::string fileName, VEC<DBL3>& data, double norm = 1.0, std::string data_type = "bin8");

	//write an OOMMF OVF2 file containing uniform scalar data
	//you can also choose the type of data output : data_type = bin4 for single precision binary, data_type = bin8 for double precision binary, or data_type = text
	BError Write_OVF2_VEC(std::string fileName, VEC<double>& data, std::string data_type = "bin8");
};