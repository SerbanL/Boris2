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

	//read an OOMMF OVF2 file containing uniform scalar data, and set the data VEC from it
	template <typename VECType>
	BError Read_OVF2_SCA(std::string fileName, VECType& data);

	//read an OOMMF OVF2 file containing uniform vector data, and set the data VEC from it
	template <typename VECType>
	BError Read_OVF2_VEC(std::string fileName, VECType& data);

	//read on OVF2 file into a VEC, which can store either a scalar or vector type : the ovf file stored data must match
	template <typename VECType>
	BError Read_OVF2(std::string fileName, VECType& data)
	{
		typedef typename std::remove_reference<VECType>::type VType;
		typedef typename contained_type<VType>::type SType;

		if (std::is_fundamental<SType>::value) {

			//scalar
			return Read_OVF2_SCA(fileName, data);
		}
		else {

			//vector
			return Read_OVF2_VEC(fileName, data);
		}
	}

	//write an OOMMF OVF2 file containing uniform scalar data
	//you can also choose the type of data output : data_type = bin4 for single precision binary, data_type = bin8 for double precision binary, or data_type = text
	template <typename VECType>
	BError Write_OVF2_SCA(std::string fileName, VECType& data, std::string data_type = "bin8");

	//write an OOMMF OVF2 file containing uniform vector data
	//you can write normalized data to norm (divide by it, no normalization by default)
	//you can also choose the type of data output : data_type = bin4 for single precision binary, data_type = bin8 for double precision binary, or data_type = text
	template <typename VECType>
	BError Write_OVF2_VEC(std::string fileName, VECType& data, std::string data_type = "bin8", double norm = 1.0);
};