#include "stdafx.h"
#include "OVF2_Handlers.h"

OVF2::OVF2(void)
{
	headers.push_back("# OOMMF OVF 2.0", OVF2_HEADER);
	headers.push_back("# meshtype: ", MESHTYPE);
	headers.push_back("# meshunit: ", MESHUNIT);
	headers.push_back("# valuedim: ", VALUEDIM);
	headers.push_back("# xmin: ", XMIN);
	headers.push_back("# ymin: ", YMIN);
	headers.push_back("# zmin: ", ZMIN);
	headers.push_back("# xmax: ", XMAX);
	headers.push_back("# ymax: ", YMAX);
	headers.push_back("# zmax: ", ZMAX);
	headers.push_back("# xnodes: ", XNODES);
	headers.push_back("# ynodes: ", YNODES);
	headers.push_back("# znodes: ", ZNODES);
	headers.push_back("# xstepsize: ", XSTEP);
	headers.push_back("# ystepsize: ", YSTEP);
	headers.push_back("# zstepsize: ", ZSTEP);
	headers.push_back("# Begin: data ", BEGIN_DATA);
	headers.push_back("# End: data ", END_DATA);
	headers.push_back("# Begin: Segment", BEGIN_SEGMENT);
	headers.push_back("# End: Segment", END_SEGMENT);
	headers.push_back("# Begin: Header", BEGIN_HEADER);
	headers.push_back("# End: Header", END_HEADER);

	data_headers.push_back("binary 4", DATA_BINARY4);
	data_headers.push_back("binary 8", DATA_BINARY8);
	data_headers.push_back("text", DATA_TEXT);
}

template BError OVF2::Read_OVF2_SCA(std::string fileName, VEC<float>& data);
template BError OVF2::Read_OVF2_SCA(std::string fileName, VEC<double>& data);
template BError OVF2::Read_OVF2_SCA(std::string fileName, VEC_VC<float>& data);
template BError OVF2::Read_OVF2_SCA(std::string fileName, VEC_VC<double>& data);

//read an OOMMF OVF2 file containing uniform scalar data, and set the data VEC from it
template <typename VECType>
BError OVF2::Read_OVF2_SCA(std::string fileName, VECType& data)
{
	BError error(__FUNCTION__);

	std::ifstream bdin;
	bdin.open(fileName.c_str(), std::ios::in | std::ios::binary);

	//need to find the following lines:
	//
	//"# meshtype: " -> this must be followed by "rectangular"; not reading any other types currently
	//"# meshunit: " -> this must be followed by "m" or "nm"
	//"# valuedim: " -> this must be followed by "1" since this is a scalar type

	//"# xmin: " -> value to follow
	//"# ymin: " -> value to follow
	//"# zmin: " -> value to follow
	//"# xmax: " -> value to follow
	//"# ymax: " -> value to follow
	//"# zmax: " -> value to follow
	//"# xnodes: " -> value to follow
	//"# ynodes: " -> value to follow
	//"# znodes : " -> value to follow
	//"# xstepsize: " -> value to follow
	//"# ystepsize: " -> value to follow
	//"# zstepsize : " -> value to follow

	//"# Begin: Data " or "# Begin: data " -> type of data to follow; Accepting "binary 4", "binary 8", "text". When this is found then stop searching for any other lines, and if conditions met start getting data.

	bool meshtype_rectangular = false;
	double meshunit = 1.0;
	int valuedim = 0;
	int data_bytes = 0;

	Rect meshRect;
	INT3 n;
	DBL3 h;

	//return 0 if something went wrong : abort loading OVF file
	//return 1 if everything fine, continue
	//return 2 if start of data header found
	auto scan_line = [&](std::string line) -> int {

		if (line.find(headers(MESHTYPE)) != std::string::npos) {

			std::string value = line.substr(headers(MESHTYPE).length());
			if (value != "rectangular") return 0;
			meshtype_rectangular = true;
		}

		else if (line.find(headers(MESHUNIT)) != std::string::npos) {

			std::string value = line.substr(headers(MESHUNIT).length());
			if (value == "m") meshunit = 1.0;
			else if (value == "nm") meshunit = 1e-9;
			else return 0;
		}

		else if (line.find(headers(VALUEDIM)) != std::string::npos) {

			std::string value = line.substr(headers(VALUEDIM).length());
			if (value == "1") valuedim = 1;
			else return 0;
		}

		else if (line.find(headers(XMIN)) != std::string::npos) {

			std::string value = line.substr(headers(XMIN).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.s.x = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(YMIN)) != std::string::npos) {

			std::string value = line.substr(headers(YMIN).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.s.y = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(ZMIN)) != std::string::npos) {

			std::string value = line.substr(headers(ZMIN).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.s.z = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(XMAX)) != std::string::npos) {

			std::string value = line.substr(headers(XMAX).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.e.x = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(YMAX)) != std::string::npos) {

			std::string value = line.substr(headers(YMAX).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.e.y = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(ZMAX)) != std::string::npos) {

			std::string value = line.substr(headers(ZMAX).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.e.z = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(XNODES)) != std::string::npos) {

			std::string value = line.substr(headers(XNODES).length());

			n.x = ToNum(value);
		}

		else if (line.find(headers(YNODES)) != std::string::npos) {

			std::string value = line.substr(headers(YNODES).length());

			n.y = ToNum(value);
		}

		else if (line.find(headers(ZNODES)) != std::string::npos) {

			std::string value = line.substr(headers(ZNODES).length());

			n.z = ToNum(value);
		}

		else if (line.find(headers(XSTEP)) != std::string::npos) {

			std::string value = line.substr(headers(XSTEP).length());

			h.x = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(YSTEP)) != std::string::npos) {

			std::string value = line.substr(headers(YSTEP).length());

			h.y = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(ZSTEP)) != std::string::npos) {

			std::string value = line.substr(headers(ZSTEP).length());

			h.z = (double)ToNum(value) * meshunit;
		}

		else if (lowercase(line).find(lowercase(headers(BEGIN_DATA))) != std::string::npos) {

			//convert to lowercase for comparison as some implementations may choose to have "Data" instead of "data"
			std::string value = lowercase(line.substr(headers(BEGIN_DATA).length()));

			//data headers already in lowercase but better use the lowercase conversion anyway to reduce possibility of bugs in the future
			if (value == lowercase(data_headers(DATA_BINARY4))) {

				data_bytes = 4;
				return 2;
			}
			else if (value == lowercase(data_headers(DATA_BINARY8))) {

				data_bytes = 8;
				return 2;
			}
			else if (value == lowercase(data_headers(DATA_TEXT))) {

				//value of 1 means it's a text input
				data_bytes = 1;
				return 2;
			}
		}

		return 1;
	};

	//read a text line from input stream
	auto read_line = [&](std::ifstream& bdin, char* line) -> bool {

		if (bdin.getline(line, FILEROWCHARS)) {

			//some file editors will append a carriage return at the end of the line, which messes std::string comparisons here
			//e.g. Python text file write methods (write, writelines) will append a carriage return in addition to a newline
			//if we find a carriage return at the end of the read line simply replace it with a std::string termination
			if (std::string(line).substr(std::string(line).length() - 1) == "\r")
				line[std::string(line).length() - 1] = '\0';

			return true;
		}
		else return false;
	};

	if (bdin.is_open()) {

		char line[FILEROWCHARS];

		//must be an OVF 2.0 file
		if (!read_line(bdin, line) || std::string(line) != headers(OVF2_HEADER)) {

			bdin.close();
			return error(BERROR_COULDNOTLOADFILE_VERSIONMISMATCH);
		}

		while (read_line(bdin, line)) {

			int check = scan_line(line);

			if (!check) {

				//something is wrong : unrecognised or incorrect format
				bdin.close();
				return error(BERROR_COULDNOTLOADFILE);
			}

			if (check == 2) {

				//start of data header found : before getting data make sure all the values are correct
				if (!meshtype_rectangular || valuedim != 1 || meshRect.IsNull() || n == INT3() || h == DBL3() || !data_bytes) {

					//something is wrong
					bdin.close();
					return error(BERROR_COULDNOTLOADFILE);
				}

				//now get data : for binary 4 data there must be the 1234567.0 value present next. for binary 8 data there must be 123456789012345.0 value present
				if (data_bytes == 4) {

					float value;
					bdin.read(reinterpret_cast<char*>(&value), sizeof(float));

					if (value != 1234567.0) {

						//not found check value : incorrect format
						bdin.close();
						return error(BERROR_COULDNOTLOADFILE);
					}
				}
				else if (data_bytes == 8) {

					double value;
					bdin.read(reinterpret_cast<char*>(&value), sizeof(double));

					if (value != 123456789012345.0) {

						//not found check value : incorrect format
						bdin.close();
						return error(BERROR_COULDNOTLOADFILE);
					}
				}
				else if (data_bytes == 1) {

					//text input
				}
				else {

					//incorrect number of bytes
					bdin.close();
					return error(BERROR_COULDNOTLOADFILE);
				}

				//resize data VEC
				if (!data.resize(h, meshRect)) {

					data.clear();
					bdin.close();
					return error(BERROR_OUTOFMEMORY_NCRIT);
				}

				if (data.n != n) {

					//dimensions don't match : something is wrong
					data.clear();
					bdin.close();
					return error(BERROR_COULDNOTLOADFILE);
				}

				for (int k = 0; k < n.z; k++) {
					for (int j = 0; j < n.y; j++) {
						for (int i = 0; i < n.x; i++) {

							if (data_bytes == 4) {

								float V;

								bdin.read(reinterpret_cast<char*>(&V), sizeof(float));

								data[INT3(i, j, k)] = V;
							}
							else if (data_bytes == 8) {

								double V;

								bdin.read(reinterpret_cast<char*>(&V), sizeof(double));

								data[INT3(i, j, k)] = V;
							}
							else if (data_bytes == 1) {

								read_line(bdin, line);

								data[INT3(i, j, k)] = ToNum(trim_leading_spaces(std::string(line)));
							}
						}
					}
				}

				break;
			}
		}

		bdin.close();
	}
	else error(BERROR_COULDNOTOPENFILE);

	return error;
}

template BError OVF2::Read_OVF2_VEC(std::string fileName, VEC<VAL3<float>>& data);
template BError OVF2::Read_OVF2_VEC(std::string fileName, VEC<VAL3<double>>& data);
template BError OVF2::Read_OVF2_VEC(std::string fileName, VEC_VC<VAL3<float>>& data);
template BError OVF2::Read_OVF2_VEC(std::string fileName, VEC_VC<VAL3<double>>& data);

//read an OOMMF OVF2 file containing uniform vector data, and set the data VEC from it
template <typename VECType>
BError OVF2::Read_OVF2_VEC(std::string fileName, VECType& data)
{
	BError error(__FUNCTION__);

	std::ifstream bdin;
	bdin.open(fileName.c_str(), std::ios::in | std::ios::binary);

	//need to find the following lines:
	//
	//"# meshtype: " -> this must be followed by "rectangular"; not reading any other types currently
	//"# meshunit: " -> this must be followed by "m" or "nm"
	//"# valuedim: " -> this must be followed by "3" since this is a vector quantity

	//"# xmin: " -> value to follow
	//"# ymin: " -> value to follow
	//"# zmin: " -> value to follow
	//"# xmax: " -> value to follow
	//"# ymax: " -> value to follow
	//"# zmax: " -> value to follow
	//"# xnodes: " -> value to follow
	//"# ynodes: " -> value to follow
	//"# znodes : " -> value to follow
	//"# xstepsize: " -> value to follow
	//"# ystepsize: " -> value to follow
	//"# zstepsize : " -> value to follow

	//"# Begin: Data " or "# Begin: data " -> type of data to follow; Accepting "binary 4", "binary 8", "text". When this is found then stop searching for any other lines, and if conditions met start getting data.

	bool meshtype_rectangular = false;
	double meshunit = 1.0;
	int valuedim = 0;
	int data_bytes = 0;

	Rect meshRect;
	INT3 n;
	DBL3 h;

	//return 0 if something went wrong : abort loading OVF file
	//return 1 if everything fine, continue
	//return 2 if start of data header found
	auto scan_line = [&](std::string line) -> int {

		if (line.find(headers(MESHTYPE)) != std::string::npos) {

			std::string value = line.substr(headers(MESHTYPE).length());
			if (value != "rectangular") return 0;
			meshtype_rectangular = true;
		}

		else if (line.find(headers(MESHUNIT)) != std::string::npos) {

			std::string value = line.substr(headers(MESHUNIT).length());
			if (value == "m") meshunit = 1.0;
			else if (value == "nm") meshunit = 1e-9;
			else return 0;
		}

		else if (line.find(headers(VALUEDIM)) != std::string::npos) {

			std::string value = line.substr(headers(VALUEDIM).length());
			if (value == "3") valuedim = 3;
			else return 0;
		}

		else if (line.find(headers(XMIN)) != std::string::npos) {

			std::string value = line.substr(headers(XMIN).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.s.x = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(YMIN)) != std::string::npos) {

			std::string value = line.substr(headers(YMIN).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.s.y = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(ZMIN)) != std::string::npos) {

			std::string value = line.substr(headers(ZMIN).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.s.z = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(XMAX)) != std::string::npos) {

			std::string value = line.substr(headers(XMAX).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.e.x = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(YMAX)) != std::string::npos) {

			std::string value = line.substr(headers(YMAX).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.e.y = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(ZMAX)) != std::string::npos) {

			std::string value = line.substr(headers(ZMAX).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.e.z = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(XNODES)) != std::string::npos) {

			std::string value = line.substr(headers(XNODES).length());

			n.x = ToNum(value);
		}

		else if (line.find(headers(YNODES)) != std::string::npos) {

			std::string value = line.substr(headers(YNODES).length());

			n.y = ToNum(value);
		}

		else if (line.find(headers(ZNODES)) != std::string::npos) {

			std::string value = line.substr(headers(ZNODES).length());

			n.z = ToNum(value);
		}

		else if (line.find(headers(XSTEP)) != std::string::npos) {

			std::string value = line.substr(headers(XSTEP).length());

			h.x = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(YSTEP)) != std::string::npos) {

			std::string value = line.substr(headers(YSTEP).length());

			h.y = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(ZSTEP)) != std::string::npos) {

			std::string value = line.substr(headers(ZSTEP).length());

			h.z = (double)ToNum(value) * meshunit;
		}

		else if (lowercase(line).find(lowercase(headers(BEGIN_DATA))) != std::string::npos) {

			//convert to lowercase for comparison as some implementations may choose to have "Data" instead of "data"
			std::string value = lowercase(line.substr(headers(BEGIN_DATA).length()));

			//data headers already in lowercase but better use the lowercase conversion anyway to reduce possibility of bugs in the future
			if (value == lowercase(data_headers(DATA_BINARY4))) {

				data_bytes = 4;
				return 2;
			}
			else if (value == lowercase(data_headers(DATA_BINARY8))) {

				data_bytes = 8;
				return 2;
			}
			else if (value == lowercase(data_headers(DATA_TEXT))) {

				//value of 1 means it's a text input
				data_bytes = 1;
				return 2;
			}
		}

		return 1;
	};

	//read a text line from input stream
	auto read_line = [&](std::ifstream& bdin, char* line) -> bool {

		if (bdin.getline(line, FILEROWCHARS)) {

			//some file editors will append a carriage return at the end of the line, which messes std::string comparisons here
			//e.g. Python text file write methods (write, writelines) will append a carriage return in addition to a newline
			//if we find a carriage return at the end of the read line simply replace it with a std::string termination
			if (std::string(line).substr(std::string(line).length() - 1) == "\r")
				line[std::string(line).length() - 1] = '\0';

			return true;
		}
		else return false;
	};

	if (bdin.is_open()) {

		char line[FILEROWCHARS];

		//must be an OVF 2.0 file
		if (!read_line(bdin, line) || std::string(line) != headers(OVF2_HEADER)) {

			bdin.close();
			return error(BERROR_COULDNOTLOADFILE_VERSIONMISMATCH);
		}

		while (read_line(bdin, line)) {

			int check = scan_line(line);

			if (!check) {

				//something is wrong : unrecognised or incorrect format
				bdin.close();
				return error(BERROR_COULDNOTLOADFILE);
			}

			if (check == 2) {

				//start of data header found : before getting data make sure all the values are correct
				if (!meshtype_rectangular || valuedim != 3 || meshRect.IsNull() || n == INT3() || h == DBL3() || !data_bytes) {

					//something is wrong
					bdin.close();
					return error(BERROR_COULDNOTLOADFILE);
				}

				//now get data : for binary 4 data there must be the 1234567.0 value present next. for binary 8 data there must be 123456789012345.0 value present
				if (data_bytes == 4) {

					float value;
					bdin.read(reinterpret_cast<char*>(&value), sizeof(float));

					if (value != 1234567.0) {

						//not found check value : incorrect format
						bdin.close();
						return error(BERROR_COULDNOTLOADFILE);
					}
				}
				else if (data_bytes == 8) {

					double value;
					bdin.read(reinterpret_cast<char*>(&value), sizeof(double));

					if (value != 123456789012345.0) {

						//not found check value : incorrect format
						bdin.close();
						return error(BERROR_COULDNOTLOADFILE);
					}
				}
				else if (data_bytes == 1) {

					//text input
				}
				else {

					//incorrect number of bytes
					bdin.close();
					return error(BERROR_COULDNOTLOADFILE);
				}

				//resize data VEC
				if (!data.resize(h, meshRect)) {

					data.clear();
					bdin.close();
					return error(BERROR_OUTOFMEMORY_NCRIT);
				}

				if (data.n != n) {

					//dimensions don't match : something is wrong
					data.clear();
					bdin.close();
					return error(BERROR_COULDNOTLOADFILE);
				}

				for (int k = 0; k < n.z; k++) {
					for (int j = 0; j < n.y; j++) {
						for (int i = 0; i < n.x; i++) {

							if (data_bytes == 4) {

								float Vx, Vy, Vz;

								bdin.read(reinterpret_cast<char*>(&Vx), sizeof(float));
								bdin.read(reinterpret_cast<char*>(&Vy), sizeof(float));
								bdin.read(reinterpret_cast<char*>(&Vz), sizeof(float));

								data[INT3(i, j, k)] = FLT3(Vx, Vy, Vz);
							}
							else if (data_bytes == 8) {

								double Vx, Vy, Vz;

								bdin.read(reinterpret_cast<char*>(&Vx), sizeof(double));
								bdin.read(reinterpret_cast<char*>(&Vy), sizeof(double));
								bdin.read(reinterpret_cast<char*>(&Vz), sizeof(double));

								data[INT3(i, j, k)] = DBL3(Vx, Vy, Vz);
							}
							else if (data_bytes == 1) {

								read_line(bdin, line);

								//typically data is space-seaprated, but allow tab-separated text data too
								std::vector<std::string> fields = split(trim_leading_spaces(std::string(line)), { " ", "\t" });

								if (fields.size() >= 3) {

									data[INT3(i, j, k)] = DBL3(ToNum(fields[0]), ToNum(fields[1]), ToNum(fields[2]));
								}
							}
						}
					}
				}

				break;
			}
		}

		bdin.close();
	}
	else error(BERROR_COULDNOTOPENFILE);

	return error;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template BError OVF2::Write_OVF2_VEC(std::string fileName, VEC<FLT3>& data, std::string data_type, double norm);
template BError OVF2::Write_OVF2_VEC(std::string fileName, VEC<DBL3>& data, std::string data_type, double norm);
template BError OVF2::Write_OVF2_VEC(std::string fileName, VEC_VC<FLT3>& data, std::string data_type, double norm);
template BError OVF2::Write_OVF2_VEC(std::string fileName, VEC_VC<DBL3>& data, std::string data_type, double norm);

//write an OOMMF OVF2 file containing uniform vector data
//you can write normalized data to norm (divide by it, no normalization by default)
//you can also choose the type of data output : data_type = bin4 for single precision binary, data_type = bin8 for double precision binary, or data_type = text
template <typename VECType>
BError OVF2::Write_OVF2_VEC(std::string fileName, VECType& data, std::string data_type, double norm)
{
	BError error(__FUNCTION__);

	if (data_type == "bin4") data_type = data_headers(DATA_BINARY4);
	else if (data_type == "bin8") data_type = data_headers(DATA_BINARY8);
	else if (data_type == "text") data_type = data_headers(DATA_TEXT);
	else return error(BERROR_INCORRECTNAME);

	std::ofstream bdout;
	bdout.open(fileName.c_str(), std::ios::out | std::ios::binary);

	ExtractFilenameDirectory(fileName);

	bdout << headers(OVF2_HEADER) << std::endl;
	bdout << "#" << std::endl;
	bdout << "# Segment count: 1" << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(BEGIN_SEGMENT) << std::endl;
	bdout << headers(BEGIN_HEADER) << std::endl;
	bdout << "#" << std::endl;
	bdout << "# Title: " << fileName << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(MESHUNIT) + "m" << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(MESHTYPE) + "rectangular" << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(XMIN) + ToString(data.rect.s.x) << std::endl;
	bdout << headers(YMIN) + ToString(data.rect.s.y) << std::endl;
	bdout << headers(ZMIN) + ToString(data.rect.s.z) << std::endl;
	bdout << headers(XMAX) + ToString(data.rect.e.x) << std::endl;
	bdout << headers(YMAX) + ToString(data.rect.e.y) << std::endl;
	bdout << headers(ZMAX) + ToString(data.rect.e.z) << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(XNODES) + ToString(data.n.x) << std::endl;
	bdout << headers(YNODES) + ToString(data.n.y) << std::endl;
	bdout << headers(ZNODES) + ToString(data.n.z) << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(XSTEP) + ToString(data.h.x) << std::endl;
	bdout << headers(YSTEP) + ToString(data.h.y) << std::endl;
	bdout << headers(ZSTEP) + ToString(data.h.z) << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(VALUEDIM) + "3" << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(END_HEADER) << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(BEGIN_DATA) + data_type << std::endl;

	if (data_type == data_headers(DATA_BINARY4)) {

		float value = 1234567.0;

		char* binary_data = reinterpret_cast<char*>(&value);
		bdout.write(binary_data, sizeof(float));
	}
	else if (data_type == data_headers(DATA_BINARY8)) {

		double value = 123456789012345.0;

		char* binary_data = reinterpret_cast<char*>(&value);
		bdout.write(binary_data, sizeof(double));
	}

	for (int k = 0; k < data.n.z; k++) {
		for (int j = 0; j < data.n.y; j++) {
			for (int i = 0; i < data.n.x; i++) {

				if (data_type == data_headers(DATA_BINARY4)) {

					FLT3 value = data[INT3(i, j, k)] / norm;

					char* binary_data = reinterpret_cast<char*>(&value);
					bdout.write(binary_data, sizeof(FLT3));
				}
				else if (data_type == data_headers(DATA_BINARY8)) {

					DBL3 value = data[INT3(i, j, k)] / norm;
						
					char* binary_data = reinterpret_cast<char*>(&value);
					bdout.write(binary_data, sizeof(DBL3));
				}
				else {

					DBL3 value = data[INT3(i, j, k)] / norm;

					bdout << value.x << " " << value.y << " " << value.z << std::endl;
				}
			}
		}
	}

	bdout << headers(END_DATA) + data_type << std::endl;
	bdout << headers(END_SEGMENT) << std::endl;

	bdout.close();

	return error;
}

template BError OVF2::Write_OVF2_SCA(std::string fileName, VEC<float>& data, std::string data_type);
template BError OVF2::Write_OVF2_SCA(std::string fileName, VEC<double>& data, std::string data_type);
template BError OVF2::Write_OVF2_SCA(std::string fileName, VEC_VC<float>& data, std::string data_type);
template BError OVF2::Write_OVF2_SCA(std::string fileName, VEC_VC<double>& data, std::string data_type);

//write an OOMMF OVF2 file containing uniform scalar data
//you can also choose the type of data output : data_type = bin4 for single precision binary, data_type = bin8 for double precision binary, or data_type = text
template <typename VECType>
BError OVF2::Write_OVF2_SCA(std::string fileName, VECType& data, std::string data_type)
{
	BError error(__FUNCTION__);

	if (data_type == "bin4") data_type = data_headers(DATA_BINARY4);
	else if (data_type == "bin8") data_type = data_headers(DATA_BINARY8);
	else if (data_type == "text") data_type = data_headers(DATA_TEXT);
	else return error(BERROR_INCORRECTNAME);

	std::ofstream bdout;
	bdout.open(fileName.c_str(), std::ios::out | std::ios::binary);

	ExtractFilenameDirectory(fileName);

	bdout << headers(OVF2_HEADER) << std::endl;
	bdout << "#" << std::endl;
	bdout << "# Segment count: 1" << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(BEGIN_SEGMENT) << std::endl;
	bdout << headers(BEGIN_HEADER) << std::endl;
	bdout << "#" << std::endl;
	bdout << "# Title: " << fileName << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(MESHUNIT) + "m" << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(MESHTYPE) + "rectangular" << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(XMIN) + ToString(data.rect.s.x) << std::endl;
	bdout << headers(YMIN) + ToString(data.rect.s.y) << std::endl;
	bdout << headers(ZMIN) + ToString(data.rect.s.z) << std::endl;
	bdout << headers(XMAX) + ToString(data.rect.e.x) << std::endl;
	bdout << headers(YMAX) + ToString(data.rect.e.y) << std::endl;
	bdout << headers(ZMAX) + ToString(data.rect.e.z) << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(XNODES) + ToString(data.n.x) << std::endl;
	bdout << headers(YNODES) + ToString(data.n.y) << std::endl;
	bdout << headers(ZNODES) + ToString(data.n.z) << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(XSTEP) + ToString(data.h.x) << std::endl;
	bdout << headers(YSTEP) + ToString(data.h.y) << std::endl;
	bdout << headers(ZSTEP) + ToString(data.h.z) << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(VALUEDIM) + "1" << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(END_HEADER) << std::endl;
	bdout << "#" << std::endl;
	bdout << headers(BEGIN_DATA) + data_type << std::endl;

	if (data_type == data_headers(DATA_BINARY4)) {

		float value = 1234567.0;

		char* binary_data = reinterpret_cast<char*>(&value);
		bdout.write(binary_data, sizeof(float));
	}
	else if (data_type == data_headers(DATA_BINARY8)) {

		double value = 123456789012345.0;

		char* binary_data = reinterpret_cast<char*>(&value);
		bdout.write(binary_data, sizeof(double));
	}

	for (int k = 0; k < data.n.z; k++) {
		for (int j = 0; j < data.n.y; j++) {
			for (int i = 0; i < data.n.x; i++) {

				if (data_type == data_headers(DATA_BINARY4)) {

					float value = data[INT3(i, j, k)];

					char* binary_data = reinterpret_cast<char*>(&value);
					bdout.write(binary_data, sizeof(float));
				}
				else if (data_type == data_headers(DATA_BINARY8)) {

					double value = data[INT3(i, j, k)];

					char* binary_data = reinterpret_cast<char*>(&value);
					bdout.write(binary_data, sizeof(double));
				}
				else {

					bdout << data[INT3(i, j, k)] << std::endl;
				}
			}
		}
	}

	bdout << headers(END_DATA) + data_type << std::endl;
	bdout << headers(END_SEGMENT) << std::endl;

	bdout.close();

	return error;
}