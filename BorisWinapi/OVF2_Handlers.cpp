#include "stdafx.h"
#include "OVF2_Handlers.h"


OVF2::OVF2(void)
{
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
	headers.push_back("# Begin: Data ", BEGIN_DATA_UC);
	headers.push_back("# Begin: data ", BEGIN_DATA_LC);
}

//read an OOMMF OVF2 file containing uniform vector data, and set the data VEC from it
BError OVF2::Read_OVF2_VEC(std::string fileName, VEC<DBL3>& data)
{
	BError error(__FUNCTION__);

	std::ifstream bdin;
	bdin.open(fileName.c_str(), std::ios::in + std::ios::binary);

	//need to find the following lines:
	//
	//"# meshtype: " -> this must be followed by "rectangular"; not reading any other types currently
	//"# meshunit: " -> this must be followed by "m" or "nm"
	//"# valuedim: " -> this must be followed by "3"; not reading any other types currently

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

	//"# Begin: Data " or "# Begin: data " -> type of data to follow; currently accepting "Binary 4" / "binary 4", "Binary 8" / "binary 8". When this is found then stop searching for any other lines, and if conditions met start getting data.

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
	auto scan_line = [&](string line) -> int {

		if (line.find(headers(MESHTYPE)) != std::string::npos) {

			string value = line.substr(headers(MESHTYPE).length());
			if (value != "rectangular") return 0;
			meshtype_rectangular = true;
		}

		else if (line.find(headers(MESHUNIT)) != std::string::npos) {

			string value = line.substr(headers(MESHUNIT).length());
			if (value == "m") meshunit = 1.0;
			else if (value == "nm") meshunit = 1e-9;
			else return 0;
		}

		else if (line.find(headers(VALUEDIM)) != std::string::npos) {

			string value = line.substr(headers(VALUEDIM).length());
			if (value == "3") valuedim = 3;
			else return 0;
		}

		else if (line.find(headers(XMIN)) != std::string::npos) {

			string value = line.substr(headers(XMIN).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.s.x = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(YMIN)) != std::string::npos) {

			string value = line.substr(headers(YMIN).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.s.y = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(ZMIN)) != std::string::npos) {

			string value = line.substr(headers(ZMIN).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.s.z = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(XMAX)) != std::string::npos) {

			string value = line.substr(headers(XMAX).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.e.x = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(YMAX)) != std::string::npos) {

			string value = line.substr(headers(YMAX).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.e.y = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(ZMAX)) != std::string::npos) {

			string value = line.substr(headers(ZMAX).length());
			//meshunit should have been read by now (value of 1.0 is the default)
			meshRect.e.z = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(XNODES)) != std::string::npos) {

			string value = line.substr(headers(XNODES).length());

			n.x = ToNum(value);
		}

		else if (line.find(headers(YNODES)) != std::string::npos) {

			string value = line.substr(headers(YNODES).length());

			n.y = ToNum(value);
		}

		else if (line.find(headers(ZNODES)) != std::string::npos) {

			string value = line.substr(headers(ZNODES).length());

			n.z = ToNum(value);
		}

		else if (line.find(headers(XSTEP)) != std::string::npos) {

			string value = line.substr(headers(XSTEP).length());

			h.x = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(YSTEP)) != std::string::npos) {

			string value = line.substr(headers(YSTEP).length());

			h.y = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(ZSTEP)) != std::string::npos) {

			string value = line.substr(headers(ZSTEP).length());

			h.z = (double)ToNum(value) * meshunit;
		}

		else if (line.find(headers(BEGIN_DATA_UC)) != std::string::npos) {

			string value = line.substr(headers(BEGIN_DATA_UC).length());

			if (value == "Binary 4" || value == "binary 4") {

				data_bytes = 4;
				return 2;
			}
			else if (value == "Binary 8" || value == "binary 8") {

				data_bytes = 8;
				return 2;
			}
		}

		else if (line.find(headers(BEGIN_DATA_LC)) != std::string::npos) {

			string value = line.substr(headers(BEGIN_DATA_LC).length());

			if (value == "Binary 4" || value == "binary 4") {

				data_bytes = 4;
				return 2;
			}
			else if (value == "Binary 8" || value == "binary 8") {

				data_bytes = 8;
				return 2;
			}
		}

		return 1;
	};

	if (bdin.is_open()) {

		char line[FILEROWCHARS];

		//must be an OVF 2.0 file
		if (!bdin.getline(line, FILEROWCHARS) || std::string(line) != "# OOMMF OVF 2.0") {

			bdin.close();
			return error(BERROR_COULDNOTLOADFILE_VERSIONMISMATCH);
		}

		while (bdin.getline(line, FILEROWCHARS)) {

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

								data[INT3(i, j, k)] = DBL3(Vx, Vy, Vz);
							}
							else if (data_bytes == 4) {

								double Vx, Vy, Vz;

								bdin.read(reinterpret_cast<char*>(&Vx), sizeof(double));
								bdin.read(reinterpret_cast<char*>(&Vy), sizeof(double));
								bdin.read(reinterpret_cast<char*>(&Vz), sizeof(double));

								data[INT3(i, j, k)] = DBL3(Vx, Vy, Vz);
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