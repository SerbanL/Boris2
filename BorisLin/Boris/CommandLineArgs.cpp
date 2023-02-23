#include "stdafx.h"

#include "CommandLineArgs.h"

enum CLARGS_ {

	CLARGS_SHOWVERSION,

	CLARGS_SHOWHELP,

	CLARGS_SERVERPORT,

	CLARGS_CUDADEVICE,

	CLARGS_PYSCRIPT,

	CLARGS_PYSCRIPT_PARALLELSPLIT,

	CLARGS_PYSCRIPT_DELETESOURCE
};

void CLArgs::process_command_line_args(int argc, char *argv[])
{
	std::vector<CLARGS_> clargs;
	std::vector<std::vector<std::string>> clargs_fields;
	std::vector<int> argnum;

	command_line_string.clear();

	progName = argv[0];
	command_line_string += progName;

	//first scan all input arguments, and extract their type (in clargs) and argument numbers (in argnum)
	for (int i = 1; i < argc; i++) {

		std::string arg(argv[i]);
		command_line_string += " " + arg;

		if (arg == "-version" || arg == "-v") {

			clargs.push_back(CLARGS_SHOWVERSION);
			argnum.push_back(i);
		}
		if (arg == "-help" || arg == "-h") {

			clargs.push_back(CLARGS_SHOWHELP);
			argnum.push_back(i);
		}
		if (arg == "-p") {

			clargs.push_back(CLARGS_SERVERPORT);
			argnum.push_back(i);
		}
		else if (arg == "-g") {

			clargs.push_back(CLARGS_CUDADEVICE);
			argnum.push_back(i);
		}
		else if (arg == "-s") {

			clargs.push_back(CLARGS_PYSCRIPT);
			argnum.push_back(i);
		}
		else if (arg == "-sp") {

			clargs.push_back(CLARGS_PYSCRIPT_PARALLELSPLIT);
			argnum.push_back(i);
		}
		else if (arg == "-sd") {

			clargs.push_back(CLARGS_PYSCRIPT_DELETESOURCE);
			argnum.push_back(i);
		}
	}

	if (!clargs.size()) { show_help = true; return; }

	//now take a second pass and extract fields for each input argument detected
	for (int idx = 0; idx < argnum.size(); idx++) {

		//argument index in argv
		int i = argnum[idx];

		//next argument index in argv
		int j = argc;
		if (idx < argnum.size() - 1) j = argnum[idx + 1];

		//thus the fields for this argument stretch from i + 1 to j - 1 inclusive
		std::vector<std::string> clarg_fields;
		for (int fidx = i + 1; fidx < j; fidx++) {

			clarg_fields.push_back(argv[fidx]);
		}

		//thus the fields for this argument stretch from i + 1 to j - 1 inclusive
		if (clargs[idx] != CLARGS_PYSCRIPT) {

			clargs_fields.push_back(clarg_fields);
		}
		else {

			std::string text;
			for (int idx = 0; idx < clarg_fields.size(); idx++) {

				text += clarg_fields[idx];
				if (idx != clarg_fields.size() - 1) text += " ";
			}

			clargs_fields.push_back({ text });
		}
	}

	//finally set globals after a final check the fields are sensible
	for (int idx = 0; idx < clargs.size(); idx++) {

		switch (clargs[idx]) {

		case CLARGS_SHOWVERSION:
			show_version = true;
			break;

		case CLARGS_SHOWHELP:
			show_help = true;
			break;

		case CLARGS_SERVERPORT:
			if (clargs_fields[idx].size() >= 1) server_port = clargs_fields[idx][0];
			if (clargs_fields[idx].size() >= 2) server_pwd = clargs_fields[idx][1];
			break;

		case CLARGS_CUDADEVICE:
			if (clargs_fields[idx].size() >= 1) {

				std::stringstream ss(clargs_fields[idx][0]);
				ss >> cudaDevice;
			}
			break;

		case CLARGS_PYSCRIPT:
			if (clargs_fields[idx].size() >= 1) python_script = clargs_fields[idx][0];
			break;

		case CLARGS_PYSCRIPT_PARALLELSPLIT:
		{
			int stride = 1;
			if (clargs_fields[idx].size() >= 1) {

				std::stringstream ss(clargs_fields[idx][0]);
				ss >> stride;
				if (stride < 1) break;
			}

			//stride from 1 upwards
			python_script_parallel.push_back(stride);

			if (clargs_fields[idx].size() >= 2) {

				for (int fidx = 1; fidx < clargs_fields[idx].size(); fidx++) {

					int offset;
					std::stringstream ss(clargs_fields[idx][fidx]);
					ss >> offset;
					if (offset < stride && offset >= 0) python_script_parallel.push_back(offset);
					else {

						//if a mistake was made in specification then disable this setting
						python_script_parallel.clear();
						break;
					}
				}
			}
			//need at least 2 fields here
			else {

				python_script_parallel.clear();
				break;
			}
		}
			break;

		case CLARGS_PYSCRIPT_DELETESOURCE:
			if (clargs_fields[idx].size() >= 1) {

				std::stringstream ss(clargs_fields[idx][0]);
				ss >> python_script_deletesource;
			}
			break;
		}
	}
}
