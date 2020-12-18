#include "stdafx.h"
#include "Simulation.h"

//save/load flags: start_check_updates, start_scriptserver, log_errors
void Simulation::Load_Startup_Flags(void)
{
	char line[FILEROWCHARS];

	std::ifstream bdin;
	bdin.open(startup_options_file.c_str(), std::ios::in);

	if (bdin.is_open()) {

		while (bdin.getline(line, FILEROWCHARS)) {

			//Check for updates?
			if (std::string(line) == STRINGIFY(start_check_updates)) {

				if (bdin.getline(line, FILEROWCHARS)) start_check_updates = ToNum(std::string(line));
			}

			//Script server start?
			if (std::string(line) == STRINGIFY(start_scriptserver)) {

				if (bdin.getline(line, FILEROWCHARS)) start_scriptserver = ToNum(std::string(line));
			}

			//Script server configuration : port
			if (std::string(line) == STRINGIFY(server_port)) {

				if (bdin.getline(line, FILEROWCHARS)) server_port = std::string(line);
			}

			//Script server configuration : password
			if (std::string(line) == STRINGIFY(server_pwd)) {

				if (bdin.getline(line, FILEROWCHARS)) server_pwd = std::string(line);
			}

			//Script server configuration : sleep
			if (std::string(line) == STRINGIFY(server_recv_sleep_ms)) {

				if (bdin.getline(line, FILEROWCHARS)) server_recv_sleep_ms = ToNum(std::string(line));
			}

			//Log errors?
			if (std::string(line) == STRINGIFY(log_errors)) {

				if (bdin.getline(line, FILEROWCHARS)) log_errors = ToNum(std::string(line));
			}

			//Set user-defined number of threads? (0 means the maximum number of available threads will be set, otherwise set indicated number)
			if (std::string(line) == STRINGIFY(OmpThreads)) {

				if (bdin.getline(line, FILEROWCHARS)) OmpThreads = ToNum(std::string(line));
				if (OmpThreads == 0 || OmpThreads > omp_get_num_procs()) OmpThreads = omp_get_num_procs();
			}
		}

		bdin.close();
	}
}

void Simulation::Save_Startup_Flags(void)
{
	std::ofstream bdout;
	bdout.open(startup_options_file.c_str(), std::ios::out);

	if (bdout.is_open()) {

		//Check for updates?
		bdout << STRINGIFY(start_check_updates) << std::endl;
		bdout << start_check_updates << std::endl;

		//Script server start?
		bdout << STRINGIFY(start_scriptserver) << std::endl;
		bdout << start_scriptserver << std::endl;

		//Script server configuration : port
		bdout << STRINGIFY(server_port) << std::endl;
		bdout << server_port << std::endl;

		//Script server configuration : password
		bdout << STRINGIFY(server_pwd) << std::endl;
		bdout << server_pwd << std::endl;

		//Script server configuration : sleep
		bdout << STRINGIFY(server_recv_sleep_ms) << std::endl;
		bdout << server_recv_sleep_ms << std::endl;

		//Log errors?
		bdout << STRINGIFY(log_errors) << std::endl;
		bdout << log_errors << std::endl;

		//Set user-defined number of threads? (0 means the maximum number of available threads will be set, otherwise set indicated number)
		bdout << STRINGIFY(OmpThreads) << std::endl;

		if (OmpThreads == omp_get_num_procs()) bdout << 0 << std::endl;
		else bdout << OmpThreads << std::endl;

		bdout.close();
	}
}