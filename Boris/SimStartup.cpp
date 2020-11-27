#include "stdafx.h"
#include "Simulation.h"

//save/load flags: start_check_updates, start_scriptserver, log_errors
void Simulation::Load_Startup_Flags(void)
{
	char line[FILEROWCHARS];

	ifstream bdin;
	bdin.open(startup_options_file.c_str(), ios::in);

	if (bdin.is_open()) {

		while (bdin.getline(line, FILEROWCHARS)) {

			//Check for updates?
			if (string(line) == STRINGIFY(start_check_updates)) {

				if (bdin.getline(line, FILEROWCHARS)) start_check_updates = ToNum(string(line));
			}

			//Script server start?
			if (string(line) == STRINGIFY(start_scriptserver)) {

				if (bdin.getline(line, FILEROWCHARS)) start_scriptserver = ToNum(string(line));
			}

			//Script server configuration : port
			if (string(line) == STRINGIFY(server_port)) {

				if (bdin.getline(line, FILEROWCHARS)) server_port = string(line);
			}

			//Script server configuration : sleep
			if (string(line) == STRINGIFY(server_recv_sleep_ms)) {

				if (bdin.getline(line, FILEROWCHARS)) server_recv_sleep_ms = ToNum(string(line));
			}

			//Log errors?
			if (string(line) == STRINGIFY(log_errors)) {

				if (bdin.getline(line, FILEROWCHARS)) log_errors = ToNum(string(line));
			}

			//Set user-defined number of threads? (0 means the maximum number of available threads will be set, otherwise set indicated number)
			if (string(line) == STRINGIFY(OmpThreads)) {

				if (bdin.getline(line, FILEROWCHARS)) OmpThreads = ToNum(string(line));
				if (OmpThreads == 0 || OmpThreads > omp_get_num_procs()) OmpThreads = omp_get_num_procs();
			}
		}

		bdin.close();
	}
}

void Simulation::Save_Startup_Flags(void)
{
	ofstream bdout;
	bdout.open(startup_options_file.c_str(), ios::out);

	if (bdout.is_open()) {

		//Check for updates?
		bdout << STRINGIFY(start_check_updates) << endl;
		bdout << start_check_updates << endl;

		//Script server start?
		bdout << STRINGIFY(start_scriptserver) << endl;
		bdout << start_scriptserver << endl;

		//Script server configuration : port
		bdout << STRINGIFY(server_port) << endl;
		bdout << server_port << endl;

		//Script server configuration : sleep
		bdout << STRINGIFY(server_recv_sleep_ms) << endl;
		bdout << server_recv_sleep_ms << endl;

		//Log errors?
		bdout << STRINGIFY(log_errors) << endl;
		bdout << log_errors << endl;

		//Set user-defined number of threads? (0 means the maximum number of available threads will be set, otherwise set indicated number)
		bdout << STRINGIFY(OmpThreads) << endl;

		if (OmpThreads == omp_get_num_procs()) bdout << 0 << endl;
		else bdout << OmpThreads << endl;

		bdout.close();
	}
}