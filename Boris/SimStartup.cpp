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

			if (string(line) == STRINGIFY(start_check_updates)) {

				if (bdin.getline(line, FILEROWCHARS)) start_check_updates = ToNum(string(line));
			}

			if (string(line) == STRINGIFY(start_scriptserver)) {

				if (bdin.getline(line, FILEROWCHARS)) start_scriptserver = ToNum(string(line));
			}

			if (string(line) == STRINGIFY(log_errors)) {

				if (bdin.getline(line, FILEROWCHARS)) log_errors = ToNum(string(line));
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

		bdout << STRINGIFY(start_check_updates) << endl;
		bdout << start_check_updates << endl;

		bdout << STRINGIFY(start_scriptserver) << endl;
		bdout << start_scriptserver << endl;

		bdout << STRINGIFY(log_errors) << endl;
		bdout << log_errors << endl;

		bdout.close();
	}
}