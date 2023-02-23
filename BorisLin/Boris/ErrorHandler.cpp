#include "stdafx.h"
#include "Simulation.h"

//-------------------------------------Error Handler methods

//intended to create a restore point to which we can go back if a critical error occurs
void Simulation::create_restore(void)
{ 
	//Not saving a restore point currently
}

//restore state after a critical error occurs
void Simulation::restore_state(void)
{
	StopSimulation();

	//restore last saved file
	BError error = LoadSimulation(directory + currentSimulationFile);

	if (error) {

		//restore.bsm failed : go back to default
		error.reset() = LoadSimulation(GetUserDocumentsPath() + boris_data_directory + boris_simulations_directory + "default");

		if (!error) {

			BD.DisplayConsoleError("Recovered from critical error.");
		}
		else {

			//It would take a bug to get here or else something extremely unlikely
			BD.DisplayConsoleError("!!! Could not recover from critical error !!! Program unsafe : " + error.info());
		}
	}
	
	return;
}

//show the error
void Simulation::show_error(BError berror, std::string error_message, bool verbose)
{ 
	if (verbose) {

		//if (berror.warning_set()) BD.DisplayConsoleWarning("WARNING : " + err_text);
		//else BD.DisplayConsoleError("ERROR : " + err_text);

		if (berror.warning_set()) BD.DisplayConsoleWarning(error_message);
		else BD.DisplayConsoleError(error_message);
	}

	if (log_errors) {

		std::ofstream bdout;
		bdout.open(errorlog_fileName.c_str(), std::ios::out | std::ios::app);
		bdout << Get_Date_Time() + error_message << std::endl;
		bdout.close();
	}

	UpdateScreen();
}