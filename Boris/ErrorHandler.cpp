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
void Simulation::show_error(BError berror, string error_message, bool verbose)
{ 
	string err_text = error_message + " Info : " + berror.info();

	if (verbose) BD.DisplayConsoleError("ERROR : " + err_text);

	if (log_errors) {

		ofstream bdout;
		bdout.open(errorlog_fileName.c_str(), ios::out | ios::app);
		bdout << Get_Date_Time() + err_text << endl;
		bdout.close();
	}

	UpdateScreen();
}