#include "stdafx.h"
#include "Boris.h"

#include "CompileFlags.h"
#if GRAPHICS == 0

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Program master object
Simulation *pSim = nullptr;

//Separate function to make the Simulation object, so it can be called on a separate std::thread.
//See https://stackoverflow.com/questions/64988525/bug-related-to-g-openmp-when-using-stdthread/65001887#65001887 for reason
void make_Simulation(void)
{
	pSim = new Simulation(Program_Version);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TextInput :
	public Threads<TextInput>
{

	string message;

private:

	void get_user_command(Simulation *pSim)
	{
		getline(cin, message);
	}

private:

	void input_loop(Simulation *pSim)
	{
		set_blocking_thread(THREAD_GETINPUT);
		single_call_launch<Simulation*>(&TextInput::get_user_command, pSim, THREAD_GETINPUT);

		if (message.length() && message != "exit") {

			pSim->NewMessage(message);
		}
	}

public:

	TextInput(Simulation *pSim)
	{
		while (message != "exit") {

			input_loop(pSim);
		}
	}
};

int main(void)
{
#if OPERATING_SYSTEM == OS_WIN
	ShowWindow(GetConsoleWindow(), SW_SHOWMAXIMIZED);
#else

#endif

	//Instantiate Simulation object and start main simulation loop thread (non-blocking)
	
	//See https://stackoverflow.com/questions/64988525/bug-related-to-g-openmp-when-using-stdthread/65001887#65001887
	std::thread simulation_instantiation_thread(&make_Simulation);
	simulation_instantiation_thread.join();

	//Get user input until exit condition reached
	TextInput textInput(pSim);

	delete pSim;

	return 0;
}

#endif