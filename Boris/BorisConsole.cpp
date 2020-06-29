#include "stdafx.h"
#include "Boris.h"

#include "CompileFlags.h"
#if GRAPHICS == 0

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

	Simulation *pSim;

	//Instantiate Simulation object and start main simulation loop thread (non-blocking)
	pSim = new Simulation(Program_Version);

	//Get user input until exit condition reached
	TextInput textInput(pSim);

	delete pSim;

	return 0;
}

#endif