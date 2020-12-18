#include "stdafx.h"
#include "Boris.h"

#include "CompileFlags.h"
#if GRAPHICS == 0

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Separate function to make the Simulation object, so it can be called on a separate std::thread.
//See https://stackoverflow.com/questions/64988525/bug-related-to-g-openmp-when-using-stdthread/65001887#65001887 for reason
void make_Simulation(void)
{
	pSim = new Simulation(Program_Version, server_port, server_pwd, cudaDevice);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TextInput :
	public Threads<TextInput>
{

	std::string message;

private:

	void get_user_command(Simulation *pSim)
	{
		getline(std::cin, message);
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

int main(int argc, char *argv[])
{
	//////////////////////
	//Arguments

	for (int i = 1; i < argc; i++) {

		//First argument: server port
		if (i == 1) {

			server_port = std::string(argv[i]);
		}

		//Second argument: cuda device
		if (i == 2) {

			cudaDevice = ToNum(std::string(argv[i]));
		}

		//Third argument: window options (front/back)
		if (i == 3) {

			window_startup_option = std::string(argv[i]);
		}

		//Fourth argument: server password
		if (i == 4) {

			server_pwd = std::string(argv[i]);
		}
	}

#if OPERATING_SYSTEM == OS_WIN
	ShowWindow(GetConsoleWindow(), SW_SHOWMAXIMIZED);
#else

#endif

	//////////////////////
	//Instantiate Simulation object and start main simulation loop thread (non-blocking)
	
	//See https://stackoverflow.com/questions/64988525/bug-related-to-g-openmp-when-using-stdthread/65001887#65001887
	std::thread simulation_instantiation_thread(&make_Simulation);
	simulation_instantiation_thread.join();

	//////////////////////
	//Get user input until exit condition reached

	TextInput textInput(pSim);

	//////////////////////
	//Finish

	delete pSim;

	return 0;
}

#endif