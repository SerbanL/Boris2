#include "stdafx.h"
#include "Boris.h"

#include "CommandLineArgs.h"

#include "CompileFlags.h"
#if GRAPHICS == 0

#if PYTHON_EMBEDDING == 1
#include "PythonScripting.h"
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Separate function to make the Simulation object, so it can be called on a separate std::thread.
//See https://stackoverflow.com/questions/64988525/bug-related-to-g-openmp-when-using-stdthread/65001887#65001887 for reason
void make_Simulation(CLArgs clargs)
{
	pSim = new Simulation(
		Program_Version, clargs.progName,
		clargs.server_port, clargs.server_pwd,
		clargs.cudaDevice,
		clargs.python_script, clargs.python_script_parallel, clargs.python_script_deletesource);

#if PYTHON_EMBEDDING == 1
	//Startup embedded Python interpreter on this thread
	PyImport_AppendInittab("BPython", &PyInit_PythonEmbeddedModule);
	Py_Initialize();

	//Append to path standard working directories (the Simulations, BorisPythonScripts, and current executable directories)
	PyObject* sysPath = PySys_GetObject((char*)"path");

	//Simulations directory in documents
	std::string pathToModuleDirectory = pSim->Get_Simulations_Directory();
	PyList_Append(sysPath, (PyUnicode_FromString(pathToModuleDirectory.c_str())));

	//the executable directory
	pathToModuleDirectory = GetExeDirectory();
	PyList_Append(sysPath, (PyUnicode_FromString(pathToModuleDirectory.c_str())));

	//BorisPythonScripts directory, found in Documents/Boris Data/ (this is where NetSocks is found)
	pathToModuleDirectory = pSim->Get_BorisPythonScripts_Directory();
	PyList_Append(sysPath, (PyUnicode_FromString(pathToModuleDirectory.c_str())));
#endif

	//Run embedded Python script specified on command line?
	if (clargs.python_script.length()) pSim->NewMessage("runscript " + clargs.python_script);
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

	CLArgs clargs;
	clargs.process_command_line_args(argc, argv);

	if (argc > 1) {

		if (clargs.show_version) {

			//show version instead of launching

			std::stringstream version;
			version << "BORIS Computational Spintronics 2022, version " << (double)Program_Version / 100.0;
			std::cout << version.str() << std::endl;

			return 0;
		}

		if (clargs.show_help) {

			//show help instead of launching

			std::cout << clargs.show_help_text() << std::endl;

			return 0;
		}
	}

#if OPERATING_SYSTEM == OS_WIN
	ShowWindow(GetConsoleWindow(), SW_SHOWMAXIMIZED);
#else

#endif

	//////////////////////
	//Instantiate Simulation object and start main simulation loop thread (non-blocking)
	
	//See https://stackoverflow.com/questions/64988525/bug-related-to-g-openmp-when-using-stdthread/65001887#65001887
	std::thread simulation_instantiation_thread(&make_Simulation, clargs);
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