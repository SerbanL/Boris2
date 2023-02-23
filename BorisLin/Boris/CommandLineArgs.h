#pragma once

#include <vector>
#include <string>
#include <sstream>

////////////////////////////////////
//Version

static int Program_Version = 380;

////////////////////////////////////
//Command-line arguments

struct CLArgs {

	std::string progName = "";
	
	//-p
	std::string server_port = "1542";
	std::string server_pwd = "";

	//-g
	int cudaDevice = -1;

	//-s
	std::string python_script = "";

	//-sp
	std::vector<int> python_script_parallel = {};

	//-sd
	int python_script_deletesource = 0;

	//the input arguments as a string
	std::string command_line_string;

	//show version instead of launching?
	bool show_version = false;

	//show help instead of launching?
	bool show_help = false;

	std::string show_help_text(void)
	{
		std::string help;

		help += "Available command line arguments : \n\n";
		
		help += "-help / -h \n";
		help += "Show help. Doesn't launch program.\n\n";

		help += "-version / -v \n";
		help += "Show program version. Doesn't launch program.\n\n";

		help += "-p serverport (password)\n";
		help += "Set script server port and optional password. Default port 1542, and no password.\n\n";

		help += "-g gpu_select\n";
		help += "Set gpu selection:\n";
		help += "gpu_select = -1: automatic gpu selection (default).\n";
		help += "gpu_select = 0, 1, 2, ...: select indicated gpu if available.\n\n";

		help += "-s python_script_file\n";
		help += "Launch program and execute indicated Python script file (doesn't have to include termination). If path not given then use the default Simulations directory. Other flags shown below modify execution behaviour of Python script.\n\n";

		help += "-sp stride offset1 (offset2, (...))\n";
		help += "When running Python script from console (specified with -s flag), specify how for loops should be executed in parallel.\n";
		help += "Embedded Python scripts will attempt to run for loops in parallel, if they are immediately preceded by the directive: #pragma parallel for\n";
		help += "Default behaviour is to split the for loop iterations evenly between all available GPUs.\n";
		help += "This setting allows specifying the stride and offsets to be executed by this program instance, in this order: stride offset1 offset2 ...\n";
		help += "stride is the total number of GPUs or machines used to execute the script, offset is any number between 0 and stride - 1.\n";
		help += "Example scenario : 2 machines are available, first with 2 GPUs and a second with 1 GPU. In this case 2/3rd of the work should be executed on the first machine, and 1/3rd on the second. Thus use a stride of 3 and give offsets 0 and 1 to the first machine, and offset 2 to the second.\n";
		help += "Thus on the first machine launch: Boris.exe -s script.py -sp 3 0 1, and on the second machine launch: Boris.exe -s script.py -sp 3 2.\n";
		help += "Note, the first machine will automatically split the two offsets between the two GPUs. Warning: do not assign more offsets than there are GPUs on any machine.\n\n";

		help += "-sd 0/1\n";
		help += "After execution of Python script from the command line on a multi-GPU machine, delete any temporary Python script files generated : 1 (default), or keep: 0.\n\n";
		
		return help;
	}

	//extract all available command line arguments and store them in globals above
	void process_command_line_args(int argc, char *argv[]);
};
