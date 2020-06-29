#include "stdafx.h"
#include "Simulation.h"

#if GRAPHICS == 1
void Simulation::NewMessage(AC_ aCode, INT2 mouse, string data)
{
	//console command received - process it and return

	if (aCode == AC_CONSOLECOMMAND || aCode == AC_CONSOLEENTRY || aCode == AC_CONSOLECOMMAND_ENTRY || aCode == AC_CONSOLECOMMAND_NOPARAMS_ENTRY) {

		//data has the command with any parameters
		//if there are no parameters then there are no spaces in data, since commands do not include spaces
		//if there are spaces then there must be parameters
		//first entry before any space is the command id number (entry in CMD_ enum) but as a string - convert it
		//all other fields after first space are command parameters, so pass them in as a string

		auto convert_data = [&](string& data) -> string {

			string command_with_parameters;

			vector<string> fields = split(data, " ");

			if (fields.size() > 1) {

				int command_index = ToNum(fields[0]);
				command_with_parameters = commands.get_key_from_index(command_index) + " " + combine(subvec(fields, 1), " ");
			}
			else {

				int command_index = ToNum(data);
				command_with_parameters = commands.get_key_from_index(command_index);
			}

			return command_with_parameters;
		};

		switch (aCode) {

		case AC_CONSOLEENTRY:
		{
			single_call_launch<string>(&Simulation::SetConsoleEntryLineText, data, THREAD_HANDLEMESSAGE2);
		}
		break;

		case AC_CONSOLECOMMAND_ENTRY:
		case AC_CONSOLECOMMAND:
		{
			string command_with_parameters = convert_data(data);

			set_blocking_thread(THREAD_HANDLEMESSAGE);
			if (single_call_launch<string>(&Simulation::HandleCommand, command_with_parameters, THREAD_HANDLEMESSAGE) != THREAD_HANDLEMESSAGE) {

				err_hndl.show_error(BERROR_BUSY);
			}
		}
		break;

		case AC_CONSOLECOMMAND_NOPARAMS_ENTRY:
		{
			string command_with_parameters = convert_data(data);

			vector<string> fields = split(command_with_parameters, " ");

			set_blocking_thread(THREAD_HANDLEMESSAGE);
			if (single_call_launch<string>(&Simulation::HandleCommand, fields[0], THREAD_HANDLEMESSAGE) != THREAD_HANDLEMESSAGE) {

				err_hndl.show_error(BERROR_BUSY);
			}
		}
		break;

		}

		if (aCode == AC_CONSOLECOMMAND_ENTRY || aCode == AC_CONSOLECOMMAND_NOPARAMS_ENTRY) {

			single_call_launch<string>(&Simulation::SetConsoleEntryLineText, convert_data(data) + " ", THREAD_HANDLEMESSAGE2);
		}

		return;
	}

	//not a console command

	//Dispatch message and check if anything needs to be done here as a result (e.g. a console command entered)
	ActionOutcome result = BD.NewMessage_ThreadSafe(aCode, mouse, data);

	//Check for special action outcomes which must be handled by the top object (Simulation)

	//dispatch command to command handler - call it on its own unique thread - will not get called if that thread is already active (i.e. still processing previous command
	if (result.IsCodeSet(AO_MESSAGERETURNED)) {

		//full command formed (enter key pressed)
		single_call_launch<string>(&Simulation::HandleCommand, result.text, THREAD_HANDLEMESSAGE);
	}

	//focus mesh but keep camera orientation
	else if (result.IsCodeSet(AO_MESHFOCUS2)) {

		single_call_launch<string>(&Simulation::HandleCommand, commands.get_key_from_index(CMD_MESHFOCUS2) + " " + result.text, THREAD_HANDLEMESSAGE);
	}

	//text entered in console
	else if (result.IsCodeSet(AO_TEXTRETURNED)) {

		//try to autocomplete after a key press (and do not allow incorrect commands)

		//only try to autocomplete the command word (not the parameters): as soon as a space is entered then command word is considered formed.
		if (result.text.find(" ") != string::npos) return;

		//get all commands which contain the returned text at the start
		string consoleLineText;
		//if first character is '?' don't include it in autocomplete
		if (result.text[0] == '?') { result.text = result.text.substr(1); consoleLineText = "?"; }

		vector<int> indexes = commands.find_keystart(result.text);

		//if no matches found then delete last character (e.g. wrong character entered) as long as it results in a correct partial command word
		if (!indexes.size()) {

			string newText = result.text.substr(0, result.text.length() - 1);

			indexes = commands.find_keystart(newText);
			if (indexes.size()) {

				consoleLineText += newText;
				SetConsoleEntryLineText(consoleLineText);
			}
		}
		else {

			//if only one match then autocomplete
			if (indexes.size() == 1) {

				consoleLineText += commands.get_key_from_index(indexes[0]);
				SetConsoleEntryLineText(consoleLineText);
			}
		}
	}

	//must update screen
	else if (result.IsCodeSet(AO_RECALCULATEMESHDISPLAY)) {

		single_call_launch<string>(&Simulation::HandleCommand, commands.get_key_from_index(CMD_UPDATESCREEN), THREAD_HANDLEMESSAGE);
	}

	//load simulation file
	else if (result.IsCodeSet(AO_FILEDROPPEDINCONSOLE)) {

		string command = commands.get_key_from_index(CMD_LOADSIM) + " " + result.text;

		single_call_launch<string>(&Simulation::HandleCommand, command, THREAD_HANDLEMESSAGE);
	}

	//load mask file
	else if (result.IsCodeSet(AO_FILEDROPPEDINMESH)) {

		string command = commands.get_key_from_index(CMD_LOADMASKFILE) + " " + result.text;

		single_call_launch<string>(&Simulation::HandleCommand, command, THREAD_HANDLEMESSAGE);
	}
}
#else

void Simulation::NewMessage(string message)
{
	//full command formed (enter key pressed)
	single_call_launch<string>(&Simulation::HandleCommand, message, THREAD_HANDLEMESSAGE);
}

#endif