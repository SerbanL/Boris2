#pragma once

#include <vector>
#include <string>

using namespace std;

//BERROR_ enum with all the possible error identifiers
enum BERROR_ 
{ 
	BERROR_NONE = 0, 
	BERROR_COMMAND_NOTRECOGNIZED,			//command not recognized
	BERROR_BUSY,							//An action could not be completed: program busy
	BERROR_OPERATIONFAILED,					//Failed in an unspecified way
	BERROR_OUTOFMEMORY_CRIT,				//Out of memory : critical, must restore program
	BERROR_OUTOFMEMORY_NCRIT,				//Out of memory : not critical, don't need restore
	BERROR_OUTOFGPUMEMORY_CRIT,				//Out of gpu memory : critical, must restore program
	BERROR_OUTOFGPUMEMORY_NCRIT,			//Out of gpu memory : not critical, program can still work
	BERROR_GPUERROR_CRIT,					//gpu error : critical, must restore program
	BERROR_INCORRECTMODCONFIG,				//Modules not configured correctly. Cannot start simulation
	BERROR_INCORRECTCONFIG,					//Incorrect configuration (e.g. eval method doesn't match ode allowed evals)
	BERROR_INCORRECTARRAYS,					//Incorrect arrays used
	BERROR_INCORRECTVALUE,					//Incorrect value used
	BERROR_INCORRECTSTRING,					//Incorrect string used
	BERROR_PARAMMISMATCH,					//Getting parameters from console command : parameters entered do not match expected parameters
	BERROR_PARAMMISMATCH_SHOW,				//not silent
	BERROR_PARAMOUTOFBOUNDS,				//Getting parameters from console command : a parameter is out of specified bounds
	BERROR_INCORRECTNAME,					//An incorrect name was given (e.g. a mesh name)
	BERROR_INCORRECTACTION,					//The action of a command is incorrect for the current state and was not executed
	BERROR_INCORRECTACTION_SILENT,			//The action of a command is incorrect for the current state and was not executed - don't show error
	BERROR_COULDNOTOPENFILE,				//A file couldn't be opened
	BERROR_COULDNOTLOADFILE,				//A file couldn't be loaded
	BERROR_COULDNOTLOADFILE_VERSIONMISMATCH,//A file couldn't be loaded because there is a mismatch in version number
	BERROR_COULDNOTLOADFILE_CRIT,			//A file couldn't be loaded - critical error (e.g. failed to load simulation file)
	BERROR_COULDNOTSAVEFILE,				//A file couldn't be saved
	BERROR_COULDNOTCONNECT,					//Could not connect to a server
	BERROR_NOTAVAILABLE,					//Requested configuration not available for this program version / workstation (e.g. no cuda-enabled gpu)
	BERROR_NOTMAGNETIC,						//Focused mesh must be magnetic.
	BERROR_NOTANTIFERROMAGNETIC,			//Focused mesh must be anti-ferromagnetic.
	BERROR_NOTATOMISTIC,					//Focused mesh must be atomistic.
	BERROR_MESHNAMEINEXISTENT,				//Mesh name doesn't exist.
	BERROR_NOTCOMPUTED,						//Not computed.
	BERROR_SPINSOLVER_FIT,					//Must be ferromagnetic mesh with transport module added and spin transport solver enabled.Must also have either Ts or Tsi computed.
	BERROR_SPINSOLVER_FIT2,					//Must be ferromagnetic mesh with transport module added and spin transport solver enabled.
	BERROR_SPINSOLVER_FIT3,					//Must be ferromagnetic mesh with transport module added and spin transport solver enabled.hm_mesh must be a metal mesh with transport module added.
	BERROR_SPINSOLVER_FIT4,					//Must give metal and ferromagnetic meshes in this order.
	BERROR_NOTDEFINED,						//Name not defined
	BERROR_ENUMSIZE
};

enum ERRLEV_ 
{ 
	ERRLEV_SILENT,							//Error occured but don't display it
	ERRLEV_NCRIT,							//More serious error so should display it
	ERRLEV_CRIT								//Critical error : must use restore point to restore program state
};

//the Error Handler should be held by a top object which implements the methods :
//
//1. void create_restore(void); -> intended to create a restore point to which we can go back if a critical error occurs
//2. void restore_state(BERROR_ berror, string error_message); -> restore state after a critical error occurs
//3. void show_error(BERROR_ berror, string error_message); -> show the error
//
//to enable error handling, a method which returns a BError object should be invoked via "ErrorHandler::call" or "ErrorHandler::qcall"
//"ErrorHandler::call" is the full method with restore point
//"ErrorHandler::qcall" doesn't create/load a restore point and is meant for methods which are guaranteed not to return a critical error
//
//Functions which return BError objects should first declare it:
//
//BError error("some local info");
//
//If an error occurs in this function set it as : error(BERROR_SOME_ERROR_CODE);
//
//At the end return the error object.
//The function may also invoke other functions which return a BError object. 
//The top function will contain the error info: error.info(); 
//and code: error.code() if set: 
//if(error) 
//{ 
//	string error_info = error.info(); 
//	BERROR_ error_code = error.code(); 
//}

//------------------------------------------------------------------------------- ERROR CONTAINER

class BError {

	//the error code
	BERROR_ err_code = BERROR_NONE;

	//info for the error code containing the whole chain of functions returning errors, each with their own info
	vector<string> err_info;

public:

	BError(void) {}

	BError(BERROR_ err_code_) 
	{
		err_code = err_code_;
	}

	BError(BERROR_ err_code_, string local_info)
	{
		err_code = err_code_;
		err_info.push_back(local_info);
	}

	//make error object with local info - further info gets added later from functions returning errors
	BError(string local_info)
	{
		err_info.push_back(local_info);
	}

	BError(const BError& copyThis)
	{
		*this = copyThis;
	}

	//get error code
	BERROR_ code(void) const { return err_code; }

	//get error info as a combined string
	string info(void) const 
	{ 
		string info_string;

		for (int idx = 0; idx < err_info.size(); idx++) {

			info_string += err_info[idx];
			if (idx != err_info.size() - 1) info_string += " > ";
		}

		return info_string;
	}

	//use to check if error code set
	operator bool() const { return err_code; }

	//set an error code
	BError& operator()(BERROR_ set_err_code) { err_code = set_err_code; return *this; }
	BError& operator()(BERROR_ set_err_code, string local_info) 
	{ 
		err_code = set_err_code; 
		err_info.push_back(local_info);
		return *this; 
	}

	//use to get error from a returning function : if error returned then copy info
	BError& operator=(const BError& rhs)
	{
		if (rhs.code()) {

			err_code = rhs.code();
			err_info.reserve(err_info.size() + rhs.err_info.size());
			err_info.insert(err_info.end(), rhs.err_info.begin(), rhs.err_info.end());
		}

		return *this;
	}

	//compare with a BERROR_ code
	bool operator==(BERROR_ code) const { return (err_code == code); }

	//reset error code and return object to receive new error, if any; e.g. if(error == BERROR_SOMETHING) error.reset() = function(...);
	BError& reset(void) { err_code = BERROR_NONE; return *this; }

	void clear(void) { err_code = BERROR_NONE; err_info.clear(); err_info.shrink_to_fit(); }
};

//------------------------------------------------------------------------------- ERROR HANDLER

template <typename Owner>
class ErrorHandler {

private:

	//the master object holding the Error Handler
	Owner * pOwner;

	//error messages
	vector<pair<string, ERRLEV_>> errors;

private:

	string get_error_string(BError error) { return errors[error.code()].first; }

	bool is_critical_error(BError error) { return (errors[error.code()].second == ERRLEV_CRIT); }

	bool is_silent_error(BError error) { return (errors[error.code()].second == ERRLEV_SILENT); }

public:

	ErrorHandler(Owner* pOwner_);

	void show_error(BError error, bool verbose = true)
	{
		if (!is_silent_error(error)) pOwner->show_error(error, get_error_string(error), verbose);
	}

	//full call, no parameters
	template <typename Object>
	bool call(BError(Object::*method)(), Object* pObject)
	{
		pOwner->create_restore();

		BError error = (pObject->*method)();

		if (error) {

			if (is_critical_error(error)) pOwner->restore_state();

			if (!is_silent_error(error)) pOwner->show_error(error, get_error_string(error), true);

			//error occured
			return true;
		}

		//no error
		return false;
	}

	//full call, with parameters
	template <typename Object, typename ... PType>
	bool call(BError(Object::*method)(PType ...), Object* pObject, PType ... params)
	{
		pOwner->create_restore();

		BError error = (pObject->*method)(params...);

		if (error) {

			if (is_critical_error(error)) pOwner->restore_state();

			if (!is_silent_error(error)) pOwner->show_error(error, get_error_string(error), true);

			//error occured
			return true;
		}

		//no error
		return false;
	}

	//quick call, no parameters
	template <typename Object>
	bool qcall(BError(Object::*method)(), Object* pObject)
	{
		BError error = (pObject->*method)();

		if (error) {

			if (!is_silent_error(error)) pOwner->show_error(error, get_error_string(error), true);

			//error occured
			return true;
		}

		//no error
		return false;
	}

	//quick call, with parameters
	template <typename Object, typename ... PType>
	bool qcall(BError(Object::*method)(PType ...), Object* pObject, PType ... params)
	{
		BError error = (pObject->*method)(params...);

		if (error) {

			if (!is_silent_error(error)) pOwner->show_error(error, get_error_string(error), true);

			//error occured
			return true;
		}

		//no error
		return false;
	}
};

//error messages and criticality defined here
template <typename Owner>
ErrorHandler<Owner>::ErrorHandler(Owner* pOwner_) :
	pOwner(pOwner_)
{
	errors.resize(BERROR_ENUMSIZE);

	errors[BERROR_NONE] = pair<string, ERRLEV_>( "", ERRLEV_SILENT );
	errors[BERROR_COMMAND_NOTRECOGNIZED] = pair<string, ERRLEV_>("Command not recognized.", ERRLEV_NCRIT);
	errors[BERROR_BUSY] = pair<string, ERRLEV_>("Busy... please wait.", ERRLEV_NCRIT);
	errors[BERROR_OPERATIONFAILED] = pair<string, ERRLEV_>("Operation failed.", ERRLEV_NCRIT);
	errors[BERROR_OUTOFMEMORY_CRIT] = pair<string, ERRLEV_>( "Out of memory.", ERRLEV_CRIT);
	errors[BERROR_OUTOFMEMORY_NCRIT] = pair<string, ERRLEV_>( "Out of memory.", ERRLEV_NCRIT);
	errors[BERROR_OUTOFGPUMEMORY_CRIT] = pair<string, ERRLEV_>( "Out of GPU memory.", ERRLEV_CRIT);
	errors[BERROR_OUTOFGPUMEMORY_NCRIT] = pair<string, ERRLEV_>("Out of GPU memory.", ERRLEV_NCRIT);
	errors[BERROR_GPUERROR_CRIT] = pair<string, ERRLEV_>("GPU error.", ERRLEV_CRIT);
	errors[BERROR_INCORRECTMODCONFIG] = pair<string, ERRLEV_>( "Incorrect modules configuration.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTCONFIG] = pair<string, ERRLEV_>( "Incorrect configuration.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTARRAYS] = pair<string, ERRLEV_>( "Incorrect arrays.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTVALUE] = pair<string, ERRLEV_>( "Incorrect value.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTSTRING] = pair<string, ERRLEV_>("Incorrect string.", ERRLEV_NCRIT);
	errors[BERROR_PARAMMISMATCH] = pair<string, ERRLEV_>( "Parameters mismatch.", ERRLEV_SILENT);
	errors[BERROR_PARAMMISMATCH_SHOW] = pair<string, ERRLEV_>("Parameters mismatch.", ERRLEV_NCRIT);
	errors[BERROR_PARAMOUTOFBOUNDS] = pair<string, ERRLEV_>( "Parameter out of bounds.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTNAME] = pair<string, ERRLEV_>( "Incorrect name.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTACTION] = pair<string, ERRLEV_>( "Incorrect action for current program state.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTACTION_SILENT] = pair<string, ERRLEV_>( "Incorrect action for current program state.", ERRLEV_SILENT);
	errors[BERROR_COULDNOTOPENFILE] = pair<string, ERRLEV_>("Could not open file.", ERRLEV_NCRIT);
	errors[BERROR_COULDNOTLOADFILE] = pair<string, ERRLEV_>( "Could not load file.", ERRLEV_NCRIT);
	errors[BERROR_COULDNOTLOADFILE_VERSIONMISMATCH] = pair<string, ERRLEV_>("Could not load file : version mismatch.", ERRLEV_NCRIT);
	errors[BERROR_COULDNOTLOADFILE_CRIT] = pair<string, ERRLEV_>( "Could not load file.", ERRLEV_CRIT);
	errors[BERROR_COULDNOTSAVEFILE] = pair<string, ERRLEV_>( "Could not save file.", ERRLEV_NCRIT);
	errors[BERROR_NOTAVAILABLE] = pair<string, ERRLEV_>( "Hardware not available or not enabled.", ERRLEV_NCRIT);
	errors[BERROR_COULDNOTCONNECT] = pair<string, ERRLEV_>("Couldn't connect.", ERRLEV_NCRIT);
	errors[BERROR_NOTMAGNETIC] = pair<string, ERRLEV_>("Focused mesh must be magnetic.", ERRLEV_NCRIT);
	errors[BERROR_NOTANTIFERROMAGNETIC] = pair<string, ERRLEV_>("Focused mesh must be antiferromagnetic (or use 2 sub-lattice model).", ERRLEV_NCRIT);
	errors[BERROR_NOTATOMISTIC] = pair<string, ERRLEV_>("Focused mesh must be atomistic.", ERRLEV_NCRIT);
	errors[BERROR_MESHNAMEINEXISTENT] = pair<string, ERRLEV_>("Mesh name doesn't exist.", ERRLEV_NCRIT);
	errors[BERROR_NOTCOMPUTED] = pair<string, ERRLEV_>("Not computed.", ERRLEV_NCRIT);
	errors[BERROR_SPINSOLVER_FIT] = pair<string, ERRLEV_>("Must be ferromagnetic mesh with transport module added and spin transport solver enabled. Must also have either Ts or Tsi computed.", ERRLEV_NCRIT);
	errors[BERROR_SPINSOLVER_FIT2] = pair<string, ERRLEV_>("Must be ferromagnetic mesh with transport module added and spin transport solver enabled.", ERRLEV_NCRIT);
	errors[BERROR_SPINSOLVER_FIT3] = pair<string, ERRLEV_>("Must be ferromagnetic mesh with transport module added and spin transport solver enabled. hm_mesh must be a metal mesh with transport module added.", ERRLEV_NCRIT);
	errors[BERROR_SPINSOLVER_FIT4] = pair<string, ERRLEV_>("Must give metal and ferromagnetic meshes in this order.", ERRLEV_NCRIT);
	errors[BERROR_NOTDEFINED] = pair<string, ERRLEV_>("Name not defined.", ERRLEV_NCRIT);
}