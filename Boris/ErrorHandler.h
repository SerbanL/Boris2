#pragma once

#include <vector>
#include <string>



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
	BERROR_CUDAVERSIONMISMATCH_NCRIT,		//Device CUDA version doesn't match program CUDA architecture
	BERROR_INCORRECTMODCONFIG,				//Modules not configured correctly. Cannot start simulation
	BERROR_INCORRECTCONFIG,					//Incorrect configuration (e.g. eval method doesn't match ode allowed evals)
	BERROR_INCORRECTARRAYS,					//Incorrect arrays used
	BERROR_INCORRECTVALUE,					//Incorrect value used
	BERROR_INCORRECTSTRING,					//Incorrect std::string used
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
	BERROR_NOHEAT,							//Heat computation not enabled
	BERROR_NOTRANSPORT,						//Transport computation not enabled
	BERROR_MESHNAMEINEXISTENT,				//Mesh name doesn't exist.
	BERROR_NOTCOMPUTED,						//Not computed.
	BERROR_SPINSOLVER_FIT,					//Must be ferromagnetic mesh with transport module added and spin transport solver enabled.Must also have either Ts or Tsi computed.
	BERROR_SPINSOLVER_FIT2,					//Must be ferromagnetic mesh with transport module added and spin transport solver enabled.
	BERROR_SPINSOLVER_FIT3,					//Must be ferromagnetic mesh with transport module added and spin transport solver enabled.hm_mesh must be a metal mesh with transport module added.
	BERROR_SPINSOLVER_FIT4,					//Must give metal and ferromagnetic meshes in this order.
	BERROR_NOTDEFINED,						//Name not defined
	BERROR_ENUMSIZE
};

//BWARNING_ enum with all the possible error identifiers : warnings are simply displayed messages, but not errors, so everything proceeds as normal
enum BWARNING_
{
	BWARNING_NONE = 0,
	BWARNING_INCORRECTCELLSIZE,				//cellsize set is incorrect
	BWARNING_NOGPUINITIALIZATION,			//could not initialize on GPU ... initialized on CPU instead
	BWARNING_ENUMSIZE
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
//2. void restore_state(BERROR_ berror, std::string error_message); -> restore state after a critical error occurs
//3. void show_error(BERROR_ berror, std::string error_message); -> show the error
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
//	std::string error_info = error.info(); 
//	BERROR_ error_code = error.code(); 
//}

//------------------------------------------------------------------------------- ERROR CONTAINER

class BError {

	//the error code
	BERROR_ err_code = BERROR_NONE;

	BWARNING_ warn_code = BWARNING_NONE;

	//accumulated info for the error or warning code containing the whole chain of functions returning errors or warnings, each with their own info
	std::vector<std::string> accum_info;

public:

	//None
	BError(void) {}

	//Errors
	BError(BERROR_ err_code_) 
	{
		err_code = err_code_;
	}

	BError(BERROR_ err_code_, std::string local_info)
	{
		err_code = err_code_;
		accum_info.push_back(local_info);
	}

	//Warnings
	BError(BWARNING_ warn_code_)
	{
		warn_code = warn_code_;
	}

	BError(BWARNING_ warn_code_, std::string local_info)
	{
		warn_code = warn_code_;
		accum_info.push_back(local_info);
	}

	//make error object with local info - further info gets added later from functions returning errors
	BError(std::string local_info)
	{
		accum_info.push_back(local_info);
	}

	//copy constructor
	BError(const BError& copyThis)
	{
		*this = copyThis;
	}

	//use to get error from a returning function : if error or warning returned then copy info
	BError& operator=(const BError& rhs)
	{
		if (rhs.code() || rhs.wcode()) {

			err_code = rhs.code();
			warn_code = rhs.wcode();
			accum_info.reserve(accum_info.size() + rhs.accum_info.size());
			accum_info.insert(accum_info.end(), rhs.accum_info.begin(), rhs.accum_info.end());
		}

		return *this;
	}

	//get error code
	BERROR_ code(void) const { return err_code; }

	//get warning code
	BWARNING_ wcode(void) const { return warn_code; }

	//get error info as a combined std::string
	std::string info(void) const 
	{ 
		std::string info_string;

		for (int idx = 0; idx < accum_info.size(); idx++) {

			info_string += accum_info[idx];
			if (idx != accum_info.size() - 1) info_string += " > ";
		}

		return info_string;
	}

	//use to check if error code set
	operator bool() const { return err_code; }

	//check if a warning is set
	bool warning_set(void) const { return warn_code; }

	//set an error code
	BError& operator()(BERROR_ set_err_code) { err_code = set_err_code; return *this; }
	BError& operator()(BERROR_ set_err_code, std::string local_info) 
	{ 
		err_code = set_err_code; 
		accum_info.push_back(local_info);
		return *this; 
	}

	//set a warning code
	BError& operator()(BWARNING_ set_warn_code) { warn_code = set_warn_code; return *this; }
	BError& operator()(BWARNING_ set_warn_code, std::string local_info)
	{
		warn_code = set_warn_code;
		accum_info.push_back(local_info);
		return *this;
	}

	//compare with a BERROR_ code
	bool operator==(BERROR_ code) const { return (err_code == code); }
	bool operator!=(BERROR_ code) const { return (err_code != code); }

	//compare with a BWARNING_ code
	bool operator==(BWARNING_ code) const { return (warn_code == code); }
	bool operator!=(BWARNING_ code) const { return (warn_code != code); }

	//reset error code and return object to receive new error, if any; e.g. if(error == BERROR_SOMETHING) error.reset() = function(...);
	BError& reset(void) 
	{ 
		err_code = BERROR_NONE;
		warn_code = BWARNING_NONE;
		return *this; 
	}

	void clear(void) 
	{ 
		err_code = BERROR_NONE;
		warn_code = BWARNING_NONE;
		accum_info.clear();
		accum_info.shrink_to_fit();
	}
};

//------------------------------------------------------------------------------- ERROR HANDLER

template <typename Owner>
class ErrorHandler {

private:

	//the master object holding the Error Handler
	Owner * pOwner;

	//error messages
	std::vector<std::pair<std::string, ERRLEV_>> errors;

	//warnings messages
	std::vector<std::string> warnings;

private:

	//std::string get_error_string(BError error) const { return errors[error.code()].first; }

	bool is_critical_error(BError error) const { return (errors[error.code()].second == ERRLEV_CRIT); }

	bool is_silent_error(BError error) const { return (errors[error.code()].second == ERRLEV_SILENT); }

	bool is_warning(BError error) const { return (error.wcode() != BWARNING_NONE); }

public:

	ErrorHandler(Owner* pOwner_);

	void show_error(BError error, bool verbose = true) const
	{
		if (!is_silent_error(error) || is_warning(error)) pOwner->show_error(error, get_error_text(error), verbose);
	}

	std::string get_error_text(BError error) const
	{
		if (is_warning(error)) return "WARNING : " + warnings[error.wcode()] + " Info : " + error.info();
		else return "ERROR : " + errors[error.code()].first + " Info : " + error.info();
	}

	//full call, no parameters
	template <typename Object>
	bool call(BError& error, BError(Object::*method)(), Object* pObject)
	{
		pOwner->create_restore();

		error = (pObject->*method)();

		if (error) {

			if (is_critical_error(error)) pOwner->restore_state();

			//error occured
			return true;
		}

		//no error
		return false;
	}

	//full call, with parameters
	template <typename Object, typename ... PType>
	bool call(BError& error, BError(Object::*method)(PType ...), Object* pObject, PType ... params)
	{
		pOwner->create_restore();

		error = (pObject->*method)(params...);

		if (error) {

			if (is_critical_error(error)) pOwner->restore_state();

			//error occured
			return true;
		}

		//no error
		return false;
	}

	//quick call, no parameters
	template <typename Object>
	bool qcall(BError& error, BError(Object::*method)(), Object* pObject)
	{
		error = (pObject->*method)();

		if (error) {

			//error occured
			return true;
		}

		//no error
		return false;
	}

	//quick call, with parameters
	template <typename Object, typename ... PType>
	bool qcall(BError& error, BError(Object::*method)(PType ...), Object* pObject, PType ... params)
	{
		error = (pObject->*method)(params...);

		if (error) {

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
	/////////////////////////////////////////////////////////////////////////////////////

	errors.resize(BERROR_ENUMSIZE);

	errors[BERROR_NONE] = std::pair<std::string, ERRLEV_>( "", ERRLEV_SILENT );
	errors[BERROR_COMMAND_NOTRECOGNIZED] = std::pair<std::string, ERRLEV_>("Command not recognized.", ERRLEV_NCRIT);
	errors[BERROR_BUSY] = std::pair<std::string, ERRLEV_>("Busy... please wait.", ERRLEV_NCRIT);
	errors[BERROR_OPERATIONFAILED] = std::pair<std::string, ERRLEV_>("Operation failed.", ERRLEV_NCRIT);
	errors[BERROR_OUTOFMEMORY_CRIT] = std::pair<std::string, ERRLEV_>( "Out of memory.", ERRLEV_CRIT);
	errors[BERROR_OUTOFMEMORY_NCRIT] = std::pair<std::string, ERRLEV_>( "Out of memory.", ERRLEV_NCRIT);
	errors[BERROR_OUTOFGPUMEMORY_CRIT] = std::pair<std::string, ERRLEV_>( "Out of GPU memory.", ERRLEV_CRIT);
	errors[BERROR_OUTOFGPUMEMORY_NCRIT] = std::pair<std::string, ERRLEV_>("Out of GPU memory.", ERRLEV_NCRIT);
	errors[BERROR_GPUERROR_CRIT] = std::pair<std::string, ERRLEV_>("GPU error.", ERRLEV_CRIT);
	errors[BERROR_CUDAVERSIONMISMATCH_NCRIT] = std::pair<std::string, ERRLEV_>("Device CUDA version mismatch. Use matching Boris CUDA version.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTMODCONFIG] = std::pair<std::string, ERRLEV_>( "Incorrect modules configuration.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTCONFIG] = std::pair<std::string, ERRLEV_>( "Incorrect configuration.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTARRAYS] = std::pair<std::string, ERRLEV_>( "Incorrect arrays.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTVALUE] = std::pair<std::string, ERRLEV_>( "Incorrect value.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTSTRING] = std::pair<std::string, ERRLEV_>("Incorrect std::string.", ERRLEV_NCRIT);
	errors[BERROR_PARAMMISMATCH] = std::pair<std::string, ERRLEV_>( "Parameters mismatch.", ERRLEV_SILENT);
	errors[BERROR_PARAMMISMATCH_SHOW] = std::pair<std::string, ERRLEV_>("Parameters mismatch.", ERRLEV_NCRIT);
	errors[BERROR_PARAMOUTOFBOUNDS] = std::pair<std::string, ERRLEV_>( "Parameter out of bounds.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTNAME] = std::pair<std::string, ERRLEV_>( "Incorrect name.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTACTION] = std::pair<std::string, ERRLEV_>( "Incorrect action for current program state.", ERRLEV_NCRIT);
	errors[BERROR_INCORRECTACTION_SILENT] = std::pair<std::string, ERRLEV_>( "Incorrect action for current program state.", ERRLEV_SILENT);
	errors[BERROR_COULDNOTOPENFILE] = std::pair<std::string, ERRLEV_>("Could not open file.", ERRLEV_NCRIT);
	errors[BERROR_COULDNOTLOADFILE] = std::pair<std::string, ERRLEV_>( "Could not load file.", ERRLEV_NCRIT);
	errors[BERROR_COULDNOTLOADFILE_VERSIONMISMATCH] = std::pair<std::string, ERRLEV_>("Could not load file : version mismatch.", ERRLEV_NCRIT);
	errors[BERROR_COULDNOTLOADFILE_CRIT] = std::pair<std::string, ERRLEV_>( "Could not load file.", ERRLEV_CRIT);
	errors[BERROR_COULDNOTSAVEFILE] = std::pair<std::string, ERRLEV_>( "Could not save file.", ERRLEV_NCRIT);
	errors[BERROR_NOTAVAILABLE] = std::pair<std::string, ERRLEV_>( "Hardware not available or not enabled.", ERRLEV_NCRIT);
	errors[BERROR_COULDNOTCONNECT] = std::pair<std::string, ERRLEV_>("Couldn't connect.", ERRLEV_NCRIT);
	errors[BERROR_NOTMAGNETIC] = std::pair<std::string, ERRLEV_>("Focused mesh must be magnetic.", ERRLEV_NCRIT);
	errors[BERROR_NOTANTIFERROMAGNETIC] = std::pair<std::string, ERRLEV_>("Focused mesh must be antiferromagnetic (or use 2 sub-lattice model).", ERRLEV_NCRIT);
	errors[BERROR_NOTATOMISTIC] = std::pair<std::string, ERRLEV_>("Focused mesh must be atomistic.", ERRLEV_NCRIT);
	errors[BERROR_NOHEAT] = std::pair<std::string, ERRLEV_>("Heat module not enabled.", ERRLEV_NCRIT);
	errors[BERROR_NOTRANSPORT] = std::pair<std::string, ERRLEV_>("Transport module not enabled.", ERRLEV_NCRIT);
	errors[BERROR_MESHNAMEINEXISTENT] = std::pair<std::string, ERRLEV_>("Mesh name doesn't exist.", ERRLEV_NCRIT);
	errors[BERROR_NOTCOMPUTED] = std::pair<std::string, ERRLEV_>("Not computed.", ERRLEV_NCRIT);
	errors[BERROR_SPINSOLVER_FIT] = std::pair<std::string, ERRLEV_>("Must be ferromagnetic mesh with transport module added and spin transport solver enabled. Must also have either Ts or Tsi computed.", ERRLEV_NCRIT);
	errors[BERROR_SPINSOLVER_FIT2] = std::pair<std::string, ERRLEV_>("Must be ferromagnetic mesh with transport module added and spin transport solver enabled.", ERRLEV_NCRIT);
	errors[BERROR_SPINSOLVER_FIT3] = std::pair<std::string, ERRLEV_>("Must be ferromagnetic mesh with transport module added and spin transport solver enabled. hm_mesh must be a metal mesh with transport module added.", ERRLEV_NCRIT);
	errors[BERROR_SPINSOLVER_FIT4] = std::pair<std::string, ERRLEV_>("Must give metal and ferromagnetic meshes in this order.", ERRLEV_NCRIT);
	errors[BERROR_NOTDEFINED] = std::pair<std::string, ERRLEV_>("Name not defined.", ERRLEV_NCRIT);

	/////////////////////////////////////////////////////////////////////////////////////

	warnings.resize(BWARNING_ENUMSIZE);

	warnings[BWARNING_INCORRECTCELLSIZE] = std::string("Working with incorrect cellsize.");
	warnings[BWARNING_NOGPUINITIALIZATION] = std::string("Could not initialize on GPU. Initialized on CPU instead.");

	/////////////////////////////////////////////////////////////////////////////////////
}