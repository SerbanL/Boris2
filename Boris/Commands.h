#pragma once

#include <string>

#include "ErrorHandler.h"

//Simulation commands enum
enum CMD_
{
	CMD_NONE = 0, CMD_RUN, CMD_STOP, CMD_RESET, CMD_COMPUTEFIELDS, CMD_ISRUNNING,
	CMD_CENTER, CMD_ITERUPDATE, CMD_CLEARSCREEN, CMD_REFRESHSCREEN, CMD_UPDATESCREEN,
	CMD_ADDFMESH, CMD_SETFMESH, CMD_ADDAFMESH, CMD_SETAFMESH, CMD_ADDDIPOLEMESH, CMD_ADDMETALMESH, CMD_ADDINSULATORMESH, CMD_ADDDIAMAGNETICMESH,
	CMD_ADDAMESHCUBIC, CMD_SETAMESHCUBIC,
	CMD_DELMESH, CMD_RENAMEMESH, CMD_MESHFOCUS, CMD_MESHFOCUS2, CMD_MESH, CMD_MESHRECT, CMD_SCALEMESHRECTS, CMD_CELLSIZE, CMD_ECELLSIZE, CMD_TCELLSIZE, CMD_MCELLSIZE, CMD_SCELLSIZE, CMD_FMSCELLSIZE, CMD_ESCELLSIZE, CMD_ATOMDMCELLSIZE,
	CMD_DELRECT, CMD_ADDRECT, CMD_RESETMESH, CMD_LOADMASKFILE, CMD_INDIVIDUALMASKSHAPE,
	CMD_SETANGLE, CMD_INVERTMAG, CMD_MIRRORMAG, CMD_SETRECT, CMD_DWALL, CMD_VORTEX, CMD_SKYRMION, CMD_SKYRMIONBLOCH, CMD_RANDOM,
	CMD_PBC,
	CMD_SETFIELD, CMD_SETSTRESS,
	CMD_MODULES, CMD_ADDMODULE, CMD_DELMODULE,
	CMD_MULTICONV, CMD_2DMULTICONV, CMD_NCOMMONSTATUS, CMD_NCOMMON, CMD_EXCLUDEMULTICONVDEMAG,
	CMD_ODE, CMD_SETODE, CMD_SETODEEVAL, CMD_SETATOMODE, CMD_SETDT, CMD_ASTEPCTRL, CMD_EVALSPEEDUP,
	CMD_SHOWDATA,
	CMD_CHDIR, CMD_SAVEDATAFILE, CMD_SAVECOMMENT, CMD_SAVEIMAGEFILE, CMD_DATASAVEFLAG, CMD_IMAGESAVEFLAG,
	CMD_DATA, CMD_ADDDATA, CMD_SETDATA, CMD_DELDATA, CMD_EDITDATA, CMD_ADDPINNEDDATA, CMD_DELPINNEDDATA,
	CMD_STAGES, CMD_ADDSTAGE, CMD_SETSTAGE, CMD_DELSTAGE, CMD_EDITSTAGE, CMD_EDITSTAGEVALUE, CMD_EDITSTAGESTOP, CMD_EDITDATASAVE,
	CMD_PARAMS, CMD_SETPARAM, CMD_PARAMSTEMP, CMD_CLEARPARAMSTEMP, CMD_SETPARAMTEMPEQUATION, CMD_SETPARAMTEMPARRAY, CMD_COPYPARAMS,
	CMD_COPYMESHDATA,
	CMD_PARAMSVAR, CMD_SETDISPLAYEDPARAMSVAR, CMD_CLEARPARAMSVAR, CMD_CLEARPARAMVAR, CMD_SETPARAMVAR,
	CMD_SAVESIM, CMD_LOADSIM, CMD_DEFAULT,
	CMD_DISPLAY, CMD_DISPLAYDETAILLEVEL, CMD_DISPLAYBACKGROUND, CMD_VECREP, CMD_SAVEMESHIMAGE, CMD_MAKEVIDEO, CMD_IMAGECROPPING, CMD_DISPLAYTRANSPARENCY, CMD_DISPLAYTHRESHOLDS, CMD_DISPLAYTHRESHOLDTRIGGER,
	CMD_MOVINGMESH, CMD_CLEARMOVINGMESH, CMD_MOVINGMESHASYM, CMD_MOVINGMESHTHRESH, CMD_PREPAREMOVINGMESH, CMD_PREPAREMOVINGBLOCHMESH, CMD_PREPAREMOVINGNEELMESH, CMD_PREPAREMOVINGSKYRMIONMESH, CMD_COUPLETODIPOLES, CMD_EXCHANGECOUPLEDMESHES,
	CMD_ADDELECTRODE, CMD_DELELECTRODE, CMD_CLEARELECTRODES, CMD_ELECTRODES, CMD_SETDEFAULTELECTRODES, CMD_SETELECTRODERECT, CMD_SETELECTRODEPOTENTIAL, CMD_DESIGNATEGROUND, CMD_SETPOTENTIAL, CMD_SETCURRENT, 
	CMD_TSOLVERCONFIG, CMD_SSOLVERCONFIG, CMD_SETSORDAMPING, CMD_STATICTRANSPORTSOLVER,
	CMD_TEMPERATURE, CMD_SETHEATDT, CMD_AMBIENTTEMPERATURE, CMD_ROBINALPHA, CMD_INSULATINGSIDES, CMD_CURIETEMPERATURE, CMD_CURIETEMPERATUREMATERIAL, CMD_ATOMICMOMENT, CMD_TAU, CMD_TMODEL,
	CMD_STOCHASTIC, CMD_LINKSTOCHASTIC, CMD_SETDTSTOCH, CMD_LINKDTSTOCHASTIC,
	CMD_SETDTSPEEDUP, CMD_LINKDTSPEEDUP,
	CMD_CUDA, CMD_MEMORY, CMD_SELCUDADEV,
	CMD_OPENMANUAL,
	CMD_REFINEROUGHNESS, CMD_CLEARROUGHNESS, CMD_ROUGHENMESH, CMD_SURFROUGHENJAGGED,
	CMD_GENERATE2DGRAINS, CMD_GENERATE3DGRAINS,
	CMD_BENCHTIME,
	CMD_MATERIALSDATABASE, CMD_ADDMATERIAL, CMD_SETMATERIAL, CMD_ADDMDBENTRY, CMD_DELMDBENTRY, CMD_REFRESHMDB, CMD_REQMDBSYNC, CMD_UPDATEMDB,
	CMD_SHOWLENGHTS, CMD_SHOWMCELLS,
	CMD_LOADOVF2MESH, CMD_LOADOVF2MAG, CMD_SAVEOVF2MAG, CMD_SAVEOVF2PARAMVAR, CMD_SAVEOVF2, CMD_LOADOVF2DISP, CMD_LOADOVF2STRAIN,
	CMD_SCRIPTSERVER, CMD_CHECKUPDATES,
	CMD_EQUATIONCONSTANTS, CMD_CLEAREQUATIONCONSTANTS, CMD_DELEQUATIONCONSTANT,
	CMD_FLUSHERRORLOG, CMD_ERRORLOG,
	CMD_STARTUPUPDATECHECK, CMD_STARTUPSCRIPTSERVER,
	CMD_THREADS,
	CMD_SERVERPORT, CMD_SERVERSLEEPMS,
	CMD_SHOWTC, CMD_SHOWMS, CMD_SHOWA, CMD_SHOWK,
	CMD_SKYPOSDMUL,
	CMD_MCSERIAL, CMD_MCCONSTRAIN,

	CMD_DP_CLEARALL, CMD_DP_CLEAR, CMD_DP_SHOWSIZES, CMD_DP_GET, CMD_DP_SET, CMD_DP_LOAD, CMD_DP_SAVE, CMD_DP_SAVEAPPEND, CMD_DP_SAVEASROW, CMD_DP_SAVEAPPENDASROW, CMD_DP_NEWFILE,
	CMD_DP_GETPROFILE, CMD_DP_GETEXACTPROFILE, CMD_DP_GETPATH, CMD_GETVALUE, CMD_AVERAGEMESHRECT, CMD_DP_TOPOCHARGE, CMD_DP_COUNTSKYRMIONS, CMD_DP_HISTOGRAM, CMD_DP_HISTOGRAM2,
	CMD_DP_APPEND, CMD_DP_SEQUENCE, CMD_DP_RAREFY, CMD_DP_EXTRACT, CMD_DP_ERASE,
	CMD_DP_ADD, CMD_DP_SUB, CMD_DP_MUL, CMD_DP_DIV, CMD_DP_DOTPROD,
	CMD_DP_ADDDP, CMD_DP_SUBDP, CMD_DP_MULDP, CMD_DP_DIVDP, CMD_DP_DOTPRODDP,
	CMD_DP_MINMAX, CMD_DP_MEAN, CMD_DP_GETAMPLI,
	CMD_DP_LINREG, CMD_DP_COERCIVITY, CMD_DP_REMANENCE, CMD_DP_COMPLETEHYSTERLOOP,
	CMD_DP_DUMPTDEP,
	CMD_DP_FITLORENTZ, CMD_DP_FITLORENTZ2, CMD_DP_FITSKYRMION, CMD_DP_FITSTT, CMD_DP_FITSOT, CMD_DP_FITSOTSTT, CMD_DP_CALCSOT, CMD_DP_FITADIABATIC, CMD_DP_FITNONADIABATIC,
	CMD_DP_CALCEXCHANGE, CMD_DP_CALCTOPOCHARGEDENSITY,
	CMD_DP_REPLACEREPEATS, CMD_DP_REMOVEOFFSET, CMD_CARTESIANTOPOLAR,
	CMD_DP_SMOOTH, CMD_DP_MONOTONIC,
	CMD_DP_CROSSINGSHISTOGRAM, CMD_DP_CROSSINGSFREQUENCY, CMD_DP_PEAKSFREQUENCY,

	//used for testing various things
	CMD_TEST
};

using namespace std;

//Describes the structure of a command, and helper methods for processing commands
struct CommandSpecifier {

public:

	//a unique command identifier
	CMD_ cmdId;

	//show usage and command description using formatted strings (for display in console).
	//Also show description of returned values for script clients.
	//If unit is specified, allow input using physical quantity units rather then e notation
	string usage, descr, return_descr;
	string unit;

	//limits for parameters specified here; the pair gives the lower and upper limit (inclusive). To ommit a limit set it as an Any(). 
	//Otherwise must set the Any with the expected parameter type, in the order expected!
	vector<pair<Any, Any>> limits;

private:

	template <typename Type, typename ... PType>
	BError check_limits(int index, Type param)
	{
		BError error;

		if (!limits.size()) return error;

		auto make_error_info = [&]() -> string {

			string limits_string = "Parameter " + ToString(index) + ", ";

			if (limits[index].first.IsNull()) limits_string += "Min : N/A to ";
			else {

				Type min = limits[index].first;

				if (!std::is_integral<Type>::value)
					limits_string += "Min : " + ToString(min, unit) + " to ";
				else limits_string += "Min : " + ToString(min) + " to ";
			}

			if (limits[index].second.IsNull()) limits_string += "Max : N/A";
			else {

				Type max = limits[index].second;

				if (!std::is_integral<Type>::value)
					limits_string += "Max : " + ToString(max, unit);
				else limits_string += "Max : " + ToString(max);
			}

			return limits_string;
		};

		Type first = limits[index].first;
		Type second = limits[index].second;

		//last parameter : return false on failed bounds check - NOTE, must use comparisons in form below so it works with VAL2 and VAL3 also
		if ((limits[index].first.IsNull() || param >= first) &&
			(limits[index].second.IsNull() || param <= second)) return error;

		return error(BERROR_PARAMOUTOFBOUNDS, make_error_info());
	}

	template <typename Type, typename ... PType>
	BError check_limits(int index, Type param, PType ... params)
	{
		BError error;

		if (!limits.size()) return error;

		auto make_error_info = [&]() -> string {

			string limits_string = "Parameter " + ToString(index) + ", ";

			if (limits[index].first.IsNull()) limits_string += "Min : N/A to ";
			else {

				Type min = limits[index].first;

				if(!std::is_integral<Type>::value)
					limits_string += "Min : " + ToString(min, unit) + " to ";
				else limits_string += "Min : " + ToString(min) + " to ";
			}

			if (limits[index].second.IsNull()) limits_string += "Max : N/A";
			else {

				Type max = limits[index].second;

				if (!std::is_integral<Type>::value)
					limits_string += "Max : " + ToString(max, unit);
				else limits_string += "Max : " + ToString(max);
			}

			return limits_string;
		};

		Type first = limits[index].first;
		Type second = limits[index].second;

		//return false on first parameter to fail bounds check - NOTE, must use comparisons in form below so it works with VAL2 and VAL3 also
		if (!((limits[index].first.IsNull() || param >= first) &&
			(limits[index].second.IsNull() || param <= second))) return error(BERROR_PARAMOUTOFBOUNDS, make_error_info());

		//if we have more parameters than bounds specified, use the last limits - this is useful if we have a command with an unspecified number of similar parameters (e.g. list of dp array indexes which all have the same limits)
		if (index < limits.size() - 1) ++index;

		//continue checking next parameter
		error = check_limits(index, params...);

		return error;
	}

public:

	CommandSpecifier(void) { cmdId = CMD_NONE; }
	CommandSpecifier(CMD_ cmdId) { this->cmdId = cmdId; }

	//get parameters contained as strings in components, into values passed in the parameter pack
	template <typename ... Type> 
	BError GetParameters(vector<string>& components, Type& ... values) 
	{
		BError error;

		//get parameters as a vector of Any objects. Construct  Any objects using pointers, so changes to vector entries result in changes to passed parameters - since these are in turn passed through reference changes are propagated to caller.
		vector<Any> a_values = make_vector(Any(&values)...);

		//component index
		int cIdx = 0;

		//go over all passed parameters
		for (int valIdx = 0; valIdx < (int)a_values.size(); valIdx++) {

			int cIdx_ = cIdx;

			//keep converting as long as there are enough components left - if not enough components then something is not right : return false
			if (cIdx < (int)components.size())
				cIdx += a_values[valIdx].convert_string(subvec(components, cIdx), unit);
			else return error(BERROR_PARAMMISMATCH);

			//couldn't convert : return false
			if (cIdx == cIdx_) return error(BERROR_PARAMMISMATCH);
		}
		
		//if last value is a string, add any remaining components to it using space separators.
		if (a_values.back().is_type(btype_info<std::string>())) {

			string text = a_values.back();

			if (cIdx < (int)components.size()) {

				text += " " + combine(subvec(components, cIdx), " ");
				a_values.back() = text;
			}
		}

		//got all values, now check limits
		error = check_limits(0, values...);
		
		return error;
	}

	template <typename ... Type>
	vector<string> PrepareReturnParameters(Type... values)
	{
		vector<string> outVector;

		for (string value_text : { ToString(values)... })
			JoinToVector(outVector, split(value_text, ", ", "; "));

		return outVector;
	}
};