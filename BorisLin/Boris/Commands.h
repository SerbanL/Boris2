#pragma once

#include <string>

#include "ErrorHandler.h"

//Simulation commands enum
enum CMD_
{
	CMD_NONE = 0, 
	
	//-------------------------------------------SIMULATION AND STATE CONTROL-------------------------------------------

	CMD_RUN, CMD_STOP, CMD_RESET, CMD_COMPUTEFIELDS, CMD_ISRUNNING,
	CMD_SAVESIM, CMD_LOADSIM, CMD_DEFAULT,

	CMD_NEXTSTAGE, CMD_SETSCHEDULESTAGE, CMD_RUNSTAGE,

	//-------------------------------------------SCRIPTING-------------------------------------------

	CMD_RUNSCRIPT, CMD_RUNSCRIPTNEWINSTANCE, CMD_STOPSCRIPT, CMD_SCRIPTPRINT,

	//-------------------------------------------DISPLAY CONTROL-------------------------------------------

	CMD_CENTER, CMD_ITERUPDATE, CMD_CLEARSCREEN, CMD_REFRESHSCREEN, CMD_UPDATESCREEN,
	CMD_MESHFOCUS, CMD_MESHFOCUS2,

	CMD_DISPLAY, CMD_DISPLAYDETAILLEVEL, CMD_DISPLAYRENDERTHRESH, CMD_DISPLAYBACKGROUND, CMD_VECREP, CMD_IMAGECROPPING, CMD_DISPLAYTRANSPARENCY, CMD_DISPLAYTHRESHOLDS, CMD_DISPLAYTHRESHOLDTRIGGER, 
	CMD_DISPLAYMODULE,
	
	CMD_ROTCAMABOUTORIGIN, CMD_ROTCAMABOUTAXIS, CMD_ADJUSTCAMDISTANCE, CMD_SHIFTCAMORIGIN,

	//-------------------------------------------MESHES-------------------------------------------
	
	//Mesh add/delete

	CMD_ADDFMESH, CMD_SETFMESH, CMD_ADDAFMESH, CMD_SETAFMESH, CMD_ADDDIPOLEMESH, CMD_SETDIPOLEMESH,
	CMD_ADDMETALMESH, CMD_SETMETALMESH, CMD_ADDINSULATORMESH, CMD_SETINSULATORMESH,
	CMD_ADDAMESHCUBIC, CMD_SETAMESHCUBIC,
	CMD_DELMESH, 
	
	//Mesh adjustment

	CMD_MESH,
	CMD_RENAMEMESH, 
	CMD_MESHRECT, CMD_SCALEMESHRECTS, 
	CMD_SHIFTDIPOLE, CMD_DIPOLEVELOCITY,

	//Celsizes

	CMD_CELLSIZE, CMD_ECELLSIZE, CMD_TCELLSIZE, CMD_MCELLSIZE, CMD_SCELLSIZE, CMD_FMSCELLSIZE, CMD_ESCELLSIZE, CMD_ATOMDMCELLSIZE,
	
	//Mesh shaping

	CMD_DELRECT, CMD_ADDRECT, CMD_RESETMESH, CMD_LOADMASKFILE, CMD_INDIVIDUALMASKSHAPE,
	CMD_GENERATE2DGRAINS, CMD_GENERATE3DGRAINS,

	//Shape objects

	CMD_SHAPEMOD_ROT, CMD_SHAPEMOD_REP, CMD_SHAPEMOD_DISP, CMD_SHAPEMOD_METHOD,
	CMD_SHAPE_DISK, CMD_SHAPE_RECT, CMD_SHAPE_TRIANGLE, CMD_SHAPE_ELLIPSOID, CMD_SHAPE_PYRAMID, CMD_SHAPE_TETRAHEDRON, CMD_SHAPE_CONE, CMD_SHAPE_TORUS, CMD_SHAPE_SET, CMD_SHAPE_GET,
	CMD_SETSHAPEANGLE, CMD_SHAPE_SETPARAM,

	//Magnetization texturing

	CMD_SETANGLE, CMD_SETOBJECTANGLE, CMD_INVERTMAG, CMD_MIRRORMAG, CMD_SETRECT,
	CMD_DWALL, CMD_VORTEX, CMD_SKYRMION, CMD_SKYRMIONBLOCH, CMD_RANDOM, CMD_RANDOMXY, CMD_FLOWER, CMD_ONION, CMD_CROSSTIE,

	//Important inputs and settings

	CMD_PBC,
	CMD_SETFIELD, CMD_SETSTRESS,

	CMD_COUPLETODIPOLES, CMD_EXCHANGECOUPLEDMESHES,

	//Various algorithms settings

	CMD_MOVINGMESH, CMD_CLEARMOVINGMESH, CMD_MOVINGMESHASYM, CMD_MOVINGMESHTHRESH, CMD_PREPAREMOVINGMESH, CMD_PREPAREMOVINGBLOCHMESH, CMD_PREPAREMOVINGNEELMESH, CMD_PREPAREMOVINGSKYRMIONMESH,
	CMD_SKYPOSDMUL, CMD_DWPOS_COMPONENT,

	//Mesh material parameters

	CMD_COPYMESHDATA,
	CMD_PARAMS, CMD_SETPARAM, CMD_PARAMSTEMP, CMD_CLEARPARAMSTEMP, CMD_SETPARAMTEMPEQUATION, CMD_SETPARAMTEMPARRAY, CMD_COPYPARAMS, CMD_SETKTENS,
	CMD_PARAMSVAR, CMD_SETDISPLAYEDPARAMSVAR, CMD_CLEARPARAMSVAR, CMD_SETPARAMVAR,

	//Mesh properties

	CMD_GETMESHTYPE,

	//-------------------------------------------MODULES-------------------------------------------

	//Module control

	CMD_MODULES,
	CMD_ADDMODULE, CMD_DELMODULE,

	//Demag computation control

	CMD_MULTICONV, CMD_2DMULTICONV, CMD_NCOMMONSTATUS, CMD_NCOMMON, CMD_EXCLUDEMULTICONVDEMAG,
	CMD_GPUKERNELS,

	//-------------------------------------------ODE-------------------------------------------

	//General

	CMD_ODE, CMD_SETODE, CMD_SETODEEVAL, CMD_SETATOMODE, CMD_SETDT, CMD_ASTEPCTRL,

	CMD_EVALSPEEDUP, CMD_SETDTSPEEDUP, CMD_LINKDTSPEEDUP,

	//Stochasticity

	CMD_STOCHASTIC, CMD_LINKSTOCHASTIC, CMD_SETDTSTOCH, CMD_LINKDTSTOCHASTIC,

	//-------------------------------------------DATA-------------------------------------------

	//Options and settings

	CMD_CHDIR, CMD_SAVEDATAFILE, CMD_SAVECOMMENT, CMD_SAVEIMAGEFILE, CMD_DATASAVEFLAG, CMD_IMAGESAVEFLAG, CMD_DATAPRECISION, CMD_DISKBUFFERLINES,

	//Output scheduled data

	CMD_SHOWDATA,
	CMD_DATA, CMD_ADDDATA, CMD_SETDATA, CMD_DELDATA, CMD_EDITDATA, CMD_ADDPINNEDDATA, CMD_DELPINNEDDATA,

	CMD_SAVEDATA,

	//Simulation schedule

	CMD_STAGES, CMD_ADDSTAGE, CMD_SETSTAGE, CMD_DELSTAGE, CMD_EDITSTAGE, CMD_EDITSTAGEVALUE, CMD_SETSTAGEVALUE, CMD_EDITSTAGESTOP, CMD_EDITDATASAVE,

	//Other output

	CMD_SAVEMESHIMAGE, CMD_SAVEIMAGE, CMD_MAKEVIDEO,

	//Command buffering (advanced data extraction)

	CMD_RUNCOMMBUFFER, CMD_CLEARCOMMBUFFER, CMD_BUFFERCOMMAND,

	//-------------------------------------------TRANSPORT SOLVER-------------------------------------------

	CMD_ADDELECTRODE, CMD_DELELECTRODE, CMD_CLEARELECTRODES, CMD_ELECTRODES, CMD_SETDEFAULTELECTRODES, CMD_SETELECTRODERECT, CMD_SETELECTRODEPOTENTIAL, CMD_DESIGNATEGROUND,
	CMD_SETPOTENTIAL, CMD_SETCURRENT, CMD_SETCURRENTDENSITY, CMD_OPENPOTENTIAL,
	CMD_SSOLVERCONFIG, CMD_SETSORDAMPING, CMD_STATICTRANSPORTSOLVER, CMD_DISABLETRANSPORTSOLVER,
	CMD_TMRTYPE,
	CMD_RAPBIAS_EQUATION, CMD_RAAPBIAS_EQUATION,
	CMD_TAMREQUATION,

	//-------------------------------------------HEAT SOLVER-------------------------------------------

	CMD_TSOLVERCONFIG,
	CMD_TEMPERATURE, CMD_SETHEATDT, CMD_AMBIENTTEMPERATURE, CMD_ROBINALPHA, CMD_INSULATINGSIDES, CMD_TMODEL,

	//Mesh related

	CMD_CURIETEMPERATURE, CMD_CURIETEMPERATUREMATERIAL, CMD_ATOMICMOMENT, CMD_TAU,

	//-------------------------------------------ELASTODYNAMICS SOLVER-------------------------------------------

	CMD_RESETELSOLVER, CMD_STRAINEQUATION, CMD_SHEARSTRAINEQUATION, CMD_CLEARSTRAINEQUATIONS, CMD_SETELDT, CMD_LINKDTELASTIC,
	CMD_SURFACEFIX, CMD_SURFACESTRESS, CMD_DELSURFACEFIX, CMD_DELSURFACESTRESS, CMD_EDITSURFACEFIX, CMD_EDITSURFACESTRESS, CMD_EDITSURFACESTRESSEQUATION,

	//-------------------------------------------CUDA-------------------------------------------

	CMD_CUDA, CMD_MEMORY, CMD_SELCUDADEV,

	//-------------------------------------------OTHERS-------------------------------------------

	CMD_OPENMANUAL,
	CMD_BENCHTIME,
	CMD_SHOWLENGHTS, CMD_SHOWMCELLS,
	CMD_SCRIPTSERVER, CMD_CHECKUPDATES,
	CMD_FLUSHERRORLOG, CMD_ERRORLOG,
	CMD_STARTUPUPDATECHECK, CMD_STARTUPSCRIPTSERVER,
	CMD_THREADS,
	CMD_SERVERPORT, CMD_SERVERPWD, CMD_SERVERSLEEPMS,
	CMD_NEWINSTANCE,

	//-------------------------------------------VERSION UPDATE-------------------------------------------

	CMD_VERSIONUPDATE,

	//-------------------------------------------ROUGHNESS-------------------------------------------

	CMD_REFINEROUGHNESS, CMD_CLEARROUGHNESS, CMD_ROUGHENMESH, CMD_SURFROUGHENJAGGED,

	//-------------------------------------------MATERIALS DATABASE-------------------------------------------
	
	CMD_MATERIALSDATABASE, CMD_ADDMATERIAL, CMD_SETMATERIAL, CMD_ADDMDBENTRY, CMD_DELMDBENTRY, CMD_REFRESHMDB, CMD_REQMDBSYNC, CMD_UPDATEMDB,

	//-------------------------------------------OVF2 COMMANDS-------------------------------------------
	
	CMD_LOADOVF2MAG, CMD_SAVEOVF2MAG, CMD_LOADOVF2FIELD, CMD_SAVEOVF2PARAMVAR, CMD_SAVEOVF2, CMD_LOADOVF2DISP, CMD_LOADOVF2STRAIN, CMD_LOADOVF2TEMP, CMD_LOADOVF2CURR,
	CMD_CLEARGLOBALFIELD, CMD_SHIFTGLOBALFIELD,

	//-------------------------------------------TEXT EQUATIONS-------------------------------------------

	CMD_EQUATIONCONSTANTS, CMD_CLEAREQUATIONCONSTANTS, CMD_DELEQUATIONCONSTANT,

	//-------------------------------------------INFO-------------------------------------------
	
	CMD_SHOWTC, CMD_SHOWMS, CMD_SHOWA, CMD_SHOWK,

	//-------------------------------------------MONTE CARLO-------------------------------------------

	CMD_MCSERIAL, CMD_MCDISABLE, CMD_MCCONSTRAIN, CMD_MCCOMPUTEFIELDS, CMD_MCCONEANGLELIMITS,

	//-------------------------------------------PRNG-------------------------------------------
	
	CMD_PRNGSEED,

	//-------------------------------------------DP COMMANDS-------------------------------------------
	
	//Basic control and input/output

	CMD_DP_CLEARALL, CMD_DP_CLEAR, CMD_DP_SHOWSIZES, CMD_DP_GET, CMD_DP_SET, CMD_DP_LOAD, CMD_DP_SAVE, CMD_DP_SAVEAPPEND, CMD_DP_SAVEASROW, CMD_DP_SAVEAPPENDASROW, CMD_DP_NEWFILE,

	//Profile and average extraction

	CMD_DP_GETPROFILE, CMD_DP_GETEXACTPROFILE, CMD_DP_GETAVERAGEDPROFILE, 
	CMD_AVERAGEMESHRECT, 
	CMD_GETVALUE,
		
	//Special data extraction

	CMD_DP_TOPOCHARGE, CMD_DP_COUNTSKYRMIONS, CMD_DP_CALCTOPOCHARGEDENSITY,
	CMD_DP_HISTOGRAM, CMD_DP_THAVHISTOGRAM, CMD_DP_ANGHISTOGRAM, CMD_DP_THAVANGHISTOGRAM, CMD_DP_HISTOGRAM2,

	//Simple dp modifiers

	CMD_DP_APPEND, CMD_DP_SEQUENCE, CMD_DP_RAREFY, CMD_DP_EXTRACT, CMD_DP_ERASE,
	CMD_DP_REPLACEREPEATS, CMD_DP_REMOVEOFFSET, CMD_CARTESIANTOPOLAR,

	//Operations

	CMD_DP_ADD, CMD_DP_SUB, CMD_DP_MUL, CMD_DP_DIV, CMD_DP_POW, CMD_DP_DOTPROD,
	CMD_DP_ADDDP, CMD_DP_SUBDP, CMD_DP_MULDP, CMD_DP_DIVDP, CMD_DP_DOTPRODDP,
	CMD_DP_MINMAX, CMD_DP_MEAN, CMD_DP_SUM, CMD_DP_CHUNKEDSTD, CMD_DP_GETAMPLI,
	CMD_DP_LINREG, CMD_DP_COERCIVITY, CMD_DP_REMANENCE, CMD_DP_COMPLETEHYSTERLOOP,
	CMD_DP_DUMPTDEP,
	
	CMD_DP_SMOOTH, CMD_DP_MONOTONIC,
	CMD_DP_CROSSINGSHISTOGRAM, CMD_DP_CROSSINGSFREQUENCY, CMD_DP_PEAKSFREQUENCY,

	//Fitting

	CMD_DP_FITLORENTZ, CMD_DP_FITLORENTZ2, CMD_DP_FITSKYRMION, CMD_DP_FITDW, CMD_DP_FITSTT, CMD_DP_FITSOT, CMD_DP_FITSOTSTT, CMD_DP_CALCSOT, CMD_DP_FITADIABATIC, CMD_DP_FITNONADIABATIC,

	//used for testing various things
	CMD_TEST
};



//Describes the structure of a command, and helper methods for processing commands
struct CommandSpecifier {

public:

	//a unique command identifier
	CMD_ cmdId;

	//show usage and command description using formatted strings (for display in console).
	//Also show description of returned values for script clients.
	//If unit is specified, allow input using physical quantity units rather then e notation
	std::string usage, descr, return_descr;
	std::string unit;

	//limits for parameters specified here; the pair gives the lower and upper limit (inclusive). To omit a limit set it as an Any(). 
	//Otherwise must set the Any with the expected parameter type, in the order expected!
	std::vector<std::pair<Any, Any>> limits;

private:

	template <typename Type, typename ... PType>
	BError check_limits(int index, Type param)
	{
		BError error;

		if (!limits.size()) return error;

		auto make_error_info = [&]() -> std::string {

			std::string limits_string = "Parameter " + ToString(index) + ", ";

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

		auto make_error_info = [&]() -> std::string {

			std::string limits_string = "Parameter " + ToString(index) + ", ";

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
	BError GetParameters(std::vector<std::string>& components, Type& ... values) 
	{
		BError error;

		//get parameters as a vector of Any objects. Construct  Any objects using pointers, so changes to vector entries result in changes to passed parameters - since these are in turn passed through reference changes are propagated to caller.
		std::vector<Any> a_values = make_vector(Any(&values)...);

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
		
		//if last value is a std::string, add any remaining components to it using space separators.
		if (a_values.back().is_type(btype_info<std::string>())) {

			std::string text = a_values.back();

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
	std::vector<std::string> PrepareReturnParameters(Type... values)
	{
		std::vector<std::string> outVector;

		for (std::string value_text : { ToString(values)... })
			JoinToVector(outVector, split(value_text, ", ", "; "));

		return outVector;
	}
};