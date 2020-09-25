#pragma once

#include "BorisLib.h"
#include "SimSharedData.h"

using namespace std;

//simulation stage/step settings -> add new values at the end to keep older simulation files compatible
enum SS_ { 
	SS_RELAX = 0, 
	SS_HFIELDXYZ, SS_HFIELDXYZSEQ, SS_HPOLARSEQ, SS_HFMR, 
	SS_V, SS_VSEQ, SS_VSIN, SS_VCOS, 
	SS_I, SS_ISEQ, SS_ISIN, SS_ICOS, 
	SS_T, SS_TSEQ, 
	SS_Q, SS_QSEQ, 
	SS_HFIELDEQUATION, SS_HFIELDEQUATIONSEQ,
	SS_VEQUATION, SS_VEQUATIONSEQ,
	SS_IEQUATION, SS_IEQUATIONSEQ,
	SS_TEQUATION, SS_TEQUATIONSEQ,
	SS_QEQUATION, SS_QEQUATIONSEQ,
	SS_HFIELDFILE, SS_VFILE, SS_IFILE, SS_TFILE, SS_QFILE,
	SS_TSIGPOLAR,

	SS_MONTECARLO
};

//simulation stage stop conditions -> add new values at the end to keep older simulation files compatible
enum STOP_ { STOP_NOSTOP = 0, STOP_ITERATIONS, STOP_MXH, STOP_TIME, STOP_DMDT };

//data save conditions -> add new values at the end to keep older simulation files compatible
enum DSAVE_ { DSAVE_NONE = 0, DSAVE_STAGE, DSAVE_STEP, DSAVE_ITER, DSAVE_TIME };

struct StageDescriptor : 
	public ProgramState<StageDescriptor, tuple<int, string, bool, bool>, tuple<> >
{

	int stageSetting;
	string unit;

	bool meshless;

	//we don't want to display some stage descriptors : set this to true in that case
	bool invisible = false;

	StageDescriptor(int stageSetting_, string unit_ = "", bool meshless_ = true, bool invisible_ = false) :
		ProgramStateNames(this, { VINFO(stageSetting), VINFO(unit) , VINFO(meshless), VINFO(invisible) }, {})
	{
		stageSetting = stageSetting_;
		unit = unit_;
		meshless = meshless_;
		invisible = invisible_;
	}

	StageDescriptor(void) :
		ProgramStateNames(this, { VINFO(stageSetting), VINFO(unit) , VINFO(meshless), VINFO(invisible) }, {})
	{
		stageSetting = SS_RELAX;
		unit = "";
		meshless = true;
		invisible = false;
	}

	StageDescriptor(const StageDescriptor& copyThis) :
		ProgramStateNames(this, { VINFO(stageSetting), VINFO(unit) , VINFO(meshless), VINFO(invisible) }, {})
	{
		*this = copyThis;
	}

	StageDescriptor& operator=(const StageDescriptor& copyThis)
	{
		stageSetting = copyThis.stageSetting;
		unit = copyThis.unit;
		meshless = copyThis.meshless;
		invisible = copyThis.invisible;
		return *this;
	}

	//implement pure virtual method from ProgramState
	void RepairObjectState(void) {}
};

struct StageStopDescriptor :
	public ProgramState<StageStopDescriptor, tuple<int, string>, tuple<> >
{

	int stopCondition;
	string unit;

	StageStopDescriptor(int stopCondition_, string unit_ = "") :
		ProgramStateNames(this, { VINFO(stopCondition), VINFO(unit) }, {})
	{
		stopCondition = stopCondition_;
		unit = unit_;
	}

	StageStopDescriptor(void) :
		ProgramStateNames(this, { VINFO(stopCondition), VINFO(unit) }, {})
	{
		stopCondition = STOP_NOSTOP;
		unit = "";
	}

	StageStopDescriptor(const StageStopDescriptor& copyThis) :
		ProgramStateNames(this, { VINFO(stopCondition), VINFO(unit) }, {})
	{
		*this = copyThis;
	}

	StageStopDescriptor& operator=(const StageStopDescriptor& copyThis)
	{
		stopCondition = copyThis.stopCondition;
		unit = copyThis.unit;
		return *this;
	}

	//implement pure virtual method from ProgramState
	void RepairObjectState(void) {}
};

struct DataSaveDescriptor :
	public ProgramState<DataSaveDescriptor, tuple<int, string>, tuple<> >
{

	int dataSaveType;
	string unit;

	DataSaveDescriptor(int dataSaveType_, string unit_ = "") :
		ProgramStateNames(this, { VINFO(dataSaveType), VINFO(unit) }, {})
	{
		dataSaveType = dataSaveType_;
		unit = unit_;
	}

	DataSaveDescriptor(void) :
		ProgramStateNames(this, { VINFO(dataSaveType), VINFO(unit) }, {})
	{	
		dataSaveType = DSAVE_NONE;
		unit = "";
	}

	DataSaveDescriptor(const DataSaveDescriptor& copyThis) :
		ProgramStateNames(this, { VINFO(dataSaveType), VINFO(unit) }, {})
	{
		*this = copyThis;
	}

	DataSaveDescriptor& operator=(const DataSaveDescriptor& copyThis)
	{
		dataSaveType = copyThis.dataSaveType;
		unit = copyThis.unit;
		return *this;
	}

	//implement pure virtual method from ProgramState
	void RepairObjectState(void) {}
};

class StageConfig :
	public SimulationSharedData,
	public ProgramState<StageConfig, tuple< StageDescriptor, Any, StageStopDescriptor, Any, string, DataSaveDescriptor, Any >, tuple<> >
{

private:

	//stage type (together with associated properties) and any particular value set for this stage type
	StageDescriptor stageDescriptor;
	Any setValue;

	//stop condition (together with associated properties) and any particular value set for this stop condition
	StageStopDescriptor stopDescriptor;
	Any stopValue;
	
	//mesh name for the stage setting, if applicable
	string meshName;

	//data saving configuration and threshold value for saving data (e.g. number of iterations, time interval etc.)
	DataSaveDescriptor dataSaveDescriptor;
	Any dsaveValue;

private:

	//----------------------------------- VALUES SETTERS : SPECIAL

	//depending on the set stage type, we may need to adjust the stage value depending on other settings
	void adjust_setValue_specials(void)
	{
		//depending on the set stage type, we may need to adjust the stage value depending on the set stop value
		if ((SS_)stageDescriptor.stageSetting == SS_HFIELDFILE) {

			FILESEQ3 fseq = setValue;
			fseq.set_filename_and_time_resolution(directory, fseq.get_fileName(), (double)stopValue);
			setValue = fseq;
		}

		else if ((SS_)stageDescriptor.stageSetting == SS_VFILE ||
			(SS_)stageDescriptor.stageSetting == SS_IFILE ||
			(SS_)stageDescriptor.stageSetting == SS_TFILE ||
			(SS_)stageDescriptor.stageSetting == SS_QFILE) {

			FILESEQ fseq = setValue;
			fseq.set_filename_and_time_resolution(directory, fseq.get_fileName(), (double)stopValue);
			setValue = fseq;
		}
	}

public:

	StageConfig(void) :
		ProgramStateNames(this, 
			{
				VINFO(stageDescriptor), VINFO(setValue),
				VINFO(stopDescriptor), VINFO(stopValue),
				VINFO(meshName),
				VINFO(dataSaveDescriptor), VINFO(dsaveValue)
			}, {})
	{}

	//make a new stage
	StageConfig(StageDescriptor &stageDescriptor_, StageStopDescriptor &stopDescriptor_, string meshName_ = "") :
		ProgramStateNames(this,
			{
				VINFO(stageDescriptor), VINFO(setValue),
				VINFO(stopDescriptor), VINFO(stopValue),
				VINFO(meshName),
				VINFO(dataSaveDescriptor), VINFO(dsaveValue)
			}, {})
	{
		stageDescriptor = stageDescriptor_;

		if(!stageDescriptor.meshless)
			meshName = meshName_;
		else meshName = "";

		stopDescriptor = stopDescriptor_;

		//default dataSaveDescriptor is to save nothing (DSAVE_NONE) - no need to set it, DataSaveDescriptor constructor defaults to this.
	}

	StageConfig(const StageConfig& copyThis) :
		ProgramStateNames(this,
			{
				VINFO(stageDescriptor), VINFO(setValue),
				VINFO(stopDescriptor), VINFO(stopValue),
				VINFO(meshName),
				VINFO(dataSaveDescriptor), VINFO(dsaveValue)
			}, {})
	{
		*this = copyThis;
	}

	StageConfig& operator=(const StageConfig& copyThis)
	{
		stageDescriptor = copyThis.stageDescriptor;
		setValue = copyThis.setValue;
		stopDescriptor = copyThis.stopDescriptor;
		stopValue = copyThis.stopValue;
		meshName = copyThis.meshName;
		dataSaveDescriptor = copyThis.dataSaveDescriptor;
		dsaveValue = copyThis.dsaveValue;
		return *this;
	}

	//like the assignment operator but only keep the general stage data, not stage value
	void copy_stage_general_data(const StageConfig& copyThis)
	{
		stopDescriptor = copyThis.stopDescriptor;
		stopValue = copyThis.stopValue;
		meshName = copyThis.meshName;
		dataSaveDescriptor = copyThis.dataSaveDescriptor;
		dsaveValue = copyThis.dsaveValue;
	}

	~StageConfig() { }

	//implement pure virtual method from ProgramState
	void RepairObjectState(void) { adjust_setValue_specials(); }

	//----------------------------------- VALUES SETTERS - Call after constructor

	template <typename Type> 
	void set_value(Type value) 
	{ 
		setValue = Any(); 
		setValue = value;
		
		adjust_setValue_specials();
	}
	
	template <typename Type> 
	void set_stopvalue(Type value)
	{ 
		stopValue = Any(); 
		stopValue = value; 

		//for certain stage types we need to adjust them when we change the stopvalue
		adjust_setValue_specials();
	}

	template <typename Type> 
	void set_dsavevalue(Type value) 
	{ 
		dsaveValue = Any(); 
		dsaveValue = value; 
	}

	void clear_stopvalue(void) { stopValue = Any(); }
	void clear_stagevalue(void) { setValue = Any(); }
	void clear_dsavevalue(void) { dsaveValue = Any(); }

	void set_stoptype(StageStopDescriptor &stopDescriptor) { this->stopDescriptor = stopDescriptor; }

	void set_dsavetype(DataSaveDescriptor &dataSaveDescriptor) { this->dataSaveDescriptor = dataSaveDescriptor; }

	//----------------------------------- GETTERS

	//get stop condition
	STOP_ stop_condition(void) { return (STOP_)stopDescriptor.stopCondition; }

	//get stage type
	SS_ stage_type(void) { return (SS_)stageDescriptor.stageSetting; }

	//get save data type
	DSAVE_ dsave_type(void) { return (DSAVE_)dataSaveDescriptor.dataSaveType; }

	string meshname(void) { return meshName; }

	//get number of steps configured for this stage (will be zero if single value set)
	int number_of_steps(void) { 
	
		//number of steps is zero (a single value set this stage) unless this is a sequence of steps : check possible types of sequences and get their number of steps.
		if(setValue.is_type(btype_info<SEQ>())) {

			SEQ sequence = setValue;
			return sequence.number_of_steps();
		}
		else if(setValue.is_type(btype_info<SEQ3>())) {

			SEQ3 sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(btype_info<SEQP>())) {

			SEQP sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(btype_info<COSSEQ>())) {

			COSSEQ sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(btype_info<COSSEQ3>())) {

			COSSEQ3 sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(btype_info<SINOSC>())) {

			SINOSC sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(btype_info<SINOSC3>())) {

			SINOSC3 sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(btype_info<COSOSC>())) {

			COSOSC sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(btype_info<COSOSC3>())) {

			COSOSC3 sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(btype_info<StringSequence>())) {

			StringSequence sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(btype_info<FILESEQ>())) {

			FILESEQ sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(btype_info<FILESEQ3>())) {

			FILESEQ3 sequence = setValue;
			return sequence.number_of_steps();
		}

		//cannot determine number of stage steps
		else return 0;
	}

	//----------------------------------- OTHER SETTERS

	void set_meshname(string meshName) { if(!stageDescriptor.meshless) this->meshName = meshName; }

	void set_stagevalue_fromcomponents(const vector<string>& value_components) 
	{ 
		if (!setValue.IsNull()) {

			setValue.convert_string(value_components, stageDescriptor.unit);
			adjust_setValue_specials();
		}
	}

	void set_stagevalue_fromstring(const string& value_string) 
	{ 
		if (!setValue.IsNull()) {

			setValue.convert_string(value_string, stageDescriptor.unit);
			adjust_setValue_specials();
		}
	}

	void set_stopvalue_fromstring(const string& value_string) 
	{ 
		if (!stopValue.IsNull()) {

			stopValue.convert_string(value_string, stopDescriptor.unit);
			
			//for certain stage types we need to adjust them when we change the stopvalue
			adjust_setValue_specials();
		}
	}

	void set_dsavevalue_fromstring(const string& value_string) { if(!dsaveValue.IsNull()) dsaveValue.convert_string(value_string, dataSaveDescriptor.unit); }

	//----------------------------------- GET STAGE SET VALUE

	bool IsValueSet(void) { return !setValue.IsNull(); }

	template <typename Type>
	Type get_value(int step = 0) { 
		
		//if a sequence then must get value for given step from within Sequence object.
		if(setValue.is_type(btype_info<SEQ>())) {

			SEQ sequence = setValue;
			double value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if(setValue.is_type(btype_info<SEQ3>())) {

			SEQ3 sequence = setValue;
			DBL3 value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(btype_info<SEQP>())) {

			SEQP sequence = setValue;
			DBL3 value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(btype_info<COSSEQ>())) {

			COSSEQ sequence = setValue;
			double value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(btype_info<COSSEQ3>())) {

			COSSEQ3 sequence = setValue;
			DBL3 value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(btype_info<SINOSC>())) {

			SINOSC sequence = setValue;
			double value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(btype_info<SINOSC3>())) {

			SINOSC3 sequence = setValue;
			DBL3 value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(btype_info<COSOSC>())) {

			COSOSC sequence = setValue;
			double value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(btype_info<COSOSC3>())) {

			COSOSC3 sequence = setValue;
			DBL3 value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(btype_info<StringSequence>())) {

			StringSequence sequence = setValue;
			std::string value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(btype_info<FILESEQ>())) {

			FILESEQ sequence = setValue;
			double value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(btype_info<FILESEQ3>())) {

			FILESEQ3 sequence = setValue;
			DBL3 value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}

		//not a sequence, so get value directly
		else return setValue;
	}

	string get_value_string(void) { return setValue.convert_to_string( stageDescriptor.unit ); }
	
	//----------------------------------- GET STOP CONDITION VALUE

	bool IsStopValueSet(void) { return !stopValue.IsNull(); }

	Any& get_stopvalue(void) { 
		
		return stopValue;
	}

	string get_stopvalue_string(void) {

		return stopValue.convert_to_string( stopDescriptor.unit );
	}

	//----------------------------------- GET DATA SAVE CONDITION VALUE

	bool IsdSaveValueSet(void) { return !dsaveValue.IsNull(); }

	Any& get_dsavevalue(void) { 
		
		return dsaveValue;
	}

	string get_dsavevalue_string(void) {

		return dsaveValue.convert_to_string( dataSaveDescriptor.unit );
	}
};
