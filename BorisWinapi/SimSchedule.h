#pragma once

#include "BorisLib.h"

using namespace std;

//simulation stage/step settings -> add new values at the end to keep older simulation files compatible
enum SS_ { SS_RELAX = 0, SS_HFIELDXYZ, SS_HFIELDXYZSEQ, SS_HPOLARSEQ, SS_HFMR, SS_V, SS_VSEQ, SS_VSIN, SS_VCOS, SS_I, SS_ISEQ, SS_ISIN, SS_ICOS, SS_T, SS_TSEQ };

//simulation stage stop conditions -> add new values at the end to keep older simulation files compatible
enum STOP_ { STOP_NOSTOP = 0, STOP_ITERATIONS, STOP_MXH, STOP_TIME, STOP_DMDT };

//data save conditions -> add new values at the end to keep older simulation files compatible
enum DSAVE_ { DSAVE_NONE = 0, DSAVE_STAGE, DSAVE_STEP, DSAVE_ITER, DSAVE_TIME };

struct StageDescriptor : 
	public ProgramState<StageDescriptor, tuple<int, string, bool>, tuple<> >
{

	int stageSetting;
	string unit;

	bool meshless;

	StageDescriptor(int stageSetting_, string unit_ = "", bool meshless_ = true) :
		ProgramStateNames(this, { VINFO(stageSetting), VINFO(unit) , VINFO(meshless) }, {})
	{
		stageSetting = stageSetting_;
		unit = unit_;
		meshless = meshless_;
	}

	StageDescriptor(void) :
		ProgramStateNames(this, { VINFO(stageSetting), VINFO(unit) , VINFO(meshless) }, {})
	{
		stageSetting = SS_RELAX;
		unit = "";
		meshless = true;
	}

	StageDescriptor(const StageDescriptor& copyThis) :
		ProgramStateNames(this, { VINFO(stageSetting), VINFO(unit) , VINFO(meshless) }, {})
	{
		*this = copyThis;
	}

	StageDescriptor& operator=(const StageDescriptor& copyThis)
	{
		stageSetting = copyThis.stageSetting;
		unit = copyThis.unit;
		meshless = copyThis.meshless;
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

	~StageConfig() { }

	//implement pure virtual method from ProgramState
	void RepairObjectState(void) {}

	//----------------------------------- VALUES SETTERS - Call after constructor

	template <typename Type> void set_value(Type value) { setValue = Any(); setValue = value; }
	template <typename Type> void set_stopvalue(Type value) { stopValue = Any(); stopValue = value; }
	template <typename Type> void set_dsavevalue(Type value) { dsaveValue = Any(); dsaveValue = value; }

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
		if(setValue.is_type(typeid(SEQ))) {

			SEQ sequence = setValue;
			return sequence.number_of_steps();
		}
		else if(setValue.is_type(typeid(SEQ3))) {

			SEQ3 sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(typeid(SEQP))) {

			SEQP sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(typeid(COSSEQ))) {

			COSSEQ sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(typeid(COSSEQ3))) {

			COSSEQ3 sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(typeid(SINOSC))) {

			SINOSC sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(typeid(SINOSC3))) {

			SINOSC3 sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(typeid(COSOSC))) {

			COSOSC sequence = setValue;
			return sequence.number_of_steps();
		}
		else if (setValue.is_type(typeid(COSOSC3))) {

			COSOSC3 sequence = setValue;
			return sequence.number_of_steps();
		}
		else return 0;
	}

	//----------------------------------- OTHER SETTERS

	void set_meshname(string meshName) { if(!stageDescriptor.meshless) this->meshName = meshName; }

	void set_stagevalue_fromcomponents(const vector<string>& value_components) { if(!setValue.IsNull()) setValue.convert_string(value_components, stageDescriptor.unit); }

	void set_stagevalue_fromstring(const string& value_string) { if(!setValue.IsNull()) setValue.convert_string(value_string, stageDescriptor.unit); }

	void set_stopvalue_fromstring(const string& value_string) { if(!stopValue.IsNull()) stopValue.convert_string(value_string, stopDescriptor.unit); }

	void set_dsavevalue_fromstring(const string& value_string) { if(!dsaveValue.IsNull()) dsaveValue.convert_string(value_string, dataSaveDescriptor.unit); }

	//----------------------------------- GET STAGE SET VALUE

	bool IsValueSet(void) { return !setValue.IsNull(); }

	template <typename Type>
	Type get_value(int step = 0) { 
		
		//if a sequence then must get value for given step from within Sequence object.
		if(setValue.is_type(typeid(SEQ))) {

			SEQ sequence = setValue;
			double value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if(setValue.is_type(typeid(SEQ3))) {

			SEQ3 sequence = setValue;
			DBL3 value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(typeid(SEQP))) {

			SEQP sequence = setValue;
			DBL3 value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(typeid(COSSEQ))) {

			COSSEQ sequence = setValue;
			double value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(typeid(COSSEQ3))) {

			COSSEQ3 sequence = setValue;
			DBL3 value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(typeid(SINOSC))) {

			SINOSC sequence = setValue;
			double value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(typeid(SINOSC3))) {

			SINOSC3 sequence = setValue;
			DBL3 value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(typeid(COSOSC))) {

			COSOSC sequence = setValue;
			double value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}
		else if (setValue.is_type(typeid(COSOSC3))) {

			COSOSC3 sequence = setValue;
			DBL3 value = sequence.value(step);
			return *reinterpret_cast<Type*>(&value);
		}

		//not a sequence, so get value directly
		else return static_cast<Type>(setValue);
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
