#include "stdafx.h"
#include "Simulation.h"

void Simulation::AddGenericStage(SS_ stageType, string meshName) 
{
	switch(stageType) {

	case SS_RELAX:
	{
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_HFIELDXYZ:
	{
		//zero field with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH), meshName);
		stageConfig.set_value(DBL3(0, 0, 0));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_HFIELDXYZSEQ:
	{
		//zero field with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH), meshName);
		stageConfig.set_value(SEQ3(DBL3(-1e5, 0, 0), DBL3(1e5, 0, 0), 100));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_HPOLARSEQ:
	{
		//zero field with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH), meshName);
		stageConfig.set_value(SEQP(DBL3(-1e5, 90, 0), DBL3(1e5, 90, 0), 100));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_HFMR:
	{
		//Bias field along y with Hrf along x. 1 GHz.
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_TIME), meshName);
		stageConfig.set_value(COSSEQ3(DBL3(0, 1e6, 0), DBL3(1e3, 0, 0), 20, 100));
		stageConfig.set_stopvalue(50e-12);

		simStages.push_back(stageConfig);
	}
	break;
	
	case SS_HFIELDEQUATION:
	{
		//zero field with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH), meshName);
		stageConfig.set_value(std::string("0, 0, 0"));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_HFIELDEQUATIONSEQ:
	{
		//zero field with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH), meshName);
		stageConfig.set_value(StringSequence(std::string("1: 0, 0, 0")));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_HFIELDFILE:
	{
		//zero field with STOP_TIME
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_TIME), meshName);
		stageConfig.set_value(FILESEQ3(directory, "file.txt", 1e-9));
		stageConfig.set_stopvalue(1e-9);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_V:
	{
		//zero potential with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH));
		stageConfig.set_value(0.0);
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_VSEQ:
	{
		//V 0.0 to 1.0 V in 10 steps with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH));
		stageConfig.set_value(SEQ(0.0, 1.0, 10));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	/*
	case SS_VSIN:
	{
		//10 mV oscillation 1 GHz for 100 cycles
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_TIME));
		stageConfig.set_value(SINOSC(10e-3, 20, 100));
		stageConfig.set_stopvalue(50e-12);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_VCOS:
	{
		//10 mV oscillation 1 GHz for 100 cycles
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_TIME));
		stageConfig.set_value(COSOSC(10e-3, 20, 100));
		stageConfig.set_stopvalue(50e-12);

		simStages.push_back(stageConfig);
	}
	break;
	*/

	case SS_VEQUATION:
	{
		//zero potential with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH));
		stageConfig.set_value(std::string("0"));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_VEQUATIONSEQ:
	{
		//zero potential with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH));
		stageConfig.set_value(StringSequence(std::string("1: 0")));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_VFILE:
	{
		//zero field with STOP_TIME
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_TIME), meshName);
		stageConfig.set_value(FILESEQ(directory, "file.txt", 1e-9));
		stageConfig.set_stopvalue(1e-9);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_I:
	{
		//zero current with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH));
		stageConfig.set_value(0.0);
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_ISEQ:
	{
		//I 0.0 to 1.0 mA in 10 steps with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH));
		stageConfig.set_value(SEQ(0.0, 1.0e-3, 10));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	/*
	case SS_ISIN:
	{
		//1 mA oscillation 1 GHz for 100 cycles
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_TIME));
		stageConfig.set_value(SINOSC(1e-3, 20, 100));
		stageConfig.set_stopvalue(50e-12);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_ICOS:
	{
		//1 mA oscillation 1 GHz for 100 cycles
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_TIME));
		stageConfig.set_value(COSOSC(1e-3, 20, 100));
		stageConfig.set_stopvalue(50e-12);

		simStages.push_back(stageConfig);
	}
	break;
	*/

	case SS_IEQUATION:
	{
		//zero current with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH));
		stageConfig.set_value(std::string("0"));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_IEQUATIONSEQ:
	{
		//zero current with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH));
		stageConfig.set_value(StringSequence(std::string("1: 0")));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_IFILE:
	{
		//zero field with STOP_TIME
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_TIME), meshName);
		stageConfig.set_value(FILESEQ(directory, "file.txt", 1e-9));
		stageConfig.set_stopvalue(1e-9);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_T:
	{
		//zero temperature with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH), meshName);
		stageConfig.set_value(0.0);
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_TSEQ:
	{
		//T 0.0 to 300K in 10 steps with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH), meshName);
		stageConfig.set_value(SEQ(0.0, 300.0, 10));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_TEQUATION:
	{
		//zero temperature with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH), meshName);
		stageConfig.set_value(std::string("0"));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;
		
	case SS_TEQUATIONSEQ:
	{
		//zero temperature with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH), meshName);
		stageConfig.set_value(StringSequence(std::string("1: 0")));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_TFILE:
	{
		//zero field with STOP_TIME
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_TIME), meshName);
		stageConfig.set_value(FILESEQ(directory, "file.txt", 1e-9));
		stageConfig.set_stopvalue(1e-9);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_Q:
	{
		//1e19 W/m3 heat source with STOP_TIME of 10ns
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_TIME), meshName);
		stageConfig.set_value(1e19);
		stageConfig.set_stopvalue(10e-9);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_QSEQ:
	{
		//Q 0.0 to 1e19 W/m3 in 10 steps with STOP_TIME of 1ns
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_TIME), meshName);
		stageConfig.set_value(SEQ(0.0, 1e19, 10));
		stageConfig.set_stopvalue(1e-9);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_QEQUATION:
	{
		//zero Q with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH), meshName);
		stageConfig.set_value(std::string("0"));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_QEQUATIONSEQ:
	{
		//zero Q with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH), meshName);
		stageConfig.set_value(StringSequence(std::string("1: 0")));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_QFILE:
	{
		//zero field with STOP_TIME
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_TIME), meshName);
		stageConfig.set_value(FILESEQ(directory, "file.txt", 1e-9));
		stageConfig.set_stopvalue(1e-9);

		simStages.push_back(stageConfig);
	}
	break;

	case SS_TSIGPOLAR:
	{
		//zero stress with STOP_MXH
		StageConfig stageConfig = StageConfig(stageDescriptors(stageType), stageStopDescriptors(STOP_MXH), meshName);
		stageConfig.set_value(DBL3(0, 0, 0));
		stageConfig.set_stopvalue(1e-4);

		simStages.push_back(stageConfig);
	}
	break;
	}
}

void Simulation::DeleteStage(int stageIndex) 
{
	simStages.erase(stageIndex);
}

void Simulation::SetGenericStopCondition(int index, STOP_ stopType) 
{
	if(!GoodIdx(simStages.last(), index)) return;

	switch(stopType) {

	case STOP_NOSTOP:
		simStages[index].set_stoptype( stageStopDescriptors(stopType) );
		simStages[index].clear_stopvalue();
		break;

	case STOP_ITERATIONS:
		simStages[index].set_stoptype( stageStopDescriptors(stopType) );
		simStages[index].set_stopvalue(1000);
		break;

	case STOP_MXH:
		simStages[index].set_stoptype( stageStopDescriptors(stopType) );
		simStages[index].set_stopvalue(1e-4);
		break;

	case STOP_DMDT:
		simStages[index].set_stoptype(stageStopDescriptors(stopType));
		simStages[index].set_stopvalue(1e-5);
		break;

	case STOP_TIME:
		simStages[index].set_stoptype( stageStopDescriptors(stopType) );
		simStages[index].set_stopvalue(10e-9);
		break;
	}
}

void Simulation::SetGenericDataSaveCondition(int index, DSAVE_ dsaveType)
{
	switch(dsaveType) {

	case DSAVE_NONE:
		simStages[index].set_dsavetype( dataSaveDescriptors(dsaveType) );
		simStages[index].clear_dsavevalue();
		break;

	case DSAVE_STAGE:
		simStages[index].set_dsavetype( dataSaveDescriptors(dsaveType) );
		simStages[index].clear_dsavevalue();
		break;

	case DSAVE_STEP:
		simStages[index].set_dsavetype( dataSaveDescriptors(dsaveType) );
		simStages[index].clear_dsavevalue();
		break;

	case DSAVE_ITER:
		simStages[index].set_dsavetype( dataSaveDescriptors(dsaveType) );
		simStages[index].set_dsavevalue(100);
		break;

	case DSAVE_TIME:
		simStages[index].set_dsavetype( dataSaveDescriptors(dsaveType) );
		simStages[index].set_dsavevalue(1e-9);
		break;
	}
}

void Simulation::EditStageType(int index, SS_ stageType, string meshName) 
{
	//if same stage type as before just change the mesh name
	if(GoodIdx(simStages.last(), index) && simStages[index].stage_type() == stageType) {
		
		simStages[index].set_meshname(meshName);
	}
	else {

		//new stage type at this index so set a generic stage to start off with
		AddGenericStage(stageType, meshName);
		simStages.move(simStages.last(), index);
		simStages.erase(index + 1);
	}
}

void Simulation::EditStageValue(int stageIndex, string value_string) 
{
	bool adjust_special = false;

	//in some special cases, changing the stage value can result in a stage type change
	if (simStages[stageIndex].stage_type() == SS_HFIELDEQUATION || simStages[stageIndex].stage_type() == SS_HFIELDEQUATIONSEQ) {

		size_t pos = value_string.find(":");
		//change from equation to equation sequence
		if (pos != std::string::npos && simStages[stageIndex].stage_type() == SS_HFIELDEQUATION) {

			adjust_special = true;

			AddGenericStage(SS_HFIELDEQUATIONSEQ, simStages[stageIndex].meshname());
		}
		//change from equation sequence to equation
		else if (pos == std::string::npos && simStages[stageIndex].stage_type() == SS_HFIELDEQUATIONSEQ) {

			adjust_special = true;

			AddGenericStage(SS_HFIELDEQUATION, simStages[stageIndex].meshname());
		}
	}

	else if (simStages[stageIndex].stage_type() == SS_VEQUATION || simStages[stageIndex].stage_type() == SS_VEQUATIONSEQ) {

		size_t pos = value_string.find(":");
		//change from equation to equation sequence
		if (pos != std::string::npos && simStages[stageIndex].stage_type() == SS_VEQUATION) {

			adjust_special = true;

			AddGenericStage(SS_VEQUATIONSEQ, simStages[stageIndex].meshname());
		}
		//change from equation sequence to equation
		else if (pos == std::string::npos && simStages[stageIndex].stage_type() == SS_VEQUATIONSEQ) {

			adjust_special = true;

			AddGenericStage(SS_VEQUATION, simStages[stageIndex].meshname());
		}
	}

	else if (simStages[stageIndex].stage_type() == SS_IEQUATION || simStages[stageIndex].stage_type() == SS_IEQUATIONSEQ) {

		size_t pos = value_string.find(":");
		//change from equation to equation sequence
		if (pos != std::string::npos && simStages[stageIndex].stage_type() == SS_IEQUATION) {

			adjust_special = true;

			AddGenericStage(SS_IEQUATIONSEQ, simStages[stageIndex].meshname());
		}
		//change from equation sequence to equation
		else if (pos == std::string::npos && simStages[stageIndex].stage_type() == SS_IEQUATIONSEQ) {

			adjust_special = true;

			AddGenericStage(SS_IEQUATION, simStages[stageIndex].meshname());
		}
	}

	else if (simStages[stageIndex].stage_type() == SS_TEQUATION || simStages[stageIndex].stage_type() == SS_TEQUATIONSEQ) {

		size_t pos = value_string.find(":");
		//change from equation to equation sequence
		if (pos != std::string::npos && simStages[stageIndex].stage_type() == SS_TEQUATION) {

			adjust_special = true;

			AddGenericStage(SS_TEQUATIONSEQ, simStages[stageIndex].meshname());
		}
		//change from equation sequence to equation
		else if (pos == std::string::npos && simStages[stageIndex].stage_type() == SS_TEQUATIONSEQ) {

			adjust_special = true;

			AddGenericStage(SS_TEQUATION, simStages[stageIndex].meshname());
		}
	}

	else if (simStages[stageIndex].stage_type() == SS_QEQUATION || simStages[stageIndex].stage_type() == SS_QEQUATIONSEQ) {

		size_t pos = value_string.find(":");
		//change from equation to equation sequence
		if (pos != std::string::npos && simStages[stageIndex].stage_type() == SS_QEQUATION) {

			adjust_special = true;

			AddGenericStage(SS_QEQUATIONSEQ, simStages[stageIndex].meshname());
		}
		//change from equation sequence to equation
		else if (pos == std::string::npos && simStages[stageIndex].stage_type() == SS_QEQUATIONSEQ) {

			adjust_special = true;

			AddGenericStage(SS_QEQUATION, simStages[stageIndex].meshname());
		}
	}

	if (!adjust_special) {

		//edit current stage
		simStages[stageIndex].set_stagevalue_fromstring(value_string);
	}
	else {

		//current stage type has been changed : we have a new generic stage of the required type at the end

		//set required value and copy current generate stage data
		simStages[simStages.last()].set_stagevalue_fromstring(value_string);
		simStages[simStages.last()].copy_stage_general_data(simStages[stageIndex]);

		//replace current stage
		simStages.move(simStages.last(), stageIndex);
		simStages.erase(stageIndex + 1);
	}
}

void Simulation::EditStageStopCondition(int index, STOP_ stopType, string stopValueString) 
{
	//if same stop condition as before just change the stop value
	if(GoodIdx(simStages.last(), index) && simStages[index].stop_condition() == stopType) {

		if(stopValueString.length()) simStages[index].set_stopvalue_fromstring(stopValueString);
	}
	else {

		SetGenericStopCondition(index, stopType);
		if(stopValueString.length()) simStages[index].set_stopvalue_fromstring(stopValueString);
	}
}

void Simulation::EditDataSaveCondition(int index, DSAVE_ dsaveType, string dsaveValueString)
{
	//if same saving condition as before just change the value
	if(GoodIdx(simStages.last(), index) && simStages[index].dsave_type() == dsaveType) {

		if(dsaveValueString.length()) simStages[index].set_dsavevalue_fromstring(dsaveValueString);
	}
	else {

		SetGenericDataSaveCondition(index, dsaveType);
		if(dsaveValueString.length()) simStages[index].set_dsavevalue_fromstring(dsaveValueString);
	}
}

void Simulation::UpdateStageMeshNames(string oldMeshName, string newMeshName) 
{
	for(int idx = 0; idx < simStages.size(); idx++) {

		if (simStages[idx].meshname() == oldMeshName) {

			simStages[idx].set_meshname(newMeshName);
		}
	}
}

INT2 Simulation::Check_and_GetStageStep()
{
	//first make sure stage value is correct - this could only happen if stages have been deleted. If incorrect just reset back to 0.
	if(stage_step.major >= simStages.size()) stage_step = INT2();

	//mak sure step value is correct - if incorrect reset back to zero.
	if(stage_step.minor > simStages[stage_step.major].number_of_steps()) stage_step.minor = 0;

	return stage_step;
}

void Simulation::CheckSimulationSchedule(void) 
{
	//if stage index exceeds number of stages then just set it to the end : stages must have been deleted whilst simulation running.
	if(stage_step.major >= simStages.size()) stage_step.major = simStages.last();

	switch( simStages[ stage_step.major ].stop_condition() ) {

	case STOP_NOSTOP:
		break;

	case STOP_ITERATIONS:

		if( SMesh.GetStageIteration() >= (int)simStages[ stage_step.major ].get_stopvalue() ) AdvanceSimulationSchedule();
		break;

	case STOP_MXH:

		if( SMesh.Get_mxh() <= (double)simStages[ stage_step.major ].get_stopvalue() ) AdvanceSimulationSchedule();
		break;

	case STOP_DMDT:

		if (SMesh.Get_dmdt() <= (double)simStages[stage_step.major].get_stopvalue()) AdvanceSimulationSchedule();
		break;

	case STOP_TIME:
	{
		//See comments for DSAVE_TIME below in Simulation::CheckSaveDataConditions() - same thing applies here.

		double time = SMesh.GetTime();
		double tstop = (double)simStages[stage_step.major].get_stopvalue();
		double dT = SMesh.GetTimeStep();

		double delta = time - floor_fixedepsilon(time / tstop) * tstop;
		if (delta < dT * 0.99) AdvanceSimulationSchedule();
	}
		break;
	}
}

void Simulation::CheckSaveDataConditions() 
{
	switch (simStages[stage_step.major].dsave_type()) {

	case DSAVE_NONE:
	case DSAVE_STAGE:
	case DSAVE_STEP:
		//step and stage save data is done in AdvanceSimulationSchedule when step or stage ending is detected
		break;

	case DSAVE_ITER:
		if (!(SMesh.GetIteration() % (int)simStages[stage_step.major].get_dsavevalue())) SaveData();
		break;

	case DSAVE_TIME:
	{
		
		double time = SMesh.GetTime();
		double tsave = (double)simStages[stage_step.major].get_dsavevalue();
		double dT = SMesh.GetTimeStep();

		//the floor_fixedepsilon is important - don't use floor!
		//the reason for this, if time / tsave ends up being very close, but slightly less, than an integer, e.g. 1.999 due to a floating point error, then floor will round it down, whereas really it should be rounded up.
		//thus with floor only you can end up not saving data points where you should be.
		//Also delta < dT * 0.99 check below is important : if using just delta < dT check, delta can be slightly smaller than dT but within a floating point error close to it - thus we end up double-saving some data points!
		//this happens especially if tsave / dT is an integer -> thus most of the time.

		//Also use floor_epsilon instead of floor_fixedepsilon:
		//the epsilon value used in floor_epsilon for small time/tsave values is too small, meaning the time/tsave will again round down when it should be rounding up -> SaveData() will never be called
		//for large time/tsave value the floor_epsilon value is too coarse, meaning you can save data when you don't want to
		//The epsilon in floor_fixedepsilon is coarse, but not too coarse, which results in correct behaviour over a wide range of time/tsave values -> covers the relevant range.

		double delta = time - floor_fixedepsilon(time / tsave) * tsave;
		if (delta < dT * 0.99) SaveData();
	}
		break;
	}
}

void Simulation::AdvanceSimulationSchedule(void) 
{
	//assume stage_step.major is correct

	//do we need to iterate the transport solver? 
	//if static_transport_solver is true then the transport solver was stopped from iterating before reaching the end of a stage or step
	if (static_transport_solver) {

		//turn off flag for now to enable iterating the transport solver
		static_transport_solver = false;

#if COMPILECUDA == 1
		if (cudaEnabled) {

			SMesh.UpdateTransportSolverCUDA();
		}
		else {

			SMesh.UpdateTransportSolver();
		}
#else
		SMesh.UpdateTransportSolver();
#endif

		//turn flag back on
		static_transport_solver = true;
	}

	//first try to increment the step number
	if(stage_step.minor < simStages[stage_step.major].number_of_steps()) {

		//save data at end of current step?
		if(simStages[stage_step.major].dsave_type() == DSAVE_STEP) SaveData();

		//next step and set value for it
		stage_step.minor++;
		SetSimulationStageValue();
	}
	else {

		//save data at end of current stage?
		if(simStages[stage_step.major].dsave_type() == DSAVE_STAGE ||
		   simStages[stage_step.major].dsave_type() == DSAVE_STEP) 
			SaveData();

		//next stage
		stage_step.major++;
		stage_step.minor = 0;

		//if not at the end then set stage value for given stage_step
		if(stage_step.major < simStages.size()) {

			SetSimulationStageValue();
		}
		else {

			//schedule reached end: stop simulation. Note, since this routine is called from Simulate routine, which runs on the THREAD_LOOP thread, cannot stop THREAD_LOOP from within it: stop it from another thread.
			single_call_launch(&Simulation::StopSimulation, THREAD_HANDLEMESSAGE);

			//set stage step to start of last stage but without resetting anything : the idea is user can edit the stopping value for the last stage, e.g. add more time to it, then run simulation some more
			//if you want a complete reset then reset command must be issued in console
			stage_step = INT2(stage_step.major - 1, 0);
		}
	}
}

void Simulation::SetSimulationStageValue(void) 
{
	SMesh.NewStageODE();

	//assume stage_step is correct (if called from AdvanceSimulationSchedule it will be. could also be called directly at the start of a simulation with stage_step reset, so it's also correct).

	switch( simStages[stage_step.major].stage_type() ) {

	case SS_RELAX:
	break;

	case SS_HFIELDXYZ:
	case SS_HFIELDXYZSEQ:
	case SS_HPOLARSEQ:
	case SS_HFMR:
	case SS_HFIELDFILE:
	{
		string meshName = simStages[stage_step.major].meshname();

		DBL3 appliedField = simStages[stage_step.major].get_value<DBL3>(stage_step.minor);

		if (SMesh.contains(meshName)) SMesh[meshName]->CallModuleMethod(&ZeemanBase::SetField, appliedField);
		else if (meshName == SMesh.superMeshHandle) {

			for (int idx = 0; idx < SMesh.size(); idx++) {
				SMesh[idx]->CallModuleMethod(&ZeemanBase::SetField, appliedField);
			}
		}
	}
	break;

	case SS_TSIGPOLAR:
	{
		string meshName = simStages[stage_step.major].meshname();

		DBL3 appliedStress = simStages[stage_step.major].get_value<DBL3>(stage_step.minor);

		if (SMesh.contains(meshName)) SMesh[meshName]->CallModuleMethod(&MElastic::SetUniformStress, appliedStress);
		else if (meshName == SMesh.superMeshHandle) {

			for (int idx = 0; idx < SMesh.size(); idx++) {
				SMesh[idx]->CallModuleMethod(&MElastic::SetUniformStress, appliedStress);
			}
		}
	}
	break;

	case SS_HFIELDEQUATION:
	{
		string meshName = simStages[stage_step.major].meshname();

		std::string equation_text = simStages[stage_step.major].get_value<std::string>(stage_step.minor);

		if (SMesh.contains(meshName)) SMesh[meshName]->CallModuleMethod(&ZeemanBase::SetFieldEquation, equation_text, 0);
		else if (meshName == SMesh.superMeshHandle) {

			for (int idx = 0; idx < SMesh.size(); idx++) {
				SMesh[idx]->CallModuleMethod(&ZeemanBase::SetFieldEquation, equation_text, 0);
			}
		}
	}
	break;

	case SS_HFIELDEQUATIONSEQ:
	{
		string meshName = simStages[stage_step.major].meshname();

		std::string equation_text = simStages[stage_step.major].get_value<std::string>(stage_step.minor);
		//for a equation sequence we have "n: actual equation", where n is the number of steps
		std::string equation_equation_text = equation_text.substr(equation_text.find_first_of(':') + 1);

		if (SMesh.contains(meshName)) SMesh[meshName]->CallModuleMethod(&ZeemanBase::SetFieldEquation, equation_equation_text, stage_step.minor);
		else if (meshName == SMesh.superMeshHandle) {

			for (int idx = 0; idx < SMesh.size(); idx++) {
				SMesh[idx]->CallModuleMethod(&ZeemanBase::SetFieldEquation, equation_equation_text, stage_step.minor);
			}
		}
	}
	break;

	case SS_VEQUATION:
	{
		std::string equation_text = simStages[stage_step.major].get_value<std::string>(stage_step.minor);

		SMesh.CallModuleMethod(&STransport::SetPotentialEquation, equation_text, 0);
	}
	break;

	case SS_VEQUATIONSEQ:
	{
		std::string equation_text = simStages[stage_step.major].get_value<std::string>(stage_step.minor);

		//for a equation sequence we have "n: actual equation", where n is the number of steps
		std::string equation_equation_text = equation_text.substr(equation_text.find_first_of(':') + 1);

		SMesh.CallModuleMethod(&STransport::SetPotentialEquation, equation_equation_text, stage_step.minor);
	}
	break;

	case SS_V:
	case SS_VSEQ:
	//case SS_VSIN:
	//case SS_VCOS:
	case SS_VFILE:
	{
		double potential = simStages[stage_step.major].get_value<double>(stage_step.minor);

		SMesh.CallModuleMethod(&STransport::SetPotential, potential, true);
	}
	break;

	case SS_I:
	case SS_ISEQ:
	//case SS_ISIN:
	//case SS_ICOS:
	case SS_IFILE:
	{
		double current = simStages[stage_step.major].get_value<double>(stage_step.minor);

		SMesh.CallModuleMethod(&STransport::SetCurrent, current, true);
	}
	break;

	case SS_IEQUATION:
	{
		std::string equation_text = simStages[stage_step.major].get_value<std::string>(stage_step.minor);

		SMesh.CallModuleMethod(&STransport::SetCurrentEquation, equation_text, 0);
	}
	break;

	case SS_IEQUATIONSEQ:
	{
		std::string equation_text = simStages[stage_step.major].get_value<std::string>(stage_step.minor);

		//for a equation sequence we have "n: actual equation", where n is the number of steps
		std::string equation_equation_text = equation_text.substr(equation_text.find_first_of(':') + 1);

		SMesh.CallModuleMethod(&STransport::SetCurrentEquation, equation_equation_text, stage_step.minor);
	}
	break;

	case SS_T:
	case SS_TSEQ:
	case SS_TFILE:
	{
		string meshName = simStages[stage_step.major].meshname();

		double temperature = simStages[stage_step.major].get_value<double>(stage_step.minor);

		if (SMesh.contains(meshName)) SMesh[meshName]->SetBaseTemperature(temperature);
		else if (meshName == SMesh.superMeshHandle) {

			//all meshes
			for (int idx = 0; idx < SMesh.size(); idx++) {

				SMesh[idx]->SetBaseTemperature(temperature);
			}
		}
	}
	break;

	case SS_TEQUATION:
	{
		string meshName = simStages[stage_step.major].meshname();

		std::string equation_text = simStages[stage_step.major].get_value<std::string>(stage_step.minor);

		if (SMesh.contains(meshName)) SMesh[meshName]->SetBaseTemperatureEquation(equation_text, 0);
		else if (meshName == SMesh.superMeshHandle) {

			//all meshes
			for (int idx = 0; idx < SMesh.size(); idx++) {

				SMesh[idx]->SetBaseTemperatureEquation(equation_text, 0);
			}
		}
	}
	break;

	case SS_TEQUATIONSEQ:
	{
		string meshName = simStages[stage_step.major].meshname();

		std::string equation_text = simStages[stage_step.major].get_value<std::string>(stage_step.minor);

		//for a equation sequence we have "n: actual equation", where n is the number of steps
		std::string equation_equation_text = equation_text.substr(equation_text.find_first_of(':') + 1);

		if (SMesh.contains(meshName)) SMesh[meshName]->SetBaseTemperatureEquation(equation_text, stage_step.minor);
		else if (meshName == SMesh.superMeshHandle) {

			//all meshes
			for (int idx = 0; idx < SMesh.size(); idx++) {

				SMesh[idx]->SetBaseTemperatureEquation(equation_text, stage_step.minor);
			}
		}
	}
	break;

	case SS_Q:
	case SS_QSEQ:
	case SS_QFILE:
	{
		string meshName = simStages[stage_step.major].meshname();

		double Qvalue = simStages[stage_step.major].get_value<double>(stage_step.minor);

		if (SMesh.contains(meshName)) SMesh[meshName]->Q = Qvalue;
		else if (meshName == SMesh.superMeshHandle) {

			//all meshes
			for (int idx = 0; idx < SMesh.size(); idx++) {

				SMesh[idx]->Q = Qvalue;
			}
		}
	}
	break;

	case SS_QEQUATION:
	{
		string meshName = simStages[stage_step.major].meshname();

		std::string equation_text = simStages[stage_step.major].get_value<std::string>(stage_step.minor);

		if (SMesh.contains(meshName)) SMesh[meshName]->CallModuleMethod(&HeatBase::SetQEquation, equation_text, 0);
		else if (meshName == SMesh.superMeshHandle) {

			for (int idx = 0; idx < SMesh.size(); idx++) {
				SMesh[idx]->CallModuleMethod(&HeatBase::SetQEquation, equation_text, 0);
			}
		}
	}
	break;

	case SS_QEQUATIONSEQ:
	{
		string meshName = simStages[stage_step.major].meshname();

		std::string equation_text = simStages[stage_step.major].get_value<std::string>(stage_step.minor);

		//for a equation sequence we have "n: actual equation", where n is the number of steps
		std::string equation_equation_text = equation_text.substr(equation_text.find_first_of(':') + 1);

		if (SMesh.contains(meshName)) SMesh[meshName]->CallModuleMethod(&HeatBase::SetQEquation, equation_text, stage_step.minor);
		else if (meshName == SMesh.superMeshHandle) {

			for (int idx = 0; idx < SMesh.size(); idx++) {
				SMesh[idx]->CallModuleMethod(&HeatBase::SetQEquation, equation_text, stage_step.minor);
			}
		}
	}
	break;
	}
}