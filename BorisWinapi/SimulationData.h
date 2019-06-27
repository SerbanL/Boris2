#pragma once

#include <string>

#include "BorisLib.h"

using namespace std;

//Simulation data available as outputs
enum DATA_ 
{ 
	DATA_NONE = 0,		 
	DATA_STAGESTEP, DATA_TIME, DATA_STAGETIME, DATA_ITERATIONS, DATA_SITERATIONS, DATA_DT, DATA_MXH,
	DATA_AVM, DATA_HA, 
	DATA_JC, DATA_JSX, DATA_JSY, DATA_JSZ, DATA_V, DATA_S, DATA_ELC, DATA_POTENTIAL, DATA_CURRENT, DATA_RESISTANCE,
	DATA_E_DEMAG, DATA_E_EXCH, DATA_E_SURFEXCH, DATA_E_ZEE, DATA_E_ANIS, DATA_E_ROUGH,
	DATA_DWSHIFT, DATA_SKYSHIFT,
	DATA_TRANSPORT_ITERSTOCONV, DATA_TRANSPORT_SITERSTOCONV, DATA_TRANSPORT_CONVERROR, DATA_TRANSPORT_ASOR,
	DATA_TEMP, DATA_HEATDT,
};

//Specifier for available output data : this is stored in a vector with lut indexing, where DATA_ values are used for the major id - the DatumSpecifier corresponds to it
//Further there's a key (a handle) associated with it (so DatumSpecifier needs to be stored in a vector_key_lut)
struct DatumSpecifier {

	//Label to display in data box (if available). If Label not set then this datum is not available for display in data box
	string Label;

	//physical unit (if any)
	string unit;

	//number of components for this data entry (e.g a vector quantity has 3 components, etc.). This info is just used for the data save file header when labelling the data saved in columns.
	int components;

	//flags: meshless means this datum is not associated with any mesh. boxless means it doesn't need a box
	bool meshless, boxless;

	DatumSpecifier(string Label_, int components_, string unit_ = "", bool meshless_ = true, bool boxless_ = true) {

		Label = Label_;
		components = components_;
		meshless = meshless_;
		boxless = boxless_;
		unit = unit_;
	}

	DatumSpecifier(void) {

		components = 1;
		meshless = true;
		boxless = true;
	}
};

//Configuration for a particular datum, which may be selected to be saved during a simulation, or display in the data box
struct DatumConfig :
	public ProgramState<DatumConfig, tuple<int, string, Rect>, tuple<>>
{

	//the particular datum this configuration is for
	int datumId;

	//mesh name this configured datum is associated with. If blank then mesh name is not applicable. For super-mesh use SuperMesh:superMeshHandle (i.e. "supermesh")
	string meshName;

	//for same data a Box may be applicable : e.g. obtain average value in given mesh only in this box (subset of mesh size)
	Rect rectangle;

	DatumConfig(int datumId_, string meshName_, Rect rectangle_) : ProgramStateNames(this, { VINFO(datumId), VINFO(meshName), VINFO(rectangle) }, {})
	{
		datumId = datumId_;
		meshName = meshName_;
		rectangle = rectangle_;
	}

	DatumConfig(int datumId_, string meshName_) : ProgramStateNames(this, { VINFO(datumId), VINFO(meshName), VINFO(rectangle) }, {})
	{
		datumId = datumId_;
		meshName = meshName_;
	}

	DatumConfig(int datumId_) : ProgramStateNames(this, { VINFO(datumId), VINFO(meshName), VINFO(rectangle) }, {})
	{
		datumId = datumId_;
	}

	DatumConfig(void) : ProgramStateNames(this, { VINFO(datumId), VINFO(meshName), VINFO(rectangle) }, {})
	{
		datumId = DATA_NONE;
	}
	 
	DatumConfig(const DatumConfig& copyThis) : ProgramStateNames(this, { VINFO(datumId), VINFO(meshName), VINFO(rectangle) }, {})
	{
		*this = copyThis;
	}

	DatumConfig& operator=(const DatumConfig& copyThis)
	{
		datumId = copyThis.datumId;
		meshName = copyThis.meshName;
		rectangle = copyThis.rectangle;
		return *this;
	}

	//implement pure virtual method from ProgramState
	void RepairObjectState(void) {}
};