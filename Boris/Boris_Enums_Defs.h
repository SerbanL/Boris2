#pragma once

#include "CompileFlags.h"

//All widely applicable enums and defines go here. Further smaller scope enums and defines are given where needed.

//when UpdateConfiguration method is called, pass a value from this enum to signal the reason it was called
enum UPDATECONFIG_ {

	////////////////////////////
	//GENERIC
	////////////////////////////

	//no reason given
	UPDATECONFIG_GENERIC,

	////////////////////////////
	//PROGRAM-WIDE SCOPE
	////////////////////////////

	//Brute force update : when this message is received all objects handling UpdateConfiguration method should update fully
	UPDATECONFIG_FORCEUPDATE,

	//Generic message issued by RepairObjectState method
	UPDATECONFIG_REPAIROBJECTSTATE,

	//Cuda state has changed
	UPDATECONFIG_SWITCHCUDASTATE,

	////////////////////////////
	//SUPERMESH
	////////////////////////////

	//A supermesh cellsize has changed
	UPDATECONFIG_SMESH_CELLSIZE,

	////////////////////////////
	//MESHES
	////////////////////////////

	//mesh shape has changed (not the rectangle but the shape inside the rectangle)
	UPDATECONFIG_MESHSHAPECHANGE,

	//The mesh has changed: cellsize, rectangle, or number of cells
	UPDATECONFIG_MESHCHANGE,

	//a new mesh has been added (all meshes added through the AddMesh method in supermesh so that method would signal this)
	UPDATECONFIG_MESHADDED,

	//a mesh has been deleted (all meshes deleted through the DelMesh method in supermesh so that method would signal this)
	UPDATECONFIG_MESHDELETED,

	////////////////////////////
	//PARAMETERS
	////////////////////////////

	//Param value changed
	UPDATECONFIG_PARAMVALUECHANGED,

	//mesh param settings changed
	UPDATECONFIG_PARAMCHANGED,

	////////////////////////////
	//ODE SOLVER
	////////////////////////////

	//Equation or evaluation method or settings changed
	UPDATECONFIG_ODE_SOLVER,

	//Moving mesh algorithm settings changed
	UPDATECONFIG_ODE_MOVEMESH,

	////////////////////////////
	//MODULES
	////////////////////////////

	//A module was added
	UPDATECONFIG_MODULEADDED,

	//A module was deleted
	UPDATECONFIG_MODULEDELETED,

	////////////////////////////
	//SPECIFIC MODULES
	////////////////////////////

	//SDemag or Demag module convolution type or settings change
	UPDATECONFIG_DEMAG_CONVCHANGE,

	//Change in roughness module
	UPDATECONFIG_ROUGHNESS_CHANGE,

	//Transport module electrode changed
	UPDATECONFIG_TRANSPORT_ELECTRODE,

	//Heat solver temperature model type changed
	UPDATECONFIG_HEAT_MODELTYPE,

	////////////////////////////
	//UpdateConfiguration_Values MESSAGES
	////////////////////////////

	//User constants changed for text equation objects
	UPDATECONFIG_TEQUATION_CONSTANTS,

	//Clear all text equations
	UPDATECONFIG_TEQUATION_CLEAR
};

namespace ucfg {

	//version without forcing flags
	template <typename Flag>
	bool __check_cfgflags(UPDATECONFIG_ cfgMessage, Flag flag)
	{
		if (cfgMessage == flag) return true;
		else return false;
	}

	//version without forcing flags
	template <typename Flag, typename ... Flags>
	bool __check_cfgflags(UPDATECONFIG_ cfgMessage, Flag flag, Flags... flags)
	{
		if (cfgMessage == flag) return true;
		else return __check_cfgflags(cfgMessage, flags...);
	}

	//version with forcing flags - use this
	template <typename Flag>
	bool check_cfgflags(UPDATECONFIG_ cfgMessage, Flag flag) 
	{ 
		if (cfgMessage == UPDATECONFIG_FORCEUPDATE ||
			cfgMessage == UPDATECONFIG_SWITCHCUDASTATE ||
			cfgMessage == UPDATECONFIG_REPAIROBJECTSTATE ||
			cfgMessage == UPDATECONFIG_MESHADDED ||
			cfgMessage == UPDATECONFIG_MESHDELETED ||
			cfgMessage == UPDATECONFIG_MODULEADDED ||
			cfgMessage == UPDATECONFIG_MODULEDELETED) return true;

		if (cfgMessage == flag) return true;
		else return false;
	}

	//version with forcing flags - use this
	template <typename Flag, typename ... Flags>
	bool check_cfgflags(UPDATECONFIG_ cfgMessage, Flag flag, Flags... flags)
	{
		if (cfgMessage == UPDATECONFIG_FORCEUPDATE ||
			cfgMessage == UPDATECONFIG_SWITCHCUDASTATE ||
			cfgMessage == UPDATECONFIG_REPAIROBJECTSTATE ||
			cfgMessage == UPDATECONFIG_MESHADDED ||
			cfgMessage == UPDATECONFIG_MESHDELETED ||
			cfgMessage == UPDATECONFIG_MODULEADDED ||
			cfgMessage == UPDATECONFIG_MODULEDELETED) return true;

		if (cfgMessage == flag) return true;
		else return __check_cfgflags(cfgMessage, flags...);
	}
};

#define CONVERSIONPRECISION 6						//Precision when converting to/from strings

#define MAXSIMSPACE		2.0							//Maximum side length of simulation space (m)
#define MINMESHSPACE	5e-11						//Minimum mesh side length (m). Also used to snap mesh rectangles to this resolution.
#define MAXFIELD		1e10						//Maximum field strength (A/m)
#define MAXSTRESS		1e15						//Maximum mechanical stress (Pa)
#define MINODERELERROR		1e-12					//Minimum relative error for ode solver that can be entered
#define MAXODERELERROR		1e-2					//Maximum relative error for ode solver that can be entered
#define MINTIMESTEP		1e-18						//Minimum time step that can be entered (s)
#define MAXTIMESTEP		1e-6						//Maximum time step that can be entered (s)

//when creating a mesh limit the number of cells. The mesh can exceed these values but user will have to adjust the cellsize manually.
#define MAXSTARTINGCELLS_X	2048
#define MAXSTARTINGCELLS_Y	2048
#define MAXSTARTINGCELLS_Z	512

//the default cellsize when creating a mesh (cubic) up to given number of maximum number of cells
#define DEFAULTCELLSIZE	5e-9

//the default atomistic cellsize when creating an atomistic mesh
#define DEFAULTATOMCELLSIZE 2e-10

//minimum and maximum damping values for fixed SOR damping algorithm
#define MINSORDAMPING	0.1
#define MAXSORDAMPING	2.0

//skyrmion definition settings for the skyrmion and skyrmionbloch commands
#define SKYRMION_RING_WIDTH 0.6
#define SKYRMION_TANH_CONST 2.0

