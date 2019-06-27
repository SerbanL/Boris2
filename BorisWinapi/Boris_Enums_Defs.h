#pragma once

#include "CompileFlags.h"

//All widely applicable enums and defines go here. Further smaller scope enums and defines are given where needed.

//when UpdateConfiguration method is called you can pass a value from this enum to signal the reason it was called - some modules need to know (e.g. SDemag needs to know if shape has changed)
enum UPDATECONFIG_ {

	//no reason given - Update Configuration without consideration for special cases
	UPDATECONFIG_GENERIC,

	//mesh shape has changed (not the rectangle but the shape inside the rectangle)
	UPDATECONFIG_MESHSHAPECHANGE,

	//use this message to force updates by calling the required UpdateConfig method directly in the required object (which must check for it)
	UPDATECONFIG_FORCEUPDATE
};

#if COMPILECUDA == 1
#if SINGLEPRECISION == 1

	#define cuReal cufftReal
	#define cuComplex cufftComplex

#else

	#define cuReal cufftDoubleReal
	#define cuComplex cufftDoubleComplex

#endif
#endif

#define CONVERSIONPRECISION 6						//Precision when converting to/from strings

#define MAXSIMSPACE		2.0							//Maximum side length of simulation space (m)
#define MINMESHSPACE	1e-11						//Minimum mesh side length (m)
#define MAXFIELD		1e8							//Maximum field strength (A/m)
#define MINTIMESTEP		1e-18						//Mnimum time step that can be entered (s)
#define MAXTIMESTEP		1e-11						//Maximum time step that can be entered (s)

//when creating a mesh limit the number of cells. The mesh can exceed these values but user will have to adjust the cellsize manually.
#define MAXSTARTINGCELLS_X	2048
#define MAXSTARTINGCELLS_Y	2048
#define MAXSTARTINGCELLS_Z	64

//the default cellsize when creating a mesh (cubic) up to given number of maximum number of cells
#define DEFAULTCELLSIZE	5e-9

//minimum and maximum damping values for fixed SOR damping algorithm
#define MINSORDAMPING	0.1
#define MAXSORDAMPING	2.0

//skyrmion definition settings for the skyrmion and skyrmionbloch commands
#define SKYRMION_RING_WIDTH 0.6
#define SKYRMION_TANH_CONST 2.0

