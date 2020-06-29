#pragma once

//heat-flow temperature model : 

//1-temperature model, 2-temperature model

//default setting is a 1-temperature model

enum TMTYPE_ {

	TMTYPE_DEFAULT = -1,
	TMTYPE_NONE = 0,
	TMTYPE_1TM,
	TMTYPE_2TM,
	TMTYPE_NUMMODELS
};
