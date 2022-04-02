#pragma once

//Simulation data available as outputs - add new entries at the end to keep older simulation files compatible
enum DATA_
{
	DATA_NONE = 0,

	//Simulation schedule related data
	DATA_STAGESTEP = 1, DATA_TIME = 2, DATA_STAGETIME = 3, DATA_ITERATIONS = 4, DATA_SITERATIONS = 5, 
	DATA_DT = 6, DATA_MXH = 7, DATA_DMDT = 34,

	//Mesh quantities output, magnetic data
	DATA_AVM = 8, DATA_AVM2 = 36, DATA_HA = 9,
	DATA_MX_MINMAX = 41, DATA_MY_MINMAX = 42, DATA_MZ_MINMAX = 43, DATA_M_MINMAX = 44,
	DATA_THAVM = 59,
	
	DATA_RESPUMP = 50, DATA_IMSPUMP = 51, DATA_RESPUMP2 = 52, DATA_IMSPUMP2 = 53, DATA_RESPUMP12 = 63, DATA_RESPUMP21 = 64,

	//Transport quantities output
	DATA_JC = 10, DATA_JSX = 11, DATA_JSY = 12, DATA_JSZ = 13, DATA_V = 14, DATA_S = 15, DATA_ELC = 16, 
	
	//Transport settings data
	DATA_POTENTIAL = 17, DATA_CURRENT = 18, DATA_RESISTANCE = 19,

	//Transport solver data
	DATA_TRANSPORT_ITERSTOCONV = 28, DATA_TRANSPORT_SITERSTOCONV = 29, DATA_TRANSPORT_CONVERROR = 30,

	//Transport special
	DATA_TMR = 60,

	//Heat module output data
	DATA_TEMP = 31, DATA_TEMP_L = 38,

	DATA_HEATDT = 32,

	//Energies
	DATA_E_DEMAG = 20, DATA_E_EXCH = 21, DATA_E_SURFEXCH = 22, DATA_E_ZEE = 23, DATA_E_STRAY = 61, DATA_E_ANIS = 24, DATA_E_ROUGH = 25, DATA_E_MOPTICAL = 49, DATA_E_MELASTIC = 37, DATA_E_TOTAL = 33,
	
	//Torques
	DATA_T_EXCH = 55, DATA_T_SURFEXCH = 54, DATA_T_ZEE = 56, DATA_T_STRAY = 62, DATA_T_ANIS = 57,

	//Algorithms data
	DATA_DWSHIFT = 26, DATA_SKYSHIFT = 27, 
	DATA_DWPOS_X = 45, DATA_DWPOS_Y = 46, DATA_DWPOS_Z = 47,
	DATA_SKYPOS = 35, DATA_Q_TOPO = 40,
	DATA_MONTECARLOPARAMS = 48,

	//Special
	DATA_COMMBUFFER = 58,

	//Previously used by DATA_E_EXCH_MAX, now deleted
	DATA_RESERVED = 39
};
//Current maximum : 64
