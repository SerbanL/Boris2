#pragma once

//Modules (MODS_ entries are super-mesh versions and not available as normal module handles
//Add new entries at the end to keep older simulation files compatible (or rather keep the integer values already defined, so can re-arrange inside enum)
//If you need to delete a module in the future you'll need to keep a dummy entry for it in this enum (can mark it as such, e.g. MOD_OBSOLETE1) although I can't see that occuring.
enum MOD_ {
	MOD_ALL = -1, MOD_ERROR = 0,

	//demag
	MOD_DEMAG_N = 1, MOD_DEMAG = 2, MODS_SDEMAG = 3, MOD_SDEMAG_DEMAG = 19,

	//exchange
	MOD_EXCHANGE = 4, MOD_DMEXCHANGE = 5, MOD_IDMEXCHANGE = 6, MOD_SURFEXCHANGE = 7,
	
	//others
	MOD_ZEEMAN = 8, MOD_ROUGHNESS = 18, MOD_MELASTIC = 20, MOD_MOPTICAL = 21,
	
	//anisotropy
	MOD_ANIUNI = 9, MOD_ANICUBI = 10, MOD_ANIBI = 24, MOD_ANITENS = 25,
	
	//transport or transport-related
	MOD_TRANSPORT = 11, MODS_STRANSPORT = 12, MODS_OERSTED = 13,
	
	//stray field
	MODS_STRAYFIELD = 14,
	
	//heat
	MOD_HEAT = 15, MODS_SHEAT = 16,
	
	//spin torque
	MOD_SOTFIELD = 17, MOD_STFIELD = 23,
	
	//atomistic dipole-dipole
	MOD_ATOM_DIPOLEDIPOLE = 22
}; //highest integer : 25
