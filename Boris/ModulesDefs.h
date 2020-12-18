#pragma once



//Modules (MODS_ entries are super-mesh versions and not available as normal module handles
//Add new entries at the end to keep older simulation files compatible
//If you need to delete a module in the future you'll need to keep a dummy entry for it in this enum (can mark it as such, e.g. MOD_OBSOLETE1) although I can't see that occuring.
enum MOD_ {
	MOD_ALL = -1,
	MOD_ERROR = 0,
	MOD_DEMAG_N, MOD_DEMAG, MODS_SDEMAG,
	MOD_EXCHANGE, MOD_DMEXCHANGE, MOD_IDMEXCHANGE, MOD_SURFEXCHANGE,
	MOD_ZEEMAN,
	MOD_ANIUNI, MOD_ANICUBI,
	MOD_TRANSPORT, MODS_STRANSPORT, MODS_OERSTED,
	MODS_STRAYFIELD,
	MOD_HEAT, MODS_SHEAT,
	MOD_SOTFIELD,
	MOD_ROUGHNESS,
	MOD_SDEMAG_DEMAG,
	MOD_MELASTIC,
	MOD_MOPTICAL,
	MOD_ATOM_DIPOLEDIPOLE,

	MOD_STFIELD
};
