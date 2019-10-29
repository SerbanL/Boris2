#pragma once

//Quantities displayable on screen. Add new values at the end to keep older simulation files compatible.
enum MESHDISPLAY_ {
	MESHDISPLAY_NONE = 0, MESHDISPLAY_SM_DEMAG, MESHDISPLAY_SM_OERSTED, MESHDISPLAY_SM_STRAYH,
	MESHDISPLAY_MAGNETIZATION, MESHDISPLAY_EFFECTIVEFIELD,
	MESHDISPLAY_CURRDENSITY, MESHDISPLAY_VOLTAGE, MESHDISPLAY_ELCOND, MESHDISPLAY_SACCUM, MESHDISPLAY_JSX, MESHDISPLAY_JSY, MESHDISPLAY_JSZ, MESHDISPLAY_TS, MESHDISPLAY_TSI,
	MESHDISPLAY_TEMPERATURE, MESHDISPLAY_PARAMVAR, MESHDISPLAY_ROUGHNESS,
	MESHDISPLAY_MAGNETIZATION2, MESHDISPLAY_MAGNETIZATION12, MESHDISPLAY_EFFECTIVEFIELD2, MESHDISPLAY_EFFECTIVEFIELD12, MESHDISPLAY_CUSTOM_VEC, MESHDISPLAY_CUSTOM_SCA
};

//mesh types
enum MESH_ { MESH_SUPERMESH = 0, MESH_FERROMAGNETIC, MESH_DIPOLE, MESH_METAL, MESH_INSULATOR, MESH_ANTIFERROMAGNETIC };
