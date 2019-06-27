#pragma once

#include "Simulation.h"

#include "ColorDefs.h"

void Simulation::MakeIOInfo(void)
{
	string versionupdate_info =
		string("[tc1,1,0,1/tc]<b>Program Version Status</b>") +
		string("\n[tc1,1,0,1/tc]click: open download page\n");

	ioInfo.push_back(versionupdate_info, IOI_PROGRAMUPDATESTATUS);

	//Data box entry, showing the label of a given entry in Simulation::dataBoxList : minorId is the minor id of elements in Simulation::dataBoxList (major id there is always 0), auxId is the number of the interactive object in the list (i.e. entry number as it appears in data box in order). textId is the mesh name (if associated with this data type)
	//Note this entry must always represent the entry in Simulation::dataBoxList with the index in auxId.
	//IOI_DATABOXFIELDLABEL

	//A set or available module for a given mesh: minorId in InteractiveObjectProperties is an entry from MOD_ enum identifying the module, auxId contains the unique mesh id number this module refers to
	//IOI_MODULE
	//super-mesh module : minor type is an entry from MOD_ enum
	//IOI_SMODULE

	string modulegeneric_info =
		string("[tc1,1,0,1/tc]<b>computational module</b>") +
		string("\n[tc1,1,0,1/tc]<i>green: on, red: off</i>") +
		string("\n[tc1,1,0,1/tc]left-click: enable") +
		string("\n[tc1,1,0,1/tc]right-click: disable\n");

	ioInfo.set(modulegeneric_info + string("<i><b>Stoner-Wohlfarth demag"), INT2(IOI_MODULE, MOD_DEMAG_N));
	ioInfo.set(modulegeneric_info + string("<i><b>Full demag field"), INT2(IOI_MODULE, MOD_DEMAG));
	ioInfo.set(modulegeneric_info + string("<i><b>Direct exchange interaction"), INT2(IOI_MODULE, MOD_EXCHANGE6NGBR));
	ioInfo.set(modulegeneric_info + string("<i><b>Dzyaloshinskii-Moriya interaction - bulk"), INT2(IOI_MODULE, MOD_DMEXCHANGE));
	ioInfo.set(modulegeneric_info + string("<i><b>Dzyaloshinskii-Moriya interaction - interfacial"), INT2(IOI_MODULE, MOD_IDMEXCHANGE));
	ioInfo.set(modulegeneric_info + string("<i><b>Surface exchange interaction\n<i><b>Couple to adjacent meshes along z\n<i><b>with surfexchange module enabled"), INT2(IOI_MODULE, MOD_SURFEXCHANGE));
	ioInfo.set(modulegeneric_info + string("<i><b>Applied field term"), INT2(IOI_MODULE, MOD_ZEEMAN));
	ioInfo.set(modulegeneric_info + string("<i><b>Magnetocrystalline anisotropy: uniaxial"), INT2(IOI_MODULE, MOD_ANIUNI));
	ioInfo.set(modulegeneric_info + string("<i><b>Magnetocrystalline anisotropy: cubic"), INT2(IOI_MODULE, MOD_ANICUBI));
	ioInfo.set(modulegeneric_info + string("<i><b>Charge and spin transport"), INT2(IOI_MODULE, MOD_TRANSPORT));
	ioInfo.set(modulegeneric_info + string("<i><b>Heat equation solver"), INT2(IOI_MODULE, MOD_HEAT));
	ioInfo.set(modulegeneric_info + string("<i><b>Spin-orbit torque field\n<i><b>Results in Slonczewski-like torques"), INT2(IOI_MODULE, MOD_SOTFIELD));
	ioInfo.set(modulegeneric_info + string("<i><b>Physical roughness\n<i><b>Demag term corrections when approximating shapes\n<i><b>from a fine mesh with a coarse mesh."), INT2(IOI_MODULE, MOD_ROUGHNESS));

	ioInfo.set(modulegeneric_info + string("<i><b>Supermesh demag field"), INT2(IOI_SMODULE, MODS_SDEMAG));
	ioInfo.set(modulegeneric_info + string("<i><b>Oersted field in electric supermesh"), INT2(IOI_SMODULE, MODS_OERSTED));
	ioInfo.set(modulegeneric_info + string("<i><b>Stray field from dipole meshes"), INT2(IOI_SMODULE, MODS_STRAYFIELD));

	//Available/set ode and evaluation method for magnetization : minorId is an entry from ODE_ (the equation), auxId is the EVAL_ entry (the evaluation method), textId is the name of the evaluation method
	//IOI_ODE

	string ode_info =
		string("[tc1,1,0,1/tc]<b>Evaluation method for</b>\n[tc1,1,0,1/tc]<b>dM/dt equation</b>") +
		string("\n[tc1,1,0,1/tc]<i>green: on, red: off</i>") +
		string("\n[tc1,1,0,1/tc]click: change state\n");

	ioInfo.push_back(ode_info, IOI_ODE);

	//Shows a mesh name : minorId is the unique mesh id number, textId is the mesh name (below are similar objects but used in different lists, so these lists need updating differently).
	//auxId is also used : value 1 means update list, value 0 means do not update list (but delete line if mesh is deleted).
	//IOI_MESH_FORPARAMS
	//IOI_MESH_FORPARAMSTEMP
	//IOI_MESH_FORPARAMSVAR
	//IOI_MESH_FORMODULES
	//IOI_MESH_FORMESHLIST
	//IOI_MESH_FORDISPLAYOPTIONS
	//IOI_MESH_FORTEMPERATURE
	//IOI_MESH_FORHEATBOUNDARIES
	//IOI_MESH_FORCURIEANDMOMENT

	string mesh_info =
		string("[tc1,1,0,1/tc]<b>mesh name</b>") +
		string("\n[tc1,1,0,1/tc]<i>green: selected, red: not selected</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit") +
		string("\n[tc1,1,0,1/tc]left-click: select") +
		string("\n[tc1,1,0,1/tc]right-click: delete");

	ioInfo.push_back(mesh_info, IOI_MESH_FORPARAMS);
	ioInfo.push_back(mesh_info, IOI_MESH_FORPARAMSTEMP);
	ioInfo.push_back(mesh_info, IOI_MESH_FORPARAMSVAR);
	ioInfo.push_back(mesh_info, IOI_MESH_FORMODULES);
	ioInfo.push_back(mesh_info, IOI_MESH_FORMESHLIST);
	ioInfo.push_back(mesh_info, IOI_MESH_FORDISPLAYOPTIONS);
	ioInfo.push_back(mesh_info, IOI_MESH_FORTEMPERATURE);
	ioInfo.push_back(mesh_info, IOI_MESH_FORHEATBOUNDARIES);
	ioInfo.push_back(mesh_info, IOI_MESH_FORCURIEANDMOMENT);

	//Shows ferromagnetic super-mesh rectangle (unit m) : textId is the mesh rectangle for the ferromagnetic super-mesh
	//IOI_FMSMESHRECTANGLE

	string fsmeshrect_info =
		string("[tc1,1,0,1/tc]<b>Ferromagnetic supermesh rectangle</b>") +
		string("\n[tc1,1,0,1/tc]<i>Smallest rectangle containing</i>\n[tc1,1,0,1/tc]<i>all ferromagnetic meshes</i>");

	ioInfo.push_back(fsmeshrect_info, IOI_FMSMESHRECTANGLE);

	//Shows electric super-mesh rectangle (unit m) : textId is the mesh rectangle for the electric super-mesh
	//IOI_ESMESHRECTANGLE

	string esmeshrect_info =
		string("[tc1,1,0,1/tc]<b>Electric supermesh rectangle</b>") +
		string("\n[tc1,1,0,1/tc]<i>Smallest rectangle containing</i>\n[tc1,1,0,1/tc]<i>all meshes with</i>\n[tc1,1,0,1/tc]<i>enabled Transport module</i>");

	ioInfo.push_back(esmeshrect_info, IOI_ESMESHRECTANGLE);

	//Shows ferromagnetic super-mesh cellsize (units m) : textId is the mesh cellsize for the ferromagnetic super-mesh
	//IOI_FMSMESHCELLSIZE

	string fsmeshcell_info =
		string("[tc1,1,0,1/tc]<b>Ferromagnetic supermesh cellsize</b>") +
		string("\n[tc1,1,0,1/tc]<i>Discretization cellsize</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(fsmeshcell_info, IOI_FMSMESHCELLSIZE);

	//Shows electric super-mesh cellsize (units m) : textId is the mesh cellsize for the electric super-mesh
	//IOI_ESMESHCELLSIZE

	string esmeshcell_info =
		string("[tc1,1,0,1/tc]<b>Electric supermesh cellsize</b>") +
		string("\n[tc1,1,0,1/tc]<i>Discretization cellsize</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(esmeshcell_info, IOI_ESMESHCELLSIZE);

	//Shows mesh rectangle (units m) : minorId is the unique mesh id number, textId is the mesh rectangle
	//IOI_MESHRECTANGLE

	string meshrect_info =
		string("[tc1,1,0,1/tc]<b>Mesh rectangle</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(meshrect_info, IOI_MESHRECTANGLE);

	//Shows magnetic mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	//IOI_MESHCELLSIZE

	string meshcell_info =
		string("[tc1,1,0,1/tc]<b>Ferromagnetic cellsize</b>") +
		string("\n[tc1,1,0,1/tc]<i>Discretization cellsize</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(meshcell_info, IOI_MESHCELLSIZE);

	//Shows electric mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	//IOI_MESHECELLSIZE

	string emeshcell_info =
		string("[tc1,1,0,1/tc]<b>Transport solver cellsize</b>") +
		string("\n[tc1,1,0,1/tc]<i>Discretization cellsize</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(emeshcell_info, IOI_MESHECELLSIZE);

	//Shows thermal mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	//IOI_MESHTCELLSIZE

	string tmeshcell_info =
		string("[tc1,1,0,1/tc]<b>Heat solver cellsize</b>") +
		string("\n[tc1,1,0,1/tc]<i>Discretization cellsize</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(tmeshcell_info, IOI_MESHTCELLSIZE);

	//Simulation output data, specifically used for showing values in console : minorId is the DATA_ id, textId is the data handle
	//IOI_SHOWDATA
	//Simulation output data, specifically used to construct output data list : minorId is the DATA_ id, textId is the data handle
	//IOI_DATA

	string showdata_info_generic =
		string("[tc1,1,0,1/tc]<b>Output data</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: show value") +
		string("\n[tc1,1,0,1/tc]right-click: pin to data box") +
		string("\n[tc1,1,0,1/tc]drag: pin to data box\n");

	ioInfo.set(showdata_info_generic + string("<i><b>Simulation stage and step</i>"), INT2(IOI_SHOWDATA, DATA_STAGESTEP));
	ioInfo.set(showdata_info_generic + string("<i><b>Simulation total time</i>"), INT2(IOI_SHOWDATA, DATA_TIME));
	ioInfo.set(showdata_info_generic + string("<i><b>Simulation stage time</i>"), INT2(IOI_SHOWDATA, DATA_STAGETIME));
	ioInfo.set(showdata_info_generic + string("<i><b>Simulation iterations</i>"), INT2(IOI_SHOWDATA, DATA_ITERATIONS));
	ioInfo.set(showdata_info_generic + string("<i><b>Simulation stage iterations</i>"), INT2(IOI_SHOWDATA, DATA_SITERATIONS));
	ioInfo.set(showdata_info_generic + string("<i><b>Magnetisation solver time step</i>"), INT2(IOI_SHOWDATA, DATA_DT));
	ioInfo.set(showdata_info_generic + string("<i><b>Magnetisation relaxation</i>"), INT2(IOI_SHOWDATA, DATA_MXH));
	ioInfo.set(showdata_info_generic + string("<i><b>Average magnetisation</i>"), INT2(IOI_SHOWDATA, DATA_AVM));
	ioInfo.set(showdata_info_generic + string("<i><b>Applied magnetic field</i>"), INT2(IOI_SHOWDATA, DATA_HA));
	ioInfo.set(showdata_info_generic + string("<i><b>Average charge current density</i>"), INT2(IOI_SHOWDATA, DATA_JC));
	ioInfo.set(showdata_info_generic + string("<i><b>Average spin x-current density</i>"), INT2(IOI_SHOWDATA, DATA_JSX));
	ioInfo.set(showdata_info_generic + string("<i><b>Average spin y-current density</i>"), INT2(IOI_SHOWDATA, DATA_JSY));
	ioInfo.set(showdata_info_generic + string("<i><b>Average spin z-current density</i>"), INT2(IOI_SHOWDATA, DATA_JSZ));
	ioInfo.set(showdata_info_generic + string("<i><b>Average charge potential</i>"), INT2(IOI_SHOWDATA, DATA_V));
	ioInfo.set(showdata_info_generic + string("<i><b>Average spin accumulation</i>"), INT2(IOI_SHOWDATA, DATA_S));
	ioInfo.set(showdata_info_generic + string("<i><b>Average electrical conductivity</i>"), INT2(IOI_SHOWDATA, DATA_ELC));
	ioInfo.set(showdata_info_generic + string("<i><b>Potential drop between electrodes</i>"), INT2(IOI_SHOWDATA, DATA_POTENTIAL));
	ioInfo.set(showdata_info_generic + string("<i><b>Charge current into\n<i><b>ground electrode and\n<i>net current (error term)</i>"), INT2(IOI_SHOWDATA, DATA_CURRENT));
	ioInfo.set(showdata_info_generic + string("<i><b>Total resistance (V/I)</i>"), INT2(IOI_SHOWDATA, DATA_RESISTANCE));
	ioInfo.set(showdata_info_generic + string("<i><b>Energy density: demag</i>"), INT2(IOI_SHOWDATA, DATA_E_DEMAG));
	ioInfo.set(showdata_info_generic + string("<i><b>Energy density: exchange</i>"), INT2(IOI_SHOWDATA, DATA_E_EXCH));
	ioInfo.set(showdata_info_generic + string("<i><b>Energy density: surface exchange</i>"), INT2(IOI_SHOWDATA, DATA_E_SURFEXCH));
	ioInfo.set(showdata_info_generic + string("<i><b>Energy density: applied H field</i>"), INT2(IOI_SHOWDATA, DATA_E_ZEE));
	ioInfo.set(showdata_info_generic + string("<i><b>Energy density: anisotropy</i>"), INT2(IOI_SHOWDATA, DATA_E_ANIS));
	ioInfo.set(showdata_info_generic + string("<i><b>Energy density: roughness</i>"), INT2(IOI_SHOWDATA, DATA_E_ROUGH));
	ioInfo.set(showdata_info_generic + string("<i><b>Domain wall shift\n<i><b>for moving mesh</i>"), INT2(IOI_SHOWDATA, DATA_DWSHIFT));
	ioInfo.set(showdata_info_generic + string("<i><b>Skyrmion shift in the xy plan\n<i><b>Only use with output save data</i>"), INT2(IOI_SHOWDATA, DATA_SKYSHIFT));
	ioInfo.set(showdata_info_generic + string("<i><b>Transport solver:\n<i><b>V iterations to convergence</i>"), INT2(IOI_SHOWDATA, DATA_TRANSPORT_ITERSTOCONV));
	ioInfo.set(showdata_info_generic + string("<i><b>Transport solver:\n<i><b>S iterations to convergence</i>"), INT2(IOI_SHOWDATA, DATA_TRANSPORT_SITERSTOCONV));
	ioInfo.set(showdata_info_generic + string("<i><b>Transport solver:\n<i><b>achieved convergence error</i>"), INT2(IOI_SHOWDATA, DATA_TRANSPORT_CONVERROR));
	ioInfo.set(showdata_info_generic + string("<i><b>Transport solver:\n<i><b>min and max aSOR damping</i>"), INT2(IOI_SHOWDATA, DATA_TRANSPORT_ASOR));
	ioInfo.set(showdata_info_generic + string("<i><b>Average temperature</i>"), INT2(IOI_SHOWDATA, DATA_TEMP));
	ioInfo.set(showdata_info_generic + string("<i><b>Heat solver time step</i>"), INT2(IOI_SHOWDATA, DATA_HEATDT));

	string data_info_generic =
		string("[tc1,1,0,1/tc]<b>Output data</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: add to output list") +
		string("\n[tc1,1,0,1/tc]right-click: pin to data box") +
		string("\n[tc1,1,0,1/tc]drag: pin to data box\n");

	ioInfo.set(data_info_generic + string("<i><b>Simulation stage and step</i>"), INT2(IOI_DATA, DATA_STAGESTEP));
	ioInfo.set(data_info_generic + string("<i><b>Simulation total time</i>"), INT2(IOI_DATA, DATA_TIME));
	ioInfo.set(data_info_generic + string("<i><b>Simulation stage time</i>"), INT2(IOI_DATA, DATA_STAGETIME));
	ioInfo.set(data_info_generic + string("<i><b>Simulation iterations</i>"), INT2(IOI_DATA, DATA_ITERATIONS));
	ioInfo.set(data_info_generic + string("<i><b>Simulation stage iterations</i>"), INT2(IOI_DATA, DATA_SITERATIONS));
	ioInfo.set(data_info_generic + string("<i><b>Magnetisation solver time step</i>"), INT2(IOI_DATA, DATA_DT));
	ioInfo.set(data_info_generic + string("<i><b>Magnetisation relaxation</i>"), INT2(IOI_DATA, DATA_MXH));
	ioInfo.set(data_info_generic + string("<i><b>Average magnetisation</i>"), INT2(IOI_DATA, DATA_AVM));
	ioInfo.set(data_info_generic + string("<i><b>Applied magnetic field</i>"), INT2(IOI_DATA, DATA_HA));
	ioInfo.set(data_info_generic + string("<i><b>Average charge current density</i>"), INT2(IOI_DATA, DATA_JC));
	ioInfo.set(data_info_generic + string("<i><b>Average spin x-current density</i>"), INT2(IOI_DATA, DATA_JSX));
	ioInfo.set(data_info_generic + string("<i><b>Average spin y-current density</i>"), INT2(IOI_DATA, DATA_JSY));
	ioInfo.set(data_info_generic + string("<i><b>Average spin z-current density</i>"), INT2(IOI_DATA, DATA_JSZ));
	ioInfo.set(data_info_generic + string("<i><b>Average charge potential</i>"), INT2(IOI_DATA, DATA_V));
	ioInfo.set(data_info_generic + string("<i><b>Average spin accumulation</i>"), INT2(IOI_DATA, DATA_S));
	ioInfo.set(data_info_generic + string("<i><b>Average electrical conductivity</i>"), INT2(IOI_DATA, DATA_ELC));
	ioInfo.set(data_info_generic + string("<i><b>Potential drop between electrodes</i>"), INT2(IOI_DATA, DATA_POTENTIAL));
	ioInfo.set(data_info_generic + string("<i><b>Charge current into\n<i><b>ground electrode and\n<i><b>net current (error term)</i>"), INT2(IOI_DATA, DATA_CURRENT));
	ioInfo.set(data_info_generic + string("<i><b>Total resistance (V/I)</i>"), INT2(IOI_DATA, DATA_RESISTANCE));
	ioInfo.set(data_info_generic + string("<i><b>Energy density: demag</i>"), INT2(IOI_DATA, DATA_E_DEMAG));
	ioInfo.set(data_info_generic + string("<i><b>Energy density: exchange</i>"), INT2(IOI_DATA, DATA_E_EXCH));
	ioInfo.set(data_info_generic + string("<i><b>Energy density: surface exchange</i>"), INT2(IOI_DATA, DATA_E_SURFEXCH));
	ioInfo.set(data_info_generic + string("<i><b>Energy density: applied H field</i>"), INT2(IOI_DATA, DATA_E_ZEE));
	ioInfo.set(data_info_generic + string("<i><b>Energy density: anisotropy</i>"), INT2(IOI_DATA, DATA_E_ANIS));
	ioInfo.set(data_info_generic + string("<i><b>Energy density: roughness</i>"), INT2(IOI_DATA, DATA_E_ROUGH));
	ioInfo.set(data_info_generic + string("<i><b>Domain wall shift\n<i><b>for moving mesh</i>"), INT2(IOI_DATA, DATA_DWSHIFT));
	ioInfo.set(data_info_generic + string("<i><b>Skyrmion shift in the xy plane\n<i><b>Rectangle must circumscribe skyrmion</i>"), INT2(IOI_DATA, DATA_SKYSHIFT));
	ioInfo.set(data_info_generic + string("<i><b>Transport solver:\n<i><b>V iterations to convergence</i>"), INT2(IOI_DATA, DATA_TRANSPORT_ITERSTOCONV));
	ioInfo.set(data_info_generic + string("<i><b>Transport solver:\n<i><b>S iterations to convergence</i>"), INT2(IOI_DATA, DATA_TRANSPORT_SITERSTOCONV));
	ioInfo.set(data_info_generic + string("<i><b>Transport solver:\n<i><b>achieved convergence error</i>"), INT2(IOI_DATA, DATA_TRANSPORT_CONVERROR));
	ioInfo.set(data_info_generic + string("<i><b>Transport solver:\n<i><b>min and max aSOR damping</i>"), INT2(IOI_DATA, DATA_TRANSPORT_ASOR));
	ioInfo.set(data_info_generic + string("<i><b>Average temperature</i>"), INT2(IOI_DATA, DATA_TEMP));
	ioInfo.set(data_info_generic + string("<i><b>Heat solver time step</i>"), INT2(IOI_DATA, DATA_HEATDT));

	//Show currently set directory : textId is the directory
	//IOI_DIRECTORY

	string dir_info =
		string("[tc1,1,0,1/tc]<b>Working directory</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(dir_info, IOI_DIRECTORY);

	//Show currently set save data file : textId is the file name
	//IOI_SAVEDATAFILE

	string savedatafile_info =
		string("[tc1,1,0,1/tc]<b>Output data file (.txt)</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(savedatafile_info, IOI_SAVEDATAFILE);

	//Show currently set image filename base : textId is the file name
	//IOI_SAVEIMAGEFILEBASE

	string saveimagefile_info =
		string("[tc1,1,0,1/tc]<b>Image save file base</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(saveimagefile_info, IOI_SAVEIMAGEFILEBASE);

	//Show flag status for data/image saving during a simulation : minorId is the flag value (boolean)
	//IOI_SAVEDATAFLAG

	string savedataflag_info =
		string("[tc1,1,0,1/tc]<b>Data saving switch</b>") +
		string("\n[tc1,1,0,1/tc]<i>green: on, red: off</i>") +
		string("\n[tc1,1,0,1/tc]click: change status");

	ioInfo.push_back(savedataflag_info, IOI_SAVEDATAFLAG);

	//IOI_SAVEIMAGEFLAG

	string saveimageflag_info =
		string("[tc1,1,0,1/tc]<b>Image saving switch</b>") +
		string("\n[tc1,1,0,1/tc]<i>green: on, red: off</i>") +
		string("\n[tc1,1,0,1/tc]click: change status");

	ioInfo.push_back(saveimageflag_info, IOI_SAVEIMAGEFLAG);

	//Show set output data : minorId is the minor id of elements in Simulation::saveDataList (major id there is always 0), auxId is the number of the interactive object in the list as it appears in the console, textId is the configured output data. 
	//Note this entry must always represent the entry in Simulation::saveDataList with the index in auxId.
	//IOI_OUTDATA

	string showoutdata_info =
		string("[tc1,1,0,1/tc]<b>Set output data entry</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit entry") +
		string("\n[tc1,1,0,1/tc]right-click: delete entry") +
		string("\n[tc1,1,0,1/tc]drag: move to change order\n");

	ioInfo.push_back(showoutdata_info, IOI_OUTDATA);

	//Shows a possible stage type, used for adding generic stages to the simulation schedule : minorId is the stage type (SS_ enum value, which is the majorId from stageDescriptors), textId is the stage setting handle
	//IOI_STAGE

	string stage_generic_info =
		string("[tc1,1,0,1/tc]<b>Simulation schedule stage type</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: add to schedule\n");

	ioInfo.set(stage_generic_info + string("<i><b>Nothing to set, just run solvers"), INT2(IOI_STAGE, SS_RELAX));
	ioInfo.set(stage_generic_info + string("<i><b>Set H field in Cartesian coordinates"), INT2(IOI_STAGE, SS_HFIELDXYZ));
	ioInfo.set(stage_generic_info + string("<i><b>Set H field sequence in Cartesian coordinates\n<i><b>Start to stop H values in a number of steps"), INT2(IOI_STAGE, SS_HFIELDXYZSEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Set H field sequence in Polar coordinates\n<i><b>Start to stop H values in a number of steps"), INT2(IOI_STAGE, SS_HPOLARSEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Set FMR H field sequence in Cartesian coordinates\n<i><b>Bias field, rf field, rf field steps, rf field cycles"), INT2(IOI_STAGE, SS_HFMR));
	ioInfo.set(stage_generic_info + string("<i><b>Set fixed voltage drop between electrodes"), INT2(IOI_STAGE, SS_V));
	ioInfo.set(stage_generic_info + string("<i><b>Set fixed voltage sequence\n<i><b>Start to stop V values in a number of steps"), INT2(IOI_STAGE, SS_VSEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Set a sinusoidal voltage sequence\n<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_STAGE, SS_VSIN));
	ioInfo.set(stage_generic_info + string("<i><b>Set a cosine voltage sequence\n<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_STAGE, SS_VCOS));
	ioInfo.set(stage_generic_info + string("<i><b>Set fixed current into ground electrode"), INT2(IOI_STAGE, SS_I));
	ioInfo.set(stage_generic_info + string("<i><b>Set fixed current sequence\n<i><b>Start to stop I values in a number of steps"), INT2(IOI_STAGE, SS_ISEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Set a sinusoidal current sequence\n<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_STAGE, SS_ISIN));
	ioInfo.set(stage_generic_info + string("<i><b>Set a cosine current sequence\n<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_STAGE, SS_ICOS));
	ioInfo.set(stage_generic_info + string("<i><b>Set base temperature value"), INT2(IOI_STAGE, SS_T));
	ioInfo.set(stage_generic_info + string("<i><b>Set base temperature sequence\n<i><b>Start to stop T values in a number of steps"), INT2(IOI_STAGE, SS_TSEQ));

	//Shows a stage added to the simulation schedule : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the configured stage text
	//Note this entry must always represent the entry in Simulation::simStages with the index in auxId.
	//IOI_SETSTAGE

	string setstage_info =
		string("[tc1,1,0,1/tc]<b>Set simulation schedule stage</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n") +
		string("\n[tc1,1,0,1/tc]right-click: delete\n") +
		string("\n[tc1,1,0,1/tc]drag: move to change order\n");

	ioInfo.push_back(setstage_info, IOI_SETSTAGE);

	//Shows the value to set for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the value as a string
	//IOI_SETSTAGEVALUE

	string stagevalue_generic_info =
		string("[tc1,1,0,1/tc]<b>Simulation schedule stage value</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.set(stagevalue_generic_info + string(""), INT2(IOI_SETSTAGEVALUE, SS_RELAX));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Hx, Hy, Hz"), INT2(IOI_SETSTAGEVALUE, SS_HFIELDXYZ));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Hstart; Hstop; Steps: Hstep = (Hstop - Hstart) / Steps"), INT2(IOI_SETSTAGEVALUE, SS_HFIELDXYZSEQ));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Hstart; Hstop; Steps: Hstep = (Hstop - Hstart) / Steps\n<i><b>H values in polar coordinates as:\n<i><b>strength value, polar angle, azimuthal angle"), INT2(IOI_SETSTAGEVALUE, SS_HPOLARSEQ));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Hbias; Hrf; rf steps; rf cycles\n<i><b>rf steps is the rf cycle discretization\n<i><b>rf cycles is the number of periods"), INT2(IOI_SETSTAGEVALUE, SS_HFMR));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Fixed voltage drop between electrodes"), INT2(IOI_SETSTAGEVALUE, SS_V));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Vstart; Vstop; Steps: Vstep = (Vstop - Vstart) / Steps"), INT2(IOI_SETSTAGEVALUE, SS_VSEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_SETSTAGEVALUE, SS_VSIN));
	ioInfo.set(stage_generic_info + string("<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_SETSTAGEVALUE, SS_VCOS));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Fixed current into ground electrode"), INT2(IOI_SETSTAGEVALUE, SS_I));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Istart; Istop; Steps: Istep = (Istop - Istart) / Steps"), INT2(IOI_SETSTAGEVALUE, SS_ISEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_SETSTAGEVALUE, SS_ISIN));
	ioInfo.set(stage_generic_info + string("<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_SETSTAGEVALUE, SS_ICOS));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Base temperature value"), INT2(IOI_SETSTAGEVALUE, SS_T));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Tstart; Tstop; Steps: Tstep = (Tstop - Tstart) / Steps"), INT2(IOI_SETSTAGEVALUE, SS_TSEQ));

	//Shows the stop condition for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the stop type and value as a string
	//IOI_STAGESTOPCONDITION

	string stagestop_generic_info =
		string("[tc1,1,0,1/tc]<b>Stage stop condition</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.set(stagestop_generic_info + string("<i><b>No stopping condition set"), INT2(IOI_STAGESTOPCONDITION, STOP_NOSTOP));
	ioInfo.set(stagestop_generic_info + string("<i><b>Stop when stage iterations value reached"), INT2(IOI_STAGESTOPCONDITION, STOP_ITERATIONS));
	ioInfo.set(stagestop_generic_info + string("<i><b>Stop when |mxh| falls below value"), INT2(IOI_STAGESTOPCONDITION, STOP_MXH));
	ioInfo.set(stagestop_generic_info + string("<i><b>Stop when stage time value reached"), INT2(IOI_STAGESTOPCONDITION, STOP_TIME));

	//Shows the saving condition for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the DSAVE_ value for this data save type, textId is the save type and value as a string
	//IOI_DSAVETYPE

	string stagesave_generic_info =
		string("[tc1,1,0,1/tc]<b>Data save condition</b>") +
		string("\n[tc1,1,0,1/tc]<i>green: on, red: off</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.set(stagesave_generic_info + string("<i><b>No saving set"), INT2(IOI_DSAVETYPE, DSAVE_NONE));
	ioInfo.set(stagesave_generic_info + string("<i><b>Save at the end of the stage"), INT2(IOI_DSAVETYPE, DSAVE_STAGE));
	ioInfo.set(stagesave_generic_info + string("<i><b>Save at the end of every step\n<i><b>in this stage"), INT2(IOI_DSAVETYPE, DSAVE_STEP));
	ioInfo.set(stagesave_generic_info + string("<i><b>Save after every given iterations"), INT2(IOI_DSAVETYPE, DSAVE_ITER));
	ioInfo.set(stagesave_generic_info + string("<i><b>Save after every given time interval"), INT2(IOI_DSAVETYPE, DSAVE_TIME));

	//Shows a stop condition, used to apply the same condition to all simulation stages : minorId is the STOP_ value, textId is the stop type handle
	//IOI_STAGESTOPCONDITIONALL

	string stagestopall_generic_info =
		string("[tc1,1,0,1/tc]<b>Stage stop condition\n[tc1,1,0,1/tc]<b>for all stages</b>") +
		string("\n[tc1,1,0,1/tc]click: set\n");

	ioInfo.set(stagestopall_generic_info + string("<i><b>No stopping condition set"), INT2(IOI_STAGESTOPCONDITIONALL, STOP_NOSTOP));
	ioInfo.set(stagestopall_generic_info + string("<i><b>Stop when stage iterations value reached"), INT2(IOI_STAGESTOPCONDITIONALL, STOP_ITERATIONS));
	ioInfo.set(stagestopall_generic_info + string("<i><b>Stop when |mxh| falls below value"), INT2(IOI_STAGESTOPCONDITIONALL, STOP_MXH));
	ioInfo.set(stagestopall_generic_info + string("<i><b>Stop when stage time value reached"), INT2(IOI_STAGESTOPCONDITIONALL, STOP_TIME));

	//Shows a data save condition, used to apply the same condition to all simulation stages : minorId is the DSAVE_ value, textId is the save type handle
	//IOI_DSAVETYPEALL

	string stagesaveall_generic_info =
		string("[tc1,1,0,1/tc]<b>Data save condition\n[tc1,1,0,1/tc]<b>for all stages</b>") +
		string("\n[tc1,1,0,1/tc]<i>green: on, red: off</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.set(stagesaveall_generic_info + string("<i><b>No saving set"), INT2(IOI_DSAVETYPEALL, DSAVE_NONE));
	ioInfo.set(stagesaveall_generic_info + string("<i><b>Save at the end of the stage"), INT2(IOI_DSAVETYPEALL, DSAVE_STAGE));
	ioInfo.set(stagesaveall_generic_info + string("<i><b>Save at the end of every step\n<i><b>in this stage"), INT2(IOI_DSAVETYPEALL, DSAVE_STEP));
	ioInfo.set(stagesaveall_generic_info + string("<i><b>Save every given iterations"), INT2(IOI_DSAVETYPEALL, DSAVE_ITER));
	ioInfo.set(stagesaveall_generic_info + string("<i><b>Stop every given time interval"), INT2(IOI_DSAVETYPEALL, DSAVE_TIME));

	//Shows parameter and value for a given mesh : minorId is the major id of elements in SimParams::simParams (i.e. an entry from PARAM_ enum), auxId is the unique mesh id number, textId is the parameter handle and value
	//IOI_MESHPARAM

	string param_generic_info =
		string("[tc1,1,0,1/tc]<b>Material parameter</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.set(param_generic_info + string("<i><b>Electron gyromagnetic ratio relative value"), INT2(IOI_MESHPARAM, PARAM_GREL));
	ioInfo.set(param_generic_info + string("<i><b>Gilbert damping"), INT2(IOI_MESHPARAM, PARAM_GDAMPING));
	ioInfo.set(param_generic_info + string("<i><b>Saturation magnetisation"), INT2(IOI_MESHPARAM, PARAM_MS));
	ioInfo.set(param_generic_info + string("<i><b>Demag factors Nx, Ny"), INT2(IOI_MESHPARAM, PARAM_DEMAGXY));
	ioInfo.set(param_generic_info + string("<i><b>Exchange stifness"), INT2(IOI_MESHPARAM, PARAM_A));
	ioInfo.set(param_generic_info + string("<i><b>Dzyaloshinskii-Moriya exchange"), INT2(IOI_MESHPARAM, PARAM_D));
	ioInfo.set(param_generic_info + string("<i><b>Bilinear surface exchange\n<i><b>Top mesh sets value"), INT2(IOI_MESHPARAM, PARAM_J1));
	ioInfo.set(param_generic_info + string("<i><b>Biquadratic surface exchange\n<i><b>Top mesh sets value"), INT2(IOI_MESHPARAM, PARAM_J2));
	ioInfo.set(param_generic_info + string("<i><b>Magnetocrystalline anisotropy"), INT2(IOI_MESHPARAM, PARAM_K1));
	ioInfo.set(param_generic_info + string("<i><b>Magnetocrystalline anisotropy (2nd order)"), INT2(IOI_MESHPARAM, PARAM_K2));
	ioInfo.set(param_generic_info + string("<i><b>Anisotropy symmetry axis - uniaxial\n<i><b>Cartesian unit vector"), INT2(IOI_MESHPARAM, PARAM_EA1));
	ioInfo.set(param_generic_info + string("<i><b>Anisotropy symmetry axis - cubic\n<i><b>Cartesian unit vector"), INT2(IOI_MESHPARAM, PARAM_EA2));
	ioInfo.set(param_generic_info + string("<i><b>Relative longitudinal\n<i><b>susceptibility for LLB"), INT2(IOI_MESHPARAM, PARAM_SUSREL));
	ioInfo.set(param_generic_info + string("<i><b>Relative transverse\n<i><b>susceptibility for LLB"), INT2(IOI_MESHPARAM, PARAM_SUSPREL));
	ioInfo.set(param_generic_info + string("<i><b>Applied field coefficient"), INT2(IOI_MESHPARAM, PARAM_HA));
	ioInfo.set(param_generic_info + string("<i><b>Base electrical conductivity"), INT2(IOI_MESHPARAM, PARAM_ELC));
	ioInfo.set(param_generic_info + string("<i><b>Anisotropic magnetoresistance"), INT2(IOI_MESHPARAM, PARAM_AMR));
	ioInfo.set(param_generic_info + string("<i><b>Current spin polarisation"), INT2(IOI_MESHPARAM, PARAM_P));
	ioInfo.set(param_generic_info + string("<i><b>Zhang-Li non-adiabaticity"), INT2(IOI_MESHPARAM, PARAM_BETA));
	ioInfo.set(param_generic_info + string("<i><b>Electron diffusion"), INT2(IOI_MESHPARAM, PARAM_DE));
	ioInfo.set(param_generic_info + string("<i><b>Diffusion spin polarisation"), INT2(IOI_MESHPARAM, PARAM_BETAD));
	ioInfo.set(param_generic_info + string("<i><b>Spin-Hall angle\n<i><b>In FM meshes used for SOTField module"), INT2(IOI_MESHPARAM, PARAM_SHA));
	ioInfo.set(param_generic_info + string("<i><b>Field-like spin torque coefficient\n<i><b>Used for SOTField module in FM meshes"), INT2(IOI_MESHPARAM, PARAM_FLSOT));
	ioInfo.set(param_generic_info + string("<i><b>Inverse spin-Hall angle"), INT2(IOI_MESHPARAM, PARAM_ISHA));
	ioInfo.set(param_generic_info + string("<i><b>Spin-flip length"), INT2(IOI_MESHPARAM, PARAM_LSF));
	ioInfo.set(param_generic_info + string("<i><b>Exchange rotation length"), INT2(IOI_MESHPARAM, PARAM_LEX));
	ioInfo.set(param_generic_info + string("<i><b>Spin dephasing length"), INT2(IOI_MESHPARAM, PARAM_LPH));
	ioInfo.set(param_generic_info + string("<i><b>Interface spin conductances (majority, minority)\n<i><b>Set to zero for continuous N-F interface\n<i><b>Top mesh sets value even if N"), INT2(IOI_MESHPARAM, PARAM_GI));
	ioInfo.set(param_generic_info + string("<i><b>Interface spin mixing conductance (real, imaginary)\n<i><b>Top mesh sets value even if N"), INT2(IOI_MESHPARAM, PARAM_GMIX));
	ioInfo.set(param_generic_info + string("<i><b>Spin torque efficiency in the bulk"), INT2(IOI_MESHPARAM, PARAM_TSEFF));
	ioInfo.set(param_generic_info + string("<i><b>Spin torque efficiency at interfaces"), INT2(IOI_MESHPARAM, PARAM_TSIEFF));
	ioInfo.set(param_generic_info + string("<i><b>Spin pumping efficiency"), INT2(IOI_MESHPARAM, PARAM_PUMPEFF));
	ioInfo.set(param_generic_info + string("<i><b>Thermal conductivity"), INT2(IOI_MESHPARAM, PARAM_THERMCOND));
	ioInfo.set(param_generic_info + string("<i><b>Mass density"), INT2(IOI_MESHPARAM, PARAM_DENSITY));
	ioInfo.set(param_generic_info + string("<i><b>Specific heat capacity"), INT2(IOI_MESHPARAM, PARAM_SHC));

	//Shows parameter temperature dependence for a given mesh : minorId is the major id of elements in SimParams::simParams (i.e. an entry from PARAM_ enum), auxId is the unique mesh id number, textId is the parameter temperature dependence setting
	//IOI_MESHPARAMTEMP

	string paramtemp_generic_info =
		string("[tc1,1,0,1/tc]<b>Material parameter\n[tc1,1,0,1/tc]<b>temperature dependence</b>") +
		string("\n[tc1,1,0,1/tc]<i>green: set, red: none</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n") +
		string("\n[tc1,1,0,1/tc]right-click: clear\n");

	ioInfo.push_back(paramtemp_generic_info, IOI_MESHPARAMTEMP);

	//Shows parameter spatial dependence for a given mesh : minorId is the major id of elements in SimParams::simParams (i.e. an entry from PARAM_ enum), auxId is the unique mesh id number, textId is the parameter spatial dependence setting
	//IOI_MESHPARAMVAR

	string paramvar_generic_info =
		string("[tc1,1,0,1/tc]<b>Material parameter spatial variation</b>") +
		string("\n[tc1,1,0,1/tc]<i>green: display selected, red: not selected</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n") +
		string("\n[tc1,1,0,1/tc]left-click: display select\n") +
		string("\n[tc1,1,0,1/tc]right-click: clear\n");

	ioInfo.push_back(paramvar_generic_info, IOI_MESHPARAMVAR);

	//Shows a possible formula name for mesh parameter temperature dependence : minorId is the MATPFORM_ enum value, textId is the formula name
	//IOI_MESHPARAMTEMPFORMULA

	string paramtempform_generic_info =
		string("[tc1,1,0,1/tc]<b>Temperature dependence\n[tc1,1,0,1/tc]<b>scaling formula") +
		string("\n[tc1,1,0,1/tc]drag: move to parameter to set\n");

	ioInfo.set(paramtempform_generic_info + string("<i><b>No temperature dependence"), INT2(IOI_MESHPARAMTEMPFORMULA, MATPFORM_NONE));
	ioInfo.set(paramtempform_generic_info + string("<i><b>Linear: y = c0 * x + 1\n<i><b>Can edit c0"), INT2(IOI_MESHPARAMTEMPFORMULA, MATPFORM_LINEAR));
	ioInfo.set(paramtempform_generic_info + string("<i><b>Parabolic: y = c0 * x^2 + c1 * x + 1\n<i><b>Can edit c0 and c1"), INT2(IOI_MESHPARAMTEMPFORMULA, MATPFORM_PARABOLIC));
	ioInfo.set(paramtempform_generic_info + string("<i><b>Inverse linear: y = 1 / (c0 * x + 1)\n<i><b>Can edit c0 >= 0"), INT2(IOI_MESHPARAMTEMPFORMULA, MATPFORM_INVERSELINEAR));

	//Shows a possible generator name for mesh parameter spatial dependence : minorId is the MATPVAR_ enum value, textId is the generator name
	//IOI_MESHPARAMVARGENERATOR 

	string paramvargen_generic_info =
		string("[tc1,1,0,1/tc]<b>Spatial dependence generator") +
		string("\n[tc1,1,0,1/tc]drag: move to parameter to set\n");

	ioInfo.set(paramvargen_generic_info + string("<i><b>Custom set from png file with grayscale.\n<i><b>offset, scale, filename\n<i><b>black = 0, white = 1"), INT2(IOI_MESHPARAMVARGENERATOR, MATPVAR_CUSTOM));
	ioInfo.set(paramvargen_generic_info + string("<i><b>Random with range (min, max) and seed"), INT2(IOI_MESHPARAMVARGENERATOR, MATPVAR_RANDOM));
	ioInfo.set(paramvargen_generic_info + string("<i><b>Jagged with range (min, max)\n<i><b>spacing (m) and seed"), INT2(IOI_MESHPARAMVARGENERATOR, MATPVAR_JAGGED));
	ioInfo.set(paramvargen_generic_info + string("<i><b>Circular defects with value range (min, max)\n<i><b>diameter range (min, max)\n<i><b>average spacing (m) and seed"), INT2(IOI_MESHPARAMVARGENERATOR, MATPVAR_DEFECTS));
	ioInfo.set(paramvargen_generic_info + string("<i><b>Line faults with value range (min, max)\n<i><b>length range (m) (min, max)\n<i><b>orientation range in degrees (min, max)\n<i><b>average spacing (m) and seed"), INT2(IOI_MESHPARAMVARGENERATOR, MATPVAR_FAULTS));
	ioInfo.set(paramvargen_generic_info + string("<i><b>Voronoi 2D cells with value range (min, max)\n<i><b>average spacing (m) and seed"), INT2(IOI_MESHPARAMVARGENERATOR, MATPVAR_VORONOI2D));
	ioInfo.set(paramvargen_generic_info + string("<i><b>Voronoi 3D cells with value range (min, max)\n<i><b>average spacing (m) and seed"), INT2(IOI_MESHPARAMVARGENERATOR, MATPVAR_VORONOI3D));
	ioInfo.set(paramvargen_generic_info + string("<i><b>Voronoi 2D cells with value range (min, max)\n<i><b>average spacing (m) and seed\n<i><b>Apply to Voronoi cell bounaries."), INT2(IOI_MESHPARAMVARGENERATOR, MATPVAR_VORONOIBND2D));
	ioInfo.set(paramvargen_generic_info + string("<i><b>Voronoi 3D cells with value range (min, max)\n<i><b>average spacing (m) and seed\n<i><b>Apply to Voronoi cell bounaries."), INT2(IOI_MESHPARAMVARGENERATOR, MATPVAR_VORONOIBND3D));
	ioInfo.set(paramvargen_generic_info + string("<i><b>Voronoi 2D cells for rotations.\n<i><b>Rotate vectorial parameters\n<i><b>through polar degrees (min, max)\n<i><b>and azimuthal degrees (min, max).\n<i><b>Average spacing (m) and seed"), INT2(IOI_MESHPARAMVARGENERATOR, MATPVAR_VORONOIROT2D));
	ioInfo.set(paramvargen_generic_info + string("<i><b>Voronoi 3D cells for rotations.\n<i><b>Rotate vectorial parameters\n<i><b>through polar degrees (min, max)\n<i><b>and azimuthal degrees (min, max).\n<i><b>Average spacing (m) and seed"), INT2(IOI_MESHPARAMVARGENERATOR, MATPVAR_VORONOIROT3D));

	//Shows mesh display option for a given mesh : minorId is the MESHDISPLAY_ value, auxId is the unique mesh id number, textId is the MESHDISPLAY_ handle
	//IOI_MESHDISPLAY

	string meshdisplay_generic_info =
		string("[tc1,1,0,1/tc]<b>Mesh quantity to display") +
		string("\n[tc1,1,0,1/tc]<i>green: on, red: off</i>") +
		string("\n[tc1,1,0,1/tc]click: change state\n");

	ioInfo.set(meshdisplay_generic_info + string("<i><b>Nothing displayed"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_NONE));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Magnetisation"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_MAGNETIZATION));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Total effective H field"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_EFFECTIVEFIELD));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Charge current density"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_CURRDENSITY));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Charge potential"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_VOLTAGE));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Electrical conductivity"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_ELCOND));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Spin accumulation"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_SACCUM));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Spin x-current density"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_JSX));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Spin y-current density"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_JSY));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Spin z-current density"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_JSZ));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Bulk spin torque"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_TS));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Interfacial spin torque"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_TSI));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Temperature"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_TEMPERATURE));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Mesh parameter spatial variation"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_PARAMVAR));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Roughness"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_ROUGHNESS));

	//Shows super-mesh display option : minorId is the MESHDISPLAY_ value, textId is the MESHDISPLAY_ handle
	//IOI_SMESHDISPLAY

	string smeshdisplay_generic_info =
		string("[tc1,1,0,1/tc]<b>Supermesh quantity to display") +
		string("\n[tc1,1,0,1/tc]<i>green: on, red: off</i>\n") +
		string("\n[tc1,1,0,1/tc]<i>When set, individual mesh\n") +
		string("\n[tc1,1,0,1/tc]<i>display not enabled</i>\n") +
		string("\n[tc1,1,0,1/tc]click: change state\n");

	ioInfo.set(smeshdisplay_generic_info + string("<i><b>Nothing displayed"), INT2(IOI_SMESHDISPLAY, MESHDISPLAY_NONE));
	ioInfo.set(smeshdisplay_generic_info + string("<i><b>Demagnetising field"), INT2(IOI_SMESHDISPLAY, MESHDISPLAY_SM_DEMAG));
	ioInfo.set(smeshdisplay_generic_info + string("<i><b>Oersted field"), INT2(IOI_SMESHDISPLAY, MESHDISPLAY_SM_OERSTED));
	ioInfo.set(smeshdisplay_generic_info + string("<i><b>Total dipole stray field"), INT2(IOI_SMESHDISPLAY, MESHDISPLAY_SM_STRAYH));

	//Shows movingmesh trigger settings : minorId is the unique mesh id number (if set), auxId is the trigger state (used or not used), textId is the mesh name (if set)
	//IOI_MOVINGMESH

	string movingmesh_info =
		string("[tc1,1,0,1/tc]<b>Moving mesh algorithm") +
		string("\n[tc1,1,0,1/tc]<i>green: on, red: off</i>\n") +
		string("\n[tc1,1,0,1/tc]<i>When set, trigger is\n") +
		string("\n[tc1,1,0,1/tc]<i>set on given mesh</i>\n") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n") +
		string("\n[tc1,1,0,1/tc]right-click: clear\n");

	ioInfo.push_back(movingmesh_info, IOI_MOVINGMESH);

	//Shows movingmesh trigger settings : minorId is the unique mesh id number (if set), auxId is the trigger state (used or not used), textId is the mesh name (if set)
	//IOI_MOVINGMESHASYM

	string movingmeshasym_info =
		string("[tc1,1,0,1/tc]<b>Moving mesh symmetry") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(movingmeshasym_info, IOI_MOVINGMESHASYM);

	//Shows movingmesh trigger settings : minorId is the unique mesh id number (if set), auxId is the trigger state (used or not used), textId is the mesh name (if set)
	//IOI_MOVINGMESHTHRESH

	string movingmeshthresh_info =
		string("[tc1,1,0,1/tc]<b>Moving mesh threshold") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(movingmeshthresh_info, IOI_MOVINGMESHTHRESH);

	//Shows electrode box. minorId is the minor Id in STransport::electrode_boxes, auxId is the number of the interactive object in the list (electrode index), textId is the electrode rect as a string
	//IOI_ELECTRODERECT

	string electroderect_info =
		string("[tc1,1,0,1/tc]<b>Electrode rectangle") +
		string("\n[tc1,1,0,1/tc]<i>Sets fixed potential at") +
		string("\n[tc1,1,0,1/tc]<i>intersections with mesh sides") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n") +
		string("\n[tc1,1,0,1/tc]right-click: delete\n");

	ioInfo.push_back(electroderect_info, IOI_ELECTRODERECT);

	//Shows electrode potential. minorId is the electrode index, textId is potential value as a string
	//IOI_ELECTRODEPOTENTIAL

	string electrodepotential_info =
		string("[tc1,1,0,1/tc]<b>Electrode potential") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n") +
		string("\n[tc1,1,0,1/tc]right-click: delete\n");

	ioInfo.push_back(electrodepotential_info, IOI_ELECTRODEPOTENTIAL);

	//Shows electrode ground setting. minorId is the electrode index, auxId is the setting (0 : not ground, 1 : ground)
	//IOI_ELECTRODEGROUND

	string electrodeground_info =
		string("[tc1,1,0,1/tc]<b>Electrode ground setting") +
		string("\n[tc1,1,0,1/tc]<i>green: ground, red: not ground</i>") +
		string("\n[tc1,1,0,1/tc]<i>Only one electrode can be a ground") +
		string("\n[tc1,1,0,1/tc]click: set ground\n");

	ioInfo.push_back(electrodeground_info, IOI_ELECTRODEGROUND);

	//Shows constant current source setting. auxId is the setting.
	//IOI_CONSTANTCURRENTSOURCE

	string powersupply_info =
		string("[tc1,1,0,1/tc]<b>Power supply mode") +
		string("\n[tc1,1,0,1/tc]<i>Constant current or constant voltage</i>") +
		string("\n[tc1,1,0,1/tc]click: switch mode\n");

	ioInfo.push_back(powersupply_info, IOI_CONSTANTCURRENTSOURCE);

	//Shows transport solver convergence error. textId is the convergence error value.
	//IOI_TSOLVERCONVERROR

	string tsolvererror_info =
		string("[tc1,1,0,1/tc]<b>Transport solver error - V") +
		string("\n[tc1,1,0,1/tc]<i>Convergence error setting</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(tsolvererror_info, IOI_TSOLVERCONVERROR);

	//Shows transport solver timeout iterations. auxId is the timeout value.
	//IOI_TSOLVERTIMEOUT

	string tsolvertout_info =
		string("[tc1,1,0,1/tc]<b>Transport solver timeout - V") +
		string("\n[tc1,1,0,1/tc]<i>Maximum number of iterations</i>") +
		string("\n[tc1,1,0,1/tc]<i>allowed to reach convergence</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(tsolvertout_info, IOI_TSOLVERTIMEOUT);

	//Shows spin transport solver convergence error. textId is the convergence error value.
	//IOI_SSOLVERCONVERROR

	string ssolvererror_info =
		string("[tc1,1,0,1/tc]<b>Transport solver error - S") +
		string("\n[tc1,1,0,1/tc]<i>Convergence error setting</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(ssolvererror_info, IOI_SSOLVERCONVERROR);

	//Shows spin transport solver timeout iterations. auxId is the timeout value.
	//IOI_SSOLVERTIMEOUT

	string ssolvertout_info =
		string("[tc1,1,0,1/tc]<b>Transport solver timeout - S") +
		string("\n[tc1,1,0,1/tc]<i>Maximum number of iterations</i>") +
		string("\n[tc1,1,0,1/tc]<i>allowed to reach convergence</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(ssolvertout_info, IOI_SSOLVERTIMEOUT);

	//Shows Poisson solver SOR damping type : true for fixed, false for adaptive. auxId is enabled (1)/disabled(0) status.
	//IOI_SORFIXEDDAMPING

	string sortype_info =
		string("[tc1,1,0,1/tc]<b>SOR damping type") +
		string("\n[tc1,1,0,1/tc]<i>SOR fixed or adaptive damping</i>") +
		string("\n[tc1,1,0,1/tc]click: switch mode\n");

	ioInfo.push_back(sortype_info, IOI_SORFIXEDDAMPING);

	//Shows SOR damping values when used in fixed damping mode. textId is the DBL2 damping value as a string. (DBL2 since we need different damping values for V and S solvers)
	//IOI_SORDAMPING

	string sordampingvalues_info =
		string("[tc1,1,0,1/tc]<b>SOR fixed damping") +
		string("\n[tc1,1,0,1/tc]<i>SOR damping for (V, S) solvers</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(sordampingvalues_info, IOI_SORDAMPING);

	//Shows mesh temperature. minorId is the unique mesh id number, textId is the temperature value
	//IOI_BASETEMPERATURE

	string basetemperature_info =
		string("[tc1,1,0,1/tc]<b>Mesh base temperature") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(basetemperature_info, IOI_BASETEMPERATURE);

	//Shows ambient temperature for heat equation Robin boundary conditions. minorId is the unique mesh id number, auxId is enabled/disabled status (Heat module must be active), textId is the temperature value
	//IOI_AMBIENT_TEMPERATURE

	string ambienttemperature_info =
		string("[tc1,1,0,1/tc]<b>Mesh ambient temperature") +
		string("\n[tc1,1,0,1/tc]<i>Air temperature outside of mesh</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(ambienttemperature_info, IOI_AMBIENT_TEMPERATURE);

	//Shows alpha value (W/m^2K) for heat equation Robin boundary conditions. minorId is the unique mesh id number, auxId is enabled/disabled status (Heat module must be active), textId is the value
	//IOI_ROBIN_ALPHA

	string robinalpha_info =
		string("[tc1,1,0,1/tc]<b>Robin heat flux coefficient") +
		string("\n[tc1,1,0,1/tc]<i>Boundary heat flux</i>") +
		string("\n[tc1,1,0,1/tc]<i>for Newton's law of cooling:</i>") +
		string("\n[tc1,1,0,1/tc]<i>flux = coeff * (T_boundary - T_ambient)</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(robinalpha_info, IOI_ROBIN_ALPHA);

	//Shows temperature insulating side setting for heat equation. minorId is the unique mesh id number, auxId is the status (Heat module must be active) : -1 disabled (gray), 0 not insulating (green), 1 insulating (red), textId represents the side : "x", "-x", "y", "-y", "z", "-z"
	//IOI_INSULATINGSIDE

	string insulatingside_info =
		string("[tc1,1,0,1/tc]<b>Insulating mesh side setting") +
		string("\n[tc1,1,0,1/tc]<i>green: not insulating, red: insulating</i>") +
		string("\n[tc1,1,0,1/tc]<i>If insulating no heat flux allowed</i>") +
		string("\n[tc1,1,0,1/tc]<i>else Newton's law of cooling applies</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(insulatingside_info, IOI_INSULATINGSIDE);

	//Shows mesh Curie temperature. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the temperature value
	//IOI_CURIETEMP

	string curietemperature_info =
		string("[tc1,1,0,1/tc]<b>Mesh Curie temperature") +
		string("\n[tc1,1,0,1/tc]<i>When set, parameter temperature</i>") +
		string("\n[tc1,1,0,1/tc]<i>dependencies automatically calculated for:</i>") +
		string("\n[tc1,1,0,1/tc]<i>damping, Ms, A, P, K1, K2, susrel</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(curietemperature_info, IOI_CURIETEMP);

	//Shows indicative material Curie temperature. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the temperature value
	//IOI_CURIETEMPMATERIAL

	curietemperature_info =
		string("[tc1,1,0,1/tc]<b>Material Curie temperature") +
		string("\n[tc1,1,0,1/tc]<i>This is the actual value for the material</i>") +
		string("\n[tc1,1,0,1/tc]<i>Not used in calculations, only indicative</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(curietemperature_info, IOI_CURIETEMPMATERIAL);

	//Shows atomic moment multiple of Bohr magneton. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the value
	//IOI_ATOMICMOMENT

	string atomicmoment_info =
		string("[tc1,1,0,1/tc]<b>Magnetic moment in Bohr magnetons") +
		string("\n[tc1,1,0,1/tc]<i>Used to calculate parameter</i>") +
		string("\n[tc1,1,0,1/tc]<i>temperature dependencies when T_Curie > 0</i>") +
		string("\n[tc1,1,0,1/tc]<i>In particular field-dependence is introduced</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(atomicmoment_info, IOI_ATOMICMOMENT);

	//Shows cuda enabled/disabled or n/a state. auxId is enabled (1)/disabled(0)/not available(-1) status.
	//IOI_CUDASTATE

	string cudastate_info =
		string("[tc1,1,0,1/tc]<b>CUDA computations state") +
		string("\n[tc1,1,0,1/tc]<i>green: set, red: not set</i>") +
		string("\n[tc1,1,0,1/tc]<i>When set, all computations done on the GPU</i>") +
		string("\n[tc1,1,0,1/tc]<i>otherwise done on the CPU</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(cudastate_info, IOI_CUDASTATE);

	//Shows gpu free memory. auxId is the value
	//IOI_GPUMEMFREE

	string gpufreemem_info =
		string("[tc1,1,0,1/tc]<b>GPU free memory");

	ioInfo.push_back(gpufreemem_info, IOI_GPUMEMFREE);

	//Shows gpu total memory. auxId is the value
	//IOI_GPUMEMTOTAL

	string gputotalmem_info =
		string("[tc1,1,0,1/tc]<b>GPU total memory");

	ioInfo.push_back(gputotalmem_info, IOI_GPUMEMTOTAL);

	//Shows cpu free memory. auxId is the value
	//IOI_CPUMEMFREE

	string cpufreemem_info =
		string("[tc1,1,0,1/tc]<b>CPU free memory");

	ioInfo.push_back(cpufreemem_info, IOI_CPUMEMFREE);

	//Shows cpu total memory. auxId is the value
	//IOI_CPUMEMTOTAL

	string cputotalmem_info =
		string("[tc1,1,0,1/tc]<b>CPU total memory");

	ioInfo.push_back(cputotalmem_info, IOI_CPUMEMTOTAL);

	//Shows scale_rects enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	//IOI_SCALERECTSSTATUS

	string scalerects_info =
		string("[tc1,1,0,1/tc]<b>Scale mesh rectangles status") +
		string("\n[tc1,1,0,1/tc]<i>green: set, red: not set</i>") +
		string("\n[tc1,1,0,1/tc]<i>When set, changing a mesh size will</i>") +
		string("\n[tc1,1,0,1/tc]<i>change all other meshes in proportion</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(scalerects_info, IOI_SCALERECTSSTATUS);

	//Shows coupled_to_dipoles enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	//IOI_COUPLEDTODIPOLESSTATUS

	string dipolecouple_info =
		string("[tc1,1,0,1/tc]<b>Dipole exchange coupling status") +
		string("\n[tc1,1,0,1/tc]<i>green: set, red: not set</i>") +
		string("\n[tc1,1,0,1/tc]<i>When set, for dipole-ferromagnetic mesh contacts</i>") +
		string("\n[tc1,1,0,1/tc]<i>mmoments at interface cells will be frozen</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(dipolecouple_info, IOI_COUPLEDTODIPOLESSTATUS);

	//Shows mesh roughness refinement value. minorId is the unique mesh id number, textId is the value
	//IOI_REFINEROUGHNESS

	string roughness_refine_info =
		string("[tc1,1,0,1/tc]<b>Mesh roughness") +
		string("\n[tc1,1,0,1/tc]<b>cells refinement multiplier\n") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(roughness_refine_info, IOI_REFINEROUGHNESS);

	//Multi-layered convolution configuration
	//IOI_MULTICONV, IOI_2DMULTICONV, IOI_NCOMMONSTATUS, IOI_NCOMMON

	string IOI_MULTICONV_info =
		string("[tc1,1,0,1/tc]<b>Multi-layered convolution status") +
		string("\n[tc1,1,0,1/tc]<i>green: set, red: not set, gray: N/A (</i>") +
		string("\n[tc1,1,0,1/tc]<i>When set, use multi-layered convolution</i>") +
		string("\n[tc1,1,0,1/tc]<i>instead of super-mesh convolution.</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(IOI_MULTICONV_info, IOI_MULTICONV);

	IOI_MULTICONV_info =
		string("[tc1,1,0,1/tc]<b>2D Multi-layered convolution status") +
		string("\n[tc1,1,0,1/tc]<i>green: set, red: not set, gray: N/A (</i>") +
		string("\n[tc1,1,0,1/tc]<i>When set, force multi-layered convolution</i>") +
		string("\n[tc1,1,0,1/tc]<i>to 2D in each layer.</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(IOI_MULTICONV_info, IOI_2DMULTICONV);

	IOI_MULTICONV_info =
		string("[tc1,1,0,1/tc]<b>Use default discretisation status") +
		string("\n[tc1,1,0,1/tc]<i>green: set, red: not set, gray: N/A (</i>") +
		string("\n[tc1,1,0,1/tc]<i>When set, use default common discretisation</i>") +
		string("\n[tc1,1,0,1/tc]<i>for multi-layered convolution.</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(IOI_MULTICONV_info, IOI_NCOMMONSTATUS);

	IOI_MULTICONV_info =
		string("[tc1,1,0,1/tc]<b>Common discretisation") +
		string("\n[tc1,1,0,1/tc]<i>Common discretisation</i>") +
		string("\n[tc1,1,0,1/tc]<i>for multi-layered convolution.</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(IOI_MULTICONV_info, IOI_NCOMMON);

	//Shows materials database in use. textId is the name of the database, including the path.
	//IOI_LOCALMDB

	string mdbfile_info =
		string("[tc1,1,0,1/tc]<b>Materials database file (.mdb)</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(mdbfile_info, IOI_LOCALMDB);
}

//---------------------------------------------------- MAKE INTERACTIVE OBJECT : Auxiliary method

template <typename ... PType>
string Simulation::MakeIO(IOI_ identifier, PType ... params)
{
	vector<string> params_str = make_vector<string>(ToString(params)...);

#if GRAPHICS == 1
	//string objectText, IOI_ majorId, int minorId = -1, int auxId = -1, string textId = "", D2D1_COLOR_F bgrndCol = MESSAGECOLOR
	auto MakeInteractiveObject = [&](string objectText, IOI_ majorId, int minorId = -1, int auxId = -1, string textId = "", D2D1_COLOR_F bgrndCol = MESSAGECOLOR) -> string {

		if (!textId.length()) textId = objectText;

		string newObject;

		string bgrndColString = ToString(bgrndCol.r) + "," + ToString(bgrndCol.g) + "," + ToString(bgrndCol.b) + "," + ToString(bgrndCol.a);

		newObject = "[io" + ToString((int)majorId) + "," + ToString(minorId) + "," + ToString((int)auxId) + "," + textId + "/io]<b>[or][tc1,1,1,1/tc][bc" + bgrndColString + "/bc] " + objectText + " </io>";

		return newObject;
	};
#else
	//string objectText, IOI_ majorId, int minorId = -1, int auxId = -1, string textId = "", D2D1_COLOR_F bgrndCol = MESSAGECOLOR
	auto MakeInteractiveObject = [&](string objectText, IOI_ majorId, int minorId = -1, int auxId = -1, string textId = "", int bgrndCol = MESSAGECOLOR) -> string {

		return objectText;
	};
#endif

	switch (identifier) {

	case IOI_PROGRAMUPDATESTATUS:
		return MakeInteractiveObject("Checking for updates...", IOI_PROGRAMUPDATESTATUS, 0, -1, "", UNAVAILABLECOLOR);
		break;

	case IOI_FMSMESHRECTANGLE:
		return MakeInteractiveObject(ToString(SMesh.GetFMSMeshRect(), "m"), IOI_FMSMESHRECTANGLE, 0, 0, ToString(SMesh.GetFMSMeshRect(), "m"));
		break;

	case IOI_FMSMESHCELLSIZE:
		return MakeInteractiveObject(ToString(SMesh.GetFMSMeshCellsize(), "m"), IOI_FMSMESHCELLSIZE, 0, 0, ToString(SMesh.GetFMSMeshCellsize(), "m"));
		break;

	case IOI_ESMESHRECTANGLE:
		return MakeInteractiveObject(ToString(SMesh.GetESMeshRect(), "m"), IOI_ESMESHRECTANGLE, 0, 0, ToString(SMesh.GetESMeshRect(), "m"));
		break;

	case IOI_ESMESHCELLSIZE:
		return MakeInteractiveObject(ToString(SMesh.GetESMeshCellsize(), "m"), IOI_ESMESHCELLSIZE, 0, 0, ToString(SMesh.GetESMeshCellsize(), "m"));
		break;

	case IOI_MESH_FORMESHLIST:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);
			string meshName = SMesh().get_key_from_index(meshIndex);

			return MakeInteractiveObject(meshName, IOI_MESH_FORMESHLIST, SMesh[meshIndex]->get_id(), 1, meshName);
		}
		break;

	case IOI_MESHRECTANGLE:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			return MakeInteractiveObject(ToString(SMesh[meshIndex]->GetMeshRect(), "m"), IOI_MESHRECTANGLE, SMesh[meshIndex]->get_id(), 0, ToString(SMesh[meshIndex]->GetMeshRect(), "m"));
		}
		break;

	case IOI_MESHCELLSIZE:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			if (SMesh[meshIndex]->MComputation_Enabled())
				return MakeInteractiveObject(ToString(SMesh[meshIndex]->GetMeshCellsize(), "m"), IOI_MESHCELLSIZE, SMesh[meshIndex]->get_id(), 1, ToString(SMesh[meshIndex]->GetMeshCellsize(), "m"), ONCOLOR);
			else
				return MakeInteractiveObject("N/A", IOI_MESHCELLSIZE, SMesh[meshIndex]->get_id(), 0, "N/A", OFFCOLOR);
		}
		break;

	case IOI_MESHECELLSIZE:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			if (SMesh[meshIndex]->EComputation_Enabled())
				return MakeInteractiveObject(ToString(SMesh[meshIndex]->GetMeshECellsize(), "m"), IOI_MESHECELLSIZE, SMesh[meshIndex]->get_id(), 1, ToString(SMesh[meshIndex]->GetMeshECellsize(), "m"), ONCOLOR);
			else
				return MakeInteractiveObject("N/A", IOI_MESHECELLSIZE, SMesh[meshIndex]->get_id(), 0, "N/A", OFFCOLOR);
		}
		break;
		
	case IOI_MESHTCELLSIZE:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			if (SMesh[meshIndex]->TComputation_Enabled())
				return MakeInteractiveObject(ToString(SMesh[meshIndex]->GetMeshTCellsize(), "m"), IOI_MESHTCELLSIZE, SMesh[meshIndex]->get_id(), 1, ToString(SMesh[meshIndex]->GetMeshTCellsize(), "m"), ONCOLOR);
			else
				return MakeInteractiveObject("N/A", IOI_MESHTCELLSIZE, SMesh[meshIndex]->get_id(), 0, "N/A", OFFCOLOR);
		}
		break;

	case IOI_SMESHDISPLAY:
		if (params_str.size() == 1) {

			int displayOption = ToNum(params_str[0]);

			return MakeInteractiveObject(displayHandles(displayOption), IOI_SMESHDISPLAY, displayOption, -1, displayHandles(displayOption));
		}
		break;

	case IOI_MESHDISPLAY:
		if (params_str.size() == 2) {

			int meshIndex = ToNum(params_str[0]);
			int displayOption = ToNum(params_str[1]);

			return MakeInteractiveObject(displayHandles(displayOption), IOI_MESHDISPLAY, displayOption, SMesh[meshIndex]->get_id(), displayHandles(displayOption));
		}
		break;

	case IOI_MESH_FORDISPLAYOPTIONS:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			string meshName = SMesh().get_key_from_index(meshIndex);

			return MakeInteractiveObject(meshName, IOI_MESH_FORDISPLAYOPTIONS, SMesh[meshIndex]->get_id(), 1, meshName);
		}
		break;

	case IOI_SMODULE:
		if (params_str.size() == 1) {

			int module = ToNum(params_str[0]);

			return MakeInteractiveObject(moduleHandles(module), IOI_SMODULE, module);
		}
		break;

	case IOI_MODULE:
		if (params_str.size() == 2) {

			int meshIndex = ToNum(params_str[0]);
			int module = ToNum(params_str[1]);

			string meshName = SMesh().get_key_from_index(meshIndex);

			return MakeInteractiveObject(moduleHandles(module), IOI_MODULE, module, SMesh[meshIndex]->get_id());
		}
		break;

	case IOI_MESH_FORMODULES:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			string meshName = SMesh().get_key_from_index(meshIndex);

			return MakeInteractiveObject(meshName, IOI_MESH_FORMODULES, SMesh[meshIndex]->get_id(), 1, meshName);
		}
		break;

	case IOI_ODE:
		if (params_str.size() == 2) {

			int odeId = ToNum(params_str[0]);
			int evalId = ToNum(params_str[1]);

			string evalHandle = odeEvalHandles(evalId);

			return MakeInteractiveObject(evalHandle, IOI_ODE, odeId, evalId, evalHandle);
		}
		break;

	case IOI_SHOWDATA:
		if (params_str.size() == 1) {

			int dataID = ToNum(params_str[0]);

			string dataName = dataDescriptor.get_key_from_ID(dataID);

			return MakeInteractiveObject(dataName, IOI_SHOWDATA, dataID, 0, dataName, UNAVAILABLECOLOR);
		}
		break;

	case IOI_DATA:
		if (params_str.size() == 1) {

			int dataID = ToNum(params_str[0]);

			string dataName = dataDescriptor.get_key_from_ID(dataID);

			return MakeInteractiveObject(dataName, IOI_DATA, dataID, 0, dataName, UNAVAILABLECOLOR);
		}
		break;

	case IOI_DIRECTORY:
		if (params_str.size() == 1) {

			string directory = params_str[0];

			return MakeInteractiveObject(directory, IOI_DIRECTORY, 0, 0, directory);
		}
		break;

	case IOI_SAVEDATAFILE:
		if (params_str.size() == 1) {

			string savedataFile = params_str[0];

			return MakeInteractiveObject(savedataFile, IOI_SAVEDATAFILE, 0, 0, savedataFile);
		}
		break;

	case IOI_SAVEDATAFLAG:
		return MakeInteractiveObject("On", IOI_SAVEDATAFLAG);
		break;

	case IOI_SAVEIMAGEFILEBASE:
		if (params_str.size() == 1) {

			string imageSaveFileBase = params_str[0];

			return MakeInteractiveObject(imageSaveFileBase, IOI_SAVEIMAGEFILEBASE, 0, 0, imageSaveFileBase);
		}
		break;

	case IOI_SAVEIMAGEFLAG:
		return MakeInteractiveObject("Off", IOI_SAVEIMAGEFLAG);
		break;

	case IOI_OUTDATA:
		if (params_str.size() == 1) {

			int index_in_list = ToNum(params_str[0]);

			string outputdata_text = Build_SetOutputData_Text(index_in_list);

			return MakeInteractiveObject(outputdata_text, IOI_OUTDATA, saveDataList.get_id_from_index(index_in_list).minor, index_in_list, outputdata_text);
		}
		break;

	case IOI_STAGE:
		if (params_str.size() == 1) {

			int stageID = ToNum(params_str[0]);

			string stagetype_text = stageDescriptors.get_key_from_ID(stageID);

			return MakeInteractiveObject(stagetype_text, IOI_STAGE, stageID, 0, stagetype_text, UNAVAILABLECOLOR);
		}
		break;

	case IOI_SETSTAGE:
		if (params_str.size() == 1) {

			int index_in_list = ToNum(params_str[0]);

			string outputdata_text = Build_SetStages_Text(index_in_list);

			return MakeInteractiveObject(outputdata_text, IOI_SETSTAGE, simStages.get_id_from_index(index_in_list).minor, index_in_list, outputdata_text);
		}
		break;

	case IOI_SETSTAGEVALUE:
		if (params_str.size() == 1) {

			int index_in_list = ToNum(params_str[0]);

			string setstage_text = simStages[index_in_list].get_value_string();
			
			return MakeInteractiveObject(setstage_text, IOI_SETSTAGEVALUE, simStages.get_id_from_index(index_in_list).minor, index_in_list, setstage_text);
		}
		break;

	case IOI_STAGESTOPCONDITION:
		if (params_str.size() == 1) {

			int index_in_list = ToNum(params_str[0]);

			string stopcondition_text = Build_SetStages_StopConditionText(index_in_list);

			return MakeInteractiveObject(stopcondition_text, IOI_STAGESTOPCONDITION, simStages.get_id_from_index(index_in_list).minor, index_in_list, stopcondition_text);
		}
		break;

	case IOI_DSAVETYPE:
		if (params_str.size() == 2) {

			int index_in_list = ToNum(params_str[0]);
			int dataSaveID = ToNum(params_str[1]);

			int dsaveIdx = dataSaveDescriptors.get_index_from_ID(dataSaveID);
			string savecondition_text = Build_SetStages_SaveConditionText(index_in_list, dsaveIdx);

			return MakeInteractiveObject(savecondition_text, IOI_DSAVETYPE, simStages.get_id_from_index(index_in_list).minor, dataSaveID, savecondition_text);
		}
		break;

	case IOI_STAGESTOPCONDITIONALL:
		if (params_str.size() == 1) {

			int stopID = ToNum(params_str[0]);

			string stopcondition_text = stageStopDescriptors.get_key_from_ID(stopID);

			return MakeInteractiveObject(stopcondition_text, IOI_STAGESTOPCONDITIONALL, stopID, 0, stopcondition_text, UNAVAILABLECOLOR);
		}
		break;

	case IOI_DSAVETYPEALL:
		if (params_str.size() == 1) {

			int dsaveID = ToNum(params_str[0]);

			string savetype_text = dataSaveDescriptors.get_key_from_ID(dsaveID);

			return MakeInteractiveObject(savetype_text, IOI_DSAVETYPEALL, dsaveID, 0, savetype_text, UNAVAILABLECOLOR);
		}
		break;

	case IOI_MESHPARAM:
		if (params_str.size() == 2) {

			int meshIndex = ToNum(params_str[0]);
			int paramId = ToNum(params_str[1]);

			string meshParam_text = Build_MeshParams_Text(meshIndex, (PARAM_)paramId);

			return MakeInteractiveObject(meshParam_text, IOI_MESHPARAM, paramId, SMesh[meshIndex]->get_id(), meshParam_text);
		}
		break;

	case IOI_MESH_FORPARAMS:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			int meshId = SMesh[meshIndex]->get_id();

			return MakeInteractiveObject(SMesh.key_from_meshId(meshId), IOI_MESH_FORPARAMS, meshId, 0, SMesh.key_from_meshId(meshId));
		}
		break;

	case IOI_MESHPARAMTEMP:
		if (params_str.size() == 2) {

			int meshIndex = ToNum(params_str[0]);
			int paramId = ToNum(params_str[1]);

			string tempDescriptor_text = SMesh[meshIndex]->get_paraminfo_string((PARAM_)paramId);

			return MakeInteractiveObject(tempDescriptor_text, IOI_MESHPARAMTEMP, paramId, SMesh[meshIndex]->get_id(), " ");
		}
		break;

	case IOI_MESHPARAMVAR:
		if (params_str.size() == 2) {

			int meshIndex = ToNum(params_str[0]);
			int paramId = ToNum(params_str[1]);

			string varDescriptor_text = SMesh[meshIndex]->get_paramvarinfo_string((PARAM_)paramId);

			return MakeInteractiveObject(varDescriptor_text, IOI_MESHPARAMVAR, paramId, SMesh[meshIndex]->get_id(), " ");
		}
		break;

	case IOI_MESHPARAMTEMPFORMULA:
		if (params_str.size() == 1) {

			int formulaID = ToNum(params_str[0]);
			
			string formulaName = formula_descriptor.get_key_from_ID(formulaID);

			return MakeInteractiveObject(formulaName, IOI_MESHPARAMTEMPFORMULA, formulaID, 0, formulaName);
		}
		break;

	case IOI_MESHPARAMVARGENERATOR:
		if (params_str.size() == 1) {

			int generatorID = ToNum(params_str[0]);

			string generatorName = vargenerator_descriptor.get_key_from_ID(generatorID);

			return MakeInteractiveObject(generatorName, IOI_MESHPARAMVARGENERATOR, generatorID, 0, generatorName);
		}
		break;

	case IOI_MESH_FORPARAMSTEMP:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			string meshName = SMesh.key_from_meshIdx(meshIndex);

			return MakeInteractiveObject(meshName, IOI_MESH_FORPARAMSTEMP, SMesh[meshIndex]->get_id(), 0, meshName);
		}
		break;

	case IOI_MESH_FORPARAMSVAR:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			string meshName = SMesh.key_from_meshIdx(meshIndex);

			return MakeInteractiveObject(meshName, IOI_MESH_FORPARAMSVAR, SMesh[meshIndex]->get_id(), 0, meshName);
		}
		break;

	case IOI_MOVINGMESH:
		if (SMesh.IsMovingMeshSet()) {

			int meshId = SMesh.GetId_of_MoveMeshTrigger();
			string meshName = SMesh.key_from_meshId(meshId);

			return MakeInteractiveObject(meshName, IOI_MOVINGMESH, meshId, true, meshName, ONCOLOR);
		}
		else {

			return MakeInteractiveObject("None", IOI_MOVINGMESH, -1, false, "", OFFCOLOR);
		}
		break;

	case IOI_MOVINGMESHASYM:
		return MakeInteractiveObject("Antisymmetric", IOI_MOVINGMESHASYM, 0, true);
		break;

	case IOI_MOVINGMESHTHRESH:
		return MakeInteractiveObject("0", IOI_MOVINGMESHTHRESH, 0, 0, "0");
		break;

	case IOI_CONSTANTCURRENTSOURCE:
		return MakeInteractiveObject(" ", IOI_CONSTANTCURRENTSOURCE, -1, !SMesh.CallModuleMethod(&STransport::UsingConstantCurrentSource));
		break;

	case IOI_ELECTRODERECT:
		if (params_str.size() == 1) {

			int el_index = ToNum(params_str[0]);

			pair<Rect, double> elInfo = SMesh.CallModuleMethod(&STransport::GetElectrodeInfo, el_index);
			int electrode_id = SMesh.CallModuleMethod(&STransport::GetElectrodeid, el_index);

			return MakeInteractiveObject(ToString(elInfo.first, "m"), IOI_ELECTRODERECT, electrode_id, el_index, ToString(elInfo.first, "m"));
		}
		break;

	case IOI_ELECTRODEPOTENTIAL:
		if (params_str.size() == 1) {

			int el_index = ToNum(params_str[0]);

			pair<Rect, double> elInfo = SMesh.CallModuleMethod(&STransport::GetElectrodeInfo, el_index);

			return MakeInteractiveObject(ToString(elInfo.second, "V"), IOI_ELECTRODEPOTENTIAL, el_index, 0, ToString(elInfo.second, "V"));
		}
		break;

	case IOI_ELECTRODEGROUND:
		if (params_str.size() == 1) {

			int el_index = ToNum(params_str[0]);

			if (SMesh.CallModuleMethod(&STransport::IsGroundElectrode, el_index)) {

				return MakeInteractiveObject("Grnd", IOI_ELECTRODEGROUND, el_index, true, "", ONCOLOR);
			}
			else {

				return MakeInteractiveObject("Grnd", IOI_ELECTRODEGROUND, el_index, false, "", OFFCOLOR);
			}
		}
		break;

	case IOI_TSOLVERCONVERROR:
	{
		double convergence_error = SMesh.CallModuleMethod(&STransport::GetConvergenceError);

		return MakeInteractiveObject(ToString(convergence_error), IOI_TSOLVERCONVERROR);
	}
	break;

	case IOI_TSOLVERTIMEOUT:
	{
		int iters_timeout = SMesh.CallModuleMethod(&STransport::GetConvergenceTimeout);

		return MakeInteractiveObject(ToString(iters_timeout), IOI_TSOLVERTIMEOUT);
	}
	break;

	case IOI_SSOLVERCONVERROR:
	{
		double convergence_error = SMesh.CallModuleMethod(&STransport::GetSConvergenceError);

		return MakeInteractiveObject(ToString(convergence_error), IOI_SSOLVERCONVERROR);
	}
	break;

	case IOI_SSOLVERTIMEOUT:
	{
		int iters_timeout = SMesh.CallModuleMethod(&STransport::GetSConvergenceTimeout);

		return MakeInteractiveObject(ToString(iters_timeout), IOI_SSOLVERTIMEOUT);
	}
	break;

	case IOI_SORFIXEDDAMPING:
	{
		bool fixed_SOR_damping = SMesh.CallModuleMethod(&STransport::IsFixedSORdamping);

		if (fixed_SOR_damping) return MakeInteractiveObject("Fixed", IOI_SORFIXEDDAMPING, -1, 1, "", ONCOLOR);
		else return MakeInteractiveObject("Adaptive", IOI_SORFIXEDDAMPING, -1, 0, "", OFFCOLOR);
	}
	break;

	case IOI_SORDAMPING:
	{
		DBL2 fixed_SOR_damping = SMesh.CallModuleMethod(&STransport::GetSORDamping);

		return MakeInteractiveObject(ToString(fixed_SOR_damping), IOI_SORDAMPING, -1, -1, ToString(fixed_SOR_damping));
	}
	break;

	case IOI_BASETEMPERATURE:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			return MakeInteractiveObject(ToString(SMesh[meshIndex]->GetAverageTemperature(), "K"), IOI_BASETEMPERATURE, SMesh[meshIndex]->get_id(), -1, ToString(SMesh[meshIndex]->GetAverageTemperature(), "K"));
		}
		break;

	case IOI_MESH_FORTEMPERATURE:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);
			string meshName = SMesh().get_key_from_index(meshIndex);

			return MakeInteractiveObject(meshName, IOI_MESH_FORTEMPERATURE, SMesh[meshIndex]->get_id(), 1, meshName);
		}
		break;

	case IOI_AMBIENT_TEMPERATURE:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			double T_ambient = SMesh[meshIndex]->CallModuleMethod(&Heat::GetAmbientTemperature);

			return MakeInteractiveObject(ToString(T_ambient, "K"), IOI_AMBIENT_TEMPERATURE, SMesh[meshIndex]->get_id(), 1, ToString(T_ambient, "K"));
		}
		break;

	case IOI_ROBIN_ALPHA:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			double alpha_boundary = SMesh[meshIndex]->CallModuleMethod(&Heat::GetAlphaBoundary);

			return MakeInteractiveObject(ToString(alpha_boundary, "W/m2K"), IOI_ROBIN_ALPHA, SMesh[meshIndex]->get_id(), 1, ToString(alpha_boundary, "W/m2K"));
		}
		break;

	case IOI_INSULATINGSIDE:
		if (params_str.size() == 2) {

			int meshIndex = ToNum(params_str[0]);
			string side_literal = params_str[1];

			return MakeInteractiveObject(side_literal + ": No", IOI_INSULATINGSIDE, SMesh[meshIndex]->get_id(), 0, side_literal);
		}
		break;

	case IOI_MESH_FORHEATBOUNDARIES:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);
			string meshName = SMesh().get_key_from_index(meshIndex);

			return MakeInteractiveObject(meshName, IOI_MESH_FORHEATBOUNDARIES, SMesh[meshName]->get_id(), 1, meshName);
		}
		break;

	case IOI_CURIETEMP:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			return MakeInteractiveObject(ToString(SMesh[meshIndex]->GetCurieTemperature(), "K"), IOI_CURIETEMP, SMesh[meshIndex]->get_id(), 1, ToString(SMesh[meshIndex]->GetCurieTemperature(), "K"));
		}
		break;

	case IOI_CURIETEMPMATERIAL:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			return MakeInteractiveObject(ToString(SMesh[meshIndex]->GetCurieTemperatureMaterial(), "K"), IOI_CURIETEMPMATERIAL, SMesh[meshIndex]->get_id(), 1, ToString(SMesh[meshIndex]->GetCurieTemperatureMaterial(), "K"));
		}
		break;

	case IOI_ATOMICMOMENT:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			return MakeInteractiveObject(ToString(SMesh[meshIndex]->GetAtomicMoment(), "uB"), IOI_ATOMICMOMENT, SMesh[meshIndex]->get_id(), 1, ToString(SMesh[meshIndex]->GetAtomicMoment(), "uB"));
		}
		break;

	case IOI_MESH_FORCURIEANDMOMENT:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);
			string meshName = SMesh().get_key_from_index(meshIndex);

			return MakeInteractiveObject(meshName, IOI_MESH_FORCURIEANDMOMENT, SMesh[meshName]->get_id(), 1, meshName);
		}
		break;

	//Shows cuda enabled/disabled or n/a state. auxId is enabled (1)/disabled(0)/not available(-1) status.
	case IOI_CUDASTATE:
		if (params_str.size() == 1) {

			int status = ToNum(params_str[0]);

			if(!cudaAvailable) return MakeInteractiveObject("N/A", IOI_CUDASTATE, -1, -1, "", UNAVAILABLECOLOR);
			else {

				if(status == 0) return MakeInteractiveObject("Off", IOI_CUDASTATE, -1, 0, "", OFFCOLOR);
				else return MakeInteractiveObject("On", IOI_CUDASTATE, -1, 1, "", ONCOLOR);
			}
		}
		break;

	//Shows gpu free memory. auxId is the value
	case IOI_GPUMEMFREE:
		if (params_str.size() == 1) {

			size_t mem_size = ToNum(params_str[0]);

			return MakeInteractiveObject(params_str[0], IOI_GPUMEMFREE, mem_size);
		}
		break;

	//Shows gpu total memory. auxId is the value
	case IOI_GPUMEMTOTAL:
		if (params_str.size() == 1) {

			size_t mem_size = ToNum(params_str[0]);

			return MakeInteractiveObject(params_str[0], IOI_GPUMEMTOTAL, mem_size);
		}
		break;

	//Shows cpu free memory. auxId is the value
	case IOI_CPUMEMFREE:
		if (params_str.size() == 1) {

			size_t mem_size = ToNum(params_str[0]);

			return MakeInteractiveObject(params_str[0], IOI_CPUMEMFREE, mem_size);
		}
		break;

	//Shows cpu total memory. auxId is the value
	case IOI_CPUMEMTOTAL:
		if (params_str.size() == 1) {

			size_t mem_size = ToNum(params_str[0]);

			return MakeInteractiveObject(params_str[0], IOI_CPUMEMTOTAL, mem_size);
		}
		break;
		
	//Shows scale_rects enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_SCALERECTSSTATUS:
		if (params_str.size() == 1) {

			int status = ToNum(params_str[0]);

			if (status == 0) return MakeInteractiveObject("Off", IOI_SCALERECTSSTATUS, -1, 0, "", OFFCOLOR);
			else return MakeInteractiveObject("On", IOI_SCALERECTSSTATUS, -1, 1, "", ONCOLOR);
		}
		break;

	//Shows coupled_to_dipoles enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_COUPLEDTODIPOLESSTATUS:
		if (params_str.size() == 1) {

			int status = ToNum(params_str[0]);

			if (status == 0) return MakeInteractiveObject("Off", IOI_COUPLEDTODIPOLESSTATUS, -1, 0, "", OFFCOLOR);
			else return MakeInteractiveObject("On", IOI_COUPLEDTODIPOLESSTATUS, -1, 1, "", ONCOLOR);
		}
		break;

	//Shows mesh roughness refinement value. minorId is the unique mesh id number, auxId is enabled (1)/disabled(0) status. textId is the value
	case IOI_REFINEROUGHNESS:
		if (params_str.size() == 1) {

			string meshName = params_str[0];

			if (SMesh[meshName]->IsModuleSet(MOD_ROUGHNESS))
				return MakeInteractiveObject(ToString(SMesh[meshName]->CallModuleMethod(&Roughness::get_refine)), IOI_REFINEROUGHNESS, SMesh[meshName]->get_id(), 1, ToString(SMesh[meshName]->CallModuleMethod(&Roughness::get_refine)), ONCOLOR);
			else 
				return MakeInteractiveObject("N/A", IOI_REFINEROUGHNESS, SMesh[meshName]->get_id(), 0, "N/A", UNAVAILABLECOLOR);
		}
		break;

	//Shows status of multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	case IOI_MULTICONV:
		if (params_str.size() == 1) {

			return MakeInteractiveObject("N/A", IOI_MULTICONV, 0, -1, "", UNAVAILABLECOLOR);
		}
		break;

	//Shows status of force 2D multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	case IOI_2DMULTICONV:
		if (params_str.size() == 1) {

			return MakeInteractiveObject("N/A", IOI_2DMULTICONV, 0, -1, "", UNAVAILABLECOLOR);
		}
		break;

	//Shows status of use default n for multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	case IOI_NCOMMONSTATUS:
		if (params_str.size() == 1) {

			return MakeInteractiveObject("N/A", IOI_NCOMMONSTATUS, 0, -1, "", UNAVAILABLECOLOR);
		}
		break;

	//Shows n_common for multi-layered convolution. auxId is the status (-1 : N/A, otherwise available). textId is the value as a SZ3.
	case IOI_NCOMMON:
		if (params_str.size() == 1) {

			return MakeInteractiveObject("N/A", IOI_NCOMMON, 0, -1, "0, 0, 0", UNAVAILABLECOLOR);
		}
		break;

	case IOI_LOCALMDB:
		if (params_str.size() == 1) {

			string mdbFile = params_str[0];

			return MakeInteractiveObject(mdbFile, IOI_LOCALMDB, 0, 0, mdbFile);
		}
		break;
	}

	return "";
}