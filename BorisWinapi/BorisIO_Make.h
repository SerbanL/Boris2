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
	ioInfo.set(modulegeneric_info + string("<i><b>Magneto-elastic term"), INT2(IOI_MODULE, MOD_MELASTIC));
	ioInfo.set(modulegeneric_info + string("<i><b>Charge and spin transport"), INT2(IOI_MODULE, MOD_TRANSPORT));
	ioInfo.set(modulegeneric_info + string("<i><b>Heat equation solver"), INT2(IOI_MODULE, MOD_HEAT));
	ioInfo.set(modulegeneric_info + string("<i><b>Spin-orbit torque field\n<i><b>Results in Slonczewski-like torques"), INT2(IOI_MODULE, MOD_SOTFIELD));
	ioInfo.set(modulegeneric_info + string("<i><b>Physical roughness\n<i><b>Demag term corrections when approximating shapes\n<i><b>from a fine mesh with a coarse mesh."), INT2(IOI_MODULE, MOD_ROUGHNESS));

	ioInfo.set(modulegeneric_info + string("<i><b>Supermesh demag field"), INT2(IOI_SMODULE, MODS_SDEMAG));
	ioInfo.set(modulegeneric_info + string("<i><b>Oersted field in electric supermesh"), INT2(IOI_SMODULE, MODS_OERSTED));
	ioInfo.set(modulegeneric_info + string("<i><b>Stray field from dipole meshes"), INT2(IOI_SMODULE, MODS_STRAYFIELD));

	//Available/set ode : minorId is an entry from ODE_ (the equation)
	//IOI_ODE

	string ode_info =
		string("[tc1,1,0,1/tc]<b>dM/dt equation</b>") +
		string("\n[tc1,1,0,1/tc]<i>green: on, red: off</i>") +
		string("\n[tc1,1,0,1/tc]click: change state\n");

	ioInfo.push_back(ode_info, IOI_ODE);

	//Set ODE time step: textId is the value
	//IOI_ODEDT

	string odedt_info =
		string("[tc1,1,0,1/tc]<b>ODE Time Step</b>") +
		string("\n[tc1,1,0,1/tc]double-click: edit\n");

	ioInfo.push_back(odedt_info, IOI_ODEDT);

	//Set heat equation time step: textId is the value
	//IOI_HEATDT

	string odeheatdt_info =
		string("[tc1,1,0,1/tc]<b>Heat Equation Time Step</b>") +
		string("\n[tc1,1,0,1/tc]double-click: edit\n");

	ioInfo.push_back(odeheatdt_info, IOI_HEATDT);

	//Available/set ode evaluation method for ode : minorId is an entry from ODE_ (the equation), auxId is the EVAL_ entry (the evaluation method), textId is the name of the evaluation method
	//IOI_ODE_EVAL

	string ode_eval_info =
		string("[tc1,1,0,1/tc]<b>Evaluation method for dM/dt equation</b>") +
		string("\n[tc1,1,0,1/tc]<i>green: on, red: off, gray: unavailable</i>") +
		string("\n[tc1,1,0,1/tc]click: change state\n");

	ioInfo.push_back(ode_eval_info, IOI_ODE_EVAL);

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
	//IOI_MESH_FORTMODEL
	//IOI_MESH_FORPBC
	//IOI_MESH_FOREXCHCOUPLING
	//IOI_MESH_FORSTOCHASTICITY

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
	ioInfo.push_back(mesh_info, IOI_MESH_FORTMODEL);
	ioInfo.push_back(mesh_info, IOI_MESH_FORPBC);
	ioInfo.push_back(mesh_info, IOI_MESH_FOREXCHCOUPLING);
	ioInfo.push_back(mesh_info, IOI_MESH_FORSTOCHASTICITY);

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

	//Shows mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	//IOI_MESHMCELLSIZE:

	string mechmeshcell_info =
		string("[tc1,1,0,1/tc]<b>Mechanical solver cellsize</b>") +
		string("\n[tc1,1,0,1/tc]<i>Discretization cellsize</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(mechmeshcell_info, IOI_MESHMCELLSIZE);

	//Shows stochastic cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	//IOI_MESHSCELLSIZE:

	string stochmeshcell_info =
		string("[tc1,1,0,1/tc]<b>Stochasticity cellsize</b>") +
		string("\n[tc1,1,0,1/tc]<i>Stochastic fields discretization</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(stochmeshcell_info, IOI_MESHSCELLSIZE);

	//Shows link stochastic flag : minorId is the unique mesh id number, auxId is the value off (0), on (1), N/A (-1)
	//IOI_LINKSTOCHASTIC

	string linkstoch_info =
		string("[tc1,1,0,1/tc]<b>Stochasticity cellsize flag</b>") +
		string("\n[tc1,1,0,1/tc]<i>Link to magnetic cellsize</i>") +
		string("\n[tc1,1,0,1/tc]click: change");

	ioInfo.push_back(linkstoch_info, IOI_LINKSTOCHASTIC);

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
	ioInfo.set(showdata_info_generic + string("<i><b>Magnetisation relaxation |mxh|</i>"), INT2(IOI_SHOWDATA, DATA_MXH));
	ioInfo.set(showdata_info_generic + string("<i><b>Magnetisation relaxation |dm/dt|</i>"), INT2(IOI_SHOWDATA, DATA_DMDT));
	ioInfo.set(showdata_info_generic + string("<i><b>Average magnetisation</i>"), INT2(IOI_SHOWDATA, DATA_AVM));
	ioInfo.set(showdata_info_generic + string("<i><b>Average magnetisation sub-lattice B</i>"), INT2(IOI_SHOWDATA, DATA_AVM2));
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
	ioInfo.set(showdata_info_generic + string("<i><b>Energy density: applied mechanical stress</i>"), INT2(IOI_SHOWDATA, DATA_E_MELASTIC));
	ioInfo.set(showdata_info_generic + string("<i><b>Energy density: anisotropy</i>"), INT2(IOI_SHOWDATA, DATA_E_ANIS));
	ioInfo.set(showdata_info_generic + string("<i><b>Energy density: roughness</i>"), INT2(IOI_SHOWDATA, DATA_E_ROUGH));
	ioInfo.set(showdata_info_generic + string("<i><b>Energy density: Total</i>"), INT2(IOI_SHOWDATA, DATA_E_TOTAL));
	ioInfo.set(showdata_info_generic + string("<i><b>Domain wall shift\n<i><b>for moving mesh</i>"), INT2(IOI_SHOWDATA, DATA_DWSHIFT));
	ioInfo.set(showdata_info_generic + string("<i><b>Skyrmion shift in the xy plane\n<i><b>Only use with output save data</i>"), INT2(IOI_SHOWDATA, DATA_SKYSHIFT));
	ioInfo.set(showdata_info_generic + string("<i><b>Skyrmion shift in the xy plane\n<i><b>Additional saving of x and y axis diameters\n<i><b>Only use with output save data</i>"), INT2(IOI_SHOWDATA, DATA_SKYPOS));
	ioInfo.set(showdata_info_generic + string("<i><b>Transport solver:\n<i><b>V iterations to convergence</i>"), INT2(IOI_SHOWDATA, DATA_TRANSPORT_ITERSTOCONV));
	ioInfo.set(showdata_info_generic + string("<i><b>Transport solver:\n<i><b>S iterations to convergence</i>"), INT2(IOI_SHOWDATA, DATA_TRANSPORT_SITERSTOCONV));
	ioInfo.set(showdata_info_generic + string("<i><b>Transport solver:\n<i><b>achieved convergence error</i>"), INT2(IOI_SHOWDATA, DATA_TRANSPORT_CONVERROR));
	ioInfo.set(showdata_info_generic + string("<i><b>Average temperature</i>"), INT2(IOI_SHOWDATA, DATA_TEMP));
	ioInfo.set(showdata_info_generic + string("<i><b>Average lattice temperature</i>"), INT2(IOI_SHOWDATA, DATA_TEMP_L));
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
	ioInfo.set(data_info_generic + string("<i><b>Magnetisation relaxation |mxh|</i>"), INT2(IOI_DATA, DATA_MXH));
	ioInfo.set(data_info_generic + string("<i><b>Magnetisation relaxation |dm/dt|</i>"), INT2(IOI_DATA, DATA_DMDT));
	ioInfo.set(data_info_generic + string("<i><b>Average magnetisation</i>"), INT2(IOI_DATA, DATA_AVM));
	ioInfo.set(data_info_generic + string("<i><b>Average magnetisation sub-lattice B</i>"), INT2(IOI_DATA, DATA_AVM2));
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
	ioInfo.set(data_info_generic + string("<i><b>Energy density: applied mechanical stress</i>"), INT2(IOI_DATA, DATA_E_MELASTIC));
	ioInfo.set(data_info_generic + string("<i><b>Energy density: anisotropy</i>"), INT2(IOI_DATA, DATA_E_ANIS));
	ioInfo.set(data_info_generic + string("<i><b>Energy density: roughness</i>"), INT2(IOI_DATA, DATA_E_ROUGH));
	ioInfo.set(data_info_generic + string("<i><b>Energy density: Total</i>"), INT2(IOI_DATA, DATA_E_TOTAL));
	ioInfo.set(data_info_generic + string("<i><b>Domain wall shift\n<i><b>for moving mesh</i>"), INT2(IOI_DATA, DATA_DWSHIFT));
	ioInfo.set(data_info_generic + string("<i><b>Skyrmion shift in the xy plane\n<i><b>Rectangle must circumscribe skyrmion</i>"), INT2(IOI_DATA, DATA_SKYSHIFT));
	ioInfo.set(data_info_generic + string("<i><b>Skyrmion shift in the xy plane\n<i><b>Also save x and y axis diameters\n<i><b>Rectangle must circumscribe skyrmion</i>"), INT2(IOI_DATA, DATA_SKYPOS));
	ioInfo.set(data_info_generic + string("<i><b>Transport solver:\n<i><b>V iterations to convergence</i>"), INT2(IOI_DATA, DATA_TRANSPORT_ITERSTOCONV));
	ioInfo.set(data_info_generic + string("<i><b>Transport solver:\n<i><b>S iterations to convergence</i>"), INT2(IOI_DATA, DATA_TRANSPORT_SITERSTOCONV));
	ioInfo.set(data_info_generic + string("<i><b>Transport solver:\n<i><b>achieved convergence error</i>"), INT2(IOI_DATA, DATA_TRANSPORT_CONVERROR));
	ioInfo.set(data_info_generic + string("<i><b>Average temperature</i>"), INT2(IOI_DATA, DATA_TEMP));
	ioInfo.set(data_info_generic + string("<i><b>Average lattice temperature</i>"), INT2(IOI_DATA, DATA_TEMP_L));
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
	ioInfo.set(stage_generic_info + string("<i><b>Set field using custom equation"), INT2(IOI_STAGE, SS_HFIELDEQUATION));
	ioInfo.set(stage_generic_info + string("<i><b>Set field sequence using custom equation"), INT2(IOI_STAGE, SS_HFIELDEQUATIONSEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Set from file in current directory\n<i><b>File must have tab-spaced columns as:\n<i><b>time and value all in S.I. units.\n<i><b>Time resolution set by stage time condition."), INT2(IOI_STAGE, SS_HFIELDFILE));
	ioInfo.set(stage_generic_info + string("<i><b>Set fixed voltage drop between electrodes"), INT2(IOI_STAGE, SS_V));
	ioInfo.set(stage_generic_info + string("<i><b>Set fixed voltage sequence\n<i><b>Start to stop V values in a number of steps"), INT2(IOI_STAGE, SS_VSEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Set a sinusoidal voltage sequence\n<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_STAGE, SS_VSIN));
	ioInfo.set(stage_generic_info + string("<i><b>Set a cosine voltage sequence\n<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_STAGE, SS_VCOS));
	ioInfo.set(stage_generic_info + string("<i><b>Set potential using custom equation"), INT2(IOI_STAGE, SS_VEQUATION));
	ioInfo.set(stage_generic_info + string("<i><b>Set potential sequence using custom equation"), INT2(IOI_STAGE, SS_VEQUATIONSEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Set from file in current directory\n<i><b>File must have tab-spaced columns as:\n<i><b>time and value all in S.I. units.\n<i><b>Time resolution set by stage time condition."), INT2(IOI_STAGE, SS_VFILE));
	ioInfo.set(stage_generic_info + string("<i><b>Set fixed current into ground electrode"), INT2(IOI_STAGE, SS_I));
	ioInfo.set(stage_generic_info + string("<i><b>Set fixed current sequence\n<i><b>Start to stop I values in a number of steps"), INT2(IOI_STAGE, SS_ISEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Set a sinusoidal current sequence\n<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_STAGE, SS_ISIN));
	ioInfo.set(stage_generic_info + string("<i><b>Set a cosine current sequence\n<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_STAGE, SS_ICOS));
	ioInfo.set(stage_generic_info + string("<i><b>Set current using custom equation"), INT2(IOI_STAGE, SS_IEQUATION));
	ioInfo.set(stage_generic_info + string("<i><b>Set current sequence using custom equation"), INT2(IOI_STAGE, SS_IEQUATIONSEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Set from file in current directory\n<i><b>File must have tab-spaced columns as:\n<i><b>time and value all in S.I. units.\n<i><b>Time resolution set by stage time condition."), INT2(IOI_STAGE, SS_IFILE));
	ioInfo.set(stage_generic_info + string("<i><b>Set base temperature value"), INT2(IOI_STAGE, SS_T));
	ioInfo.set(stage_generic_info + string("<i><b>Set base temperature sequence\n<i><b>Start to stop T values in a number of steps"), INT2(IOI_STAGE, SS_TSEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Set base temperature using custom equation"), INT2(IOI_STAGE, SS_TEQUATION));
	ioInfo.set(stage_generic_info + string("<i><b>Set base temperature sequence using custom equation"), INT2(IOI_STAGE, SS_TEQUATIONSEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Set from file in current directory\n<i><b>File must have tab-spaced columns as:\n<i><b>time and value all in S.I. units.\n<i><b>Time resolution set by stage time condition."), INT2(IOI_STAGE, SS_TFILE));
	ioInfo.set(stage_generic_info + string("<i><b>Set heat source value"), INT2(IOI_STAGE, SS_Q));
	ioInfo.set(stage_generic_info + string("<i><b>Set heat source sequence\n<i><b>Start to stop Q values in a number of steps"), INT2(IOI_STAGE, SS_QSEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Set heat source using custom equation"), INT2(IOI_STAGE, SS_QEQUATION));
	ioInfo.set(stage_generic_info + string("<i><b>Set heat source sequence using custom equation"), INT2(IOI_STAGE, SS_QEQUATIONSEQ));
	ioInfo.set(stage_generic_info + string("<i><b>Set from file in current directory\n<i><b>File must have tab-spaced columns as:\n<i><b>time and value all in S.I. units.\n<i><b>Time resolution set by stage time condition."), INT2(IOI_STAGE, SS_QFILE));
	ioInfo.set(stage_generic_info + string("<i><b>Set uniform stress in polar coordinates"), INT2(IOI_STAGE, SS_TSIGPOLAR));

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
	ioInfo.set(stagevalue_generic_info + string("<i><b>Text vector equation"), INT2(IOI_SETSTAGEVALUE, SS_HFIELDEQUATION));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Steps: Text vector equation"), INT2(IOI_SETSTAGEVALUE, SS_HFIELDEQUATIONSEQ));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Set values from file."), INT2(IOI_SETSTAGEVALUE, SS_HFIELDFILE));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Fixed voltage drop between electrodes"), INT2(IOI_SETSTAGEVALUE, SS_V));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Vstart; Vstop; Steps: Vstep = (Vstop - Vstart) / Steps"), INT2(IOI_SETSTAGEVALUE, SS_VSEQ));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_SETSTAGEVALUE, SS_VSIN));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_SETSTAGEVALUE, SS_VCOS));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Text equation"), INT2(IOI_SETSTAGEVALUE, SS_VEQUATION));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Steps: Text equation"), INT2(IOI_SETSTAGEVALUE, SS_VEQUATIONSEQ));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Set values from file."), INT2(IOI_SETSTAGEVALUE, SS_VFILE));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Fixed current into ground electrode"), INT2(IOI_SETSTAGEVALUE, SS_I));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Istart; Istop; Steps: Istep = (Istop - Istart) / Steps"), INT2(IOI_SETSTAGEVALUE, SS_ISEQ));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_SETSTAGEVALUE, SS_ISIN));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Amplitude, steps per cycle, cycles"), INT2(IOI_SETSTAGEVALUE, SS_ICOS));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Text equation"), INT2(IOI_SETSTAGEVALUE, SS_IEQUATION));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Steps: Text equation"), INT2(IOI_SETSTAGEVALUE, SS_IEQUATIONSEQ));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Set values from file."), INT2(IOI_SETSTAGEVALUE, SS_IFILE));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Base temperature value"), INT2(IOI_SETSTAGEVALUE, SS_T));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Tstart; Tstop; Steps: Tstep = (Tstop - Tstart) / Steps"), INT2(IOI_SETSTAGEVALUE, SS_TSEQ));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Text equation"), INT2(IOI_SETSTAGEVALUE, SS_TEQUATION));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Steps: Text equation"), INT2(IOI_SETSTAGEVALUE, SS_TEQUATIONSEQ));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Set values from file."), INT2(IOI_SETSTAGEVALUE, SS_TFILE));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Heat source value"), INT2(IOI_SETSTAGEVALUE, SS_Q));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Qstart; Qstop; Steps: Qstep = (Qstop - Qstart) / Steps"), INT2(IOI_SETSTAGEVALUE, SS_QSEQ));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Text equation"), INT2(IOI_SETSTAGEVALUE, SS_QEQUATION));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Steps: Text equation"), INT2(IOI_SETSTAGEVALUE, SS_QEQUATIONSEQ));
	ioInfo.set(stagevalue_generic_info + string("<i><b>Set values from file."), INT2(IOI_SETSTAGEVALUE, SS_QFILE));
	ioInfo.set(stagevalue_generic_info + string("<i><b>magnitude, theta, polar"), INT2(IOI_SETSTAGEVALUE, SS_TSIGPOLAR));

	//Shows the stop condition for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the stop type and value as a string
	//IOI_STAGESTOPCONDITION

	string stagestop_generic_info =
		string("[tc1,1,0,1/tc]<b>Stage stop condition</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.set(stagestop_generic_info + string("<i><b>No stopping condition set"), INT2(IOI_STAGESTOPCONDITION, STOP_NOSTOP));
	ioInfo.set(stagestop_generic_info + string("<i><b>Stop when stage iterations value reached"), INT2(IOI_STAGESTOPCONDITION, STOP_ITERATIONS));
	ioInfo.set(stagestop_generic_info + string("<i><b>Stop when |mxh| falls below value"), INT2(IOI_STAGESTOPCONDITION, STOP_MXH));
	ioInfo.set(stagestop_generic_info + string("<i><b>Stop when |dm/dt| falls below value"), INT2(IOI_STAGESTOPCONDITION, STOP_DMDT));
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
	ioInfo.set(stagestopall_generic_info + string("<i><b>Stop when |dm/dt| falls below value"), INT2(IOI_STAGESTOPCONDITIONALL, STOP_DMDT));
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
	ioInfo.set(param_generic_info + string("<i><b>Electron gyromagnetic ratio relative value"), INT2(IOI_MESHPARAM, PARAM_GREL_AFM));
	ioInfo.set(param_generic_info + string("<i><b>Gilbert damping"), INT2(IOI_MESHPARAM, PARAM_GDAMPING));
	ioInfo.set(param_generic_info + string("<i><b>Gilbert damping"), INT2(IOI_MESHPARAM, PARAM_GDAMPING_AFM));
	ioInfo.set(param_generic_info + string("<i><b>Saturation magnetisation"), INT2(IOI_MESHPARAM, PARAM_MS));
	ioInfo.set(param_generic_info + string("<i><b>Saturation magnetisation"), INT2(IOI_MESHPARAM, PARAM_MS_AFM));
	ioInfo.set(param_generic_info + string("<i><b>Demag factors Nx, Ny"), INT2(IOI_MESHPARAM, PARAM_DEMAGXY));
	ioInfo.set(param_generic_info + string("<i><b>Exchange stiffness"), INT2(IOI_MESHPARAM, PARAM_A));
	ioInfo.set(param_generic_info + string("<i><b>Exchange stiffness"), INT2(IOI_MESHPARAM, PARAM_A_AFM));
	ioInfo.set(param_generic_info + string("<i><b>Homogeneous AFM coupling"), INT2(IOI_MESHPARAM, PARAM_A_AFH));
	ioInfo.set(param_generic_info + string("<i><b>Nonhomogeneous AFM coupling"), INT2(IOI_MESHPARAM, PARAM_A_AFNH));
	ioInfo.set(param_generic_info + string("<i><b>Exchange parameter to critical temperature ratio\n<i><b>Used with 2-sublattice model."), INT2(IOI_MESHPARAM, PARAM_AFTAU));
	ioInfo.set(param_generic_info + string("<i><b>Exchange parameter to critical temperature ratio\n<i><b>Used with 2-sublattice model, cross-coupling"), INT2(IOI_MESHPARAM, PARAM_AFTAUCROSS));
	ioInfo.set(param_generic_info + string("<i><b>Dzyaloshinskii-Moriya exchange"), INT2(IOI_MESHPARAM, PARAM_D));
	ioInfo.set(param_generic_info + string("<i><b>Dzyaloshinskii-Moriya exchange"), INT2(IOI_MESHPARAM, PARAM_D_AFM));
	ioInfo.set(param_generic_info + string("<i><b>Bilinear surface exchange\n<i><b>Top mesh sets value"), INT2(IOI_MESHPARAM, PARAM_J1));
	ioInfo.set(param_generic_info + string("<i><b>Biquadratic surface exchange\n<i><b>Top mesh sets value"), INT2(IOI_MESHPARAM, PARAM_J2));
	ioInfo.set(param_generic_info + string("<i><b>Surface exchange from diamagnet\n<i><b>Diamagnetic mesh sets value"), INT2(IOI_MESHPARAM, PARAM_NETADIA));
	ioInfo.set(param_generic_info + string("<i><b>Magnetocrystalline anisotropy"), INT2(IOI_MESHPARAM, PARAM_K1));
	ioInfo.set(param_generic_info + string("<i><b>Magnetocrystalline anisotropy (2nd order)"), INT2(IOI_MESHPARAM, PARAM_K2));
	ioInfo.set(param_generic_info + string("<i><b>Magnetocrystalline anisotropy 2-sublattice"), INT2(IOI_MESHPARAM, PARAM_K1_AFM));
	ioInfo.set(param_generic_info + string("<i><b>Magnetocrystalline anisotropy (2nd order) 2-sublattice"), INT2(IOI_MESHPARAM, PARAM_K2_AFM));
	ioInfo.set(param_generic_info + string("<i><b>Anisotropy symmetry axis - uniaxial\n<i><b>Cartesian unit vector"), INT2(IOI_MESHPARAM, PARAM_EA1));
	ioInfo.set(param_generic_info + string("<i><b>Anisotropy symmetry axis - cubic\n<i><b>Cartesian unit vector"), INT2(IOI_MESHPARAM, PARAM_EA2));
	ioInfo.set(param_generic_info + string("<i><b>Relative longitudinal\n<i><b>susceptibility for LLB"), INT2(IOI_MESHPARAM, PARAM_SUSREL));
	ioInfo.set(param_generic_info + string("<i><b>Relative longitudinal\n<i><b>susceptibility for LLB 2-sublattice"), INT2(IOI_MESHPARAM, PARAM_SUSREL_AFM));
	ioInfo.set(param_generic_info + string("<i><b>Relative transverse\n<i><b>susceptibility for LLB"), INT2(IOI_MESHPARAM, PARAM_SUSPREL));
	ioInfo.set(param_generic_info + string("<i><b>Applied field coefficient"), INT2(IOI_MESHPARAM, PARAM_HA));
	ioInfo.set(param_generic_info + string("<i><b>Set temperature coefficient"), INT2(IOI_MESHPARAM, PARAM_T));
	ioInfo.set(param_generic_info + string("<i><b>Heat source"), INT2(IOI_MESHPARAM, PARAM_Q));
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
	ioInfo.set(param_generic_info + string("<i><b>Interface spin conductances (majority, minority)\n<i><b>Top mesh sets value even if N"), INT2(IOI_MESHPARAM, PARAM_GI));
	ioInfo.set(param_generic_info + string("<i><b>Interface spin mixing conductance (real, imaginary)\n<i><b>Set to zero for continuous N-F interface\n<i><b>Top mesh sets value even if N"), INT2(IOI_MESHPARAM, PARAM_GMIX));
	ioInfo.set(param_generic_info + string("<i><b>Spin torque efficiency in the bulk"), INT2(IOI_MESHPARAM, PARAM_TSEFF));
	ioInfo.set(param_generic_info + string("<i><b>Spin torque efficiency at interfaces"), INT2(IOI_MESHPARAM, PARAM_TSIEFF));
	ioInfo.set(param_generic_info + string("<i><b>Spin pumping efficiency"), INT2(IOI_MESHPARAM, PARAM_PUMPEFF));
	ioInfo.set(param_generic_info + string("<i><b>Charge pumping efficiency"), INT2(IOI_MESHPARAM, PARAM_CPUMP_EFF));
	ioInfo.set(param_generic_info + string("<i><b>Topological Hall efficiency"), INT2(IOI_MESHPARAM, PARAM_THE_EFF));
	ioInfo.set(param_generic_info + string("<i><b>Carrier density"), INT2(IOI_MESHPARAM, PARAM_NDENSITY));
	ioInfo.set(param_generic_info + string("<i><b>Thermal conductivity"), INT2(IOI_MESHPARAM, PARAM_THERMCOND));
	ioInfo.set(param_generic_info + string("<i><b>Mass density"), INT2(IOI_MESHPARAM, PARAM_DENSITY));
	ioInfo.set(param_generic_info + string("<i><b>Magnetoelastic coefficients (B1, B2)"), INT2(IOI_MESHPARAM, PARAM_MECOEFF));
	ioInfo.set(param_generic_info + string("<i><b>Young's modulus"), INT2(IOI_MESHPARAM, PARAM_YOUNGSMOD));
	ioInfo.set(param_generic_info + string("<i><b>Poisson's ratio"), INT2(IOI_MESHPARAM, PARAM_POISSONRATIO));
	ioInfo.set(param_generic_info + string("<i><b>Specific heat capacity"), INT2(IOI_MESHPARAM, PARAM_SHC));
	ioInfo.set(param_generic_info + string("<i><b>Electronic specific heat capacity"), INT2(IOI_MESHPARAM, PARAM_SHC_E));
	ioInfo.set(param_generic_info + string("<i><b>Electron coupling constant\n<i><b>2TM : electron-lattice"), INT2(IOI_MESHPARAM, PARAM_G_E));


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

	//Shows a possible temperature dependence : minorId is the type (entry from MATPTDEP_ enum)
	//IOI_MESHPARAMTEMPTYPE

	string paramtemptype_generic_info =
		string("[tc1,1,0,1/tc]<b>Temperature dependence type\n") +
		string("\n[tc1,1,0,1/tc]drag: move to parameter to set\n");

	ioInfo.set(paramtemptype_generic_info + string("<i><b>No temperature dependence"), INT2(IOI_MESHPARAMTEMPTYPE, MATPTDEP_NONE));
	ioInfo.set(paramtemptype_generic_info + string("<i><b>Array"), INT2(IOI_MESHPARAMTEMPTYPE, MATPTDEP_ARRAY));
	ioInfo.set(paramtemptype_generic_info + string("<i><b>Custom equation"), INT2(IOI_MESHPARAMTEMPTYPE, MATPTDEP_EQUATION));

	//Shows a possible generator name for mesh parameter spatial dependence : minorId is the MATPVAR_ enum value, textId is the generator name
	//IOI_MESHPARAMVARGENERATOR 

	string paramvargen_generic_info =
		string("[tc1,1,0,1/tc]<b>Spatial dependence generator") +
		string("\n[tc1,1,0,1/tc]drag: move to parameter to set\n");

	ioInfo.set(paramvargen_generic_info + string("<i><b>Set from png mask file with grayscale.\n<i><b>offset, scale, filename\n<i><b>black = 0, white = 1"), INT2(IOI_MESHPARAMVARGENERATOR, MATPVAR_MASK));
	ioInfo.set(paramvargen_generic_info + string("<i><b>Set from ovf2 file (mapped to current dimensions)."), INT2(IOI_MESHPARAMVARGENERATOR, MATPVAR_OVF2));
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
	ioInfo.set(paramvargen_generic_info + string("<i><b>User text equation.\n<i><b>Scalar or vector equation:\n<i><b>define vector equation for rotations only.\n<i><b>Rotation specified using unit vector:\n<i><b>use direction cosines."), INT2(IOI_MESHPARAMVARGENERATOR, MATPVAR_EQUATION));

	//Shows mesh display option for a given mesh : minorId is the MESHDISPLAY_ value, auxId is the unique mesh id number, textId is the MESHDISPLAY_ handle
	//IOI_MESHDISPLAY

	string meshdisplay_generic_info =
		string("[tc1,1,0,1/tc]<b>Mesh quantity to display") +
		string("\n[tc1,1,0,1/tc]<i>green (foreground): on, red: off</i>") +
		string("\n[tc1,1,0,1/tc]<i>orange (background): on, red: off</i>") +
		string("\n[tc1,1,0,1/tc]left-click: change foreground state\n") +
		string("\n[tc1,1,0,1/tc]right-click: change background state\n");

	ioInfo.set(meshdisplay_generic_info + string("<i><b>Nothing displayed"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_NONE));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Magnetisation"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_MAGNETIZATION));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Magnetisation sub-lattice 2"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_MAGNETIZATION2));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>AF Magnetisation"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_MAGNETIZATION12));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Total effective H field"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_EFFECTIVEFIELD));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Total effective H field sub-lattice 2"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_EFFECTIVEFIELD2));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>AF Total effective H field"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_EFFECTIVEFIELD12));
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
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Mechanical Displacement"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_UDISP));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Strain : xx, yy, zz"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_STRAINDIAG));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Strain : yz, xz, xy"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_STRAINODIAG));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Mesh parameter spatial variation"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_PARAMVAR));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Roughness"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_ROUGHNESS));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Custom, Vectorial"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_CUSTOM_VEC));
	ioInfo.set(meshdisplay_generic_info + string("<i><b>Custom, Scalar"), INT2(IOI_MESHDISPLAY, MESHDISPLAY_CUSTOM_SCA));

	//Shows dual mesh display transparency values : textId is the DBL2 value as a string
	//IOI_MESHDISPLAYTRANSPARENCY,

	string meshdisplay_transparency_info =
		string("[tc1,1,0,1/tc]<b>Transparency values for dual mesh display") +
		string("\n[tc1,1,0,1/tc]<i>Foreground, Background</i>") +
		string("\n[tc1,1,0,1/tc]<i>0: transparent, 1: opaque</i>") +
		string("\n[tc1,1,0,1/tc]double-click: edit\n");

	ioInfo.push_back(meshdisplay_transparency_info, IOI_MESHDISPLAYTRANSPARENCY);

	//Shows mesh display threshold values : textId is the DBL2 value as a string
	//IOI_MESHDISPLAYTHRESHOLDS,

	string meshdisplay_thresholds_info =
		string("[tc1,1,0,1/tc]<b>Threshold values for display") +
		string("\n[tc1,1,0,1/tc]<i>Minimum, Maximum</i>") +
		string("\n[tc1,1,0,1/tc]<i>Set both to 0 to disable.</i>") +
		string("\n[tc1,1,0,1/tc]double-click: edit\n") +
		string("\n[tc1,1,0,1/tc]right-click: clear\n");

	ioInfo.push_back(meshdisplay_thresholds_info, IOI_MESHDISPLAYTHRESHOLDS);

	//Shows mesh display threshold trigger type : auxId is the trigger option
	//IOI_MESHDISPLAYTHRESHOLDTRIG

	string meshdisplay_thresholdtrig_info =
		string("[tc1,1,0,1/tc]<b>Threshold trigger component") +
		string("\n[tc1,1,0,1/tc]<i>X, Y, Z, magnitude</i>") +
		string("\n[tc1,1,0,1/tc]left-click: change state forward") +
		string("\n[tc1,1,0,1/tc]right-click: change state backward\n");

	ioInfo.push_back(meshdisplay_thresholdtrig_info, IOI_MESHDISPLAYTHRESHOLDTRIG);

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

	//Shows mesh vectorial quantity display option : minorId is the unique mesh id number, auxId is the display option
	//IOI_MESHVECREP
	//IOI_SMESHVECREP

	string meshdisplay_option_info =
		string("[tc1,1,0,1/tc]<b>Vectorial quantity display option") +
		string("\n[tc1,1,0,1/tc]<i>full, X, Y, Z, direction, magnitude</i>") +
		string("\n[tc1,1,0,1/tc]left-click: change state forward") +
		string("\n[tc1,1,0,1/tc]right-click: change state backward\n");

	ioInfo.push_back(meshdisplay_option_info, IOI_MESHVECREP);
	ioInfo.push_back(meshdisplay_option_info, IOI_SMESHVECREP);

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


	//Shows SOR damping values when used in fixed damping mode. textId is the DBL2 damping value as a string. (DBL2 since we need different damping values for V and S solvers)
	//IOI_SORDAMPING

	string sordampingvalues_info =
		string("[tc1,1,0,1/tc]<b>SOR fixed damping") +
		string("\n[tc1,1,0,1/tc]<i>SOR damping for (V, S) solvers</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(sordampingvalues_info, IOI_SORDAMPING);

	//Static transport solver state. auxId is the value (0/1)
	//IOI_STATICTRANSPORT:

	string statictransport_info =
		string("[tc1,1,0,1/tc]<b>Static transport solver") +
		string("\n[tc1,1,0,1/tc]<i>If set, transport solver iterated only</i>") +
		string("\n[tc1,1,0,1/tc]<i>at end of a step or stage.</i>") +
		string("\n[tc1,1,0,1/tc]<i>You should also set a high timeout.</i>") +
		string("\n[tc1,1,0,1/tc]click: switch mode\n");

	ioInfo.push_back(statictransport_info, IOI_STATICTRANSPORT);

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

	//Shows atomic moment multiple of Bohr magneton for AF meshes. minorId is the unique mesh id number, auxId is available/not available status (must be antiferromagnetic mesh), textId is the value
	//IOI_ATOMICMOMENT_AFM

	string atomicmoment_afm_info =
		string("[tc1,1,0,1/tc]<b>Magnetic moment in Bohr magnetons") +
		string("\n[tc1,1,0,1/tc]<i>Used to calculate parameter</i>") +
		string("\n[tc1,1,0,1/tc]<i>temperature dependencies when T_Neel > 0</i>") +
		string("\n[tc1,1,0,1/tc]<i>In particular field-dependence is introduced</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(atomicmoment_afm_info, IOI_ATOMICMOMENT_AFM);

	//Shows Tc tau couplings. minorId is the unique mesh id number, auxId is available/not available status (must be antiferromagnetic mesh), textId is the value
	//IOI_TAU

	string tau_info =
		string("[tc1,1,0,1/tc]<b>Tc tau couplings") +
		string("\n[tc1,1,0,1/tc]<i>Used to calculate parameter</i>") +
		string("\n[tc1,1,0,1/tc]<i>temperature dependencies when T_Neel > 0</i>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit\n");

	ioInfo.push_back(tau_info, IOI_TAU);

	//Shows temperature model type for mesh. minorId is the unique mesh id number, auxId is the model identifier (entry from TMTYPE_ enum)
	//IOI_TMODEL

	string tmodel_info =
		string("[tc1,1,0,1/tc]<b>Temperature model type") +
		string("\n[tc1,1,0,1/tc]click: change\n");

	ioInfo.push_back(tmodel_info, IOI_TMODEL);

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
		string("\n[tc1,1,0,1/tc]<i>moments at interface cells will be frozen</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(dipolecouple_info, IOI_COUPLEDTODIPOLESSTATUS);

	//Shows log_errors enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	//IOI_ERRORLOGSTATUS:

	string logerrors_info =
		string("[tc1,1,0,1/tc]<b>Log errors status") +
		string("\n[tc1,1,0,1/tc]<i>green: set, red: not set</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(logerrors_info, IOI_ERRORLOGSTATUS);

	//Shows start_check_updates enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	//IOI_UPDATESTATUSCHECKSTARTUP

	string start_check_updates_info =
		string("[tc1,1,0,1/tc]<b>Check for updates on startup") +
		string("\n[tc1,1,0,1/tc]<i>green: set, red: not set</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(start_check_updates_info, IOI_UPDATESTATUSCHECKSTARTUP);

	//Shows start_scriptserver enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	//IOI_SCRIPTSERVERSTARTUP

	string start_scriptserver_info =
		string("[tc1,1,0,1/tc]<b>Start script server on startup") +
		string("\n[tc1,1,0,1/tc]<i>green: set, red: not set</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(start_scriptserver_info, IOI_SCRIPTSERVERSTARTUP);

	//IOI_MESHEXCHCOUPLING

	string ec_info =
		string("[tc1,1,0,1/tc]<b>Neighboring mesh exchange coupling status") +
		string("\n[tc1,1,0,1/tc]<i>green: set, red: not set</i>") +
		string("\n[tc1,1,0,1/tc]<i>When set, neighboring ferromagnetic meshes</i>") +
		string("\n[tc1,1,0,1/tc]<i>in contact with this one will be exchange coupled.</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(ec_info, IOI_MESHEXCHCOUPLING);

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

	//Adaptive time step control
	//IOI_ODERELERRFAIL, IOI_ODERELERRHIGH, IOI_ODERELERRLOW, IOI_ODEDTINCR, IOI_ODEDTMIN, IOI_ODEDTMAX

	string astep_ctrl_info =
		string("[tc1,1,0,1/tc]<b>Adaptive time step control</b>") +
		string("\n[tc1,1,0,1/tc]Fail above this error.") + 
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(astep_ctrl_info, IOI_ODERELERRFAIL);

	astep_ctrl_info =
		string("[tc1,1,0,1/tc]<b>Adaptive time step control</b>") +
		string("\n[tc1,1,0,1/tc]Decrease dT above this error.") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(astep_ctrl_info, IOI_ODERELERRHIGH);

	astep_ctrl_info =
		string("[tc1,1,0,1/tc]<b>Adaptive time step control</b>") +
		string("\n[tc1,1,0,1/tc]Increase dT below this error.") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(astep_ctrl_info, IOI_ODERELERRLOW);

	astep_ctrl_info =
		string("[tc1,1,0,1/tc]<b>Adaptive time step control</b>") +
		string("\n[tc1,1,0,1/tc]dT increase factor.") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(astep_ctrl_info, IOI_ODEDTINCR);

	astep_ctrl_info =
		string("[tc1,1,0,1/tc]<b>Adaptive time step control</b>") +
		string("\n[tc1,1,0,1/tc]Minimum dT.") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(astep_ctrl_info, IOI_ODEDTMIN);

	astep_ctrl_info =
		string("[tc1,1,0,1/tc]<b>Adaptive time step control</b>") +
		string("\n[tc1,1,0,1/tc]Maximum dT.") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit");

	ioInfo.push_back(astep_ctrl_info, IOI_ODEDTMAX);

	//Shows PBC setting for individual demag modules. minorId is the unique mesh id number, auxId is the pbc images number (0 disables pbc; -1 means setting is not available) (must be ferromagnetic mesh);
	//IOI_PBC_X,
	//IOI_PBC_Y,
	//IOI_PBC_Z

	string pbc_info =
		string("[tc1,1,0,1/tc]<b>Periodic Boundary Conditions") +
		string("\n[tc1,1,0,1/tc]<i>Applicable for magnetization</i>") +
		string("\n[tc1,1,0,1/tc]left-click: set") +
		string("\n[tc1,1,0,1/tc]right-click: clear") + 
		string("\n[tc1,1,0,1/tc]double-click: edit\n");

	ioInfo.push_back(pbc_info, IOI_PBC_X);
	ioInfo.push_back(pbc_info, IOI_PBC_Y);
	ioInfo.push_back(pbc_info, IOI_PBC_Z);

	//Shows PBC setting for supermesh/multilayered demag. auxId is the pbc images number (0 disables pbc; -1 means setting is not available)
	//IOI_SPBC_X
	//IOI_SPBC_Y
	//IOI_SPBC_Z

	ioInfo.push_back(pbc_info, IOI_SPBC_X);
	ioInfo.push_back(pbc_info, IOI_SPBC_Y);
	ioInfo.push_back(pbc_info, IOI_SPBC_Z);

	//Shows individual shape control flag. auxId is the value (0/1)
	//IOI_INDIVIDUALSHAPE

	string IOI_INDIVIDUALSHAPE_info =
		string("[tc1,1,0,1/tc]<b>Individual shape control status flag") +
		string("\n[tc1,1,0,1/tc]<i>When On, shapes are applied only</i>") +
		string("\n[tc1,1,0,1/tc]<i>to primary displayed quantities.</i>") +
		string("\n[tc1,1,0,1/tc]<i>When Off, all primary quantities</i>") +
		string("\n[tc1,1,0,1/tc]<i>are modified.</i>") +
		string("\n[tc1,1,0,1/tc]click: change status\n");

	ioInfo.push_back(IOI_INDIVIDUALSHAPE_info, IOI_INDIVIDUALSHAPE);

	//Shows image cropping settings : textId has the DBL4 value as text
	//IOI_IMAGECROPPING

	string IOI_IMAGECROPPING_info =
		string("[tc1,1,0,1/tc]<b>Image save cropping, normalized.") +
		string("\n[tc1,1,0,1/tc]double-click: edit\n");

	ioInfo.push_back(IOI_IMAGECROPPING_info, IOI_IMAGECROPPING);

	//Show user constant for text equations : minorId is the index in Simulation::userConstants, auxId is the number of the interactive object in the list as it appears in the console, textId is the constant name and value string 
	//Note this entry must always represent the entry in Simulation::userConstants with the index in auxId.
	//IOI_USERCONSTANT

	string showuserconstant_info =
		string("[tc1,1,0,1/tc]<b>User constant</b>") +
		string("\n[tc1,1,0,1/tc]dbl-click: edit entry") +
		string("\n[tc1,1,0,1/tc]right-click: delete entry\n");

	ioInfo.push_back(showuserconstant_info, IOI_USERCONSTANT);
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
		
	case IOI_MESHMCELLSIZE:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			if (SMesh[meshIndex]->MechComputation_Enabled())
				return MakeInteractiveObject(ToString(SMesh[meshIndex]->GetMeshMCellsize(), "m"), IOI_MESHMCELLSIZE, SMesh[meshIndex]->get_id(), 1, ToString(SMesh[meshIndex]->GetMeshMCellsize(), "m"), ONCOLOR);
			else
				return MakeInteractiveObject("N/A", IOI_MESHMCELLSIZE, SMesh[meshIndex]->get_id(), 0, "N/A", OFFCOLOR);
		}
		break;

	case IOI_MESHSCELLSIZE:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			if (SMesh[meshIndex]->MComputation_Enabled())
				return MakeInteractiveObject(ToString(SMesh[meshIndex]->GetMeshSCellsize(), "m"), IOI_MESHSCELLSIZE, SMesh[meshIndex]->get_id(), 1, ToString(SMesh[meshIndex]->GetMeshSCellsize(), "m"), ONCOLOR);
			else
				return MakeInteractiveObject("N/A", IOI_MESHSCELLSIZE, SMesh[meshIndex]->get_id(), 0, "N/A", OFFCOLOR);
		}
		break;

	//Shows link stochastic flag : minorId is the unique mesh id number, auxId is the value off (0), on (1), N/A (-1)
	case IOI_LINKSTOCHASTIC:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			if (SMesh[meshIndex]->MComputation_Enabled())
				return MakeInteractiveObject("On", IOI_LINKSTOCHASTIC, SMesh[meshIndex]->get_id(), 1);
			else
				return MakeInteractiveObject("N/A", IOI_LINKSTOCHASTIC, SMesh[meshIndex]->get_id(), -1, "N/A", OFFCOLOR);
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

	case IOI_MESHDISPLAYTRANSPARENCY:
		if (params_str.size() == 1) {

			return MakeInteractiveObject(params_str[0], IOI_MESHDISPLAYTRANSPARENCY, -1, -1, params_str[0]);
		}
		break;

	case IOI_MESHDISPLAYTHRESHOLDS:
		if (params_str.size() == 1) {

			return MakeInteractiveObject(params_str[0], IOI_MESHDISPLAYTHRESHOLDS, -1, -1, params_str[0]);
		}
		break;

	case IOI_MESHDISPLAYTHRESHOLDTRIG:
		return MakeInteractiveObject("Z", IOI_MESHDISPLAYTHRESHOLDTRIG, -1, (int)VEC3REP_Z);
		break;

	case IOI_MESHVECREP:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			return MakeInteractiveObject("full", IOI_MESHVECREP, SMesh[meshIndex]->get_id(), 0);
		}
		break;

	case IOI_SMESHVECREP:
		return MakeInteractiveObject("full", IOI_SMESHVECREP, 0, 0);
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
		if (params_str.size() == 1) {

			int odeId = ToNum(params_str[0]);

			string odeHandle = odeHandles(odeId);

			return MakeInteractiveObject(odeHandle, IOI_ODE, odeId);
		}
		break;

	case IOI_ODEDT:
		return MakeInteractiveObject("0", IOI_ODEDT, 0, 0, "0");
		break;

	case IOI_HEATDT:
		return MakeInteractiveObject("0", IOI_HEATDT, 0, 0, "0");
		break;

	case IOI_ODE_EVAL:
		if (params_str.size() == 1) {

			int evalId = ToNum(params_str[0]);

			string evalHandle = odeEvalHandles(evalId);

			//mark the set ode with ODE_ERROR, so the state handler will be forced to update the console object with correct state and color
			return MakeInteractiveObject(evalHandle, IOI_ODE_EVAL, ODE_ERROR, evalId, evalHandle);
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

	case IOI_MESHPARAMTEMPTYPE:
		if (params_str.size() == 1) {

			int dependenceID = ToNum(params_str[0]);
			
			string type = temperature_dependence_type(dependenceID);

			return MakeInteractiveObject(type, IOI_MESHPARAMTEMPTYPE, dependenceID, 0, type);
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

	case IOI_MESH_FORTEMPERATURE:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);
			string meshName = SMesh().get_key_from_index(meshIndex);

			return MakeInteractiveObject(meshName, IOI_MESH_FORTEMPERATURE, SMesh[meshIndex]->get_id(), 1, meshName);
		}
		break;

	case IOI_MESH_FORCURIEANDMOMENT:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);
			string meshName = SMesh().get_key_from_index(meshIndex);

			return MakeInteractiveObject(meshName, IOI_MESH_FORCURIEANDMOMENT, SMesh[meshName]->get_id(), 1, meshName);
		}
		break;

	case IOI_MESH_FORTMODEL:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);
			string meshName = SMesh().get_key_from_index(meshIndex);

			return MakeInteractiveObject(meshName, IOI_MESH_FORTMODEL, SMesh[meshName]->get_id(), 1, meshName);
		}
		break;

	case IOI_MESH_FORPBC:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);
			string meshName = SMesh().get_key_from_index(meshIndex);

			return MakeInteractiveObject(meshName, IOI_MESH_FORPBC, SMesh[meshName]->get_id(), 1, meshName);
		}
		break;

	case IOI_MESH_FOREXCHCOUPLING:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);
			string meshName = SMesh().get_key_from_index(meshIndex);

			return MakeInteractiveObject(meshName, IOI_MESH_FOREXCHCOUPLING, SMesh[meshName]->get_id(), 1, meshName);
		}
		break;

	case IOI_MESH_FORSTOCHASTICITY:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);
			string meshName = SMesh().get_key_from_index(meshIndex);

			return MakeInteractiveObject(meshName, IOI_MESH_FORSTOCHASTICITY, SMesh[meshName]->get_id(), 1, meshName);
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

	case IOI_SORDAMPING:
	{
		DBL2 fixed_SOR_damping = SMesh.CallModuleMethod(&STransport::GetSORDamping);

		return MakeInteractiveObject(ToString(fixed_SOR_damping), IOI_SORDAMPING, -1, -1, ToString(fixed_SOR_damping));
	}
	break;

	case IOI_STATICTRANSPORT:
	{
		if (static_transport_solver) return MakeInteractiveObject("On", IOI_STATICTRANSPORT, 0, 1, "", ONCOLOR);
		else return MakeInteractiveObject("Off", IOI_STATICTRANSPORT, 0, 0, "", OFFCOLOR);
	}
	break;

	case IOI_BASETEMPERATURE:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			return MakeInteractiveObject(ToString(SMesh[meshIndex]->GetAverageTemperature(), "K"), IOI_BASETEMPERATURE, SMesh[meshIndex]->get_id(), -1, ToString(SMesh[meshIndex]->GetAverageTemperature(), "K"));
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

	case IOI_ATOMICMOMENT_AFM:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			return MakeInteractiveObject(ToString(SMesh[meshIndex]->GetAtomicMoment_AFM(), "uB"), IOI_ATOMICMOMENT_AFM, SMesh[meshIndex]->get_id(), 1, ToString(SMesh[meshIndex]->GetAtomicMoment_AFM(), "uB"));
		}
		break;

	case IOI_TAU:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			return MakeInteractiveObject(ToString(SMesh[meshIndex]->GetTcCoupling()), IOI_TAU, SMesh[meshIndex]->get_id(), 1, ToString(SMesh[meshIndex]->GetTcCoupling()));
		}
		break;

	case IOI_TMODEL:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			return MakeInteractiveObject(" 1TM ", IOI_TMODEL, SMesh[meshIndex]->get_id(), 1);
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

	//Shows log_errors enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_ERRORLOGSTATUS:
		if (params_str.size() == 1) {

			int status = ToNum(params_str[0]);

			if (status == 0) return MakeInteractiveObject("Off", IOI_ERRORLOGSTATUS, -1, 0, "", OFFCOLOR);
			else return MakeInteractiveObject("On", IOI_ERRORLOGSTATUS, -1, 1, "", ONCOLOR);
		}
		break;

	//Shows start_check_updates enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_UPDATESTATUSCHECKSTARTUP:
		if (params_str.size() == 1) {

			int status = ToNum(params_str[0]);

			if (status == 0) return MakeInteractiveObject("Off", IOI_UPDATESTATUSCHECKSTARTUP, -1, 0, "", OFFCOLOR);
			else return MakeInteractiveObject("On", IOI_UPDATESTATUSCHECKSTARTUP, -1, 1, "", ONCOLOR);
		}
		break;

	//Shows start_scriptserver enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	case IOI_SCRIPTSERVERSTARTUP:
		if (params_str.size() == 1) {

			int status = ToNum(params_str[0]);

			if (status == 0) return MakeInteractiveObject("Off", IOI_SCRIPTSERVERSTARTUP, -1, 0, "", OFFCOLOR);
			else return MakeInteractiveObject("On", IOI_SCRIPTSERVERSTARTUP, -1, 1, "", ONCOLOR);
		}
		break;

	case IOI_MESHEXCHCOUPLING:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			if (SMesh[meshIndex]->MComputation_Enabled()) {

				return MakeInteractiveObject("Off", IOI_MESHEXCHCOUPLING, SMesh[meshIndex]->get_id(), 0, "", OFFCOLOR);
			}
			else {

				return MakeInteractiveObject("N/A", IOI_MESHEXCHCOUPLING, SMesh[meshIndex]->get_id(), -1, "", UNAVAILABLECOLOR);
			}
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

	case IOI_ODERELERRFAIL:
		if (params_str.size() == 1) {

			string value = params_str[0];

			return MakeInteractiveObject(value, IOI_ODERELERRFAIL, 0, 0, value);
		}
		break;

	case IOI_ODERELERRHIGH:
		if (params_str.size() == 1) {

			string value = params_str[0];

			return MakeInteractiveObject(value, IOI_ODERELERRHIGH, 0, 0, value);
		}
		break;

	case IOI_ODERELERRLOW:
		if (params_str.size() == 1) {

			string value = params_str[0];

			return MakeInteractiveObject(value, IOI_ODERELERRLOW, 0, 0, value);
		}
		break;

	case IOI_ODEDTINCR:
		if (params_str.size() == 1) {

			string value = params_str[0];

			return MakeInteractiveObject(value, IOI_ODEDTINCR, 0, 0, value);
		}
		break;

	case IOI_ODEDTMIN:
		if (params_str.size() == 1) {

			string value = params_str[0];

			return MakeInteractiveObject(value, IOI_ODEDTMIN, 0, 0, value);
		}
		break;

	case IOI_ODEDTMAX:
		if (params_str.size() == 1) {

			string value = params_str[0];

			return MakeInteractiveObject(value, IOI_ODEDTMAX, 0, 0, value);
		}
		break;

	case IOI_PBC_X:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			return MakeInteractiveObject("N/A", IOI_PBC_X, SMesh[meshIndex]->get_id(), -1, "", UNAVAILABLECOLOR);
		}
		break;

	case IOI_PBC_Y:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			return MakeInteractiveObject("N/A", IOI_PBC_Y, SMesh[meshIndex]->get_id(), -1, "", UNAVAILABLECOLOR);
		}
		break;

	case IOI_PBC_Z:
		if (params_str.size() == 1) {

			int meshIndex = ToNum(params_str[0]);

			return MakeInteractiveObject("N/A", IOI_PBC_Z, SMesh[meshIndex]->get_id(), -1, "", UNAVAILABLECOLOR);
		}
		break;

	case IOI_SPBC_X:
		return MakeInteractiveObject("N/A", IOI_SPBC_X, 0, -1, "", UNAVAILABLECOLOR);
		break;

	case IOI_SPBC_Y:
		return MakeInteractiveObject("N/A", IOI_SPBC_Y, 0, -1, "", UNAVAILABLECOLOR);
		break;

	case IOI_SPBC_Z:
		return MakeInteractiveObject("N/A", IOI_SPBC_Z, 0, -1, "", UNAVAILABLECOLOR);
		break;

	case IOI_INDIVIDUALSHAPE:
		if (params_str.size() == 1) {

			int status = ToNum(params_str[0]);

			if (status) {

				return MakeInteractiveObject("On", IOI_INDIVIDUALSHAPE, 0, 1, "", ONCOLOR);
			}
			else {

				return MakeInteractiveObject("Off", IOI_INDIVIDUALSHAPE, 0, 0, "", OFFCOLOR);
			}
		}
		break;

	case IOI_IMAGECROPPING:
		return MakeInteractiveObject("0, 0, 1, 1", IOI_IMAGECROPPING, 0, 0, "0, 0, 1, 1");
		break;
	
	case IOI_USERCONSTANT:
		if (params_str.size() == 3) {

			string constant_name = params_str[0];
			double value = ToNum(params_str[1]);
			int index_in_list = ToNum(params_str[2]);

			string userConstant_text = Build_EquationConstants_Text(index_in_list);

			return MakeInteractiveObject(userConstant_text, IOI_USERCONSTANT, index_in_list, index_in_list, userConstant_text);
		}
		break;

	}

	return "";
}
