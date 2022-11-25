#pragma once



//Interactive console object identifiers - this is the majorId in InteractiveObjectProperties
enum IOI_ 
{ 	
	//Shows program version update status : auxIdis the status as -1: attempting to connect, 0: connection failure, 1: program up to date, 2: update available
	IOI_PROGRAMUPDATESTATUS,

	//Shows start_check_updates enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	IOI_UPDATESTATUSCHECKSTARTUP,

	//Shows start_scriptserver enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	IOI_SCRIPTSERVERSTARTUP,

	//Data box entry, showing the label of a given entry in Simulation::dataBoxList : minorId is the minor id of elements in Simulation::dataBoxList (major id there is always 0), auxId is the number of the interactive object in the list (i.e. entry number as it appears in data box in order). 
	//textId is the mesh name (if associated with this data type)
	//Note this entry must always represent the entry in Simulation::dataBoxList with the index in auxId.
	IOI_DATABOXFIELDLABEL,

	//A set or available module for a given mesh: minorId in InteractiveObjectProperties is an entry from MOD_ enum identifying the module, auxId contains the unique mesh id number this module refers to
	IOI_MODULE,

	//super-mesh module : minor type is an entry from MOD_ enum
	IOI_SMODULE,

	//Module used for effective field display for a given mesh: minorId in InteractiveObjectProperties is an entry from MOD_ enum identifying the module, auxId contains the unique mesh id number this module refers to, textId is the MOD_ value used for display
	IOI_DISPLAYMODULE,

	//Available/set ode : minorId is an entry from ODE_ (the equation)
	IOI_ODE,

	//Available/set ode for atomistic meshes: minorId is an entry from ODE_ (the equation)
	IOI_ATOMODE,

	//Set ODE time step: textId is the value
	IOI_ODEDT,

	//Set heat equation time step: textId is the value
	IOI_HEATDT,

	//Set elastodynamics equation time step: textId is the value
	IOI_ELDT,

	//Link elastodynamics time-step to ODE dT flag : auxId is the value
	IOI_LINKELDT,

	//Set stochastic time-step: textId is the value
	IOI_STOCHDT,

	//Link stochastic time-step to ODE dT flag : auxId is the value
	IOI_LINKSTOCHDT,

	//Set evaluation speedup time-step: textId is the value
	IOI_SPEEDUPDT,

	//Link evaluation speedup time-step to ODE dT flag : auxId is the value
	IOI_LINKSPEEDUPDT,

	//Available/set evaluation method for ode : minorId is an entry from ODE_ as : micromagnetic equation value + 100 * atomistic equation value, auxId is the EVAL_ entry (the evaluation method), textId is the name of the evaluation method
	IOI_ODE_EVAL,

	//Shows a mesh name : minorId is the unique mesh id number, textId is the mesh name (below are similar objects but used in different lists, so these lists need updating differently).
	//auxId is also used : value 1 means update list, value 0 means do not update list (but delete line if mesh is deleted).
	IOI_MESH_FORPARAMS,
	IOI_MESH_FORPARAMSTEMP,
	IOI_MESH_FORPARAMSVAR,
	IOI_MESH_FORMODULES,
	IOI_MESH_FORDISPLAYMODULES,
	IOI_MESH_FORMESHLIST,
	IOI_MESH_FORDISPLAYOPTIONS,
	IOI_MESH_FORTEMPERATURE,
	IOI_MESH_FORHEATBOUNDARIES,
	IOI_MESH_FORCURIEANDMOMENT,
	IOI_MESH_FORTMODEL,
	IOI_MESH_FORPBC,
	IOI_MESH_FOREXCHCOUPLING, 
	IOI_MESH_FORSTOCHASTICITY,
	IOI_MESH_FORELASTICITY,
	IOI_MESH_FORSPEEDUP,
	IOI_MESH_FORSKYPOSDMUL,
	IOI_MESH_FORMC,
	IOI_MESH_FORDIPOLESHIFT,
	IOI_MESH_FORTMR,

	//Shows ferromagnetic super-mesh rectangle (unit m) : textId is the mesh rectangle for the ferromagnetic super-mesh
	IOI_FMSMESHRECTANGLE,
	//Shows electric super-mesh rectangle (unit m) : textId is the mesh rectangle for the electric super-mesh
	IOI_ESMESHRECTANGLE,
	
	//Shows ferromagnetic super-mesh cellsize (units m) : textId is the mesh cellsize for the ferromagnetic super-mesh
	IOI_FMSMESHCELLSIZE,
	//Shows electric super-mesh cellsize (units m) : textId is the mesh cellsize for the electric super-mesh
	IOI_ESMESHCELLSIZE,
	
	//Shows mesh rectangle (units m) : minorId is the unique mesh id number, textId is the mesh rectangle
	IOI_MESHRECTANGLE,
	
	//Shows magnetic mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	IOI_MESHCELLSIZE,
	//Shows electric mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	IOI_MESHECELLSIZE,
	//Shows thermal mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	IOI_MESHTCELLSIZE,
	//Shows mechanical mesh cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	IOI_MESHMCELLSIZE,
	//Shows stochastic cellsize (units m) : minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	IOI_MESHSCELLSIZE,
	//Shows link stochastic flag : minorId is the unique mesh id number, auxId is the value off (0), on (1), N/A (-1)
	IOI_LINKSTOCHASTIC,

	//Shows macrocell size (units m) for atomistic meshes: minorId is the unique mesh id number, auxId is enabled/disabled status, textId is the mesh cellsize
	IOI_MESHDMCELLSIZE,

	//Shows evaluation speedup type: auxId is the type value.
	IOI_SPEEDUPMODE,

	//Simulation output data, specifically used for showing values in console : minorId is the DATA_ id, textId is the data handle
	IOI_SHOWDATA,

	//Simulation output data, specifically used to construct output data list : minorId is the DATA_ id, textId is the data handle
	IOI_DATA,

	//Show currently set directory : textId is the directory
	IOI_DIRECTORY,

	//Show currently set save data file : textId is the file name
	IOI_SAVEDATAFILE,

	//Show currently set image filename base : textId is the file name
	IOI_SAVEIMAGEFILEBASE,

	//Show flag status for data/image saving during a simulation : minorId is the flag value (boolean)
	IOI_SAVEDATAFLAG,
	IOI_SAVEIMAGEFLAG,

	//Show set output data : minorId is the minor id of elements in Simulation::saveDataList (major id there is always 0), auxId is the number of the interactive object in the list as it appears in the console, textId is the configured output data. 
	//Note this entry must always represent the entry in Simulation::saveDataList with the index in auxId.
	IOI_OUTDATA,

	//Shows a possible stage type, used for adding generic stages to the simulation schedule : minorId is the stage type (SS_ enum value, which is the majorId from stageDescriptors), textId is the stage setting handle
	IOI_STAGE,

	//Shows a stage added to the simulation schedule : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the configured stage text
	//Note this entry must always represent the entry in Simulation::simStages with the index in auxId.
	IOI_SETSTAGE,

	//Shows the value to set for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the value as a std::string
	IOI_SETSTAGEVALUE,

	//Shows the stop condition for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the number of the interactive object in the list, textId is the stop type and value as a std::string
	IOI_STAGESTOPCONDITION,

	//Shows the saving condition for the simulation schedule stage : minorId is the minor id of elements in Simulation::simStages (major id there is always 0), auxId is the DSAVE_ value for this data save type, textId is the save type and value as a std::string
	IOI_DSAVETYPE,

	//Shows a stop condition, used to apply the same condition to all simulation stages : minorId is the STOP_ value, textId is the stop type handle
	IOI_STAGESTOPCONDITIONALL,

	//Shows a data save condition, used to apply the same condition to all simulation stages : minorId is the DSAVE_ value, textId is the save type handle
	IOI_DSAVETYPEALL,

	//Shows parameter and value for a given mesh : minorId is the major id of elements in SimParams::simParams (i.e. an entry from PARAM_ enum), auxId is the unique mesh id number, textId is the parameter handle and value
	IOI_MESHPARAM,

	//Shows parameter temperature dependence for a given mesh : minorId is the major id of elements in SimParams::simParams (i.e. an entry from PARAM_ enum), auxId is the unique mesh id number, textId is the parameter temperature dependence setting
	IOI_MESHPARAMTEMP,

	//Shows parameter spatial dependence for a given mesh : minorId is the major id of elements in SimParams::simParams (i.e. an entry from PARAM_ enum), auxId is the unique mesh id number, textId is the parameter spatial dependence setting
	IOI_MESHPARAMVAR,

	//Shows a possible temperature dependence : minorId is the type (entry from MATPTDEP_ enum)
	IOI_MESHPARAMTEMPTYPE,

	//Shows a possible generator name for mesh parameter spatial dependence : minorId is the MATPVAR_ enum value, textId is the generator name
	IOI_MESHPARAMVARGENERATOR,

	//Shows mesh display option for a given mesh : minorId is the MESHDISPLAY_ value, auxId is the unique mesh id number, textId is the MESHDISPLAY_ handle
	IOI_MESHDISPLAY,
	//Shows super-mesh display option : minorId is the MESHDISPLAY_ value, textId is the MESHDISPLAY_ handle
	IOI_SMESHDISPLAY,

	//Shows dual mesh display transparency values : textId is the DBL2 value as a std::string
	IOI_MESHDISPLAYTRANSPARENCY,
	
	//Shows mesh display threshold values : textId is the DBL2 value as a std::string
	IOI_MESHDISPLAYTHRESHOLDS,

	//Shows mesh display threshold trigger type : auxId is the trigger option
	IOI_MESHDISPLAYTHRESHOLDTRIG,

	//Shows mesh vectorial quantity display option : minorId is the unique mesh id number, auxId is the display option
	IOI_MESHVECREP,
	//Shows supermesh vectorial quantity display option : auxId is the display option
	IOI_SMESHVECREP,

	//Shows movingmesh trigger settings : minorId is the unique mesh id number (if set), auxId is the trigger state (used or not used), textId is the mesh name (if set)
	IOI_MOVINGMESH,

	//Shows movingmesh symmetry : auxId is the asymmetry status (1: asymmetric, 0: symmetric)
	IOI_MOVINGMESHASYM,

	//Shows movingmesh threshold : textId is the threshold value as a std::string
	IOI_MOVINGMESHTHRESH,

	//Shows electrode box. minorId is the minor Id in STransport::electrode_boxes, auxId is the number of the interactive object in the list (electrode index), textId is the electrode rect as a std::string
	IOI_ELECTRODERECT,

	//Shows electrode potential. minorId is the electrode index, textId is potential value as a std::string
	IOI_ELECTRODEPOTENTIAL,

	//Shows electrode ground setting. minorId is the electrode index, auxId is the setting (0 : not ground, 1 : ground)
	IOI_ELECTRODEGROUND,

	//Shows constant current source setting. auxId is the setting.
	IOI_CONSTANTCURRENTSOURCE,

	//Shows transport solver convergence error. textId is the convergence error value.
	IOI_TSOLVERCONVERROR,

	//Shows transport solver timeout iterations. auxId is the timeout value.
	IOI_TSOLVERTIMEOUT,

	//Shows spin transport solver convergence error. textId is the convergence error value.
	IOI_SSOLVERCONVERROR,

	//Shows spin transport solver timeout iterations. auxId is the timeout value.
	IOI_SSOLVERTIMEOUT,

	//Shows SOR damping values when used in fixed damping mode. textId is the DBL2 damping value as a std::string. (DBL2 since we need different damping values for V and S solvers)
	IOI_SORDAMPING,

	//Shows tmr type setting. minorId is the unique mesh id number, auxId is the value.
	IOI_TMRTYPE,

	//Shows mesh temperature. minorId is the unique mesh id number, textId is the temperature value
	IOI_BASETEMPERATURE,

	//Shows ambient temperature for heat equation Robin boundary conditions. minorId is the unique mesh id number, auxId is enabled/disabled status (Heat module must be active), textId is the temperature value
	IOI_AMBIENT_TEMPERATURE,

	//Shows alpha value (W/m^2K) for heat equation Robin boundary conditions. minorId is the unique mesh id number, auxId is enabled/disabled status (Heat module must be active), textId is the value
	IOI_ROBIN_ALPHA,

	//Shows temperature insulating side setting for heat equation. minorId is the unique mesh id number, auxId is the status (Heat module must be active) : -1 disabled (gray), 0 not insulating (green), 1 insulating (red), textId represents the side : "x", "-x", "y", "-y", "z", "-z"
	IOI_INSULATINGSIDE,

	//Shows mesh Curie temperature. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the temperature value
	IOI_CURIETEMP,

	//Shows indicative material Curie temperature. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the temperature value
	IOI_CURIETEMPMATERIAL,

	//Shows atomic moment multiple of Bohr magneton. minorId is the unique mesh id number, auxId is available/not available status (must be ferromagnetic mesh), textId is the value
	IOI_ATOMICMOMENT,

	//Shows atomic moment multiple of Bohr magneton for AF meshes. minorId is the unique mesh id number, auxId is available/not available status (must be antiferromagnetic mesh), textId is the value
	IOI_ATOMICMOMENT_AFM,

	//Shows Tc tau couplings. minorId is the unique mesh id number, auxId is available/not available status (must be antiferromagnetic mesh), textId is the value
	IOI_TAU,

	//Shows temperature model type for mesh. minorId is the unique mesh id number, auxId is the model identifier (entry from TMTYPE_ enum)
	IOI_TMODEL,

	//Shows cuda enabled/disabled or n/a state. auxId is enabled (1)/disabled(0)/not available(-1) status.
	IOI_CUDASTATE,

	//Shows CUDA device information and state. minorId is the device number (from 1 up), auxId is enabled (1)/disabled(0)/not available(-1) status. 
	IOI_CUDADEVICE,

	//Shows gpu free memory. auxId is the value
	IOI_GPUMEMFREE,

	//Shows gpu total memory. auxId is the value
	IOI_GPUMEMTOTAL,

	//Shows cpu free memory. auxId is the value
	IOI_CPUMEMFREE,

	//Shows cpu total memory. auxId is the value
	IOI_CPUMEMTOTAL,

	//Shows scale_rects enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	IOI_SCALERECTSSTATUS,
		
	//Shows coupled_to_dipoles enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	IOI_COUPLEDTODIPOLESSTATUS,

	//Shows dipole velocity value. minorId is the unique mesh id number. textId is the value
	IOI_DIPOLEVELOCITY,
	//Shows dipole shift clipping value. minorId is the unique mesh id number. textId is the value
	IOI_DIPOLESHIFTCLIP,

	//Shows diagonal strain set equation. minorId is the unique mesh id number. textId is the equation. auxId is enabled(1)/disabled(0) status.
	IOI_STRAINEQUATION,
	//Shows shear strain set equation. minorId is the unique mesh id number. textId is the equation. auxId is enabled(1)/disabled(0) status.
	IOI_SHEARSTRAINEQUATION,

	//Shows fixed surface rectangle. minorId is the minor Id in SMElastic::fixed_u_surfaces, auxId is the number of the interactive object in the list (electrode index), textId is the surface rect as a std::string
	IOI_SURFACEFIX,

	//Shows stress surface rectangle. minorId is the minor Id in SMElastic::stress_surfaces_rect, auxId is the number of the interactive object in the list (electrode index), textId is the surface rect as a std::string
	IOI_SURFACESTRESS,

	//Shows stress surface equation. minorId is the index in SMElastic::stress_surfaces_equations, auxId is the number of the interactive object in the list (electrode index), textId is the equation
	IOI_SURFACESTRESSEQ,

	//Shows log_errors enabled/disabled state. auxId is enabled (1)/disabled(0) status.
	IOI_ERRORLOGSTATUS,

	//Shows number of threads. auxId is the value.
	IOI_THREADS,

	//Shows server port. auxId is the value.
	IOI_SERVERPORT,
		
	//Shows server password. textId is the password.
	IOI_SERVERPWD,

	//Shows server sleep time in ms. auxId is the value.
	IOI_SERVERSLEEPMS,

	//Shows neighboring meshes exchange coupling setting for this mesh. minorId is the unique mesh id number, auxId is the status (1/0 : on/off, -1 : not available: must be ferromagnetic mesh)
	IOI_MESHEXCHCOUPLING,

	//Shows mesh roughness refinement value. minorId is the unique mesh id number, auxId is enabled (1)/disabled(0) status. textId is the value
	IOI_REFINEROUGHNESS,

	//Shows status of multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	IOI_MULTICONV,

	//Shows status of force 2D multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	IOI_2DMULTICONV,

	//Shows status of use default n for multi-layered convolution. auxId is the status (-1 : N/A, 0 : Off, 1 : On)
	IOI_NCOMMONSTATUS,

	//Shows n_common for multi-layered convolution. auxId is the status (-1 : N/A, otherwise available). textId is the value as a SZ3.
	IOI_NCOMMON,

	//Shows status of gpu kernels demag initialization. auxId is the status (0 : Off, 1 : On)
	IOI_GPUKERNELS,

	//Shows materials database in use. textId is the name of the database, including the path.
	IOI_LOCALMDB,

	//Shows relative error fail threshold for ode eval. textId is the value.
	IOI_ODERELERRFAIL,

	//Shows dT increase factor. textId is the value.
	IOI_ODEDTINCR,

	//Shows minimum dT value. textId is the value.
	IOI_ODEDTMIN,

	//Shows maximum dT value. textId is the value.
	IOI_ODEDTMAX,

	//Shows PBC setting for individual demag modules. minorId is the unique mesh id number, auxId is the pbc images number (0 disables pbc) (must be ferromagnetic mesh);
	IOI_PBC_X,
	IOI_PBC_Y,
	IOI_PBC_Z,

	//Shows PBC setting for supermesh/multilayered demag. auxId is the pbc images number (0 disables pbc)
	IOI_SPBC_X,
	IOI_SPBC_Y,
	IOI_SPBC_Z,

	//Shows individual shape control flag. auxId is the value (0/1)
	IOI_INDIVIDUALSHAPE,

	//Static transport solver state. auxId is the value (0/1)
	IOI_STATICTRANSPORT,

	//Disabled transport solver state. auxId is the value (0/1)
	IOI_DISABLEDTRANSPORT,

	//Shows image cropping settings : textId has the DBL4 value as text
	IOI_IMAGECROPPING,

	//Show user constant for text equations : minorId is the index in Simulation::userConstants, auxId is the number of the interactive object in the list as it appears in the console, textId is the constant name and value std::string 
	//Note this entry must always represent the entry in Simulation::userConstants with the index in auxId.
	IOI_USERCONSTANT,

	//Show skypos diameter multiplier : minorId is the unique mesh id number, textId is the multiplier as a std::string
	IOI_SKYPOSDMUL,

	//Shows dwpos fitting component. auxId is the value (-1, 0, 1, 2)
	IOI_DWPOSCOMPONENT,

	//Shows Monte-Carlo computation type (serial/parallel) : minorId is the unique mesh id number, auxId is the status (0 : parallel, 1 : serial, -1 : N/A)
	IOI_MCCOMPUTATION,
	
	//Shows Monte-Carlo algorithm type : minorId is the unique mesh id number, auxId is the type (0 : classical, 1 : constrained, -1 : N/A), textId is the constrained DBL3 direction.
	IOI_MCTYPE,

	//Shows Monte-Carlo disabled/enabled status : minorId is the unique mesh id number, auxId is the status (1 : disabled, 0 : enabled).
	IOI_MCDISABLED,

	//Shows Monte-Carlo computefields state flag : auxId is the state (0: disabled, 1: enabled)
	IOI_MCCOMPUTEFIELDS,

	//Shows shape rotation setting: textId is the value as text (DBL3)
	IOI_SHAPEROT,

	//Shows shape repetition setting: textId is the value as text (INT3)
	IOI_SHAPEREP,

	//Shows shape displacement setting: textId is the value as text (DBL3)
	IOI_SHAPEDSP,

	//Shows shape method setting: textId is the value as text (method)
	IOI_SHAPEMET,

	//Shows display render detail level: textId is the value
	IOI_DISPRENDER_DETAIL,

	//Shows display render threshold 1: auxId is the value
	IOI_DISPRENDER_THRESH1,

	//Shows display render threshold 2: auxId is the value
	IOI_DISPRENDER_THRESH2,

	//Shows display render threshold 3: auxId is the value
	IOI_DISPRENDER_THRESH3
};