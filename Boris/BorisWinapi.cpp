#include "stdafx.h"
#include "Boris.h"

#include "CompileFlags.h"
#if GRAPHICS == 1

#include <Commdlg.h>

#define MAX_LOADSTRING 2000

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Separate function to make the Simulation object, so it can be called on a separate std::thread.
//See https://stackoverflow.com/questions/64988525/bug-related-to-g-openmp-when-using-stdthread/65001887#65001887 for reason
void make_Simulation(HWND hWnd)
{
	pSim = new Simulation(hWnd, Program_Version, server_port, server_pwd, cudaDevice);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Global Variables:
HINSTANCE hInst;								// current instance
WCHAR szTitle[MAX_LOADSTRING];	 				// The title bar text
WCHAR szWindowClass[MAX_LOADSTRING];			// the main window class name

HMENU hMenubar;

// Forward declarations of functions included in this code module:
ATOM				MyRegisterClass(HINSTANCE hInstance);
BOOL				InitInstance(HINSTANCE, int);
LRESULT CALLBACK	WndProc(HWND, UINT, WPARAM, LPARAM);
INT_PTR CALLBACK	About(HWND, UINT, WPARAM, LPARAM);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//----------------------MENU ITEMS

enum IDM_ {
	//FILE
	IDM_FILE_DEFAULT, IDM_FILE_LOAD, IDM_FILE_LOADEXAMPLE, IDM_FILE_SAVE, IDM_FILE_SAVEAS, IDM_FILE_SAVEOVF2, IDM_FILE_SAVEIMAGE, IDM_FILE_MAKEVIDEO, IDM_FILE_IMAGECROPPING, IDM_FILE_CHDIR, IDM_FILE_QUIT,
	//SIMULATION
	IDM_SIM_RUN, IDM_SIM_COMPUTEFIELDS, IDM_SIM_STOP, IDM_SIM_RESET, IDM_SIM_CUDA, IDM_SIM_SERVER, IDM_SIM_SCHEDULE, IDM_SIM_DATA, IDM_SIM_SHOWDATA, IDM_SIM_BENCHTIME,
	//MESH
	IDM_MESH_SHOW, IDM_MESH_LOADMASK, IDM_MESH_MASKALL, IDM_MESH_SCALERECTS, IDM_MESH_RESET, IDM_MESH_ADDMATERIAL, IDM_MESH_ADDFERROMAGNET, IDM_MESH_ADDANTIFERROMAGNET, IDM_MESH_ADDMETAL, IDM_MESH_ADDINSULATOR, IDM_MESH_ADDDIPOLE, IDM_MESH_ATOMSIMPLECUBIC, IDM_MESH_EXCHANGECOUPLED, IDM_MESH_COUPLETODIPOLES, IDM_MESH_COPY,
	//SHAPES
	IDM_SHAPES_MODIFIERS, IDM_SHAPES_DELRECT, IDM_SHAPES_RECT, IDM_SHAPES_DISK, IDM_SHAPES_TRIANGLE, IDM_SHAPES_TETRAHEDRON, IDM_SHAPES_PYRAMID, IDM_SHAPES_CONE, IDM_SHAPES_ELLIPSOID, IDM_SHAPES_TORUS,
	//CONFIGURATION
	IDM_CFG_MODULES, IDM_CFG_ODE, IDM_CFG_TIMESTEP, IDM_CFG_ASTEP, IDM_CFG_MULTICONV, IDM_CFG_PBC,
	//PARAMETERS
	IDM_PARAM_SHOW, IDM_PARAM_TDEP, IDM_PARAM_SDEP, IDM_PARAM_COPY,
	//MAGNETIZATION
	IDM_MAG_SET, IDM_MAG_SETOBJECT, IDM_MAG_FIELD, IDM_MAG_DWALL, IDM_MAG_VORTEX, IDM_MAG_SKYRMION, IDM_MAG_BLOCHSKYRMION, IDM_MAG_RANDOM, IDM_MAG_2DGRAINS, IDM_MAG_3DGRAINS, IDM_MAG_LOADOVF2MAG, IDM_MAG_SAVEOVF2MAG, IDM_MAG_ADDRECT, IDM_MAG_DELRECT, IDM_MAG_SETRECT, IDM_MAG_INVERT, IDM_MAG_MIRROR,
	//ROUGHNESS
	IDM_ROUGH_SHOW, IDM_ROUGH_EDIT, IDM_ROUGH_ROUGHEN, IDM_ROUGH_SURFROUGHEN, IDM_ROUGH_LOAD, IDM_ROUGH_CLEARROUGH,
	//TEMPERATURE
	IDM_TEMPERATURE_SET, IDM_TEMPERATURE_CURIE, IDM_TEMPERATURE_AMBIENT, IDM_TEMPERATURE_LOADOVF2, IDM_TEMPERATURE_TIMESTEP,
	//TRANSPORT
	IDM_TRANS_DISABLE, IDM_TRANS_DEFAULTELECTRODES, IDM_TRANS_ELECTRODES, IDM_TRANS_ADDELECTRODE, IDM_TRANS_CLEARELECTRODES, IDM_TRANS_SETV, IDM_TRANS_SETI, IDM_TRANS_SETJC, IDM_TRANS_LOADOVF2, IDM_TRANS_CONFIG,
	//MELASTIC
	IDM_MELASTIC_SETSTRESS, IDM_MELASTIC_LOADOVF2DISP, IDM_MELASTIC_LOADOVF2STRAIN,
	//VIEW
	IDM_VIEW_DISPLAY, IDM_VIEW_UPDATEFREQUENCY, IDM_VIEW_RENDER, IDM_VIEW_CENTER, IDM_VIEW_CLEARCONSOLE,
	//ALGORITHMS
	IDM_ALGO_MOVINGMESH, IDM_ALGO_BLOCHMOVINGMESH, IDM_ALGO_NEELMOVINGMESH, IDM_ALGO_SKYMOVINGMESH, IDM_ALGO_CLEARMOVINGMESH, IDM_ALGO_MOVINGMESHENABLED,
	//DATA PROCESSING
	IDM_DP_CLEARALL, IDM_DP_SHOWSIZES, IDM_DP_LOAD, IDM_DP_SAVE, IDM_DP_GETPROFILE, IDM_AVERAGEMESHRECT, IDM_DP_TOPOCHARGE, IDM_DP_SKYRMIONCOUNT, IDM_DP_MUL, IDM_DP_DOTPROD, IDM_DP_ADDDP, IDM_DP_SUBDP, IDM_DP_MINMAX, IDM_DP_MEAN, IDM_DP_LINREG, IDM_DP_COERCIVITY, IDM_DP_REMANENCE, IDM_DP_COMPLETEHYSTER, IDM_DP_DUMPTDEP,
	IDM_DP_FITLORENTZ, IDM_DP_FITSKYRMION, IDM_DP_FITDW, IDM_DP_FITSTT, IDM_DP_FITSOT, IDM_DP_FITSOTSTT, IDM_DP_FITADIA, IDM_DP_FITNONADIA, IDM_DP_REMOVEOFFSET, IDM_DP_CARTESIANTOPOLAR, IDM_DP_SMOOTH,
	//ABOUT
	IDM_ABOUT_MANUAL, IDM_ABOUT_MDB, IDM_ABOUT_UPDATEMDB, IDM_ABOUT_CHECKUPDATES, IDM_ABOUT_STARTUPOPTIONS, IDM_ABOUT_ABOUT
};

void AddMenus(HWND hwnd)
{
	hMenubar = CreateMenu();

	//FILE
	HMENU hMenu_File = CreateMenu();

	//---
	AppendMenuW(hMenu_File, MF_STRING, IDM_FILE_DEFAULT,       L"&Default\tCtrl-N");
	AppendMenuW(hMenu_File, MF_STRING, IDM_FILE_CHDIR,         L"&Set Directory\tCtrl-D");
	AppendMenuW(hMenu_File, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_File, MF_STRING, IDM_FILE_LOAD,          L"&Load\tCtrl-O");
	AppendMenuW(hMenu_File, MF_STRING, IDM_FILE_LOADEXAMPLE, L"&Load Example");
	AppendMenuW(hMenu_File, MF_STRING, IDM_FILE_SAVE,          L"&Save\tCtrl-S");
	AppendMenuW(hMenu_File, MF_STRING, IDM_FILE_SAVEAS,        L"&Save As\tCtrl-A");
	AppendMenuW(hMenu_File, MF_STRING, IDM_FILE_SAVEOVF2, L"&Save OVF2");
	AppendMenuW(hMenu_File, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_File, MF_STRING, IDM_FILE_SAVEIMAGE,     L"&Save Image\tCtrl-I");
	AppendMenuW(hMenu_File, MF_STRING, IDM_FILE_MAKEVIDEO,     L"&Make Video");
	AppendMenuW(hMenu_File, MF_STRING, IDM_FILE_IMAGECROPPING, L"&Set Cropping");
	AppendMenuW(hMenu_File, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_File, MF_STRING, IDM_FILE_QUIT,          L"&Quit\tAlt-F4");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_File, L"&File");

	//SIMULATION
	HMENU hMenu_Sim = CreateMenu();

	//---
	AppendMenuW(hMenu_Sim, MF_STRING, IDM_SIM_RUN, L"&Run\tF5");
	AppendMenuW(hMenu_Sim, MF_STRING, IDM_SIM_COMPUTEFIELDS, L"&Compute Fields");
	AppendMenuW(hMenu_Sim, MF_STRING, IDM_SIM_STOP, L"&Stop\tF6");
	AppendMenuW(hMenu_Sim, MF_STRING, IDM_SIM_RESET, L"&Reset\tF7");
	AppendMenuW(hMenu_Sim, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Sim, MF_STRING, IDM_SIM_CUDA, L"&CUDA");
	if (pSim) {

		int status = pSim->GetCUDAStatus();

		if (status == -1) EnableMenuItem(hMenu_Sim, IDM_SIM_CUDA, MF_DISABLED);
		else if (status == 1) CheckMenuItem(hMenu_Sim, IDM_SIM_CUDA, MF_CHECKED);
	}
	AppendMenuW(hMenu_Sim, MF_STRING, IDM_SIM_SERVER, L"&Script Server");
	AppendMenuW(hMenu_Sim, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Sim, MF_STRING, IDM_SIM_SCHEDULE, L"&Schedule\tAlt-S");
	AppendMenuW(hMenu_Sim, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Sim, MF_STRING, IDM_SIM_DATA, L"&Output Data\tAlt-D");
	AppendMenuW(hMenu_Sim, MF_STRING, IDM_SIM_SHOWDATA, L"&Show Data");
	AppendMenuW(hMenu_Sim, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Sim, MF_STRING, IDM_SIM_BENCHTIME, L"&Bench Time");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_Sim, L"&Simulation");

	//MESH
	HMENU hMenu_Mesh = CreateMenu();

	//---
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_SHOW, L"&Show Meshes\tAlt-M");
	AppendMenuW(hMenu_Mesh, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_LOADMASK, L"&Load Mask");
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_RESET, L"&Reset Mesh\tCtrl-R");
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_MASKALL, L"&Mask All");
	CheckMenuItem(hMenu_Mesh, IDM_MESH_MASKALL, MF_CHECKED);
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_SCALERECTS, L"&Scale All");
	AppendMenuW(hMenu_Mesh, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_ADDMATERIAL, L"&Add Material");
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_ADDFERROMAGNET, L"&Add Ferromagnet");
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_ADDANTIFERROMAGNET, L"&Add Antiferromagnet");
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_ADDDIPOLE, L"&Add Dipole");
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_ADDMETAL, L"&Add Metal");
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_ADDINSULATOR, L"&Add Insulator");
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_ATOMSIMPLECUBIC, L"&Add Atomistic Simple Cubic");
	AppendMenuW(hMenu_Mesh, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_EXCHANGECOUPLED, L"&Exchange Coupling");
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_COUPLETODIPOLES, L"&Couple to Dipoles");
	AppendMenuW(hMenu_Mesh, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Mesh, MF_STRING, IDM_MESH_COPY, L"&Copy Data");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_Mesh, L"&Mesh");

	//SHAPES
	HMENU hMenu_Shapes = CreateMenu();

	AppendMenuW(hMenu_Shapes, MF_STRING, IDM_SHAPES_MODIFIERS, L"&Modifiers");
	AppendMenuW(hMenu_Shapes, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Shapes, MF_STRING, IDM_SHAPES_DELRECT, L"&Clear All");
	AppendMenuW(hMenu_Shapes, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Shapes, MF_STRING, IDM_SHAPES_RECT, L"&Rectangle");
	AppendMenuW(hMenu_Shapes, MF_STRING, IDM_SHAPES_DISK, L"&Disk");
	AppendMenuW(hMenu_Shapes, MF_STRING, IDM_SHAPES_TRIANGLE, L"&Triangle");
	AppendMenuW(hMenu_Shapes, MF_STRING, IDM_SHAPES_TETRAHEDRON, L"&Tetrahedron");
	AppendMenuW(hMenu_Shapes, MF_STRING, IDM_SHAPES_PYRAMID, L"&Pyramid");
	AppendMenuW(hMenu_Shapes, MF_STRING, IDM_SHAPES_CONE, L"&Cone");
	AppendMenuW(hMenu_Shapes, MF_STRING, IDM_SHAPES_ELLIPSOID, L"&Ellipsoid");
	AppendMenuW(hMenu_Shapes, MF_STRING, IDM_SHAPES_TORUS, L"&Torus");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_Shapes, L"&Shapes");

	//CONFIGURATION
	HMENU hMenu_Cfg = CreateMenu();

	AppendMenuW(hMenu_Cfg, MF_STRING, IDM_CFG_MODULES, L"&Modules\tAlt-F");
	AppendMenuW(hMenu_Cfg, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Cfg, MF_STRING, IDM_CFG_ODE, L"&ODE\tAlt-E");
	AppendMenuW(hMenu_Cfg, MF_STRING, IDM_CFG_TIMESTEP, L"&Time Step");
	AppendMenuW(hMenu_Cfg, MF_STRING, IDM_CFG_ASTEP, L"&Adaptive Step");
	AppendMenuW(hMenu_Cfg, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Cfg, MF_STRING, IDM_CFG_MULTICONV, L"&Multilayered Convolution");
	AppendMenuW(hMenu_Cfg, MF_STRING, IDM_CFG_PBC, L"&PBC");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_Cfg, L"&Configuration");

	//PARAMETERS
	HMENU hMenu_Param = CreateMenu();

	AppendMenuW(hMenu_Param, MF_STRING, IDM_PARAM_SHOW, L"&Values\tAlt-P");
	AppendMenuW(hMenu_Param, MF_STRING, IDM_PARAM_TDEP, L"&Temperature Dependence\tAlt-T");
	AppendMenuW(hMenu_Param, MF_STRING, IDM_PARAM_SDEP, L"&Spatial Variation\tAlt-V");
	AppendMenuW(hMenu_Param, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Param, MF_STRING, IDM_PARAM_COPY, L"&Copy Data");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_Param, L"&Parameters");

	//MAGNETIZATION
	HMENU hMenu_Mag = CreateMenu();

	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_SET, L"&Uniform");
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_SETOBJECT, L"&Object");
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_DWALL, L"&Domain Wall");
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_VORTEX, L"&Vortex");
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_SKYRMION, L"&Skyrmion");
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_BLOCHSKYRMION, L"&Bloch Skyrmion");
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_RANDOM, L"&Random");
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_INVERT, L"&Invert");
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_MIRROR, L"&Mirror");
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_SETRECT, L"&Set Rectangle");
	AppendMenuW(hMenu_Mag, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_LOADOVF2MAG, L"&Load OVF2 File");
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_SAVEOVF2MAG, L"&Save OVF2 File");
	AppendMenuW(hMenu_Mag, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_2DGRAINS, L"&2D Grains");
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_3DGRAINS, L"&3D Grains");
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_ADDRECT, L"&Add Rectangle");
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_DELRECT, L"&Del Rectangle");
	AppendMenuW(hMenu_Mag, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Mag, MF_STRING, IDM_MAG_FIELD, L"&Field");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_Mag, L"&Magnetization");

	//ROUGHNESS
	HMENU hMenu_Rough = CreateMenu();

	AppendMenuW(hMenu_Rough, MF_STRING, IDM_ROUGH_SHOW, L"&Show Roughness");
	AppendMenuW(hMenu_Rough, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Rough, MF_STRING, IDM_ROUGH_EDIT, L"&Edit Refinement");
	AppendMenuW(hMenu_Rough, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Rough, MF_STRING, IDM_ROUGH_ROUGHEN, L"&Generate");
	AppendMenuW(hMenu_Rough, MF_STRING, IDM_ROUGH_SURFROUGHEN, L"&Surface Generate");
	AppendMenuW(hMenu_Rough, MF_STRING, IDM_ROUGH_LOAD, L"&Load Grayscale File");
	AppendMenuW(hMenu_Rough, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Rough, MF_STRING, IDM_ROUGH_CLEARROUGH, L"&Clear");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_Rough, L"&Roughness");

	//TEMPERATURE
	HMENU hMenu_Temp = CreateMenu();

	AppendMenuW(hMenu_Temp, MF_STRING, IDM_TEMPERATURE_SET, L"&Set");
	AppendMenuW(hMenu_Temp, MF_STRING, IDM_TEMPERATURE_CURIE, L"&Curie");
	AppendMenuW(hMenu_Temp, MF_STRING, IDM_TEMPERATURE_AMBIENT, L"&Ambient");
	AppendMenuW(hMenu_Temp, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Temp, MF_STRING, IDM_TEMPERATURE_LOADOVF2, L"&Load OVF2 Temperature");
	AppendMenuW(hMenu_Temp, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Temp, MF_STRING, IDM_TEMPERATURE_TIMESTEP, L"&Time Step");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_Temp, L"&Temperature");

	//TRANSPORT
	HMENU hMenu_Trans = CreateMenu();

	AppendMenuW(hMenu_Trans, MF_STRING, IDM_TRANS_DISABLE, L"&Disable Iteration");
	AppendMenuW(hMenu_Trans, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Trans, MF_STRING, IDM_TRANS_DEFAULTELECTRODES, L"&Set Default Electrodes");
	AppendMenuW(hMenu_Trans, MF_STRING, IDM_TRANS_ELECTRODES, L"&Show Electrodes");
	AppendMenuW(hMenu_Trans, MF_STRING, IDM_TRANS_ADDELECTRODE, L"&Add Electrode");
	AppendMenuW(hMenu_Trans, MF_STRING, IDM_TRANS_CLEARELECTRODES, L"&Clear Electrodes");
	AppendMenuW(hMenu_Trans, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Trans, MF_STRING, IDM_TRANS_SETV, L"&Set Potential");
	AppendMenuW(hMenu_Trans, MF_STRING, IDM_TRANS_SETI, L"&Set Current");
	AppendMenuW(hMenu_Trans, MF_STRING, IDM_TRANS_SETJC, L"&Set Current Density");
	AppendMenuW(hMenu_Trans, MF_STRING, IDM_TRANS_LOADOVF2, L"&Load OVF2 Current Density");
	AppendMenuW(hMenu_Trans, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Trans, MF_STRING, IDM_TRANS_CONFIG, L"&Configuration");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_Trans, L"&Transport");

	//MELASTIC
	HMENU hMenu_Melastic = CreateMenu();

	AppendMenuW(hMenu_Melastic, MF_STRING, IDM_MELASTIC_SETSTRESS, L"&Set Stress");
	AppendMenuW(hMenu_Melastic, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Melastic, MF_STRING, IDM_MELASTIC_LOADOVF2DISP, L"&Load OVF2 Displacement");
	AppendMenuW(hMenu_Melastic, MF_STRING, IDM_MELASTIC_LOADOVF2STRAIN, L"&Load OVF2 Strain");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_Melastic, L"&Mechanical");

	//VIEW
	HMENU hMenu_View = CreateMenu();
	
	AppendMenuW(hMenu_View, MF_STRING, IDM_VIEW_DISPLAY, L"&Display\tCtrl-Alt-D");
	AppendMenuW(hMenu_View, MF_STRING, IDM_VIEW_CENTER, L"&Center Mesh\tCtrl-Alt-C");
	AppendMenuW(hMenu_View, MF_STRING, IDM_VIEW_UPDATEFREQUENCY, L"&Update Frequency\tCtrl-Alt-I");
	AppendMenuW(hMenu_View, MF_STRING, IDM_VIEW_RENDER, L"&Render Options");
	AppendMenuW(hMenu_View, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_View, MF_STRING, IDM_VIEW_CLEARCONSOLE, L"&Clear Console");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_View, L"&View");

	//ALGORITHMS
	HMENU hMenu_Algo = CreateMenu();

	AppendMenuW(hMenu_Algo, MF_STRING, IDM_ALGO_MOVINGMESH, L"&Set Moving Mesh");
	AppendMenuW(hMenu_Algo, MF_STRING, IDM_ALGO_BLOCHMOVINGMESH, L"&Set Bloch Moving Mesh");
	AppendMenuW(hMenu_Algo, MF_STRING, IDM_ALGO_NEELMOVINGMESH, L"&Set Neel Moving Mesh");
	AppendMenuW(hMenu_Algo, MF_STRING, IDM_ALGO_SKYMOVINGMESH, L"&Set Skyrmion Moving Mesh");
	AppendMenuW(hMenu_Algo, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Algo, MF_STRING, IDM_ALGO_CLEARMOVINGMESH, L"&Undo Moving Mesh");
	AppendMenuW(hMenu_Algo, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_Algo, MF_STRING, IDM_ALGO_MOVINGMESHENABLED, L"&Moving Mesh Enabled");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_Algo, L"&Algorithms");

	//DATA PROCESSING
	HMENU hMenu_DP = CreateMenu();

	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_LOAD, L"&Load Data");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_SAVE, L"&Save Data");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_CLEARALL, L"&Clear All");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_SHOWSIZES, L"&Show Arrays");
	AppendMenuW(hMenu_DP, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_MUL, L"&Multiply");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_DOTPROD, L"&Dot Product");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_ADDDP, L"&Add Arrays");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_SUBDP, L"&Subtract Arrays");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_REMOVEOFFSET, L"&Remove Offset");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_SMOOTH, L"&Smooth");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_CARTESIANTOPOLAR, L"&Cartesian To Polar");
	AppendMenuW(hMenu_DP, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_GETPROFILE, L"&Get Profile");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_AVERAGEMESHRECT, L"&Average Mesh");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_DUMPTDEP, L"&Get Temperature Dependence");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_TOPOCHARGE, L"&Topological Charge");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_SKYRMIONCOUNT, L"&Count skyrmions");
	AppendMenuW(hMenu_DP, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_MINMAX, L"&Find Min-Max");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_MEAN, L"&Find Mean");
	AppendMenuW(hMenu_DP, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_LINREG, L"&Linear Regression");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_FITLORENTZ, L"&Fit Lorentzian");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_FITSKYRMION, L"&Fit Skyrmion");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_FITDW, L"&Fit Domain Wall");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_FITSTT, L"&Fit STT");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_FITSOT, L"&Fit SOT");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_FITSOTSTT, L"&Fit SOT and STT");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_FITADIA, L"&Fit Adiabatic");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_FITNONADIA, L"&Fit Nonadiabatic");
	AppendMenuW(hMenu_DP, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_COERCIVITY, L"&Coercivity");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_REMANENCE, L"&Remanence");
	AppendMenuW(hMenu_DP, MF_STRING, IDM_DP_COMPLETEHYSTER, L"&Complete Hysteresis");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_DP, L"&Data Processing");

	//ABOUT
	HMENU hMenu_About = CreateMenu();

	AppendMenuW(hMenu_About, MF_STRING, IDM_ABOUT_MANUAL, L"&Manual");
	AppendMenuW(hMenu_About, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_About, MF_STRING, IDM_ABOUT_MDB, L"&Materials Database");
	AppendMenuW(hMenu_About, MF_STRING, IDM_ABOUT_UPDATEMDB, L"&Update Database");
	AppendMenuW(hMenu_About, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_About, MF_STRING, IDM_ABOUT_CHECKUPDATES, L"&Check Updates");
	AppendMenuW(hMenu_About, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_About, MF_STRING, IDM_ABOUT_STARTUPOPTIONS, L"&Startup Options");
	AppendMenuW(hMenu_About, MF_SEPARATOR, 0, NULL);
	AppendMenuW(hMenu_About, MF_STRING, IDM_ABOUT_ABOUT, L"&About");

	AppendMenuW(hMenubar, MF_POPUP, (UINT_PTR)hMenu_About, L"&About");

	SetMenu(hwnd, hMenubar);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

int APIENTRY _tWinMain(_In_ HINSTANCE hInstance,
					   _In_opt_ HINSTANCE hPrevInstance,
					   _In_ LPTSTR    lpCmdLine,
					   _In_ int       nCmdShow)
{
	//////////////////////

	UNREFERENCED_PARAMETER(hPrevInstance);
	UNREFERENCED_PARAMETER(lpCmdLine);

 	// TODO: Place code here.
	MSG msg;

	// Initialize global strings
	LoadString(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
	LoadString(hInstance, IDC_BORISWINAPI, szWindowClass, MAX_LOADSTRING);
	
	//////////////////////
	//overwrite the title bar

	std::string sTitle = std::string("Boris v") + ToString((double)Program_Version / 100.0);

#if COMPILECUDA == 1
	sTitle += std::string(" CUDA ") + ToString((double)__CUDA_ARCH__ / 100);
#if SINGLEPRECISION == 1
	sTitle += std::string(" SP");
#else
	sTitle += std::string(" DP");
#endif
#else
	sTitle += std::string(" Lite.");
#endif

	copy(sTitle.begin(), sTitle.end(), szTitle);
	szTitle[sTitle.size()] = 0;

	//////////////////////
	//Arguments

	LPWSTR *szArgList;
	int argc;

	szArgList = CommandLineToArgvW(GetCommandLine(), &argc);

	for (int i = 1; i < argc; i++)
	{
		//First argument: server port
		if (i == 1) {

			server_port = WideStringtoString(szArgList[i]);
		}

		//Second argument: cuda device
		if (i == 2) {

			cudaDevice = ToNum(WideStringtoString(szArgList[i]));
		}

		//Third argument: window options (front/back)
		if (i == 3) {

			window_startup_option = WideStringtoString(szArgList[i]);
		}

		//Fourth argument: server password
		if (i == 4) {

			server_pwd = WideStringtoString(szArgList[i]);
		}
	}

	LocalFree(szArgList);

	//////////////////////
	//Startup

	MyRegisterClass(hInstance);

	// Perform application initialization: (use SW_MAXIMIZE to have the base window maximized)
	if (!InitInstance (hInstance, nCmdShow))
	{
		return FALSE;
	}

	HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_BORISWINAPI));
	
	//Windows message loop
	while(GetMessage(&msg, nullptr, 0, 0)) {

		if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg)) {

			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
	}
	
	return (int) msg.wParam;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//
//  FUNCTION: MyRegisterClass()
//
//  PURPOSE: Registers the window class.
//
ATOM MyRegisterClass(HINSTANCE hInstance)
{
	WNDCLASSEX wcex;

	wcex.cbSize = sizeof(WNDCLASSEX);

	wcex.style			= CS_HREDRAW | CS_VREDRAW;
	wcex.lpfnWndProc	= WndProc;
	wcex.cbClsExtra		= 0;
	wcex.cbWndExtra		= 0;
	wcex.hInstance		= hInstance;
	wcex.hIcon			= LoadIcon(hInstance, MAKEINTRESOURCE(IDI_BORISWINAPI));
	wcex.hCursor		= LoadCursor(nullptr, IDC_ARROW);
	wcex.hbrBackground	= (HBRUSH)(COLOR_WINDOW+1);
	wcex.lpszMenuName	= MAKEINTRESOURCE(IDC_BORISWINAPI);
	wcex.lpszClassName	= szWindowClass;
	wcex.hIconSm		= LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

	return RegisterClassEx(&wcex);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//
//   FUNCTION: InitInstance(HINSTANCE, int)
//
//   PURPOSE: Saves instance handle and creates main window
//
//   COMMENTS:
//
//        In this function, we save the instance handle in a global variable and
//        create and display the main program window.
//
BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
   HWND hWnd;

   hInst = hInstance; // Store instance handle in our global variable

   hWnd = CreateWindow(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW | WS_EX_ACCEPTFILES,
	   CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);

   if (!hWnd)
   {
      return FALSE;
   }

   ShowWindow(hWnd, SW_SHOWMAXIMIZED);

   if (window_startup_option == "back") {

	   SetWindowPos(hWnd, HWND_BOTTOM, 0, 0, 0, 0, SWP_NOACTIVATE | SWP_NOMOVE | SWP_NOSIZE);
   }

   UpdateWindow(hWnd);
   
   //Stop the "Program has stopped working" pop-up error box after closing program. I really cannot figure out why it appears sometimes!!!
   //!!! NOTE !!! If this is not enabled and the program closes with the "Program has stopped working" error, then the exe file is put on some windows watch list
   //(some kind of hook inserted in the running code to check function calls?), which will slow down execution dramatically. 
   //It seems this is done by the DiagTrack service (Diagnostics Tracking Service) - stop DPS service (Diagnostic Policy Service) then DiagTrack.
   SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);

   //Set normal arrow cursor - we are answering WM_SETCURSOR message and doing nothing there so need to set cursor separately
   SetCursor(LoadCursor(nullptr, IDC_ARROW));
   
   //Set to accept dropped files - file drop handled by WM_DROPFILES
   DragAcceptFiles(hWnd, true);

   //Allow WM_DROPFILES through UIPI filter for Windows 10 : if running in Administrator mode this message will be filtered since drag and drop from an unelevated source to an elevated recipient is blocked by default
   //This is what you have to do to force Windows 10 to allow drag and drop (marvellous isn't it!):
   ChangeWindowMessageFilterEx(hWnd, WM_DROPFILES, MSGFLT_ALLOW, nullptr);
   ChangeWindowMessageFilterEx(hWnd, WM_COPYDATA, MSGFLT_ALLOW, nullptr);
   unsigned int WM_COPYGLOBALDATA = 0x0049;
   ChangeWindowMessageFilterEx(hWnd, WM_COPYGLOBALDATA, MSGFLT_ALLOW, nullptr);
   
   //Instantiate Simulation object and start main simulation loop thread (non-blocking)
   
   //See https://stackoverflow.com/questions/64988525/bug-related-to-g-openmp-when-using-stdthread/65001887#65001887
   std::thread simulation_instantiation_thread(&make_Simulation, hWnd);
   simulation_instantiation_thread.join();

   //Add menus - this needs to be after creating Simulation as we need to get status flags to disable some menu items (e.g. CUDA available)
   AddMenus(hWnd);

   return TRUE;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//------------------------- SET DIRECTORY DIALOG

static int CALLBACK BrowseCallbackProc(HWND hwnd, UINT uMsg, LPARAM lParam, LPARAM lpData)
{
	if (uMsg == BFFM_INITIALIZED) {

		WCHAR* path = reinterpret_cast<WCHAR*>(lpData);
		
		::SendMessage(hwnd, BFFM_SETSELECTION, true, (LPARAM)path);
	}

	return 0;
}

std::string BrowseFolder(void)
{
	WCHAR path[MAX_PATH];

	BROWSEINFO bi = { 0 };
	bi.lpszTitle = L"Set default directory";
	bi.ulFlags = BIF_RETURNONLYFSDIRS | BIF_NEWDIALOGSTYLE;
	bi.lpfn = BrowseCallbackProc;

	LPITEMIDLIST pidl = SHBrowseForFolder(&bi);

	if (pidl != 0)
	{
		//get the name of the folder and put it in path
		SHGetPathFromIDList(pidl, path);

		//free memory used
		IMalloc * imalloc = 0;
		if (SUCCEEDED(SHGetMalloc(&imalloc)))
		{
			imalloc->Free(pidl);
			imalloc->Release();
		}

		return WideStringtoString(std::wstring(path));
	}

	return "";
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//------------------------- LOAD FILE DIALOG

std::string OpenFileDialog(HWND hWnd, std::string initial_dir = "")
{
	OPENFILENAME ofn;
	WCHAR szFile[260] = { 0 };

	// Initialize OPENFILENAME
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.hwndOwner = hWnd;
	ofn.lpstrFile = szFile;
	ofn.nMaxFile = sizeof(szFile);
	ofn.lpstrFilter = L"All\0*.*\0Text\0*.TXT\0";
	ofn.nFilterIndex = 1;
	ofn.lpstrFileTitle = nullptr;
	ofn.nMaxFileTitle = 0;
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

	if (initial_dir.length()) ofn.lpstrInitialDir = StringtoWCHARPointer(initial_dir);
	else ofn.lpstrInitialDir = nullptr;

	if (GetOpenFileName(&ofn) == TRUE) return WideStringtoString(std::wstring(szFile));
	else return "";
}

//------------------------- SAVE FILE DIALOG

std::string SaveFileDialog(HWND hWnd, std::string initial_dir = "")
{
	OPENFILENAME ofn;
	WCHAR szFile[260] = { 0 };

	// Initialize OPENFILENAME
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.hwndOwner = hWnd;
	ofn.lpstrFile = szFile;
	ofn.nMaxFile = sizeof(szFile);
	ofn.lpstrFilter = L"All\0*.*\0Text\0*.TXT\0";
	ofn.nFilterIndex = 1;
	ofn.lpstrFileTitle = nullptr;
	ofn.nMaxFileTitle = 0;
	ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

	if (initial_dir.length()) ofn.lpstrInitialDir = StringtoWCHARPointer(initial_dir);
	else ofn.lpstrInitialDir = nullptr;

	if (GetSaveFileName(&ofn) == TRUE) return WideStringtoString(std::wstring(szFile));
	else return "";
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//------------------------- WNDPROC

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	PAINTSTRUCT ps;
	HDC hdc = 0;
	
	switch (message) {
	
	case WM_CREATE:
		break;

	case WM_PAINT:
		hdc = BeginPaint(hWnd, &ps);
		if(pSim) pSim->RefreshScreen();
		EndPaint(hWnd, &ps);
		break;

	case WM_CLOSE:
		{
			//clean up and exit
			if(pSim) delete pSim;
			pSim = nullptr;

			//Destroy window : WM_DESTROY will be issued
			DestroyWindow(hWnd);
		}
		break;

	case WM_DESTROY:
	{
		//Quit program : WM_QUIT will be issued
		PostQuitMessage(0);
	}
		break;

	case WM_MENUSELECT:
	{
		//Make sure all menu items wich can be checked/unchecked reflect the simulation state correctly
		if (pSim) {

			CheckMenuItem(hMenubar, IDM_SIM_SERVER, (pSim->GetScriptServerStatus() ? MF_CHECKED : MF_UNCHECKED));
			
			CheckMenuItem(hMenubar, IDM_MESH_MASKALL, (pSim->GetIndividualShapeStatus() ? MF_UNCHECKED : MF_CHECKED));
			
			CheckMenuItem(hMenubar, IDM_MESH_SCALERECTS, (pSim->GetScaleRectsStatus() ? MF_CHECKED : MF_UNCHECKED));
			
			CheckMenuItem(hMenubar, IDM_MESH_COUPLETODIPOLES, (pSim->GetCoupleToDipolesStatus() ? MF_CHECKED : MF_UNCHECKED));
			
			if (pSim->GetCUDAStatus() != -1) CheckMenuItem(hMenubar, IDM_SIM_CUDA, (pSim->GetCUDAStatus() ? MF_CHECKED : MF_UNCHECKED));
			
			CheckMenuItem(hMenubar, IDM_ALGO_MOVINGMESHENABLED, (pSim->GetMovingMeshStatus() ? MF_CHECKED : MF_UNCHECKED));

			CheckMenuItem(hMenubar, IDM_TRANS_DISABLE, (pSim->GetDisabledTransportStatus() ? MF_CHECKED : MF_UNCHECKED));
		}
	}
		break;

	case WM_COMMAND:
		switch (LOWORD(wParam)) {
			
			//FILE
		case IDM_FILE_DEFAULT_HK:
		case IDM_FILE_DEFAULT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DEFAULT));
			break;

		case IDM_FILE_LOAD_HK:
		case IDM_FILE_LOAD:
		{
			std::string fileName = OpenFileDialog(hWnd);
			if (pSim && fileName.length()) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_LOADSIM) + " " + fileName);
		}
			break;

		case IDM_FILE_LOADEXAMPLE:
		{
			std::string fileName = OpenFileDialog(hWnd, GetUserDocumentsPath() + std::string("Boris Data/Examples"));
			if (pSim && fileName.length()) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_LOADSIM) + " " + fileName);
		}
		break;
			
		case IDM_FILE_SAVE_HK:
		case IDM_FILE_SAVE:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SAVESIM));
			break;

		case IDM_FILE_SAVEAS_HK:
		case IDM_FILE_SAVEAS:
		{
			std::string fileName = SaveFileDialog(hWnd);
			if (pSim && fileName.length()) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SAVESIM) + " " + fileName);
			
		}
			break;

		case IDM_FILE_SAVEOVF2:
		{
			std::string fileName = SaveFileDialog(hWnd);
			if (pSim && fileName.length()) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SAVEOVF2) + " " + fileName);
		}
			break;

		case IDM_FILE_SAVEIMAGE_HK:
		case IDM_FILE_SAVEIMAGE:
		{
			std::string fileName = SaveFileDialog(hWnd);
			if (pSim && fileName.length()) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SAVEMESHIMAGE) + " " + fileName);
		}
			break;

		case IDM_FILE_MAKEVIDEO:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_MAKEVIDEO));
			break;

		case IDM_FILE_IMAGECROPPING:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_IMAGECROPPING));
			break;

		case IDM_FILE_CHDIR_HK:
		case IDM_FILE_CHDIR:
		{
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_CHDIR));
			std::string directory = BrowseFolder();
			if (pSim && directory.length()) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_CHDIR) + " " + directory);
		}
			break;

		case IDM_FILE_QUIT_HK:
		case IDM_FILE_QUIT:

			SendMessage(hWnd, WM_CLOSE, 0, 0);
			break;

			//MESH
		case IDM_MESH_SHOW_HK:
		case IDM_MESH_SHOW:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_MESH));
			break;

		case IDM_MESH_LOADMASK:
		{
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_LOADMASKFILE));
			std::string fileName = OpenFileDialog(hWnd);
			if (pSim && fileName.length()) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_LOADMASKFILE) + " " + fileName);
		}
			break;

		case IDM_MESH_MASKALL:
		{
			UINT state = GetMenuState(hMenubar, IDM_MESH_MASKALL, MF_BYCOMMAND);

			if (state == MF_CHECKED) {

				CheckMenuItem(hMenubar, IDM_MESH_MASKALL, MF_UNCHECKED);

				if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_INDIVIDUALMASKSHAPE) + " 1");
			}
			else {

				CheckMenuItem(hMenubar, IDM_MESH_MASKALL, MF_CHECKED);

				if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_INDIVIDUALMASKSHAPE) + " 0");
			}
		}
			break;

		case IDM_MESH_SCALERECTS:
		{
			UINT state = GetMenuState(hMenubar, IDM_MESH_SCALERECTS, MF_BYCOMMAND);

			if (state == MF_CHECKED) {

				CheckMenuItem(hMenubar, IDM_MESH_SCALERECTS, MF_UNCHECKED);

				if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SCALEMESHRECTS) + " 0");
			}
			else {

				CheckMenuItem(hMenubar, IDM_MESH_SCALERECTS, MF_CHECKED);

				if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SCALEMESHRECTS) + " 1");
			}
		}
			break;

		case IDM_MESH_RESET_HK:
		case IDM_MESH_RESET:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_RESETMESH));
			break;

		case IDM_MESH_ADDMATERIAL:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ADDMATERIAL));
			break;

		case IDM_MESH_ADDFERROMAGNET:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ADDFMESH));
			break;

		case IDM_MESH_ADDANTIFERROMAGNET:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ADDAFMESH));
			break;

		case IDM_MESH_ADDMETAL:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ADDMETALMESH));
			break;

		case IDM_MESH_ADDINSULATOR:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ADDINSULATORMESH));
			break;

		case IDM_MESH_ADDDIPOLE:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ADDDIPOLEMESH));
			break;

		case IDM_MESH_ATOMSIMPLECUBIC:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ADDAMESHCUBIC));
			break;

		case IDM_MESH_EXCHANGECOUPLED:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_EXCHANGECOUPLEDMESHES));
			break;

		case IDM_MESH_COUPLETODIPOLES:
		{
			UINT state = GetMenuState(hMenubar, IDM_MESH_COUPLETODIPOLES, MF_BYCOMMAND);

			if (state == MF_CHECKED) {

				CheckMenuItem(hMenubar, IDM_MESH_COUPLETODIPOLES, MF_UNCHECKED);

				if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_COUPLETODIPOLES) + " 0");
			}
			else {

				CheckMenuItem(hMenubar, IDM_MESH_COUPLETODIPOLES, MF_CHECKED);

				if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_COUPLETODIPOLES) + " 1");
			}
		}
			break;

		case IDM_MESH_COPY:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_COPYMESHDATA));
			break;

			//SHAPES
		case IDM_SHAPES_MODIFIERS:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SHAPEMOD_ROT));
			break;

		case IDM_SHAPES_DELRECT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DELRECT));
			break;

		case IDM_SHAPES_RECT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SHAPE_RECT));
			break;

		case IDM_SHAPES_DISK:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SHAPE_DISK));
			break;

		case IDM_SHAPES_TRIANGLE:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SHAPE_TRIANGLE));
			break;

		case IDM_SHAPES_TETRAHEDRON:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SHAPE_TETRAHEDRON));
			break;

		case IDM_SHAPES_PYRAMID:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SHAPE_PYRAMID));
			break;

		case IDM_SHAPES_CONE:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SHAPE_CONE));
			break;

		case IDM_SHAPES_ELLIPSOID:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SHAPE_ELLIPSOID));
			break;

		case IDM_SHAPES_TORUS:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SHAPE_TORUS));
			break;

			//SIMULATION
		case IDM_SIM_RUN_HK:
		case IDM_SIM_RUN:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_RUN));
			break;

		case IDM_SIM_COMPUTEFIELDS:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_COMPUTEFIELDS));
			break;

		case IDM_SIM_STOP_HK:
		case IDM_SIM_STOP:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_STOP));
			break;

		case IDM_SIM_RESET_HK:
		case IDM_SIM_RESET:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_RESET));
			break;

		case IDM_SIM_CUDA:
			if (pSim) {

				if (pSim->GetCUDAStatus() != -1) {

					UINT state = GetMenuState(hMenubar, IDM_SIM_CUDA, MF_BYCOMMAND);

					if (state == MF_CHECKED) {

						CheckMenuItem(hMenubar, IDM_SIM_CUDA, MF_UNCHECKED);

						if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_CUDA) + " 0");
					}
					else {

						CheckMenuItem(hMenubar, IDM_SIM_CUDA, MF_CHECKED);

						if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_CUDA) + " 1");
					}
				}
			}
			break;

		case IDM_SIM_SERVER:
		{
			UINT state = GetMenuState(hMenubar, IDM_SIM_SERVER, MF_BYCOMMAND);

			if (state == MF_CHECKED) {

				CheckMenuItem(hMenubar, IDM_SIM_SERVER, MF_UNCHECKED);

				if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SCRIPTSERVER) + " 0");
			}
			else {

				CheckMenuItem(hMenubar, IDM_SIM_SERVER, MF_CHECKED);

				if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SCRIPTSERVER) + " 1");
			}
		}
			break;

		case IDM_SIM_SCHEDULE_HK:
		case IDM_SIM_SCHEDULE:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_STAGES));
			break;

		case IDM_SIM_DATA_HK:
		case IDM_SIM_DATA:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DATA));
			break;
			
		case IDM_SIM_SHOWDATA:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SHOWDATA));
			break;

		case IDM_SIM_BENCHTIME:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_BENCHTIME));
			break;

		//CONFIGURATION
		case IDM_CFG_MODULES_HK:
		case IDM_CFG_MODULES:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_MODULES));
			break;

		case IDM_CFG_ODE_HK:
		case IDM_CFG_ODE:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ODE));
			break;
			
		case IDM_CFG_TIMESTEP:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SETDT));
			break;

		case IDM_CFG_ASTEP:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ASTEPCTRL));
			break;

		case IDM_CFG_MULTICONV:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_MULTICONV));
			break;

		case IDM_CFG_PBC:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_PBC));
			break;

		//PARAMETERS
		case IDM_PARAM_SHOW_HK:
		case IDM_PARAM_SHOW:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_PARAMS));
			break;

		case IDM_PARAM_TDEP_HK:
		case IDM_PARAM_TDEP:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_PARAMSTEMP));
			break;

		case IDM_PARAM_SDEP_HK:
		case IDM_PARAM_SDEP:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_PARAMSVAR));
			break;

		case IDM_PARAM_COPY:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_COPYPARAMS));
			break;

		//MAGNETIZATION
		case IDM_MAG_SET:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_NOPARAMS_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SETANGLE) + " 90 0");
			break;

		case IDM_MAG_SETOBJECT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SETOBJECTANGLE));
			break;

		case IDM_MAG_FIELD:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_NOPARAMS_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SETFIELD) + " 1e3 90 0");
			break;

		case IDM_MAG_DWALL:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_NOPARAMS_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DWALL) + " x y 80nm 0");
			break;

		case IDM_MAG_VORTEX:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_NOPARAMS_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_VORTEX) + " 1 1 1");
			break;

		case IDM_MAG_SKYRMION:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_NOPARAMS_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SKYRMION) + " -1 -1 40nm 40nm 40nm");
			break;

		case IDM_MAG_BLOCHSKYRMION:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_NOPARAMS_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SKYRMIONBLOCH) + " -1 -1 40nm 40nm 40nm");
			break;

		case IDM_MAG_2DGRAINS:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_NOPARAMS_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_GENERATE2DGRAINS) + " 40nm");
			break;

		case IDM_MAG_3DGRAINS:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_NOPARAMS_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_GENERATE3DGRAINS) + " 40nm");
			break;

		case IDM_MAG_RANDOM:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_RANDOM));
			break;

		case IDM_MAG_LOADOVF2MAG:
		{
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_LOADOVF2MAG));
			std::string fileName = OpenFileDialog(hWnd);
			if (pSim && fileName.length()) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_LOADOVF2MAG) + " " + fileName);
		}
			break;

		case IDM_MAG_SAVEOVF2MAG:
		{
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SAVEOVF2MAG));
			std::string fileName = SaveFileDialog(hWnd);
			if (pSim && fileName.length()) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SAVEOVF2MAG) + " " + fileName);
		}
			break;

		case IDM_MAG_ADDRECT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ADDRECT));
			break;

		case IDM_MAG_DELRECT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DELRECT));
			break;

		case IDM_MAG_SETRECT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SETRECT));
			break;

		case IDM_MAG_INVERT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_INVERTMAG));
			break;

		case IDM_MAG_MIRROR:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_MIRRORMAG));
			break;

		case IDM_ROUGH_SHOW:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DISPLAY) + " Roughness");
			break;

		case IDM_ROUGH_EDIT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_REFINEROUGHNESS));
			break;

		case IDM_ROUGH_ROUGHEN:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_NOPARAMS_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ROUGHENMESH) + " 4nm");
			break;

		case IDM_ROUGH_SURFROUGHEN:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_NOPARAMS_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SURFROUGHENJAGGED) + " 4nm 40nm");
			break;

		case IDM_ROUGH_LOAD:
		{
			std::string fileName = OpenFileDialog(hWnd);
			if (pSim && fileName.length()) pSim->NewMessage(AC_CONSOLECOMMAND_NOPARAMS_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_LOADMASKFILE) + " 4nm " + fileName);
		}
			break;

		case IDM_ROUGH_CLEARROUGH:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_CLEARROUGHNESS));
			break;

		//TEMPERATURE
		case IDM_TEMPERATURE_SET:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_TEMPERATURE));
			break;

		case IDM_TEMPERATURE_CURIE:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_CURIETEMPERATURE));
			break;

		case IDM_TEMPERATURE_AMBIENT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_AMBIENTTEMPERATURE));
			break;

		case IDM_TEMPERATURE_LOADOVF2:
		{
			std::string fileName = OpenFileDialog(hWnd);
			if (pSim && fileName.length()) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_LOADOVF2TEMP) + " " + fileName);
		}
			break;

		case IDM_TEMPERATURE_TIMESTEP:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SETHEATDT));
			break;

		//TRANSPORT
		case IDM_TRANS_DISABLE:
		{
			UINT state = GetMenuState(hMenubar, IDM_TRANS_DISABLE, MF_BYCOMMAND);

			if (state == MF_CHECKED) {

				if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DISABLETRANSPORTSOLVER) + " 0");
				if (!pSim->GetDisabledTransportStatus()) CheckMenuItem(hMenubar, IDM_TRANS_DISABLE, MF_UNCHECKED);
			}
			else {

				if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DISABLETRANSPORTSOLVER) + " 1");
				if (pSim->GetDisabledTransportStatus()) CheckMenuItem(hMenubar, IDM_TRANS_DISABLE, MF_CHECKED);
			}
		}
			break;

		case IDM_TRANS_DEFAULTELECTRODES:
			if (pSim) {

				pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SETDEFAULTELECTRODES));
				pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ELECTRODES));
			}
			break;

		case IDM_TRANS_ELECTRODES:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ELECTRODES));
			break;

		case IDM_TRANS_ADDELECTRODE:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ADDELECTRODE));
			break;

		case IDM_TRANS_CLEARELECTRODES:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_CLEARELECTRODES));
			break;

		case IDM_TRANS_SETV:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SETPOTENTIAL));
			break;

		case IDM_TRANS_SETI:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SETCURRENT));
			break;

		case IDM_TRANS_SETJC:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SETCURRENTDENSITY));
			break;

		case IDM_TRANS_LOADOVF2:
		{
			std::string fileName = OpenFileDialog(hWnd);
			if (pSim && fileName.length()) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_LOADOVF2CURR) + " " + fileName);
		}
			break;

		case IDM_TRANS_CONFIG:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_TSOLVERCONFIG));
			break;

		//MELASTIC
		case IDM_MELASTIC_SETSTRESS:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_SETSTRESS));
			break;

		case IDM_MELASTIC_LOADOVF2DISP:
		{
			std::string fileName = OpenFileDialog(hWnd);
			if (pSim && fileName.length()) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_LOADOVF2DISP) + " " + fileName);
		}
			break;

		case IDM_MELASTIC_LOADOVF2STRAIN:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_LOADOVF2STRAIN));
			break;

		//VIEW
		case IDM_VIEW_DISPLAY_HK:
		case IDM_VIEW_DISPLAY:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DISPLAY));
			break;

		case IDM_VIEW_CENTER_HK:
		case IDM_VIEW_CENTER:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_CENTER));
			break;

		case IDM_VIEW_UPDATEFREQUENCY_HK:
		case IDM_VIEW_UPDATEFREQUENCY:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_ITERUPDATE));
			break;

		case IDM_VIEW_RENDER:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DISPLAYDETAILLEVEL));
			break;

		case IDM_VIEW_CLEARCONSOLE:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_CLEARSCREEN));
			break;

		//ALGORITHMS
		case IDM_ALGO_MOVINGMESH:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_PREPAREMOVINGMESH));
			break;

		case IDM_ALGO_BLOCHMOVINGMESH:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_PREPAREMOVINGBLOCHMESH));
			break;

		case IDM_ALGO_NEELMOVINGMESH:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_PREPAREMOVINGNEELMESH));
			break;

		case IDM_ALGO_SKYMOVINGMESH:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_PREPAREMOVINGSKYRMIONMESH));
			break;

		case IDM_ALGO_CLEARMOVINGMESH:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_CLEARMOVINGMESH));
			break;

		case IDM_ALGO_MOVINGMESHENABLED:
		{ 
			UINT state = GetMenuState(hMenubar, IDM_ALGO_MOVINGMESHENABLED, MF_BYCOMMAND);

			if (state == MF_CHECKED) {

				CheckMenuItem(hMenubar, IDM_ALGO_MOVINGMESHENABLED, MF_UNCHECKED);

				if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_MOVINGMESH) + " 0");
			}
			else {

				CheckMenuItem(hMenubar, IDM_ALGO_MOVINGMESHENABLED, MF_CHECKED);

				if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_MOVINGMESH) + " 1");
			}
		}
			break;

		//DATA PROCESSING
		case IDM_DP_CLEARALL:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_CLEARALL));
			break;

		case IDM_DP_SHOWSIZES:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_SHOWSIZES));
			break;

		case IDM_DP_LOAD:
		{
			std::string fileName = OpenFileDialog(hWnd);
			if (pSim && fileName.length()) pSim->NewMessage(AC_CONSOLECOMMAND_NOPARAMS_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_LOAD) + " " + fileName + " ");
		}
			break;

		case IDM_DP_SAVE:
		{
			std::string fileName = SaveFileDialog(hWnd);
			if (pSim && fileName.length()) pSim->NewMessage(AC_CONSOLECOMMAND_NOPARAMS_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_SAVE) + " " + fileName + " ");
		}
			break;

		case IDM_DP_GETPROFILE:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_GETEXACTPROFILE));
			break;

		case IDM_AVERAGEMESHRECT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_AVERAGEMESHRECT));
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_AVERAGEMESHRECT));
			break;

		case IDM_DP_TOPOCHARGE:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_TOPOCHARGE));
			break;
		
		case IDM_DP_SKYRMIONCOUNT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_COUNTSKYRMIONS));
			break;
		
		case IDM_DP_MUL:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_MUL));
			break;

		case IDM_DP_DOTPROD:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_DOTPROD));
			break;

		case IDM_DP_ADDDP:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_ADDDP));
			break;

		case IDM_DP_SUBDP:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_SUBDP));
			break;

		case IDM_DP_MINMAX:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_MINMAX));
			break;

		case IDM_DP_MEAN:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_MEAN));
			break;

		case IDM_DP_LINREG:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_LINREG));
			break;

		case IDM_DP_COERCIVITY:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_COERCIVITY));
			break;

		case IDM_DP_REMANENCE:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_REMANENCE));
			break;

		case IDM_DP_COMPLETEHYSTER:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_COMPLETEHYSTERLOOP));
			break;

		case IDM_DP_DUMPTDEP:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_DUMPTDEP));
			break;

		case IDM_DP_FITLORENTZ:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_FITLORENTZ));
			break;

		case IDM_DP_FITSKYRMION:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_FITSKYRMION));
			break;

		case IDM_DP_FITDW:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_FITDW));
			break;

		case IDM_DP_FITSTT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_FITSTT));
			break;

		case IDM_DP_FITSOT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_FITSOT));
			break;

		case IDM_DP_FITSOTSTT:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_FITSOTSTT));
			break;

		case IDM_DP_FITADIA:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_FITADIABATIC));
			break;

		case IDM_DP_FITNONADIA:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_FITNONADIABATIC));
			break;

		case IDM_DP_REMOVEOFFSET:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_REMOVEOFFSET));
			break;

		case IDM_DP_CARTESIANTOPOLAR:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_CARTESIANTOPOLAR));
			break;

		case IDM_DP_SMOOTH:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND_ENTRY, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_DP_SMOOTH));
			break;

		//ABOUT
		case IDM_ABOUT_MANUAL:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_OPENMANUAL));
			break;

		case IDM_ABOUT_MDB:
			open_web_page("https://www.boris-spintronics.uk/online-materials-database/");
			break;

		case IDM_ABOUT_UPDATEMDB:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_UPDATEMDB));
			break;

		case IDM_ABOUT_CHECKUPDATES:
			if (pSim) pSim->NewMessage(AC_CONSOLECOMMAND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(CMD_CHECKUPDATES));
			break;

		case IDM_ABOUT_STARTUPOPTIONS:
			if (pSim) pSim->Show_StartupOptions();
			break;

		case IDM_ABOUT_ABOUT:
			open_web_page("https://www.boris-spintronics.uk");
			break;
		}
		break;
		
	//keyboard key down - treat some special character codes here
	case WM_KEYDOWN: 
		{
			
			//intercept CTRL+ something combinations - NOTE (!!!), if that something is an ascii character, it will also raise the WM_CHAR case, so you must check there that CTRL is not on
			if(GetKeyState(VK_CONTROL) & 0x8000) {
				
				switch(wParam) {

					case 0x56: //Ctrl+V
						{
							if(pSim) pSim->NewMessage(AC_CTRL_V, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), GetClipboardText());
						}
						break;

					case 0x43: //Ctrl+C
						{
							if(pSim) pSim->NewMessage(AC_CTRL_C, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
						}
						break;

					default:
						break;
				}
			}
			else {
				
				switch (wParam) { 
			
				case VK_LEFT:   // LEFT ARROW
					if(pSim) {

						if(!(GetKeyState(VK_SHIFT) & 0x8000)) pSim->NewMessage(AC_KBDLEFT, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
						else pSim->NewMessage(AC_KBDSHIFTLEFT, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					}
					break;
			
				case VK_RIGHT:  // RIGHT ARROW
					if(pSim) {

						if(!(GetKeyState(VK_SHIFT) & 0x8000)) pSim->NewMessage(AC_KBDRIGHT, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
						else pSim->NewMessage(AC_KBDSHIFTRIGHT, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					}
					break; 

				case VK_UP:     // UP ARROW
					if(pSim) pSim->NewMessage(AC_KBDUP, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;
			
				case VK_DOWN:   // DOWN ARROW 
					if(pSim) pSim->NewMessage(AC_KBDDN, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;
			
				case VK_HOME:   // HOME
					if(pSim) {

						if(!(GetKeyState(VK_SHIFT) & 0x8000)) pSim->NewMessage(AC_KBDHOME, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
						else pSim->NewMessage(AC_KBDSHIFTHOME, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					}
					break;
			
				case VK_END:    // END
					if(pSim) {

						if(!(GetKeyState(VK_SHIFT) & 0x8000)) pSim->NewMessage(AC_KBDEND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
						else pSim->NewMessage(AC_KBDSHIFTEND, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					}
					break;
			
				case VK_INSERT: // INS
					if(pSim) pSim->NewMessage(AC_KBDINS, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;
			
				case VK_DELETE: // DEL
					if(pSim) pSim->NewMessage(AC_KBDDEL, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;

				case VK_PRIOR:	//PAGE UP
					if(pSim) pSim->NewMessage(AC_KBDPGUP, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;

				case VK_NEXT:	//PAGE DOWN
					if(pSim) pSim->NewMessage(AC_KBDPGDN, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;

				case VK_ESCAPE:	//ESC
					if(pSim) pSim->NewMessage(AC_KBDESC, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;

				case VK_BACK:	//BACKSPACE
					if(pSim) pSim->NewMessage(AC_KBDBACK, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;

				case VK_RETURN: //ENTER KEY
					if(pSim) pSim->NewMessage(AC_KBDENTER, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
					break;

				default:
					break;
				}
			}
		}
		break;
			
	//treat ascii character codes
	case WM_CHAR:
		{

			//if an ascii character key is pressed, first make sure we don't have CTRL key on as CTRL+something combinations are handled elsewhere. If you don't check here, the something key will also be handled here - not good!
			if(!(GetKeyState(VK_CONTROL) & 0x8000)) {

				switch (wParam) {
			
				case 0x08: // Process a backspace. 
					break;
			
				case 0x0A: // Process a linefeed. 
					break;
			
				case 0x1B: // Process an escape. 
					break; 
			
				case 0x09: // Process a tab. 
					break; 
			
				case 0x0D: // Process a carriage return.
					break;
			
				default: // Process displayable characters. 			
					if(pSim) pSim->NewMessage(AC_KEYBOARD, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString((char)wParam));
					break; 
				}
			}
		}
		break;
			
	case WM_LBUTTONDBLCLK: //mouse left button double click
		//NOTE: Need CS_DBLCLKS style set when creating Instance, for it to issue WM_LBUTTONDBLCLK.
		if(pSim) pSim->NewMessage(AC_DOUBLECLICK, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		break;

	case WM_LBUTTONDOWN: //left mouse button down
		if (pSim) {
			
			if ((GetKeyState(VK_SHIFT) & 0x8000) == 0x8000) pSim->NewMessage(AC_SHIFT_MOUSELEFTDOWN, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
			else pSim->NewMessage(AC_MOUSELEFTDOWN, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		}
		break;

	case WM_LBUTTONUP:  //left mouse button up
		if(pSim) pSim->NewMessage(AC_MOUSELEFTUP, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		break;

	case WM_MBUTTONDOWN:	//middle mouse button down
		if(pSim) pSim->NewMessage(AC_MOUSEMIDDLEDOWN, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		break;

	case WM_MBUTTONUP:	//middle mouse button up
		if(pSim) pSim->NewMessage(AC_MOUSEMIDDLEUP, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		break;
		
	case WM_RBUTTONDOWN: //right mouse button down
		if(pSim) pSim->NewMessage(AC_MOUSERIGHTDOWN, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		break;

	case WM_RBUTTONUP: //right mouse button up
		if(pSim) pSim->NewMessage(AC_MOUSERIGHTUP, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		break;

	case WM_NCMOUSELEAVE:
	case WM_MOUSELEAVE:
		//User area left : signal this to display and set cursor back to arrow default
		if (pSim) pSim->NewMessage(AC_ALLWINDOWSLEFT, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
		break;

	case WM_MOUSEMOVE:
	{
		// track the mouse so we know when it leaves our window
		TRACKMOUSEEVENT tme;
		tme.cbSize = sizeof(tme);
		tme.hwndTrack = hWnd;
		tme.dwFlags = TME_LEAVE;
		tme.dwHoverTime = 1;
		TrackMouseEvent(&tme);

		if (pSim) pSim->NewMessage(AC_MOUSEMOVE, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)));
	}
		break;

	case WM_MOUSEWHEEL: //mouse wheel used
		if(pSim) pSim->NewMessage(AC_MOUSEWHEEL, INT2(GET_X_LPARAM(lParam), GET_Y_LPARAM(lParam)), ToString(GET_WHEEL_DELTA_WPARAM(wParam)/120));
		break;
		
	case WM_SETCURSOR:
		//answering this message to stop cursor from reverting to the default IDC_ARROW setting
		break;
		
	case WM_DROPFILES:
		{
			WCHAR lpszFile[MAX_LOADSTRING] = { 0 };
			UINT uFile = 0;
			HDROP hDrop = (HDROP)wParam;

			uFile = DragQueryFile(hDrop, 0xFFFFFFFF, NULL, NULL);
			if (uFile != 1) {
			
				//need only one dropped file
				DragFinish(hDrop);
				break;
			}

			lpszFile[0] = '\0';
		
			if (DragQueryFile(hDrop, 0, lpszFile, MAX_LOADSTRING)) {
			
				if (pSim) {

					POINT dropPoint;
					DragQueryPoint(hDrop, &dropPoint);

					std::string fileName = WideStringtoString(std::wstring(lpszFile));

					pSim->NewMessage(AC_DROPFILES, INT2((int)dropPoint.x, (int)dropPoint.y), fileName);
				}
			}

			DragFinish(hDrop);
		}
		break;

	default:
		return DefWindowProc(hWnd, message, wParam, lParam);
		break;
	}
	
	return 0;
}

#endif

