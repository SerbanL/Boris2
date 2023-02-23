#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "ErrorHandler.h"
#include "MeshCUDA.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC

#include "BorisCUDALib.h"

class FMesh;
class ManagedDiffEqFMCUDA;

//Store Mesh quantities as cu_obj managed cuda VECs
class FMeshCUDA :
	public MeshCUDA
{

private:

	//pointer to cpu version of this mesh
	FMesh *pFMesh;
	
public:

	//make this object by copying data from the Mesh holding this object
	FMeshCUDA(FMesh* pMesh);

	~FMeshCUDA();

	//----------------------------------- IMPORTANT CONTROL METHODS
	
	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}
	
	//Take a Monte Carlo Metropolis step in this mesh
	cuBReal Iterate_MonteCarloCUDA_Classic(cuBReal mc_cone_angledeg, double target_acceptance_rate);

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//Check if mesh needs to be moved (using the MoveMesh method) - return amount of movement required (i.e. parameter to use when calling MoveMesh).
	cuBReal CheckMoveMesh(bool antisymmetric, double threshold);

	//----------------------------------- ENABLED MESH PROPERTIES CHECKERS

	//get exchange_couple_to_meshes status flag from the cpu version
	bool GetMeshExchangeCoupling(void);

	//----------------------------------- ODE METHODS IN (ANTI)FERROMAGNETIC MESH : Mesh_FerromagneticCUDA.cu, Mesh_FerromagneticCUDA.cpp

	//return average dm/dt in the given avRect (relative rect). Here m is the direction vector.
	cuReal3 Average_dmdt(cuBox avBox);

	//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
	cuReal3 Average_mxdmdt(cuBox avBox);

	//Save current magnetization in sM VECs (e.g. useful to reset dM / dt calculation)
	void SaveMagnetization(void);

	//----------------------------------- CALCULATION METHODS

	//calculate thermodynamic average of magnetization
	cuReal3 GetThermodynamicAverageMagnetization(cuRect rectangle);

	//As for Get_Histogram, but use thermal averaging in each macrocell
	bool Get_ThAvHistogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, cuINT3 macrocell_dims);

	//As for Get_AngHistogram, but use thermal averaging in each macrocell
	bool Get_ThAvAngHistogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, cuINT3 macrocell_dims, cuReal3 ndir);

	//-----------------------------------OBJECT GETTERS

	cu_obj<ManagedDiffEqFMCUDA>& Get_ManagedDiffEqCUDA(void);
};

#else

class FMeshCUDA :
	public MeshCUDA
{

private:

public:

	//make this object by copying data from the Mesh holding this object
	FMeshCUDA(Mesh* pMesh) :
		MeshCUDA(pMesh)
	{}

	~FMeshCUDA() {}

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//----------------------------------- ODE METHODS IN (ANTI)FERROMAGNETIC MESH : Mesh_FerromagneticCUDA.cu

	//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
	cuReal3 Average_dmdt(cuBox avBox) { return DBL3(); }

	//return average m x dm/dt in the given avRect (relative rect). Here m is the direction vector.
	cuReal3 Average_mxdmdt(cuBox avBox) { return DBL3(); }
};

#endif
#endif