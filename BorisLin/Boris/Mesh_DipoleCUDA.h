#pragma once

#include "CompileFlags.h"
#if COMPILECUDA == 1

#include "ErrorHandler.h"
#include "MeshCUDA.h"

#ifdef MESH_COMPILATION_DIPOLE

#include "BorisCUDALib.h"

class DipoleMesh;

class DipoleMeshCUDA :
	public MeshCUDA
{

private:

	//pointer to cpu version of this mesh
	DipoleMesh* pDipoleMesh;

	//dipoles are used to generate stray fields, to be used in ferromagnetic meshes. 
	//When a dipole value changes (direction or strength) then the stray field generated must be recalculated - this flag is checked by the StrayField module where the recalculation is done.
	bool& recalculateStrayField;
	bool& strayField_recalculated;

public:


public:

	//make this object by copying data from the Mesh holding this object
	DipoleMeshCUDA(DipoleMesh* pDipoleMesh_);

	~DipoleMeshCUDA();

	//----------------------------------- IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//----------------------------------- VARIOUS SET METHODS

	//set magnitude for Mdipole
	void Reset_Mdipole(void);

	void SetMagAngle(cuBReal polar, cuBReal azim);

	//shift dipole mesh rectangle by given amount
	void Shift_Dipole(void);

	//----------------------------------- VARIOUS GET METHODS

	//check if stray field needs to be recalculated depending on current settings, and prepare Mdipole for stray field recalculation (set to required value)
	bool CheckRecalculateStrayField(void);
};

#else

class DipoleMeshCUDA :
	public MeshCUDA
{

private:

public:


public:

	//make this object by copying data from the Mesh holding this object
	DipoleMeshCUDA(Mesh* pMesh) :
		MeshCUDA(pMesh)
	{}

	~DipoleMeshCUDA() {}

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//call when a configuration change has occurred - some objects might need to be updated accordingly
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//----------------------------------- VARIOUS SET METHODS

	//set magnitude for Mdipole
	void Reset_Mdipole(void) {}

	//----------------------------------- VARIOUS GET METHODS

	//check if stray field needs to be recalculated depending on current settings, and prepare Mdipole for stray field recalculation (set to required value)
	bool CheckRecalculateStrayField(void) { return false; }
};

#endif
#endif

