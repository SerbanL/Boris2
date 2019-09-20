#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "MeshCUDA.h"

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

public:


public:

	//make this object by copying data from the Mesh holding this object
	DipoleMeshCUDA(DipoleMesh* pDipoleMesh_);

	~DipoleMeshCUDA();

	//----------------------------------- OTHER IMPORTANT CONTROL METHODS

	//call when the mesh dimensions have changed - sets every quantity to the right dimensions
	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	//----------------------------------- VARIOUS SET METHODS

	//set magnitude for Mdipole
	void Reset_Mdipole(void);

	void SetMagnetisationAngle(cuBReal polar, cuBReal azim);

	void Reset_recalculateStrayField(void) { recalculateStrayField = false; }

	//----------------------------------- VARIOUS GET METHODS

	//check if stray field needs to be recalculated depending on current settings, and prepare Mdipole for stray field recalculation (set to required value)
	bool Check_recalculateStrayField(void);
};

#endif

