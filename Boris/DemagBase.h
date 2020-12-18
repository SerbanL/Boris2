#pragma once

#include "BorisLib.h"
#include "Boris_Enums_Defs.h"
#include "ErrorHandler.h"



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	Abstract base class for Demag modules (there are more than one: micromagnetic, atomistic)
//
//	All modules are implementations of Modules ABC, so this is an additional ABC (so multiple inheritance).
//	The main reason for this, it allows calling a Demag-type method (defined here, implemented in the specific Demag type module) using a Modules base pointer.
//	This is very convenient since Meshes hold Modules base pointers, so e.g. if you want to call the SetPBC method on a mesh, this will only take effect if a Demag-type module is active in it.
//  Thus this works irrespective of the particular mesh type (micromagnetic or atomistic), which have different implementations of DemagBase, and is accomplished by checking return of 
//	dynamic_cast from Modules* to DemagBase*.
//	Also defines common data types and common interface, as these modules have the same structure but do things in a different way.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class DemagBase {

protected:

	//number of pbc images in each dimension (set to zero to disable).
	//There is also a copy of this in ConvolutionData inherited from Convolution - we need another copy here to detect changes
	//these pbc images are applicable in individual demag modules only
	INT3 demag_pbc_images = INT3();

public:

	DemagBase(void)
	{
	}

	virtual ~DemagBase() {}

	//-------------------Setters

	//Set PBC
	virtual BError Set_PBC(INT3 demag_pbc_images_) = 0;

	//-------------------Getters

	//Get PBC images
	INT3 Get_PBC(void) { return demag_pbc_images; }
};

