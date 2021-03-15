#pragma once

#include "BorisLib.h"

#include "DemagTFunc_Defs.h"

#include <math.h>

/////////////////////////////////////////// EXACT NEWELL EQUATIONS FOR SHORT DISTANCES

class DipoleDipoleTFunc {

private:

	//get self demag coefficients for a prism with dimensions hRatios (normalize h to largest side) - use these for macrocell tensors
	//use Newell formulas to compute these (usual N = +1/3 for a cube, but macrocell size is not restricted to a cube)
	double SelfDemag(double hx, double hy, double hz);

public:

	DipoleDipoleTFunc(void);

	//---------------------ZERO SHIFT VERSION (FOR INTERNAL FIELD) :  : DipoleDipoleTFunnc.cpp

	//3D

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. 
	//Need mesh dimensions n, convolution mesh dimensions N (this must be a power of 2 and N/2 smallest integer >= n for all dimensions -> from n to N/2 pad with zeroes).
	//Note unlike the demag tensor, h are the actual cellsize values in m, cannot be normalized. However since the moment is used in units of muB, need to multiply tensor by muB so correct result in obtained after convolution.
	//This works out well since r^3 is a very small number, but muB/r^3 has workable order of magnitude (maximum ~10^5 to 10^6 neighboring dipoles)
	//There's also the option of including a self demag term : normally the 0, 0, 0 point in the tensor is zero, but if a macrocell is considered instead it must include a self demag term.
	//If using a macrocell, bear in mind this must be the same number of atomic unit cells in each dimension (i.e. cubic if unit cell is cubic), otherwise the result in incorrect.
	bool CalcDiagTens3D(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 h, bool include_self_demag);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens3D(VEC<DBL3> &Dodiag, INT3 n, INT3 N, DBL3 h);

	//2D

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
	bool CalcDiagTens2D(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 h, bool include_self_demag);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy only) which has sizes given by N. 
	bool CalcOffDiagTens2D(std::vector<double> &Dodiag, INT3 n, INT3 N, DBL3 h);

	//---------------------ZERO SHIFT VERSION (FOR INTERNAL FIELD) WITH PBC : DipoleDipoleTFunnc_PBC.cpp

	//3D

	//Logical mesh dimensions N (this must be a power of 2).
	//Need cellsize hRatios -> these can be normalized (e.g. to have largest value 1)
	//pbc conditions are applied for directions with non-zero number of images below; if number of images are zero in a dimension, the tensor coefficients are calculated without pbc
	//These versions reduce to the non-pbc versions (almost) if all images are set to zero, but less readable so prefer to keep separate versions.

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
	bool CalcDiagTens3D_PBC(VEC<DBL3> &Ddiag, INT3 N, DBL3 h, bool include_self_demag,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens3D_PBC(VEC<DBL3> &Dodiag, INT3 N, DBL3 h,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//2D

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
	bool CalcDiagTens2D_PBC(VEC<DBL3> &Ddiag, INT3 N, DBL3 h, bool include_self_demag,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens2D_PBC(std::vector<double> &Dodiag, INT3 N, DBL3 h,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//---------------------SELF DEMAG (for macrocell) : DipoleDipoleTFunnc_Self.cpp

	DBL3 SelfDemag(DBL3 hRatios);

	DBL3 SelfDemag_PBC(DBL3 h, DBL3 n, INT3 demag_pbc_images);
};
