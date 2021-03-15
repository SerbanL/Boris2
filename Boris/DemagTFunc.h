#pragma once

#include "BorisLib.h"

#include "DemagTFunc_Defs.h"
#include "DemagTFunc_Asympt.h"

#include <math.h>

/////////////////////////////////////////// EXACT NEWELL EQUATIONS FOR SHORT DISTANCES

class DemagTFunc {

private:

	std::vector<std::vector<double>> main_sum, f_sum, g_sum, irregular_sum;

	//calculated f and g values in the entire mesh, including boundary points
	//results from f and g functions stored here, then demag tensor elements built from these values as required - this avoids unnecessary repeated calls to f and g
	//sizes (n.x + 2) * (n.y + 2) * (n.z + 2) for self demag fields
	//sizes (2 * n.x + 1) * (2 * n.y + 1) * (2 * n.z + 1) for stray fields without use of any symmetries
	VEC<double> f_vals_xx, f_vals_yy, f_vals_zz, g_vals_xy, g_vals_xz, g_vals_yz;

	//additional spaces for irregular tensor elements calculations
	VEC<double> f_vals_xx_del, f_vals_yy_del, f_vals_zz_del, g_vals_xy_del, g_vals_xz_del, g_vals_yz_del;

private:

	//---------------------f and g VECTORS COMPUTATION

	//calculate f and g values and set them in f_vals, g_vals vectors (also allocate memory)
	
	//regular versions for self-demag
	bool fill_f_vals(INT3 n, DBL3 hRatios, int asymptotic_distance);
	bool fill_g_vals(INT3 n, DBL3 hRatios, int asymptotic_distance);
	bool fill_g2D_vals(INT3 n, DBL3 hRatios, int asymptotic_distance);

	//shifted versions for stray field
	bool fill_f_vals_shifted(INT3 n, DBL3 hRatios, DBL3 shift);
	bool fill_g_vals_shifted(INT3 n, DBL3 hRatios, DBL3 shift);

	//shifted versions for stray field -> in particular for z shift only
	bool fill_f_vals_zshifted(INT3 n, DBL3 hRatios, DBL3 shift, int asymptotic_distance);
	bool fill_g_vals_zshifted(INT3 n, DBL3 hRatios, DBL3 shift, int asymptotic_distance);

	//shifted and irregular versions for stray field (h_src and h_dst only allowed to differ in z values, giving del = sz - dz)
	bool fill_f_vals_shifted_irregular(INT3 n, DBL3 h_src, DBL3 h_dst, DBL3 shift);
	bool fill_g_vals_shifted_irregular(INT3 n, DBL3 h_src, DBL3 h_dst, DBL3 shift);

	//shifted and irregular versions for stray field (h_src and h_dst only allowed to differ in z values, giving del = sz - dz) -> in particular for z shift only
	bool fill_f_vals_zshifted_irregular(INT3 n, DBL3 h_src, DBL3 h_dst, DBL3 shift, int asymptotic_distance);
	bool fill_g_vals_zshifted_irregular(INT3 n, DBL3 h_src, DBL3 h_dst, DBL3 shift, int asymptotic_distance);

	//f and g calculation functions -> basis functions
	double f(double x, double y, double z);
	double g(double x, double y, double z);

	//---------------------TENSOR ELEMENTS USING PRECOMPUTED f and g

	//REGULAR VERSIONS FOR INTERNAL FIELD

	//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
	double Ldia(int i, int j, int k, double hx, double hy, double hz, VEC<double>& f_vals);
	double Lodia(int i, int j, int k, double hx, double hy, double hz, VEC<double>& g_vals);

	//SHIFTED VERSIONS FOR STRAY FIELD

	double Ldia_shifted(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, VEC<double>& f_vals);
	double Lodia_shifted(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, VEC<double>& g_vals);

	//special versions for z shift only
	double Ldia_zshifted_xx_yy(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, VEC<double>& f_vals);
	double Ldia_zshifted_zz(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, VEC<double>& f_vals);
	double Lodia_xy_zshifted(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, VEC<double>& g_vals);
	double Lodia_xz_yz_zshifted(int i, int j, int k, int nx, int ny, int nz, double hx, double hy, double hz, VEC<double>& g_vals);

	//IRREGULAR AND SHIFTED VERSIONS FOR STRAY FIELD

	//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xx and yy components only
	double Ldia_shifted_irregular_xx_yy(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& f_vals, VEC<double>& f_vals_del);
	//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : zz components only
	double Ldia_shifted_irregular_zz(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& f_vals, VEC<double>& f_vals_del);

	//off-diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xy component only
	double Lodia_shifted_irregular_xy(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& g_vals, VEC<double>& g_vals_del);
	//off-diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xz and yz components only
	double Lodia_shifted_irregular_xz_yz(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& g_vals, VEC<double>& g_vals_del);

	//special versions for z shift only

	//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xx and yy components only
	double Ldia_zshifted_irregular_xx_yy(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& f_vals, VEC<double>& f_vals_del);
	//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : zz components only
	double Ldia_zshifted_irregular_zz(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& f_vals, VEC<double>& f_vals_del);

	//off-diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xy component only
	double Lodia_zshifted_irregular_xy(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& g_vals, VEC<double>& g_vals_del);
	//off-diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xz and yz components only
	double Lodia_zshifted_irregular_xz_yz(int i, int j, int k, int nx, int ny, int nz, double tau, VEC<double>& g_vals, VEC<double>& g_vals_del);

	//---------------------SINGLE VALUE COMPUTE

	//get diagonal component at given position for given cellsizes
	double Ldia_single(double x, double y, double z, double hx, double hy, double hz);

	//get off-diagonal component at given position for given cellsizes
	double Lodia_single(double x, double y, double z, double hx, double hy, double hz);

	//self demag coefficient only for given cellsizes
	double SelfDemag(double hx, double hy, double hz);

	//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xx and yy components only
	double Ldia_shifted_irregular_xx_yy_single(double x, double y, double z, double hx, double hy, double sz, double dz);
	//diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : zz components only
	double Ldia_shifted_irregular_zz_single(double x, double y, double z, double sx, double dx, double hy, double hz);

	//off-diagonal component for irregular tensor, where source and destination cells can differ in z cellsize : xy component only
	double Lodia_shifted_irregular_xy_single(double x, double y, double z, double hx, double hy, double sz, double dz);
	//off-diagonal components for irregular tensor, where source and destination cells can differ in z cellsize : xz and yz components only
	double Lodia_shifted_irregular_xz_yz_single(double x, double y, double z, double hx, double sy, double dy, double hz);

public:

	DemagTFunc(void);

	//---------------------ZERO SHIFT VERSION (FOR INTERNAL FIELD)

	//3D

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. 
	//Need mesh dimensions n, convolution mesh dimensions N (this must be a power of 2 and N/2 smallest integer >= n for all dimensions -> from n to N/2 pad with zeroes).
	//Need cellsize hRatios -> these can be normalized (e.g. to have largest value 1)
	//You can opt to set self-demag coefficients to zero (include_self_demag = false);
	bool CalcDiagTens3D(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 hRatios, bool include_self_demag = true, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);
	
	//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens3D(VEC<DBL3> &Dodiag, INT3 n, INT3 N, DBL3 hRatios, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//2D

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
	bool CalcDiagTens2D(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 hRatios, bool include_self_demag = true, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy only) which has sizes given by N. 
	bool CalcOffDiagTens2D(std::vector<double> &Dodiag, INT3 n, INT3 N, DBL3 hRatios, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//---------------------ZERO SHIFT VERSION (FOR INTERNAL FIELD) WITH PBC

	//3D

	//Logical mesh dimensions N (this must be a power of 2).
	//Need cellsize hRatios -> these can be normalized (e.g. to have largest value 1)
	//pbc conditions are applied for directions with non-zero number of images below; if number of images are zero in a dimension, the tensor coefficients are calculated without pbc
	//These versions reduce to the non-pbc versions (almost) if all images are set to zero, but less readable so prefer to keep separate versions.

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
	bool CalcDiagTens3D_PBC(VEC<DBL3> &Ddiag, INT3 N, DBL3 hRatios, 
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE, 
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens3D_PBC(VEC<DBL3> &Dodiag, INT3 N, DBL3 hRatios, 
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE, 
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//2D

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
	bool CalcDiagTens2D_PBC(VEC<DBL3> &Ddiag, INT3 N, DBL3 hRatios, 
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE, 
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens2D_PBC(std::vector<double> &Dodiag, INT3 N, DBL3 hRatios, 
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE, 
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//---------------------SHIFTED VERSIONS

	//3D

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. 
	bool CalcDiagTens3D_Shifted(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 hRatios, DBL3 shift, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens3D_Shifted(VEC<DBL3> &Dodiag, INT3 n, INT3 N, DBL3 hRatios, DBL3 shift, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//2D

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
	bool CalcDiagTens2D_Shifted(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 hRatios, DBL3 shift, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens2D_Shifted(VEC<DBL3> &Dodiag, INT3 n, INT3 N, DBL3 hRatios, DBL3 shift, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//---------------------SHIFTED VERSIONS with PBC

	//pbc conditions are applied for directions with non-zero number of images below; if number of images are zero in a dimension, the tensor coefficients are calculated without pbc
	//These versions reduce to the non-pbc versions (almost) if all images are set to zero, but less readable so prefer to keep separate versions.

	//3D
	
	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. 
	bool CalcDiagTens3D_Shifted_PBC(VEC<DBL3> &Ddiag, INT3 N, DBL3 hRatios, DBL3 shift, 
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE, 
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens3D_Shifted_PBC(VEC<DBL3> &Dodiag, INT3 N, DBL3 hRatios, DBL3 shift, 
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE, 
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//2D

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
	bool CalcDiagTens2D_Shifted_PBC(VEC<DBL3> &Ddiag, INT3 N, DBL3 hRatios, DBL3 shift, 
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens2D_Shifted_PBC(VEC<DBL3> &Dodiag, INT3 N, DBL3 hRatios, DBL3 shift, 
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);
		
	//---------------------SHIFTED AND IRREGULAR VERSIONS

	//2D

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. This applies for irregular cells, specifically for 2D with s.z and d.z allowed to differ; s.x, d.x resp s.y, d.y must be the same.
	bool CalcDiagTens2D_Shifted_Irregular(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 s, DBL3 d, DBL3 shift, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. This applies for irregular cells, specifically for 2D with s.z and d.z allowed to differ; s.x, d.x resp s.y, d.y must be the same.
	bool CalcOffDiagTens2D_Shifted_Irregular(VEC<DBL3> &Dodiag, INT3 n, INT3 N, DBL3 s, DBL3 d, DBL3 shift, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//---------------------SHIFTED AND IRREGULAR VERSIONS with PBC

	//pbc conditions are applied for directions with non-zero number of images below; if number of images are zero in a dimension, the tensor coefficients are calculated without pbc
	//These versions reduce to the non-pbc versions (almost) if all images are set to zero, but less readable so prefer to keep separate versions.

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. This applies for irregular cells, specifically for 2D with s.z and d.z allowed to differ; s.x, d.x resp s.y, d.y must be the same.
	bool CalcDiagTens2D_Shifted_Irregular_PBC(VEC<DBL3> &Ddiag, INT3 N, DBL3 s, DBL3 d, DBL3 shift,
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. This applies for irregular cells, specifically for 2D with s.z and d.z allowed to differ; s.x, d.x resp s.y, d.y must be the same.
	bool CalcOffDiagTens2D_Shifted_Irregular_PBC(VEC<DBL3> &Dodiag, INT3 N, DBL3 s, DBL3 d, DBL3 shift,
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//---------------------SINGLE VALUE COMPUTE (TESTING ONLY)

	//all diagonal components - single value computation version
	DBL3 Ldia_single(DBL3 dist, DBL3 h, bool minus = true);

	//all off-diagonal components - single value computation version
	DBL3 Lodia_single(DBL3 dist, DBL3 h, bool minus = true);

	//---------------------SINGLE VALUE COMPUTE (SELF DEMAG METHODS)

	//Self demag coefficients only (Dxx, Dyy, Dzz) for the given cellsize - you can calculate these separately e.g. if you set include_self_demag = false in the above methods (useful for super-mesh demag).
	DBL3 SelfDemag(DBL3 h, bool minus = true);

	//As above but also add PBC contribution to self demag
	DBL3 SelfDemag_PBC(DBL3 h, DBL3 n, INT3 demag_pbc_images, int asymptotic_distance = ASYMPTOTIC_DISTANCE, bool minus = true);
};