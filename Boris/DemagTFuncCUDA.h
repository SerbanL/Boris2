#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "DemagTFunc_Defs.h"
#include "DemagTFunc_AsymptCUDA.h"

/////////////////////////////////////////// EXACT NEWELL EQUATIONS FOR SHORT DISTANCES

class DemagTFuncCUDA {

private:

	//calculated f and g values in the entire mesh, including boundary points
	//results from f and g functions stored here, then demag tensor elements built from these values as required - this avoids unnecessary repeated calls to f and g
	//sizes (n.x + 2) * (n.y + 2) * (n.z + 2) for self demag fields
	//sizes (2 * n.x + 1) * (2 * n.y + 1) * (2 * n.z + 1) for stray fields without use of any symmetries
	cu_obj<cuVEC<double>> f_vals_xx, f_vals_yy, f_vals_zz, g_vals_xy, g_vals_xz, g_vals_yz;

	//additional spaces for irregular tensor elements calculations
	cu_obj<cuVEC<double>> f_vals_xx_del, f_vals_yy_del, f_vals_zz_del, g_vals_xy_del, g_vals_xz_del, g_vals_yz_del;

	cu_obj<DemagAsymptoticDiagCUDA> demagAsymptoticDiag_xx, demagAsymptoticDiag_yy, demagAsymptoticDiag_zz;
	cu_obj<DemagAsymptoticOffDiagCUDA> demagAsymptoticOffDiag_xy, demagAsymptoticOffDiag_xz, demagAsymptoticOffDiag_yz;

private:

	//---------------------f and g VECTORS COMPUTATION

	//calculate f and g values and set them in f_vals, g_vals vectors (also allocate memory)

	//regular versions for self-demag
	bool fill_f_vals(cuINT3 n, cuDBL3 hRatios, int asymptotic_distance);
	bool fill_g_vals(cuINT3 n, cuDBL3 hRatios, int asymptotic_distance);
	bool fill_g2D_vals(cuINT3 n, cuDBL3 hRatios, int asymptotic_distance);

	//shifted versions for stray field
	bool fill_f_vals_shifted(cuINT3 n, cuDBL3 hRatios, cuDBL3 shift);
	bool fill_g_vals_shifted(cuINT3 n, cuDBL3 hRatios, cuDBL3 shift);

	//shifted versions for stray field -> in particular for z shift only
	bool fill_f_vals_zshifted(cuINT3 n, cuDBL3 hRatios, cuDBL3 shift, int asymptotic_distance);
	bool fill_g_vals_zshifted(cuINT3 n, cuDBL3 hRatios, cuDBL3 shift, int asymptotic_distance);

	//shifted and irregular versions for stray field (h_src and h_dst only allowed to differ in z values, giving del = sz - dz)
	bool fill_f_vals_shifted_irregular(cuINT3 n, cuDBL3 h_src, cuDBL3 h_dst, cuDBL3 shift);
	bool fill_g_vals_shifted_irregular(cuINT3 n, cuDBL3 h_src, cuDBL3 h_dst, cuDBL3 shift);

	//shifted and irregular versions for stray field (h_src and h_dst only allowed to differ in z values, giving del = sz - dz) -> in particular for z shift only
	bool fill_f_vals_zshifted_irregular(cuINT3 n, cuDBL3 h_src, cuDBL3 h_dst, cuDBL3 shift, int asymptotic_distance);
	bool fill_g_vals_zshifted_irregular(cuINT3 n, cuDBL3 h_src, cuDBL3 h_dst, cuDBL3 shift, int asymptotic_distance);

	//---------------------TENSOR ELEMENTS USING PRECOMPUTED f and g

	//REGULAR VERSIONS FOR INTERNAL FIELD

	//diagonal and off-diagonal elements calculation functions - f_vals and g_vals vectors must be computed before using fill_f_vals, fill_g_vals (or fill_g2D_vals)
	//2D
	void CalcTens2D_Ldia(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33, 
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance);
	void CalcTens2D_Lodia(
		cu_arr<double>& D12,
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance);

	//3D
	void CalcTens3D_Ldia(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance);
	void CalcTens3D_Lodia(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance);

	//2D with PBC
	void CalcTens2D_Ldia_PBC(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
		int x_images, int y_images, int z_images);
	void CalcTens2D_Lodia_PBC(
		cu_arr<double>& D12,
		cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
		int x_images, int y_images, int z_images);

	//3D with PBC
	void CalcTens3D_Ldia_PBC(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
		int x_images, int y_images, int z_images);
	void CalcTens3D_Lodia_PBC(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		cuINT3 N, cuDBL3 hRatios, int sign, int asymptotic_distance,
		int x_images, int y_images, int z_images);

	//SHIFTED VERSIONS FOR STRAY FIELD

	//2D
	void CalcTens2D_Shifted_Ldia(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance);

	void CalcTens2D_Shifted_Lodia(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance);

	//2D with PBC
	void CalcTens2D_Shifted_Ldia_PBC(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
		int x_images, int y_images, int z_images);

	void CalcTens2D_Shifted_Lodia_PBC(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
		int x_images, int y_images, int z_images);

	//3D
	void CalcTens3D_Shifted_Ldia(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance);

	void CalcTens3D_Shifted_Lodia(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance);

	//3D with PBC
	void CalcTens3D_Shifted_Ldia_PBC(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
		int x_images, int y_images, int z_images);

	void CalcTens3D_Shifted_Lodia_PBC(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, int sign, int asymptotic_distance,
		int x_images, int y_images, int z_images);

	//IRREGULAR AND SHIFTED VERSIONS FOR STRAY FIELD

	//2D
	void CalcTens2D_Shifted_Irregular_Ldia(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		cuINT3 n, cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance);

	void CalcTens2D_Shifted_Irregular_Lodia(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		cuINT3 n, cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance);

	//2D with PBC
	void CalcTens2D_Shifted_Irregular_Ldia_PBC(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance,
		int x_images, int y_images, int z_images);

	void CalcTens2D_Shifted_Irregular_Lodia_PBC(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, int sign, int asymptotic_distance,
		int x_images, int y_images, int z_images);

public:

	DemagTFuncCUDA(void);

	//---------------------ZERO SHIFT VERSION (FOR INTERNAL FIELD)

	//3D

	//Compute diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. 
	//Need mesh dimensions n, convolution mesh dimensions N (this must be a power of 2 and N/2 smallest integer >= n for all dimensions -> from n to N/2 pad with zeroes).
	//Need cellsize hRatios -> these can be normalized (e.g. to have largest value 1)
	//You can opt to set self-demag coefficients to zero (include_self_demag = false);
	bool CalcDiagTens3D(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33, 
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, bool include_self_demag = true, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//Compute the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens3D(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23, 
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//2D

	//Compute diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
	bool CalcDiagTens2D(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33, 
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, bool include_self_demag = true, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//Compute the off-diagonal tensor elements (Dxy only) which has sizes given by N. 
	bool CalcOffDiagTens2D(
		cu_arr<double>& D12,
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	
	//---------------------ZERO SHIFT VERSION (FOR INTERNAL FIELD) WITH PBC

	//3D

	//Logical mesh dimensions N (this must be a power of 2).
	//Need cellsize hRatios -> these can be normalized (e.g. to have largest value 1)
	//pbc conditions are applied for directions with non-zero number of images below; if number of images are zero in a dimension, the tensor coefficients are calculated without pbc
	//These versions reduce to the non-pbc versions (almost) if all images are set to zero, but less readable so prefer to keep separate versions.

	//Compute the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
	bool CalcDiagTens3D_PBC(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		cuINT3 N, cuDBL3 hRatios,
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//Compute the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens3D_PBC(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		cuINT3 N, cuDBL3 hRatios,
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//2D

	//Compute the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
	bool CalcDiagTens2D_PBC(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33, 
		cuINT3 N, cuDBL3 hRatios, 
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//Compute the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens2D_PBC(
		cu_arr<double>& D12,
		cuINT3 N, cuDBL3 hRatios,
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//---------------------SHIFTED VERSIONS

	//3D

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. 
	bool CalcDiagTens3D_Shifted(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens3D_Shifted(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);
	
	//2D

	//Compute the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
	bool CalcDiagTens2D_Shifted(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//Compute the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens2D_Shifted(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		cuINT3 n, cuINT3 N, cuDBL3 hRatios, cuDBL3 shift, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);
	
	//---------------------SHIFTED VERSIONS with PBC

	//pbc conditions are applied for directions with non-zero number of images below; if number of images are zero in a dimension, the tensor coefficients are calculated without pbc
	//These versions reduce to the non-pbc versions (almost) if all images are set to zero, but less readable so prefer to keep separate versions.

	//3D

	//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. 
	bool CalcDiagTens3D_Shifted_PBC(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		cuINT3 N, cuDBL3 hRatios, cuDBL3 shift,
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens3D_Shifted_PBC(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		cuINT3 N, cuDBL3 hRatios, cuDBL3 shift,
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//2D

	//Compute diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
	bool CalcDiagTens2D_Shifted_PBC(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		INT3 N, DBL3 hRatios, DBL3 shift,
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//Compute off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
	bool CalcOffDiagTens2D_Shifted_PBC(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		INT3 N, DBL3 hRatios, DBL3 shift,
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);
		
	//---------------------SHIFTED AND IRREGULAR VERSIONS

	//2D

	//Compute the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. This applies for irregular cells, specifically for 2D with s.z and d.z allowed to differ; s.x, d.x resp s.y, d.y must be the same.
	bool CalcDiagTens2D_Shifted_Irregular(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33, 
		cuINT3 n, cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);

	//Compute the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. This applies for irregular cells, specifically for 2D with s.z and d.z allowed to differ; s.x, d.x resp s.y, d.y must be the same.
	bool CalcOffDiagTens2D_Shifted_Irregular(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		cuINT3 n, cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift, bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE);
	
	//---------------------SHIFTED AND IRREGULAR VERSIONS with PBC

	//pbc conditions are applied for directions with non-zero number of images below; if number of images are zero in a dimension, the tensor coefficients are calculated without pbc
	//These versions reduce to the non-pbc versions (almost) if all images are set to zero, but less readable so prefer to keep separate versions.

	//Compute the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N. This applies for irregular cells, specifically for 2D with s.z and d.z allowed to differ; s.x, d.x resp s.y, d.y must be the same.
	bool CalcDiagTens2D_Shifted_Irregular_PBC(
		cu_arr<double>& D11, cu_arr<double>& D22, cu_arr<double>& D33,
		cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift,
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);

	//Compute the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. This applies for irregular cells, specifically for 2D with s.z and d.z allowed to differ; s.x, d.x resp s.y, d.y must be the same.
	bool CalcOffDiagTens2D_Shifted_Irregular_PBC(
		cu_arr<double>& D12, cu_arr<double>& D13, cu_arr<double>& D23,
		cuINT3 N, cuDBL3 s, cuDBL3 d, cuDBL3 shift,
		bool minus = true, int asymptotic_distance = ASYMPTOTIC_DISTANCE,
		int x_images = PBC_X_IMAGES, int y_images = PBC_Y_IMAGES, int z_images = PBC_Z_IMAGES);
};

#endif