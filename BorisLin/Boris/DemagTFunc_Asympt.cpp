#include "stdafx.h"
#include "DemagTFunc_Asympt.h"

//I've adapted the asymptotic equations from OOMMF code, not my work - credit goes to M.J. Donahue.
//I've checked them, no point reinventing the wheel especially since it's quite a lot of work to write them.
//I did code the rest from scratch though (exact Newell equations implementation).

DemagAsymptoticDiag::DemagAsymptoticDiag(double hx, double hy, double hz)
{
	double hx2 = hx * hx;
	double hy2 = hy * hy;
	double hz2 = hz * hz;

	double hx4 = hx2 * hx2;
	double hy4 = hy2 * hy2;
	double hz4 = hz2 * hz2;

	double hx6 = hx4 * hx2;
	double hy6 = hy4 * hy2;
	double hz6 = hz4 * hz2;

	lead_weight = (-hx * hy * hz / (4 * PI));

	//Initialize coefficients for 1/R^5 term
	if (hx2 != hy2 || hx2 != hz2 || hy2 != hz2) {

		//Non-cubic cell
		cubic_cell = false;

		a1 = a2 = a3 = a4 = a5 = a6 = lead_weight / 4.0;
		a1 *= 8 * hx2 - 4 * hy2 - 4 * hz2;
		a2 *= -24 * hx2 + 27 * hy2 - 3 * hz2;
		a3 *= -24 * hx2 - 3 * hy2 + 27 * hz2;
		a4 *= 3 * hx2 - 4 * hy2 + 1 * hz2;
		a5 *= 6 * hx2 - 3 * hy2 - 3 * hz2;
		a6 *= 3 * hx2 + 1 * hy2 - 4 * hz2;
	}
	else { 

		//cubic cell
		cubic_cell = true;

		a1 = a2 = a3 = a4 = a5 = a6 = 0.0;
	}

	//Initialize coefficients for 1/R^7 term
	b1 = b2 = b3 = b4 = b5 = b6 = b7 = b8 = b9 = b10 = lead_weight / 16.0;

	if (cubic_cell) {

		b1 *= -14 * hx4;
		b2 *= 105 * hx4;
		b3 *= 105 * hx4;
		b4 *= -105 * hx4;
		b6 *= -105 * hx4;
		b7 *= 7 * hx4;
		b10 *= 7 * hx4;
		b5 = b8 = b9 = 0;
	}
	else {

		b1 *= 32 * hx4 - 40 * hx2*hy2 - 40 * hx2*hz2 + 12 * hy4 + 10 * hy2*hz2 + 12 * hz4;
		b2 *= -240 * hx4 + 580 * hx2*hy2 + 20 * hx2*hz2 - 202 * hy4 - 75 * hy2*hz2 + 22 * hz4;
		b3 *= -240 * hx4 + 20 * hx2*hy2 + 580 * hx2*hz2 + 22 * hy4 - 75 * hy2*hz2 - 202 * hz4;
		b4 *= 180 * hx4 - 505 * hx2*hy2 + 55 * hx2*hz2 + 232 * hy4 - 75 * hy2*hz2 + 8 * hz4;
		b5 *= 360 * hx4 - 450 * hx2*hy2 - 450 * hx2*hz2 - 180 * hy4 + 900 * hy2*hz2 - 180 * hz4;
		b6 *= 180 * hx4 + 55 * hx2*hy2 - 505 * hx2*hz2 + 8 * hy4 - 75 * hy2*hz2 + 232 * hz4;
		b7 *= -10 * hx4 + 30 * hx2*hy2 - 5 * hx2*hz2 - 16 * hy4 + 10 * hy2*hz2 - 2 * hz4;
		b8 *= -30 * hx4 + 55 * hx2*hy2 + 20 * hx2*hz2 + 8 * hy4 - 75 * hy2*hz2 + 22 * hz4;
		b9 *= -30 * hx4 + 20 * hx2*hy2 + 55 * hx2*hz2 + 22 * hy4 - 75 * hy2*hz2 + 8 * hz4;
		b10 *= -10 * hx4 - 5 * hx2*hy2 + 30 * hx2*hz2 - 2 * hy4 + 10 * hy2*hz2 - 16 * hz4;
	}

	//Initialize coefficients for 1/R^9 term
	c1 = c2 = c3 = c4 = c5 = c6 = c7 = c8 = c9 = c10 = c11 = c12 = c13 = c14 = c15 = lead_weight / 192.0;
	
	if (cubic_cell) {

		c1 *= 32 * hx6;
		c2 *= -448 * hx6;
		c3 *= -448 * hx6;
		c4 *= -150 * hx6;
		c5 *= 7620 * hx6;
		c6 *= -150 * hx6;
		c7 *= 314 * hx6;
		c8 *= -3810 * hx6;
		c9 *= -3810 * hx6;
		c10 *= 314 * hx6;
		c11 *= -16 * hx6;
		c12 *= 134 * hx6;
		c13 *= 300 * hx6;
		c14 *= 134 * hx6;
		c15 *= -16 * hx6;
	}
	else {

		c1 *= 384 * hx6 + -896 * hx4*hy2 + -896 * hx4*hz2 + 672 * hx2*hy4 + 560 * hx2*hy2*hz2 + 672 * hx2*hz4 + -120 * hy6 + -112 * hy4*hz2 + -112 * hy2*hz4 + -120 * hz6;
		c2 *= -5376 * hx6 + 22624 * hx4*hy2 + 2464 * hx4*hz2 + -19488 * hx2*hy4 + -7840 * hx2*hy2*hz2 + 672 * hx2*hz4 + 3705 * hy6 + 2198 * hy4*hz2 + 938 * hy2*hz4 + -345 * hz6;
		c3 *= -5376 * hx6 + 2464 * hx4*hy2 + 22624 * hx4*hz2 + 672 * hx2*hy4 + -7840 * hx2*hy2*hz2 + -19488 * hx2*hz4 + -345 * hy6 + 938 * hy4*hz2 + 2198 * hy2*hz4 + 3705 * hz6;
		c4 *= 10080 * hx6 + -48720 * hx4*hy2 + 1680 * hx4*hz2 + 49770 * hx2*hy4 + -2625 * hx2*hy2*hz2 + -630 * hx2*hz4 + -10440 * hy6 + -1050 * hy4*hz2 + 2100 * hy2*hz4 + -315 * hz6;
		c5 *= 20160 * hx6 + -47040 * hx4*hy2 + -47040 * hx4*hz2 + -6300 * hx2*hy4 + 133350 * hx2*hy2*hz2 + -6300 * hx2*hz4 + 7065 * hy6 + -26670 * hy4*hz2 + -26670 * hy2*hz4 + 7065 * hz6;
		c6 *= 10080 * hx6 + 1680 * hx4*hy2 + -48720 * hx4*hz2 + -630 * hx2*hy4 + -2625 * hx2*hy2*hz2 + 49770 * hx2*hz4 + -315 * hy6 + 2100 * hy4*hz2 + -1050 * hy2*hz4 + -10440 * hz6;
		c7 *= -3360 * hx6 + 17290 * hx4*hy2 + -1610 * hx4*hz2 + -19488 * hx2*hy4 + 5495 * hx2*hy2*hz2 + -588 * hx2*hz4 + 4848 * hy6 + -3136 * hy4*hz2 + 938 * hy2*hz4 + -75 * hz6;
		c8 *= -10080 * hx6 + 32970 * hx4*hy2 + 14070 * hx4*hz2 + -6300 * hx2*hy4 + -66675 * hx2*hy2*hz2 + 12600 * hx2*hz4 + -10080 * hy6 + 53340 * hy4*hz2 + -26670 * hy2*hz4 + 3015 * hz6;
		c9 *= -10080 * hx6 + 14070 * hx4*hy2 + 32970 * hx4*hz2 + 12600 * hx2*hy4 + -66675 * hx2*hy2*hz2 + -6300 * hx2*hz4 + 3015 * hy6 + -26670 * hy4*hz2 + 53340 * hy2*hz4 + -10080 * hz6;
		c10 *= -3360 * hx6 + -1610 * hx4*hy2 + 17290 * hx4*hz2 + -588 * hx2*hy4 + 5495 * hx2*hy2*hz2 + -19488 * hx2*hz4 + -75 * hy6 + 938 * hy4*hz2 + -3136 * hy2*hz4 + 4848 * hz6;
		c11 *= 105 * hx6 + -560 * hx4*hy2 + 70 * hx4*hz2 + 672 * hx2*hy4 + -280 * hx2*hy2*hz2 + 42 * hx2*hz4 + -192 * hy6 + 224 * hy4*hz2 + -112 * hy2*hz4 + 15 * hz6;
		c12 *= 420 * hx6 + -1610 * hx4*hy2 + -350 * hx4*hz2 + 672 * hx2*hy4 + 2345 * hx2*hy2*hz2 + -588 * hx2*hz4 + 528 * hy6 + -3136 * hy4*hz2 + 2198 * hy2*hz4 + -345 * hz6;
		c13 *= 630 * hx6 + -1470 * hx4*hy2 + -1470 * hx4*hz2 + -630 * hx2*hy4 + 5250 * hx2*hy2*hz2 + -630 * hx2*hz4 + 360 * hy6 + -1050 * hy4*hz2 + -1050 * hy2*hz4 + 360 * hz6;
		c14 *= 420 * hx6 + -350 * hx4*hy2 + -1610 * hx4*hz2 + -588 * hx2*hy4 + 2345 * hx2*hy2*hz2 + 672 * hx2*hz4 + -345 * hy6 + 2198 * hy4*hz2 + -3136 * hy2*hz4 + 528 * hz6;
		c15 *= 105 * hx6 + 70 * hx4*hy2 + -560 * hx4*hz2 + 42 * hx2*hy4 + -280 * hx2*hy2*hz2 + 672 * hx2*hz4 + 15 * hy6 + -112 * hy4*hz2 + 224 * hy2*hz4 + -192 * hz6;
	}
}

double DemagAsymptoticDiag::AsymptoticLdia(double x, double y, double z)
{
	double tx2 = 0, ty2 = 0, tz2 = 0;
	double R = 0, iR = 0;
	double R2 = 0, iR2 = 0;

	R2 = x * x + y * y + z * z;
	R = sqrt(R2);

	if (R) {

		tx2 = x * x / (R2*R2);
		ty2 = y * y / (R2*R2);
		tz2 = z * z / (R2*R2);

		iR = 1 / R;
		iR2 = 1 / R2;
	}

	if (iR2 <= 0.0) {

		//Asymptotic expansion doesn't apply for R==0. Don't use!
		return 0.0;
	}

	double tz4 = tz2 * tz2;
	double tz6 = tz4 * tz2;
	double term3 = (2 * tx2 - ty2 - tz2) * lead_weight;

	double term5 = 0.0;
	double term7 = 0.0;

	if (cubic_cell) {

		double ty4 = ty2 * ty2;

		term7 = ((b1 * tx2
			+ (b2*ty2 + b3 * tz2)) * tx2
			+ (b4*ty4 + b6 * tz4)) * tx2
			+ b7 * ty4*ty2 + b10 * tz6;
	}
	else {

		term5 = (a1*tx2 + (a2*ty2 + a3 * tz2)) * tx2
			+ (a4*ty2 + a5 * tz2) * ty2 + a6 * tz4;
		
		term7 = ((b1 * tx2
			+ (b2*ty2 + b3 * tz2)) * tx2
			+ ((b4*ty2 + b5 * tz2)*ty2 + b6 * tz4)) * tx2
			+ ((b7*ty2 + b8 * tz2)*ty2 + b9 * tz4) * ty2
			+ b10 * tz6;
	}

	double term9 = (((c1*tx2
		+ (c2*ty2 + c3 * tz2)) * tx2
		+ ((c4*ty2 + c5 * tz2) * ty2 + c6 * tz4)) * tx2
		+ (((c7*ty2 + c8 * tz2)*ty2 + c9 * tz4)*ty2 + c10 * tz6)) * tx2
		+ (((c11*ty2 + c12 * tz2)*ty2 + c13 * tz4)*ty2 + c14 * tz6) * ty2
		+ c15*tz4*tz4;

	//Error should be of order 1/R^11
	return (term9 + term7 + term5 + term3)*iR;
}


//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------

DemagAsymptoticOffDiag::DemagAsymptoticOffDiag(double hx, double hy, double hz)
{
	double hx2 = hx * hx;
	double hy2 = hy * hy;
	double hz2 = hz * hz;

	double hx4 = hx2 * hx2;
	double hy4 = hy2 * hy2;
	double hz4 = hz2 * hz2;

	double hx6 = hx4 * hx2;
	double hy6 = hy4 * hy2;
	double hz6 = hz4 * hz2;

	lead_weight = (-hx * hy*hz / (4 * PI));

	// Initialize coefficients for 1/R^5 term
	if (hx2 != hy2 || hx2 != hz2 || hy2 != hz2) { 
		
		//Non-cubic cell
		cubic_cell = false;

		a1 = a2 = a3 = (lead_weight*5.0) / 4.0;
		a1 *= 4 * hx2 - 3 * hy2 - 1 * hz2;
		a2 *= -3 * hx2 + 4 * hy2 - 1 * hz2;
		a3 *= -3 * hx2 - 3 * hy2 + 6 * hz2;
	}
	else { 

		//Cubic cell
		cubic_cell = true;

		a1 = a2 = a3 = 0.0;
	}

	// Initialize coefficients for 1/R^7 term
	b1 = b2 = b3 = b4 = b5 = b6 = (lead_weight*7.0) / 16.0;

	if (cubic_cell) {

		b1 *= -7 * hx4;
		b2 *= 19 * hx4;
		b3 *= 13 * hx4;
		b4 *= -7 * hx4;
		b5 *= 13 * hx4;
		b6 *= -13 * hx4;
	}
	else {

		b1 *= 16 * hx4 - 30 * hx2*hy2 - 10 * hx2*hz2 + 10 * hy4 + 5 * hy2*hz2 + 2 * hz4;
		b2 *= -40 * hx4 + 105 * hx2*hy2 - 5 * hx2*hz2 - 40 * hy4 - 5 * hy2*hz2 + 4 * hz4;
		b3 *= -40 * hx4 - 15 * hx2*hy2 + 115 * hx2*hz2 + 20 * hy4 - 35 * hy2*hz2 - 32 * hz4;
		b4 *= 10 * hx4 - 30 * hx2*hy2 + 5 * hx2*hz2 + 16 * hy4 - 10 * hy2*hz2 + 2 * hz4;
		b5 *= 20 * hx4 - 15 * hx2*hy2 - 35 * hx2*hz2 - 40 * hy4 + 115 * hy2*hz2 - 32 * hz4;
		b6 *= 10 * hx4 + 15 * hx2*hy2 - 40 * hx2*hz2 + 10 * hy4 - 40 * hy2*hz2 + 32 * hz4;
	}

	// Initialize coefficients for 1/R^9 term
	c1 = c2 = c3 = c4 = c5 = c6 = c7 = c8 = c9 = c10 = lead_weight / 64.0;

	if (cubic_cell) {

		c1 *= 48 * hx6;
		c2 *= -142 * hx6;
		c3 *= -582 * hx6;
		c4 *= -142 * hx6;
		c5 *= 2840 * hx6;
		c6 *= -450 * hx6;
		c7 *= 48 * hx6;
		c8 *= -582 * hx6;
		c9 *= -450 * hx6;
		c10 *= 180 * hx6;
	}
	else {

		c1 *= 576 * hx6 + -2016 * hx4*hy2 + -672 * hx4*hz2 + 1680 * hx2*hy4 + 840 * hx2*hy2*hz2 + 336 * hx2*hz4 + -315 * hy6 + -210 * hy4*hz2 + -126 * hy2*hz4 + -45 * hz6;
		c2 *= -3024 * hx6 + 13664 * hx4*hy2 + 448 * hx4*hz2 + -12670 * hx2*hy4 + -2485 * hx2*hy2*hz2 + 546 * hx2*hz4 + 2520 * hy6 + 910 * hy4*hz2 + 84 * hy2*hz4 + -135 * hz6;
		c3 *= -3024 * hx6 + 1344 * hx4*hy2 + 12768 * hx4*hz2 + 2730 * hx2*hy4 + -10185 * hx2*hy2*hz2 + -8694 * hx2*hz4 + -945 * hy6 + 1680 * hy4*hz2 + 2394 * hy2*hz4 + 1350 * hz6;
		c4 *= 2520 * hx6 + -12670 * hx4*hy2 + 910 * hx4*hz2 + 13664 * hx2*hy4 + -2485 * hx2*hy2*hz2 + 84 * hx2*hz4 + -3024 * hy6 + 448 * hy4*hz2 + 546 * hy2*hz4 + -135 * hz6;
		c5 *= 5040 * hx6 + -9940 * hx4*hy2 + -13580 * hx4*hz2 + -9940 * hx2*hy4 + 49700 * hx2*hy2*hz2 + -6300 * hx2*hz4 + 5040 * hy6 + -13580 * hy4*hz2 + -6300 * hy2*hz4 + 2700 * hz6;
		c6 *= 2520 * hx6 + 2730 * hx4*hy2 + -14490 * hx4*hz2 + 420 * hx2*hy4 + -7875 * hx2*hy2*hz2 + 17640 * hx2*hz4 + -945 * hy6 + 3990 * hy4*hz2 + -840 * hy2*hz4 + -3600 * hz6;
		c7 *= -315 * hx6 + 1680 * hx4*hy2 + -210 * hx4*hz2 + -2016 * hx2*hy4 + 840 * hx2*hy2*hz2 + -126 * hx2*hz4 + 576 * hy6 + -672 * hy4*hz2 + 336 * hy2*hz4 + -45 * hz6;
		c8 *= -945 * hx6 + 2730 * hx4*hy2 + 1680 * hx4*hz2 + 1344 * hx2*hy4 + -10185 * hx2*hy2*hz2 + 2394 * hx2*hz4 + -3024 * hy6 + 12768 * hy4*hz2 + -8694 * hy2*hz4 + 1350 * hz6;
		c9 *= -945 * hx6 + 420 * hx4*hy2 + 3990 * hx4*hz2 + 2730 * hx2*hy4 + -7875 * hx2*hy2*hz2 + -840 * hx2*hz4 + 2520 * hy6 + -14490 * hy4*hz2 + 17640 * hy2*hz4 + -3600 * hz6;
		c10 *= -315 * hx6 + -630 * hx4*hy2 + 2100 * hx4*hz2 + -630 * hx2*hy4 + 3150 * hx2*hy2*hz2 + -3360 * hx2*hz4 + -315 * hy6 + 2100 * hy4*hz2 + -3360 * hy2*hz4 + 1440 * hz6;
	}
}

double DemagAsymptoticOffDiag::AsymptoticLodia(double x, double y, double z)
{
	double tx2 = 0, ty2 = 0, tz2 = 0;
	double R = 0, iR = 0;
	double R2 = 0, iR2 = 0;

	R2 = x * x + y * y + z * z;
	R = sqrt(R2);

	if (R) {

		tx2 = x * x / (R2*R2);
		ty2 = y * y / (R2*R2);
		tz2 = z * z / (R2*R2);

		iR = 1 / R;
		iR2 = 1 / R2;
	}

	if (R2 <= 0.0) {

		// Asymptotic expansion doesn't apply for R==0. Don't use!
		return 0.0;
	}

	double term3 = 3 * lead_weight;

	double term5 = 0.0;

	if (!cubic_cell) {

		term5 = a1 * tx2 + a2 * ty2 + a3 * tz2;
	}

	double tz4 = tz2 * tz2;

	double term7 = (b1*tx2 + (b2*ty2 + b3 * tz2))*tx2 + (b4*ty2 + b5 * tz2)*ty2 + b6 * tz4;

	double term9 = ((c1*tx2 
		+ (c2*ty2 + c3 * tz2)) * tx2 
		+ ((c4*ty2 + c5 * tz2) * ty2 + c6 * tz4)) * tx2
		+ ((c7*ty2 + c8 * tz2)*ty2 + c9 * tz4) * ty2
		+ c10 * tz4*tz2;

	double iR5 = iR2 * iR2*iR;

	// Error should be of order 1/R^11
	return (term9 + term7 + term5 + term3) * iR5*x*y;
}