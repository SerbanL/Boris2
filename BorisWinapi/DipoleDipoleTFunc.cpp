#include "stdafx.h"
#include "DipoleDipoleTFunc.h"

DipoleDipoleTFunc::DipoleDipoleTFunc(void)
{
}

//---------------------ZERO SHIFT VERSIONS (FOR INTERNAL FIELD)

bool DipoleDipoleTFunc::CalcDiagTens3D(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 h, bool include_self_demag)
{
	//zero the tensor first
	Ddiag.set(DBL3());

	//Calculate tensor diagonal elements in corect position ready for FFT
	//Use the following symmetries to speed-up calculation (exactly the same as for demag tensor):

	//1. D22(n, m, k) = D11(m, n, k) and D33(n, m, k) = D11(k, m, n)
	//2. D11(n, m, k) = D11(|n|, |m|, |k|)

	int I = N.x;
	int J = N.y;
	int K = N.z;

	// in the first octant only need to generate values up to n, not N / 2.
	//n may not be a power of 2, in which case any points between n and N/2 will be padded with zeroes, so we don't want the tensor coefficients from n to N/2 - keep them zero.
#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {
		for (int k = 0; k < n.z; k++) {
			for (int i = 0; i < n.x; i++) {

				if (i == 0 && j == 0 && k == 0) {

					if (include_self_demag) Ddiag[0] = SelfDemag(h / h.maxdim()) * (-MUB / h.dim());
					continue;
				}

				//displacement vector
				DBL3 r = DBL3(i * h.x, j * h.y, k * h.z);

				//length of displacement vector
				double r_norm = r.norm();

				//prefactor
				double c = MUB / (4 * PI * r_norm * r_norm * r_norm);

				//unit displacement vector
				DBL3 r_dir = r / r_norm;

				//D11, D22, D33
				DBL3 val = DBL3(3 * r_dir.x*r_dir.x - 1.0, 3 * r_dir.y*r_dir.y - 1.0, 3 * r_dir.z*r_dir.z - 1.0) * c;

				Ddiag[((i + I) % I) + ((j + J) % J)*I + ((k + K) % K)*I*J] = val;
				Ddiag[((-i + I) % I) + ((j + J) % J)*I + ((k + K) % K)*I*J] = val;
				Ddiag[((i + I) % I) + ((-j + J) % J)*I + ((k + K) % K)*I*J] = val;
				Ddiag[((i + I) % I) + ((j + J) % J)*I + ((-k + K) % K)*I*J] = val;
				Ddiag[((-i + I) % I) + ((-j + J) % J)*I + ((k + K) % K)*I*J] = val;
				Ddiag[((-i + I) % I) + ((j + J) % J)*I + ((-k + K) % K)*I*J] = val;
				Ddiag[((i + I) % I) + ((-j + J) % J)*I + ((-k + K) % K)*I*J] = val;
				Ddiag[((-i + I) % I) + ((-j + J) % J)*I + ((-k + K) % K)*I*J] = val;
			}
		}
	}

	return true;
}

bool DipoleDipoleTFunc::CalcOffDiagTens3D(VEC<DBL3> &Dodiag, INT3 n, INT3 N, DBL3 h)
{
	//zero the tensor first
	Dodiag.set(DBL3());

	//Calculate demag tensor elements in corect position ready for FFT
	//Use the following symmetries to speed-up calculation (exactly the same as for demag tensor):

	//1. D13(n, m, k) = D12(n, k, m) and D23(n, m, k) = D12(m, k, n)
	//2. D12(n, m, k) = sign(n) * sign(m) D12(|n|, |m|, |k|)
	//thus:
	//3. D13(n, m, k) = sign(n) * sign(k) D12(|n|, |k|, |m|)
	//4. D23(n, m, k) = sign(m) * sign(k) D12(|m|, |k|, |n|)

	int I = N.x;
	int J = N.y;
	int K = N.z;

	//in the first octant only need to generate values up to n, not N/2.
	//n may not be a power of 2, in which case any points between n and N/2 will be padded with zeroes, so we don't want the tensor coefficients from n to N/2 - keep them zero.
#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {
		for (int k = 0; k < n.z; k++) {
			for (int i = 0; i < n.x; i++) {

				if (i == 0 && j == 0 && k == 0) continue;

				//displacement vector
				DBL3 r = DBL3(i * h.x, j * h.y, k * h.z);

				//length of displacement vector
				double r_norm = r.norm();

				//prefactor
				double c = MUB / (4 * PI * r_norm * r_norm * r_norm);

				//unit displacement vector
				DBL3 r_dir = r / r_norm;

				//D12, D13, D23
				DBL3 val = DBL3(3 * r_dir.x*r_dir.y, 3 * r_dir.x*r_dir.z, 3 * r_dir.y*r_dir.z) * c;

				Dodiag[((i + I) % I) + ((j + J) % J)*I + ((k + K) % K)*I*J] = val & DBL3(+1, +1, +1);
				Dodiag[((-i + I) % I) + ((j + J) % J)*I + ((k + K) % K)*I*J] = val & DBL3(-1, -1, +1);
				Dodiag[((i + I) % I) + ((-j + J) % J)*I + ((k + K) % K)*I*J] = val & DBL3(-1, +1, -1);
				Dodiag[((i + I) % I) + ((j + J) % J)*I + ((-k + K) % K)*I*J] = val & DBL3(+1, -1, -1);
				Dodiag[((-i + I) % I) + ((-j + J) % J)*I + ((k + K) % K)*I*J] = val & DBL3(+1, -1, -1);
				Dodiag[((-i + I) % I) + ((j + J) % J)*I + ((-k + K) % K)*I*J] = val & DBL3(-1, +1, -1);
				Dodiag[((i + I) % I) + ((-j + J) % J)*I + ((-k + K) % K)*I*J] = val & DBL3(-1, -1, +1);
				Dodiag[((-i + I) % I) + ((-j + J) % J)*I + ((-k + K) % K)*I*J] = val & DBL3(+1, +1, +1);
			}
		}
	}

	return true;
}

bool DipoleDipoleTFunc::CalcDiagTens2D(VEC<DBL3> &Ddiag, INT3 n, INT3 N, DBL3 h, bool include_self_demag)
{
	//zero the tensor first
	Ddiag.set(DBL3());

	//Calculate demag tensor diagonal elements in corect position ready for FFT
	//Use the following symmetries to speed-up calculation:

	//1. D22(n, m, k) = D11(m, n, k) and D33(n, m, k) = D11(k, m, n)
	//2. D11(n, m, k) = D11(|n|, |m|, |k|)

	int I = N.x;
	int J = N.y;

#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {
		for (int i = 0; i < n.x; i++) {

			if (i == 0 && j == 0) {

				if (include_self_demag) Ddiag[0] = SelfDemag(h / h.maxdim()) * (-MUB / h.dim());
				continue;
			}

			//displacement vector
			DBL3 r = DBL3(i * h.x, j * h.y, 0.0);

			//length of displacement vector
			double r_norm = r.norm();

			//prefactor
			double c = MUB / (4 * PI * r_norm * r_norm * r_norm);

			//unit displacement vector
			DBL3 r_dir = r / r_norm;

			//D11, D22, D33
			DBL3 val = DBL3(3 * r_dir.x*r_dir.x - 1.0, 3 * r_dir.y*r_dir.y - 1.0, -1.0) * c;

			Ddiag[((i + I) % I) + ((j + J) % J)*I] = val;
			Ddiag[((-i + I) % I) + ((j + J) % J)*I] = val;
			Ddiag[((i + I) % I) + ((-j + J) % J)*I] = val;
			Ddiag[((-i + I) % I) + ((-j + J) % J)*I] = val;
		}
	}

	return true;
}

bool DipoleDipoleTFunc::CalcOffDiagTens2D(std::vector<double> &Dodiag, INT3 n, INT3 N, DBL3 h)
{
	//zero the tensor first
	std::fill(Dodiag.begin(), Dodiag.end(), 0.0);

	//Calculate demag tensor elements in corect position ready for FFT
	//Use the following symmetries to speed-up calculation:

	//1. D13(n, m, k) = D12(n, k, m) and D23(n, m, k) = D12(m, k, n)
	//2. D12(n, m, k) = sign(n) * sign(m) D12(|n|, |m|, |k|)
	//thus:
	//3. D13(n, m, k) = sign(n) * sign(k) D12(|n|, |k|, |m|)
	//4. D23(n, m, k) = sign(m) * sign(k) D12(|m|, |k|, |n|)

	int I = N.x;
	int J = N.y;

#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {
		for (int i = 0; i < n.x; i++) {

			if (i == 0 && j == 0) continue;

			//displacement vector
			DBL3 r = DBL3(i * h.x, j * h.y, 0.0);

			//length of displacement vector
			double r_norm = r.norm();

			//prefactor
			double c = MUB / (4 * PI * r_norm * r_norm * r_norm);

			//unit displacement vector
			DBL3 r_dir = r / r_norm;

			//D12
			double val = 3 * r_dir.x*r_dir.y * c;

			Dodiag[((i + I) % I) + ((j + J) % J)*I] = val;
			Dodiag[((-i + I) % I) + ((j + J) % J)*I] = -1 * val;
			Dodiag[((i + I) % I) + ((-j + J) % J)*I] = -1 * val;
			Dodiag[((-i + I) % I) + ((-j + J) % J)*I] = val;
		}
	}

	return true;
}
