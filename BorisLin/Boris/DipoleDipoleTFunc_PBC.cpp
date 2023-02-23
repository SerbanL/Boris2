#include "stdafx.h"
#include "DipoleDipoleTFunc.h"

//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
bool DipoleDipoleTFunc::CalcDiagTens2D_PBC(
	VEC<DBL3> &Ddiag, INT3 N, DBL3 h, 
	bool include_self_demag, 
	int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but here we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	Ddiag.set(DBL3());

	//1. D22(n, m, k) = D11(m, n, k) and D33(n, m, k) = D11(k, m, n)
	//2. D11(n, m, k) = D11(|n|, |m|, |k|)

#pragma omp parallel for
	for (int j0 = 0; j0 < (y_images ? N.y : N.y / 2); j0++) {
		for (int i0 = 0; i0 < (x_images ? N.x : N.x / 2); i0++) {

			DBL3 val = DBL3();

			for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
				for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
					for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

						int i = i_img * N.x + i0;
						int j = j_img * N.y + j0;
						int k = k_img;

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
						val += DBL3(3 * r_dir.x*r_dir.x - 1.0, 3 * r_dir.y*r_dir.y - 1.0, 3 * r_dir.z*r_dir.z - 1.0) * c;
					}
				}
			}

			//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
			Ddiag[INT3(i0, j0, 0)] = val;

			if (!x_images) {

				Ddiag[INT3((N.x - i0) % N.x, j0, 0)] = val;

				if (!y_images) {

					Ddiag[INT3((N.x - i0) % N.x, (N.y - j0) % N.y, 0)] = val;
				}
			}

			if (!y_images) {

				Ddiag[INT3(i0, (N.y - j0) % N.y, 0)] = val;
			}
		}
	}

	return true;
}

//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
bool DipoleDipoleTFunc::CalcOffDiagTens2D_PBC(std::vector<double> &Dodiag, INT3 N, DBL3 h, int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but here we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	std::fill(Dodiag.begin(), Dodiag.end(), 0.0);

	//1. D13(n, m, k) = D12(n, k, m) and D23(n, m, k) = D12(m, k, n)
	//2. D12(n, m, k) = sign(n) * sign(m) D12(|n|, |m|, |k|)
	//thus:
	//3. D13(n, m, k) = sign(n) * sign(k) D12(|n|, |k|, |m|)
	//4. D23(n, m, k) = sign(m) * sign(k) D12(|m|, |k|, |n|)

#pragma omp parallel for
	for (int j0 = 0; j0 < (y_images ? N.y : N.y / 2); j0++) {
		for (int i0 = 0; i0 < (x_images ? N.x : N.x / 2); i0++) {

			double val = 0.0;

			for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
				for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
					for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

						int i = i_img * N.x + i0;
						int j = j_img * N.y + j0;
						int k = k_img;

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
						val += 3 * r_dir.x*r_dir.y * c;
					}
				}
			}

			//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
			Dodiag[i0 + j0 * N.x] = val;

			if (!x_images) {

				Dodiag[((N.x - i0) % N.x) + j0 * N.x] = -val;

				if (!y_images) {

					Dodiag[((N.x - i0) % N.x) + ((N.y - j0) % N.y) * N.x] = val;
				}
			}

			if (!y_images) {

				Dodiag[i0 + ((N.y - j0) % N.y) * N.x] = -val;
			}
		}
	}

	return true;
}

//Compute in Ddiag the diagonal tensor elements (Dxx, Dyy, Dzz) which has sizes given by N.
bool DipoleDipoleTFunc::CalcDiagTens3D_PBC(VEC<DBL3> &Ddiag, INT3 N, DBL3 h, bool include_self_demag, int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but here we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	Ddiag.set(DBL3());

	//1. D22(n, m, k) = D11(m, n, k) and D33(n, m, k) = D11(k, m, n)
	//2. D11(n, m, k) = D11(|n|, |m|, |k|)

#pragma omp parallel for
	for (int j0 = 0; j0 < (y_images ? N.y : N.y / 2); j0++) {
		for (int k0 = 0; k0 < (z_images ? N.z : N.z / 2); k0++) {
			for (int i0 = 0; i0 < (x_images ? N.x : N.x / 2); i0++) {

				DBL3 val = DBL3();

				for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
					for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
						for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

							int i = i_img * N.x + i0;
							int j = j_img * N.y + j0;
							int k = k_img * N.z + k0;

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
							val += DBL3(3 * r_dir.x*r_dir.x - 1.0, 3 * r_dir.y*r_dir.y - 1.0, 3 * r_dir.z*r_dir.z - 1.0) * c;
						}
					}
				}

				//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
				Ddiag[INT3(i0, j0, k0)] = val;

				if (!x_images) {

					Ddiag[INT3((N.x - i0) % N.x, j0, k0)] = val;

					if (!y_images) {

						Ddiag[INT3((N.x - i0) % N.x, (N.y - j0) % N.y, k0)] = val;

						if (!z_images)
							Ddiag[INT3((N.x - i0) % N.x, (N.y - j0) % N.y, (N.z - k0) % N.z)] = val;
					}

					if (!z_images) {

						Ddiag[INT3((N.x - i0) % N.x, j0, (N.z - k0) % N.z)] = val;
					}
				}

				if (!y_images) {

					Ddiag[INT3(i0, (N.y - j0) % N.y, k0)] = val;

					if (!z_images)
						Ddiag[INT3(i0, (N.y - j0) % N.y, (N.z - k0) % N.z)] = val;
				}

				if (!z_images) {

					Ddiag[INT3(i0, j0, (N.z - k0) % N.z)] = val;
				}
			}
		}
	}

	return true;
}

//Compute in Dodiag the off-diagonal tensor elements (Dxy, Dxz, Dyz) which has sizes given by N. 
bool DipoleDipoleTFunc::CalcOffDiagTens3D_PBC(VEC<DBL3> &Dodiag, INT3 N, DBL3 h, int x_images, int y_images, int z_images)
{
	//caller can have these negative (which means use inverse pbc for differential operators, but here we need them positive)
	x_images = abs(x_images);
	y_images = abs(y_images);
	z_images = abs(z_images);

	//zero the tensor first
	Dodiag.set(DBL3());

	//1. D13(n, m, k) = D12(n, k, m) and D23(n, m, k) = D12(m, k, n)
	//2. D12(n, m, k) = sign(n) * sign(m) D12(|n|, |m|, |k|)
	//thus:
	//3. D13(n, m, k) = sign(n) * sign(k) D12(|n|, |k|, |m|)
	//4. D23(n, m, k) = sign(m) * sign(k) D12(|m|, |k|, |n|)

#pragma omp parallel for
	for (int j0 = 0; j0 < (y_images ? N.y : N.y / 2); j0++) {
		for (int k0 = 0; k0 < (z_images ? N.z : N.z / 2); k0++) {
			for (int i0 = 0; i0 < (x_images ? N.x : N.x / 2); i0++) {

				DBL3 val = DBL3();

				for (int i_img = -x_images; i_img < x_images + 1; i_img++) {
					for (int j_img = -y_images; j_img < y_images + 1; j_img++) {
						for (int k_img = -z_images; k_img < z_images + 1; k_img++) {

							int i = i_img * N.x + i0;
							int j = j_img * N.y + j0;
							int k = k_img * N.z + k0;

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
							val += DBL3(3 * r_dir.x*r_dir.y, 3 * r_dir.x*r_dir.z, 3 * r_dir.y*r_dir.z) * c;
						}
					}
				}

				//For axes not using pbc we can use symmetries to fill remaining coefficients : the simple scheme below covers all possible cases.
				Dodiag[INT3(i0, j0, k0)] = val;

				if (!x_images) {

					Dodiag[INT3((N.x - i0) % N.x, j0, k0)] = val & DBL3(-1, -1, +1);

					if (!y_images) {

						Dodiag[INT3((N.x - i0) % N.x, (N.y - j0) % N.y, k0)] = val & DBL3(+1, -1, -1);

						if (!z_images)
							Dodiag[INT3((N.x - i0) % N.x, (N.y - j0) % N.y, (N.z - k0) % N.z)] = val;
					}

					if (!z_images) {

						Dodiag[INT3((N.x - i0) % N.x, j0, (N.z - k0) % N.z)] = val & DBL3(-1, +1, -1);
					}
				}

				if (!y_images) {

					Dodiag[INT3(i0, (N.y - j0) % N.y, k0)] = val & DBL3(-1, +1, -1);

					if (!z_images)
						Dodiag[INT3(i0, (N.y - j0) % N.y, (N.z - k0) % N.z)] = val & DBL3(-1, -1, +1);
				}

				if (!z_images) {

					Dodiag[INT3(i0, j0, (N.z - k0) % N.z)] = val & DBL3(+1, -1, -1);
				}
			}
		}
	}

	return true;
}