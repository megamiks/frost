#include <fstream>
#include <iostream>
#include <string> 
#include <sstream>
#include <cmath>

#define xlength 1
#define ylength 1
#define zlength 1
#define tau 0.01
#define nx 10
#define ny 10
#define nz 10

const double hx = xlength / (double(nx));
const double hy = ylength / (double(ny));
const double hz = zlength / (double(nz));
# define M_PI           3.14159265358979323846  /* pi */

/**
 * n - число уравнений (строк матрицы)
 * b - диагональ, лежащая над главной (нумеруется: [0;n-2])
 * c - главная диагональ матрицы A (нумеруется: [0;n-1])
 * a - диагональ, лежащая под главной (нумеруется: [1;n-1])
 * f - правая часть (столбец)
 * x - решение, массив x будет содержать ответ
 */
void solveMatrix(int n, double* a, double* c, double* b, double* f, double* x)
{
	double m;
	for (int i = 1; i < n; i++)
	{
		m = a[i] / c[i - 1];
		c[i] = c[i] - m * b[i - 1];
		f[i] = f[i] - m * f[i - 1];
	}

	x[n - 1] = f[n - 1] / c[n - 1];

	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = (f[i] - b[i] * x[i + 1]) / c[i];
	}
}


double f1(double x, double y, double z, double t) {
	return (3 * pow(M_PI, 2) - 1) * exp(-t) * sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z);
}

template <typename T>
std::string toString(T val)
{
	std::ostringstream oss;
	oss << val;
	return oss.str();
}

int count = 100;

int main() {


	double*** Tn;
	Tn= new double**[nx];
	for (int i = 0;i < nx;i++) 
		Tn[i] = new double* [ny];
	for (int i = 0;i < nx;i++) 
		for (int j = 0;j < ny;j++)
			Tn[i][j] = new double[nz];
	for (int i = 0;i < nx;i++)
		for (int j = 0;j < ny;j++)
			for (int k = 0;k < nz;k++)
				Tn[i][j][k] = 0;

	double*** Tne;
	Tne = new double** [nx];
	for (int i = 0;i < nx;i++)
		Tne[i] = new double* [ny];
	for (int i = 0;i < nx;i++)
		for (int j = 0;j < ny;j++)
			Tne[i][j] = new double[nz];
	for (int i = 0;i < nx;i++)
		for (int j = 0;j < ny;j++)
			for (int k = 0;k < nz;k++)
				Tne[i][j][k] = 0;

	double*** Tnew;
	Tnew = new double** [nx];
	for (int i = 0;i < nx;i++)
		Tnew[i] = new double* [ny];
	for (int i = 0;i < nx;i++)
		for (int j = 0;j < ny;j++)
			Tnew[i][j] = new double[nz];
	for (int i = 0;i < nx;i++)
		for (int j = 0;j < ny;j++)
			for (int k = 0;k < nz;k++)
				Tnew[i][j][k] = 0;
	int time = 1;
	std::ofstream fp;

	while (count > time) {
		if (time < 50)
			for (int i = 0;i < nz;i++)
				Tn[5][5][i] = 10.0;
			//предиктор по i
			double* Ti = new double[nx];
			double* am1 = new double[nx];
			double* bm1 = new double[nx];
			double* cm1 = new double[nx];
			double* solution1 = new double[nx];
			for(int j=0;j<ny;j++)
				for (int k = 0;k < nz;k++) {
					for (int i = 1;i < nx-1;i++) {
						Ti[i] = (Tn[i][j][k] / (3 * tau) + 5 * (Tn[i - 1][j][k] - 2 * Tn[i][j][k] + Tn[i + 1][j][k]) / (6 * hx * hx));
						am1[i] = (-1 / (6 * hx * hx));
						bm1[i] = (-1 / (6 * hx * hx));
						cm1[i] = (1 / (3 * tau) + 1 / (3 * hx * hx));
					}
					am1[0] = 0;
					bm1[nx - 1] = 0;
					solveMatrix(nx,am1, cm1, bm1, Ti, solution1);
					for (int i = 0;i < nx;i++)
						Tne[i][j][k] = (solution1[i] + 5 * Tn[i][j][k]) / 6;
				}
			delete [] Ti;
			delete [] am1;
			delete [] bm1;
			delete [] cm1;
			delete [] solution1;
			//предиктор по j
			double* Tj = new double[ny];
			double* am2 = new double[ny];
			double* bm2 = new double[ny];
			double* cm2 = new double[ny];
			double* solution2 = new double[ny];
			for (int i = 0;i < ny;i++)
				for (int k = 0;k < nz;k++) {
					for (int j = 1;j < ny-1;j++) {
						Tj[j] = ((6 * Tne[i][j][k] - 4 * Tn[i][j][k]) / (3 * tau) + 2 * (Tn[i][j - 1][k] - 2 * Tn[i][j][k] + Tn[i][j + 1][k]) / (3 * hy * hy));
						am2[j] = (-1 / (3 * hy * hy));
						bm2[j]=(-1 / (3 * hy * hy));
						cm2[j]=(2 / (3 * tau) + 2 / (3 * hy * hy));
					}
					am2[0] = 0;
					bm2[ny - 1] = 0;
					solveMatrix(ny, am2, cm2, bm2, Tj, solution2);
					for (int j = 0;j < ny;j++)
						Tne[i][j][k] = (solution2[j] + 2 * Tn[i][j][k]) / 3;
				}
			delete[] Tj;
			delete[] am2;
			delete[] bm2;
			delete[] cm2;
			delete[] solution2;
			//предиктор по k
			double* Tk = new double[nz];
			double* am3 = new double[nz];
			double* bm3 = new double[nz];
			double* cm3 = new double[nz];
			double* solution3 = new double[nz];
			for(int i=0;i<nx;i++)
				for (int j = 0;j < ny;j++) {
					for (int k = 1;k < nz-1;k++) {
						Tk[k] = ((2 * Tne[i][j][k] - Tn[i][j][k]) / tau + (Tn[i][j][k - 1] - 2 * Tn[i][j][k] + Tn[i][j][k + 1]) / (2 * hz * hz));
						am3[k] = (-1 / (2 * hz * hz));
						bm3[k] = (-1 / (2 * hz * hz));
						cm3[k] = (1 / tau + 1 / (hz * hz));
					}
					am3[0] = 0;
					bm3[nz - 1] = 0;
					solveMatrix(nz, am3, cm3, bm3, Tk, solution3);
					for (int k = 0;k < nz;k++)
						Tne[i][j][k] = (solution3[k] + Tn[i][j][k]) * 0.5;
				}
			delete[] Tk;
			delete[] am3;
			delete[] bm3;
			delete[] cm3;
			delete[] solution3;
			//корректор
			for (int i = 1;i < nx-1;i++)
				for (int j = 1;j < ny-1;j++)
					for (int k = 1;k < nz-1;k++)
						Tnew[i][j][k] = Tn[i][j][k] + tau * ((Tne[i - 1][j][k] - 2 * Tne[i][j][k] + Tne[i + 1][j][k]) / (hx * hx)
							+ (Tne[i][j - 1][k] - 2 * Tne[i][j][k] + Tne[i][j + 1][k]) / (hy * hy)
							+ (Tne[i][j][k - 1] - 2 * Tne[i][j][k] + Tne[i][j][k + 1]) / (hz * hz) + (f1(i * hx, j * hy, k * hz, tau * time) + f1(i * hx, j * hy, k * hz, tau * (time + 1))) * 0.5);

			for (int i = 0;i < nx;i++)
				for (int j = 0;j < ny;j++)
					for (int k = 0;k < nz;k++)
						Tn[i][j][k] = Tnew[i][j][k];
			
			//char* ttostr = itoa(time);
			std::string str = toString(time);
			std::string out = "out" + str + ".vtk";
			fp.open(out);
			fp << "# vtk DataFile Version 2.0\n";
			fp << "Mesh data from BSM\n";
			fp << "ASCII\n";
			fp << "DATASET STRUCTURED_GRID\n";
			fp << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
			fp << "POINTS " << nx * ny * nz << " float\n";
			for (int i = 0;i < nx;i++)
				for (int j = 0;j < ny;j++)
					for (int k = 0;k < nz;k++)
						fp << i * hx << " " << j * hy << " " << k * hz << "\n";
			fp << "POINT_DATA " << nx * ny * nz << "\n";
			fp << "SCALARS T float 1\n";
			fp << "LOOKUP_TABLE default\n";
			for (int i = 0;i < nx;i++)
				for (int j = 0;j < ny;j++)
					for (int k = 0;k < nz;k++)
						fp << Tn[i][j][k] << "\n";
			fp.close();
			time++;
	}
	for (int i = 0;i < nx;i++)
		for (int j = 0;j < ny;j++)
			delete[]Tn[i][j];
	for (int i = 0;i < nx;i++)
			delete[]Tn[i];
	delete[]Tn;

	for (int i = 0;i < nx;i++)
		for (int j = 0;j < ny;j++)
			delete[]Tne[i][j];
	for (int i = 0;i < nx;i++)
		delete[]Tne[i];
	delete[]Tne;

	for (int i = 0;i < nx;i++)
		for (int j = 0;j < ny;j++)
			delete[]Tnew[i][j];
	for (int i = 0;i < nx;i++)
		delete[]Tnew[i];
	delete[]Tnew;

	return 0;
}
