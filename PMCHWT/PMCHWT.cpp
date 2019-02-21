// PMCHWT.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <fstream>
#include <vector>
#include "Mesh.h"
#include "Points.h"
#include "math.h"
#include "Trianinfo.h"
#include "Products.h"
#include <complex>
#include <ppl.h>
#include <omp.h>
#include "MFIEop.h"
#include "EFIEop.h"
#include <Eigen/Dense>
#include "mpi.h"

#define COMPLEX complex<double>

#define PI           double(3.14159265358979323846)  /* pi */
#define I             COMPLEX (0,1)
#define eta0          COMPLEX(119.9169832*PI) 
#define eps         double(2.2204e-16)
#define lam         double(1) // lambda equals to 1m
#define k0           COMPLEX ((2 * PI) / lam)

using namespace std;
using namespace Eigen;

typedef Matrix<COMPLEX, Dynamic, Dynamic> MatrixXCPL;

int main(int args, char *argv[]) {
	int my_rank, numprocs;



	// MPI initialization
	MPI_Init(&args, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	vector<int>sizes(3);
	vector<double> coord;
	vector<int> topol;
	vector<int> trian;

	double num1 = 0;
	int num2 = 0;
	int num3 = 0;
	if (my_rank == 0) {

		ifstream file1("coord.txt");
		ifstream file2("topol.txt");
		ifstream file3("trian.txt");

		while (file1 >> num1) {
			coord.emplace_back(num1);
		}

		while (file2 >> num2) {
			topol.emplace_back(num2);
		}

		while (file3 >> num3) {
			trian.emplace_back(num3);
		}

		int sizeC = coord.size();
		int sizeTO = topol.size();
		int sizeTR = trian.size();

		sizes = { sizeC,sizeTO,sizeTR };

	}

	MPI_Bcast(&sizes[0], sizes.size(), MPI_INT, 0, MPI_COMM_WORLD);

	coord.resize(sizes[0]);
	topol.resize(sizes[1]);
	trian.resize(sizes[2]);



	MPI_Bcast(&coord[0], coord.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&topol[0], topol.size(), MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&trian[0], trian.size(), MPI_INT, 0, MPI_COMM_WORLD);


	for (int i = 0; i < topol.size(); ++i) {
		topol[i] = topol[i] - 1;
	}


	Mesh mesh(coord, topol, trian);

	unsigned int Nt = (topol.size()) / 3;  //number of triangles

	vector<Trianinfo> Triangles;

	int maxele = 1;

	for (int i = 0; i < trian.size(); ++i) {

		if (abs(trian[i]) >= maxele)
			maxele = abs(trian[i]);
	}


	for (int i = 0; i < Nt; ++i) {

		int* vertind = mesh.getNOvertex(i);

		double* vertcoord1 = mesh.getCoord(*vertind);
		double* vertcoord2 = mesh.getCoord(*(vertind + 1));
		double* vertcoord3 = mesh.getCoord(*(vertind + 2));

		Trianinfo trian(vertcoord1, vertcoord2, vertcoord3);
		Triangles.push_back(trian);

	}

	//DEFINICIJE PARAMETARA ZA DIELEKTRIK

	COMPLEX epsr2 = -2.0-3.0*I;
	COMPLEX mur2 = 1.0;

	COMPLEX k2 = k0*sqrt(epsr2*mur2);
	COMPLEX eta2 = eta0*sqrt(mur2 / epsr2);

	int SIZE = maxele*maxele;

	vector<COMPLEX> A1Eg(SIZE);
	vector<COMPLEX> A2Eg(SIZE); //globalne matrice
	vector<COMPLEX> A1Mg(SIZE);
	vector<COMPLEX> A2Mg(SIZE);

	vector<COMPLEX> A1El(SIZE);
	vector<COMPLEX> A2El(SIZE); //lokalne matrice
	vector<COMPLEX> A1Ml(SIZE);
	vector<COMPLEX> A2Ml(SIZE);

	fill(A1Eg.begin(), A1Eg.end(), COMPLEX(0));
	fill(A2Eg.begin(), A2Eg.end(), COMPLEX(0));
	fill(A1Mg.begin(), A1Mg.end(), COMPLEX(0));
	fill(A2Mg.begin(), A2Mg.end(), COMPLEX(0));

	fill(A1El.begin(), A1El.end(), COMPLEX(0));
	fill(A2El.begin(), A2El.end(), COMPLEX(0));
	fill(A1Ml.begin(), A1Ml.end(), COMPLEX(0));
	fill(A2Ml.begin(), A2Ml.end(), COMPLEX(0));

	MatrixXCPL A11(maxele,maxele);
	MatrixXCPL A12(maxele,maxele);
	MatrixXCPL A21(maxele,maxele);
	MatrixXCPL A22(maxele,maxele);

	MatrixXCPL A(2*maxele, 2*maxele);

	Points points;
	cout << "Calculating .." << endl;

	auto begin = std::chrono::high_resolution_clock::now();

	assemble_system_matrixEFIE(A1E, mesh, Triangles, points, Nt, maxele, k0, eta0);
	assemble_system_matrixEFIE(A2E, mesh, Triangles, points, Nt, maxele, k2, eta2);

	assemble_system_matrixMFIE(A1M, mesh, Triangles, points, Nt, maxele, k0);
	assemble_system_matrixMFIE(A2M, mesh, Triangles, points, Nt, maxele, k2);

	auto end = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < maxele; ++i) {
		for (int j = 0; j < maxele; ++j) {

			A11(i,j) = A1E[i*maxele + j] + A2E[i*maxele + j];
			A12(i,j) = -A1M[i*maxele + j] - A2M[i*maxele + j];
			A21(i,j) = A1M[i*maxele + j] + A2M[i*maxele + j];
			A22(i,j) = (COMPLEX(1)/ pow(eta0,2))*A1E[i*maxele + j] + (COMPLEX(1) / pow(eta2, 2))*A2E[i*maxele + j];
		}
	}

	A.block(0, 0, maxele, maxele) = A11;
	A.block(0, maxele, maxele, maxele) = eta0*A12;
	A.block(maxele, 0, maxele, maxele) = eta0*A21;
	A.block(maxele, maxele, maxele, maxele) =pow(eta0,2)*A22;

	cout << A(2, 0) << endl;

	std::cout << 1.0*std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9 << "s" << std::endl;

	MPI_Finalize();

    return 0;
}

