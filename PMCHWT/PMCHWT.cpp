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


void send_data(vector<COMPLEX> &local_data, int SIZE, int numprocs, int my_rank) {

	MPI_Send(&local_data[0], SIZE, MPI_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD);
}

void receive_data(vector<COMPLEX> &A11, vector<COMPLEX> &A12, vector<COMPLEX> &A21, vector<COMPLEX> &A22, int SIZE, int numprocs) {

	int n = sqrt(SIZE);

	vector<COMPLEX> temp(SIZE);

	for (int rank = 1; rank < numprocs; ++rank) {

		MPI_Recv(&temp[0], SIZE, MPI_DOUBLE_COMPLEX, rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		for (int j1 = 0; j1 < n; ++j1) {


			for (int j2 = 0; j2 != n; ++j2) {

				if ((rank == 1) || (rank == 2))
				A11[j1*n + j2] += temp[j1*n + j2];

				if ((rank == 3) || (rank == 4))
				A12[j1*n + j2] += temp[j1*n + j2];

				if ((rank == 5) || (rank == 6))
				A21[j1*n + j2] += temp[j1*n + j2];

				if ((rank == 7) || (rank == 8))
				A22[j1*n + j2] += temp[j1*n + j2];


			}

		}

	}
}




int main(int args, char *argv[]) {

	int my_rank, numprocs;

	//program ce lokalno raditi na 9 procesa, 1 master i 8 slavea, svaka matrica ce se racunati sa dva procesa

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



	vector<COMPLEX> A1El(SIZE);
	vector<COMPLEX> A2El(SIZE); //lokalne matrice
	vector<COMPLEX> A1Ml(SIZE);
	vector<COMPLEX> A2Ml(SIZE);


	fill(A1El.begin(), A1El.end(), COMPLEX(0));
	fill(A2El.begin(), A2El.end(), COMPLEX(0));
	fill(A1Ml.begin(), A1Ml.end(), COMPLEX(0));
	fill(A2Ml.begin(), A2Ml.end(), COMPLEX(0));

	
	Points points;


    if((my_rank==1)|| (my_rank == 2)){
	assemble_system_matrixEFIE(A1El, mesh, Triangles, points, Nt, maxele, k0, eta0, my_rank);
	    send_data(A1El, SIZE, numprocs, my_rank);
    }
	else if ((my_rank == 3) || (my_rank == 4)){
		assemble_system_matrixEFIE(A2El, mesh, Triangles, points, Nt, maxele, k2, eta2, my_rank); 
		send_data(A2El, SIZE, numprocs, my_rank);
	}

	else if ((my_rank == 5) || (my_rank == 6)){
		assemble_system_matrixMFIE(A1Ml, mesh, Triangles, points, Nt, maxele, k0, my_rank); 
		send_data(A1Ml, SIZE, numprocs, my_rank);
	}

	else if ((my_rank == 7) || (my_rank == 8)){
	assemble_system_matrixMFIE(A2Ml, mesh, Triangles, points, Nt, maxele, k2, my_rank);
       send_data(A2Ml, SIZE, numprocs, my_rank);
}
		

	if (my_rank == 0) {

		vector<COMPLEX> A1Eg(SIZE);
		vector<COMPLEX> A2Eg(SIZE); //globalne matrice
		vector<COMPLEX> A1Mg(SIZE);
		vector<COMPLEX> A2Mg(SIZE);

		fill(A1Eg.begin(), A1Eg.end(), COMPLEX(0));
		fill(A2Eg.begin(), A2Eg.end(), COMPLEX(0));
		fill(A1Mg.begin(), A1Mg.end(), COMPLEX(0));
		fill(A2Mg.begin(), A2Mg.end(), COMPLEX(0));

		MatrixXCPL A11(maxele, maxele);
		MatrixXCPL A12(maxele, maxele);
		MatrixXCPL A21(maxele, maxele);
		MatrixXCPL A22(maxele, maxele);

		MatrixXCPL A(2 * maxele, 2 * maxele);

       auto begin = std::chrono::high_resolution_clock::now();

		receive_data(A1Eg,A2Eg,A1Mg,A2Mg,maxele, numprocs);

	   auto end = std::chrono::high_resolution_clock::now();

	for (int i = 0; i < maxele; ++i) {
		for (int j = 0; j < maxele; ++j) {

			A11(i,j) = A1Eg[i*maxele + j] + A2Eg[i*maxele + j];
			A12(i,j) = -A1Mg[i*maxele + j] - A2Mg[i*maxele + j];
			A21(i,j) = A1Mg[i*maxele + j] + A2Mg[i*maxele + j];
			A22(i,j) = (COMPLEX(1)/ pow(eta0,2))*A1Eg[i*maxele + j] + (COMPLEX(1) / pow(eta2, 2))*A2Eg[i*maxele + j];
		}
	}

	A.block(0, 0, maxele, maxele) = A11;
	A.block(0, maxele, maxele, maxele) = eta0*A12;
	A.block(maxele, 0, maxele, maxele) = eta0*A21;
	A.block(maxele, maxele, maxele, maxele) =pow(eta0,2)*A22;


cout << A(2, 0) << endl;

std::cout << 1.0*std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1e9 << "s" << std::endl;

}


	MPI_Finalize();

    return 0;
}

