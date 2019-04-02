// PMCHWT.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "Mesh.h"
#include "Points.h"
#include "math.h"
#include "Trianinfo.h"
#include "Products.h"
#include <complex>
#include "MFIEop.h"
#include "EFIEop.h"
#include "excvecH.h"
#include "excvecE.h"
#include <time.h>
#include <Eigen/Dense>
#include <omp.h>

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
typedef Matrix<COMPLEX,Dynamic, 1>VectorXCPL;


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

	//program ce lokalno raditi na 9 procesa, 1 master i 8 slavea, svaka matrica ce se racunati sa dva procesa - tako je napisan i mora biti 
	// taj broj procesa

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
	
	clock_t start, end;
    double cpu_time_used;
    
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

    file1.clear();
	file1.seekg(0, ios::beg);
	
	file2.clear();
	file2.seekg(0, ios::beg);  // go to the beginning of files

	file3.clear();
	file3.seekg(0, ios::beg);

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

	COMPLEX epsr2 = -3.0-2.0*I;
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
		EFIE::assemble_system_matrixEFIE(A1El, mesh, Triangles, points, Nt, maxele, k0, eta0, my_rank);
		send_data(A1El, SIZE, numprocs, my_rank);
	}
	else if ((my_rank == 3) || (my_rank == 4)){
		EFIE::assemble_system_matrixEFIE(A2El, mesh, Triangles, points, Nt, maxele, k2, eta2, my_rank); 
		send_data(A2El, SIZE, numprocs, my_rank);
	}

	else if ((my_rank == 5) || (my_rank == 6)){
		MFIE::assemble_system_matrixMFIE(A1Ml, mesh, Triangles, points, Nt, maxele, k0, my_rank); 
		send_data(A1Ml, SIZE, numprocs, my_rank);
	}

	else if ((my_rank == 7) || (my_rank == 8)){
		MFIE::assemble_system_matrixMFIE(A2Ml, mesh, Triangles, points, Nt, maxele, k2, my_rank);
		send_data(A2Ml, SIZE, numprocs, my_rank);
	}
		

	if (my_rank == 0) {
	
	
		vector<COMPLEX> A1Eg(SIZE);
		vector<COMPLEX> A2Eg(SIZE); //globalne matrice
		vector<COMPLEX> A1Mg(SIZE);
		vector<COMPLEX> A2Mg(SIZE);
		
			vector<COMPLEX> H(maxele);
	        vector<COMPLEX> E(maxele);

		fill(A1Eg.begin(), A1Eg.end(), COMPLEX(0));
		fill(A2Eg.begin(), A2Eg.end(), COMPLEX(0));
		fill(A1Mg.begin(), A1Mg.end(), COMPLEX(0));
		fill(A2Mg.begin(), A2Mg.end(), COMPLEX(0));
		
			fill(H.begin(), H.end(), COMPLEX(0));
			fill(E.begin(), E.end(), COMPLEX(0));

		MatrixXCPL A11(maxele, maxele);
		MatrixXCPL A12(maxele, maxele);
		MatrixXCPL A21(maxele, maxele);
		MatrixXCPL A22(maxele, maxele);

		MatrixXCPL A(2 * maxele, 2 * maxele);
		VectorXCPL C(2*maxele);
		VectorXCPL B(2*maxele);
		
		 start = omp_get_wtime();

		receive_data(A1Eg, A2Eg, A1Mg, A2Mg, SIZE, numprocs);
		
		
		EFIE::excEFIE::assemble_exic_vector(E, mesh, Triangles, points, Nt, k0);
	    MFIE::excMFIE::assemble_exic_vector(H, mesh, Triangles, points, Nt, k0, eta0);
	    
	    end = omp_get_wtime();
	
		cpu_time_used = ((double) (end - start));
	
		cout << "System matrix assembling," << cpu_time_used << endl; 


		for (int i = 0; i < maxele; ++i) {
		
		 C(i)=E[i];
	     C(i+maxele)=eta0*H[i];
	     
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
		

	start = omp_get_wtime();
	
    B = A.partialPivLu().solve(C);//B = A.colPivHouseholderQr().solve(C);
     
	end = omp_get_wtime();
	

	
	cpu_time_used = ((double) (end - start));
	
	cout << "Linear system solving," << cpu_time_used << endl;
	
	
	B.tail(maxele)=eta0*B.tail(maxele);
	
	ofstream out1("solution.txt");
	ofstream out2("parameters.txt");
	
   if (out1.is_open())
 
  {
    
    for (int i = 0; i < 2*maxele; ++i){ 
    out1 << real(B(i));
    out1 <<',';
    out1 << imag(B(i));
    out1 << endl;
    }
    
    out1.close();
    
  }
  
   if (out2.is_open())
   {
   
   out2 << maxele;
   
   out2.close();
   
   }
  

	COMPLEX zbroj = B.sum();

	
    cout << zbroj<< endl ;

	}


	MPI_Finalize();

    return 0;
}

