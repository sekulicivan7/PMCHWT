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
#include "MFIEop.h"
#include "EFIEop.h"
#include "excvecH.h"
#include "excvecE.h"
#include <time.h>
#include <Eigen/Dense>

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

int main()
{

	ifstream file1("coord.txt");
	ifstream file2("topol.txt");
	ifstream file3("trian.txt");

	vector<double> coord;
	vector<int> topol;
	vector<int> trian;

	double num1 = 0;
	int num2 = 0;
	int num3 = 0;
	
	clock_t start, end;
    double cpu_time_used;

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

	COMPLEX epsr2 = -3.0 - 2.0*I;
	COMPLEX mur2 = 1.0;

	COMPLEX k2 = k0*sqrt(epsr2*mur2);
	COMPLEX eta2 = eta0*sqrt(mur2 / epsr2);

    
	vector<COMPLEX> A1E(maxele*maxele);
	vector<COMPLEX> A2E(maxele*maxele);
	vector<COMPLEX> A1M(maxele*maxele);
	vector<COMPLEX> A2M(maxele*maxele); 
	
	vector<COMPLEX> H(maxele);
	vector<COMPLEX> E(maxele);

	fill(A1E.begin(), A1E.end(), COMPLEX(0));
	fill(A2E.begin(), A2E.end(), COMPLEX(0));
	fill(A1M.begin(), A1M.end(), COMPLEX(0));
	fill(A2M.begin(), A2M.end(), COMPLEX(0));
	
	fill(H.begin(), H.end(), COMPLEX(0));
	fill(E.begin(), E.end(), COMPLEX(0));

	MatrixXCPL A11(maxele,maxele);
	MatrixXCPL A12(maxele,maxele);
	MatrixXCPL A21(maxele,maxele);
	MatrixXCPL A22(maxele,maxele);

	MatrixXCPL A(2*maxele, 2*maxele);
	VectorXCPL C(2*maxele);
	VectorXCPL B(2*maxele);

	Points points;
	cout << "Calculating .." << endl;
	
    start = clock();
    
	EFIE::assemble_system_matrixEFIE(A1E, mesh, Triangles, points, Nt, maxele, k0, eta0);
	EFIE::assemble_system_matrixEFIE(A2E, mesh, Triangles, points, Nt, maxele, k2, eta2);

	MFIE::assemble_system_matrixMFIE(A1M, mesh, Triangles, points, Nt, maxele, k0);
	MFIE::assemble_system_matrixMFIE(A2M, mesh, Triangles, points, Nt, maxele, k2);
	
	EFIE::excEFIE::assemble_exic_vector(E, mesh, Triangles, points, Nt, k0);
	MFIE::excMFIE::assemble_exic_vector(H, mesh, Triangles, points, Nt, k0, eta0);
	
	end = clock();
	
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	
	cout << "System matrix assembling," << cpu_time_used << endl; 

	for (int i = 0; i < maxele; ++i) {
	 
	 C(i)=E[i];
	 C(i+maxele)=eta0*H[i];
	
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
	
	start = clock();
	B = A.colPivHouseholderQr().solve(C);
	end = clock();
	
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	
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


    return 0;
}

