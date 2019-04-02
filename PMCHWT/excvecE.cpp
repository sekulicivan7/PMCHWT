#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include "Mesh.h"
#include "Points.h"
#include "math.h"
#include "Trianinfo.h"
#include "Products.h"
#include "excvecH.h"
#include <omp.h>

#define COMPLEX complex<double>

#define PI           double(3.14159265358979323846)  /* pi */
#define I             COMPLEX (0,1)
#define eps         double(2.2204e-16)
#define lam         double(1) // lambda equals to 1m


using namespace std;


namespace EFIE{

namespace excEFIE{


void assemble_exic_vector(vector<COMPLEX> &C, Mesh &mesh, vector<Trianinfo> &Triangles, Points &points, int Nt, COMPLEX k)
{
			
	
	double k_i[3];

	double e_i[3];

	const vector<vector<double>> PointsNS = points.getPointsNS();

	const vector<double> WeightsNS = points.getWeightsNS();


	k_i[0] = 0.0;
	k_i[1] = 0.0;
	k_i[2] = 1.0;

	e_i[0] = 1.0;
	e_i[1] = 0.0;
	e_i[2] = 0.0;
	
    #pragma omp parallel for
	for (unsigned int ele1 = 0; ele1 < Nt; ++ele1) {
	
	double RWGf[3];

	double pomocvF[3];

	vector<COMPLEX> alok1exvec(3);
	
	vector<double> fielpoin(3);
	
       const int* n1 = mesh.getNOvertex(ele1);
		
	   const double*p1 = mesh.getCoord(n1[0]);
		
	   const double* p2 = mesh.getCoord(n1[1]);
		
	   const double* p3 = mesh.getCoord(n1[2]);
		
		vector<int> rwg1 = mesh.getRWG(ele1);

		Trianinfo trianF = Triangles[ele1];

		vector<double> nvecF = trianF.getnorm();

		double det1 = trianF.getDeter();

		double AR1 = det1 / 2.0;
		
		fill(alok1exvec.begin(), alok1exvec.end(), COMPLEX(0));
		
		for (unsigned int nf = 0; nf != 7; ++nf) {

					double N1 = PointsNS[nf][0];
					double N2 = PointsNS[nf][1];
					double N0 = 1.0 - N1 - N2;

					double wf = WeightsNS[nf];

					fielpoin[0] = p1[0]*N1 + p2[0]*N2 + p3[0]*N0;
					fielpoin[1] = p1[1]*N1 + p2[1]*N2 + p3[1]*N0;
					fielpoin[2] = p1[2]*N1 + p2[2]*N2 + p3[2]*N0;

			for (unsigned int i = 0; i != 3; ++i) {

					double* p = mesh.getCoord(n1[i]);

					subtract(pomocvF, &fielpoin[0], p);
					
					double konstF = (1.0 / (2.0 * AR1));

					multconst(RWGf, pomocvF, konstF);
				
				    alok1exvec[i] = alok1exvec[i] + wf*det1*dot(RWGf, e_i)*exp(-I*k*dot(k_i, &fielpoin[0]));
			}
		}


		for (unsigned int i1 = 0; i1 != 3; ++i1) {
		
			const int tmp1 = rwg1[i1];
			
			const int s1 = (tmp1 < 0) ? -1 : 1;
			
			const int j1 = int(tmp1 * s1 - 1);
			
			#pragma omp critical
			C[j1] += double(s1) * alok1exvec[i1];
			
		}

	}
}

}//end of namespace

}//end of namespace
