#include <iostream>
#include <fstream>
#include <vector>
#include "Mesh.h"
#include "Points.h"
#include "math.h"
#include "Trianinfo.h"
#include "Products.h"
#include <complex>
#include "EFIEop.h"

#define COMPLEX complex<double>

#define PI           double(3.14159265358979323846)  /* pi */
#define I             COMPLEX (0,1)
#define eps         double(2.2204e-16)
#define lam         double(1) // lambda equals to 1m


using namespace std;

namespace EFIE{

 double* rhoE = new double[3];
double* rho0NE= new double[3];
double* rho1E = new double[3];
double* rho2E = new double[3];
double* rho3E = new double[3];
double* l1E = new double[3];
double* l2E = new double[3];
double* l3E = new double[3];
double* RWGfE = new double[3];
double* IkonE = new double[3];
double* POME = new double[3];
double CONST1E;
double CONST2E;
double IscaE;
double LOGplE;
double LOGmnE;
double ATANplE;
double ATANmnE;
double KE;
double* KKE = new double[3];
double* npom1E = new double[3];
double* npom2E = new double[3];
double* npom3E = new double[3];
double* POM1E = new double[3];
double* POM2E = new double[3];

vector<double> fielpoinE(3);
vector<double> sourcepoinE(3);
double* u1E= new double[3];
double* u2E = new double[3];
double* u3E = new double[3];
vector<double*> PnulvecE(3);
vector<double*> UE(3);

vector<double> LplE(3);
vector<double> LmnE(3);
vector<double> PnulE(3);
vector<double> PplE(3);
vector<double> PmnE(3);
vector<double> DE(3);
vector<double> RnulE(3);
vector<double> RplE(3);
vector<double> RmnE(3);

double* rhotmpE = new double[3];
vector<COMPLEX> alok1E(9), alok2E(9);

double* vecpomE = new double[3];
double *p1nulvecE = new double[3];
double *p2nulvecE = new double[3];
double *p3nulvecE = new double[3];
double* pomVECE = new double[3];
double* pomocvFE = new double[3];
double* pomNvecE = new double[3];
double* pomVec1E = new double[3];
double* pomVec2E = new double[3];

double *IvecE = new double[3];




bool is_regular(const Trianinfo &trianF, const Trianinfo &trianS)
{
	vector<double> cp1 = trianF.getcp();
	vector<double> cp2 = trianS.getcp();
	double len = trianF.getLen() + trianS.getLen();
	if (norm2(&cp1[0], &cp2[0]) < 0.25 * len * len) return false;
	return true;
}



void singularityEFIE(double &det1, double &AR1, double &AR2, double* p1, double* p2, double* p3,
	double* q1, double* q2, double* q3, vector<double> &nvec2,
	vector<vector<double>> &PointsS, vector<double> &WeightsS, int* n1, int* n2,
	Mesh &mesh, COMPLEX k, COMPLEX eta, vector<COMPLEX> &alok2) {


	KE = dot(&nvec2[0], &q1[0]);

	multconst(npom1E, &nvec2[0], KE);

	subtract(rho1E, &q1[0], npom1E);

	KE = dot(&nvec2[0], &q2[0]);
	multconst(npom2E, &nvec2[0], KE);
	subtract(rho2E, &q2[0], npom2E);

	KE = dot(&nvec2[0], &q3[0]);

	multconst(npom3E, &nvec2[0], KE);

	subtract(rho3E, &q3[0], npom3E);

	subtract(KKE, rho2E, rho1E);
	double norm1 = norm(KKE);

	subtract(rhotmpE, rho2E, rho1E);

	dividecon(l1E, rhotmpE, norm1);

	subtract(rhotmpE, rho3E, rho2E);
	norm1 = norm(rhotmpE);

	subtract(rhotmpE, rho3E, rho2E);
	dividecon(l2E, rhotmpE, norm1);

	subtract(rhotmpE, rho1E, rho3E);
	norm1 = norm(rhotmpE);
	dividecon(l3E, rhotmpE, norm1);



	cross(u1E, l1E, &nvec2[0]);
	cross(u2E, l2E, &nvec2[0]);
	cross(u3E, l3E, &nvec2[0]);

	UE[0] = u1E;
	UE[1] = u2E;
	UE[2] = u3E;


	for (int nf = 0; nf != 12; ++nf) {

		double N1 = PointsS[nf][0];
		double N2 = PointsS[nf][1];
		double wf = WeightsS[nf];

		double N0 = 1.0 - N1 - N2;

		fielpoinE[0] = *p1*N1 + *p2*N2 + *p3*N0;
		fielpoinE[1] = *(p1 + 1)*N1 + *(p2 + 1)*N2 + *(p3 + 1)*N0;
		fielpoinE[2] = *(p1 + 2)*N1 + *(p2 + 2)*N2 + *(p3 + 2)*N0;

		KE = dot(&nvec2[0], &fielpoinE[0]);

		multconst(vecpomE, &nvec2[0], KE);

		subtract(rhoE, &fielpoinE[0], vecpomE);

		subtract(rhotmpE, rho2E, rhoE);
		double l1pl = dot(rhotmpE, l1E);
		subtract(rhotmpE, rho1E, rhoE);
		double l1mn = dot(rhotmpE, l1E);

		subtract(rhotmpE, rho3E, rhoE);
		double l2pl = dot(rhotmpE, l2E);
		subtract(rhotmpE, rho2E, rhoE);
		double l2mn = dot(rhotmpE, l2E);

		subtract(rhotmpE, rho1E, rhoE);
		double l3pl = dot(rhotmpE, l3E);
		subtract(rhotmpE, rho3E, rhoE);
		double l3mn = dot(rhotmpE, l3E);

		LplE[0] = l1pl;
		LplE[1] = l2pl;
		LplE[2] = l3pl;


		LmnE[0] = l1mn;
		LmnE[1] = l2mn;
		LmnE[2] = l3mn;

		subtract(rhotmpE, rho2E, rhoE);
		double p1nul = abs(dot(rhotmpE, u1E));
		subtract(rhotmpE, rho3E, rhoE);
		double p2nul = abs(dot(rhotmpE, u2E));
		subtract(rhotmpE, rho1E, rhoE);
		double p3nul = abs(dot(rhotmpE, u3E));

		PnulE[0] = p1nul;
		PnulE[1] = p2nul;
		PnulE[2] = p3nul;

		subtract(rhotmpE, rho2E, rhoE);
		double p1pl = norm(rhotmpE);
		subtract(rhotmpE, rho1E, rhoE);
		double p1mn = norm(rhotmpE);

		subtract(rhotmpE, rho3E, rhoE);
		double p2pl = norm(rhotmpE);
		subtract(rhotmpE, rho2E, rhoE);
		double p2mn = norm(rhotmpE);

		subtract(rhotmpE, rho1E, rhoE);
		double p3pl = norm(rhotmpE);
		subtract(rhotmpE, rho3E, rhoE);
		double p3mn = norm(rhotmpE);

		PplE[0] = p1pl;
		PplE[1] = p2pl;
		PplE[2] = p3pl;
		PmnE[0] = p1mn;
		PmnE[1] = p2mn;
		PmnE[2] = p3mn;

		if (p1nul < 50 * eps) {

			p1nulvecE[0] = 0.0;
			p1nulvecE[1] = 0.0;
			p1nulvecE[2] = 0.0;
		}

		else
		{
			subtract(POM1E, rho2E, rhoE);
			multconst(POM2E, l1E, l1pl);
			subtract(POME, POM1E, POM2E);

			dividecon(p1nulvecE, POME, p1nul);

		}

		if (p2nul < 50 * eps) {

			p2nulvecE[0] = 0.0;
			p2nulvecE[1] = 0.0;
			p2nulvecE[2] = 0.0;
		}

		else
		{
			subtract(POM1E, rho3E, rhoE);
			multconst(POM2E, l2E, l2pl);
			subtract(POME, POM1E, POM2E);

			dividecon(p2nulvecE, POME, p2nul);

		}

		if (p3nul < 50 * eps)
		{
			p3nulvecE[0] = 0.0;
			p3nulvecE[1] = 0.0;
			p3nulvecE[2] = 0.0;
		}
		else
		{
			subtract(POM1E, rho1E, rhoE);
			multconst(POM2E, l3E, l3pl);
			subtract(POME, POM1E, POM2E);
			dividecon(p3nulvecE, POME, p3nul);

		}


		PnulvecE[0] = p1nulvecE;
		PnulvecE[1] = p2nulvecE;
		PnulvecE[2] = p3nulvecE;


		subtract(rhotmpE, &fielpoinE[0], q1);
		double d1 = dot(&nvec2[0], rhotmpE);
		subtract(rhotmpE, &fielpoinE[0], q2);
		double d2 = dot(&nvec2[0], rhotmpE);
		subtract(rhotmpE, &fielpoinE[0], q3);
		double d3 = dot(&nvec2[0], rhotmpE);

		DE[0] = d1;
		DE[1] = d2;
		DE[2] = d3;

		double R1nul = sqrt(pow(p1nul, 2) + pow(d1, 2));
		double R2nul = sqrt(pow(p2nul, 2) + pow(d2, 2));
		double R3nul = sqrt(pow(p3nul, 2) + pow(d3, 2));

		RnulE[0] = R1nul;
		RnulE[1] = R2nul;
		RnulE[2] = R3nul;

		double R1pl = sqrt(pow(p1pl, 2) + pow(d1, 2));
		double R1mn = sqrt(pow(p1mn, 2) + pow(d1, 2));

		double R2pl = sqrt(pow(p2pl, 2) + pow(d2, 2));
		double R2mn = sqrt(pow(p2mn, 2) + pow(d2, 2));

		double R3pl = sqrt(pow(p3pl, 2) + pow(d3, 2));
		double R3mn = sqrt(pow(p3mn, 2) + pow(d3, 2));

		RplE[0] = R1pl;
		RplE[1] = R2pl;
		RplE[2] = R3pl;

		RmnE[0] = R1mn;
		RmnE[1] = R2mn;
		RmnE[2] = R3mn;



		IvecE[0] = 0.0;
		IvecE[1] = 0.0;
		IvecE[2] = 0.0;

		IscaE = 0.0; // scalar part of singular integral

					///// CALCULATION OF INTEGRALSS
					// vec part
		for (unsigned int iv = 0; iv != 3; ++iv) {

			if ((RplE[iv] + LplE[iv]) < 50 * eps)
				LOGplE = 0.0;
			else
				LOGplE = log(RplE[iv] + LplE[iv]);


			if ((RmnE[iv] + LmnE[iv]) < 50 * eps)
				LOGmnE = 0.0;
			else
				LOGmnE = log(RmnE[iv] + LmnE[iv]);


			CONST1E = 0.5*(pow(RnulE[iv], 2)*(LOGplE - LOGmnE) + RplE[iv] * LplE[iv] - RmnE[iv] * LmnE[iv]);


			multconst(pomVECE, UE[iv], CONST1E);

			add(IvecE, IvecE, pomVECE);


		}



		// scalar part
		for (unsigned int is = 0; is != 3; ++is) {

			if ((RplE[is] + LplE[is]) < 50 * eps)
				LOGplE = 0.0;
			else
				LOGplE = log(RplE[is] + LplE[is]);

			if ((RmnE[is] + LmnE[is]) < 50 * eps)
				LOGmnE = 0.0;
			else
				LOGmnE = log(RmnE[is] + LmnE[is]);

			if ((pow(RnulE[is], 2) + abs(DE[is])*RplE[is]) < 50 * eps)
				ATANplE = 0.0;
			else
				ATANplE = atan((PnulE[is] * LplE[is]) / (pow(RnulE[is], 2) + abs(DE[is])*RplE[is]));

			if ((pow(RnulE[is], 2) + abs(DE[is])*RmnE[is]) < 50 * eps)
				ATANmnE= 0.0;
			else
				ATANmnE = atan((PnulE[is] * LmnE[is]) / (pow(RnulE[is], 2) + abs(DE[is])*RmnE[is]));


			IscaE = IscaE + dot(PnulvecE[is], UE[is])*(PnulE[is] * (LOGplE - LOGmnE) - abs(DE[is])*(ATANplE - ATANmnE));

		}


		// testing Galerkin
		for (unsigned int i = 0; i != 3; ++i) {


			//vector<double> p(3);
			double* p = mesh.getCoord(n1[i]);

			subtract(pomocvFE, &fielpoinE[0], p);
			double konstF = (1 / (2 * AR1));

			multconst(RWGfE, pomocvFE, konstF);



			for (unsigned int j = 0; j != 3; ++j) {

				//vector<double> q(3);
				double* q = mesh.getCoord(n2[j]);

				double K1 = dot(&nvec2[0], q);

				multconst(pomNvecE, &nvec2[0], K1);

				subtract(rho0NE, q, pomNvecE);

				double K2 = double(1 / (2 * AR2));


				subtract(rhotmpE, rhoE, rho0NE);
				multconst(pomVec1E, rhotmpE, IscaE);

				add(pomVec2E, IvecE, pomVec1E);

				multconst(IkonE, pomVec2E, K2);

				alok2[i * 3 + j] = alok2[i * 3 + j] + wf * det1*(I*k*eta*dot(RWGfE, IkonE) - ((I*eta) / k)*(double(1 / AR1))*(double(1 / AR2))*(IscaE));
			}
		}
	}
}








void assemble_system_matrixEFIE(vector<COMPLEX> &A, Mesh &mesh, vector<Trianinfo> &Triangles, Points &points, unsigned int Nt, unsigned int maxele, COMPLEX k, COMPLEX eta)


{

	vector<vector<double>> PointsNS = points.getPointsNS();

	vector<double> WeightsNS = points.getWeightsNS();

	vector<vector<double>> PointsS = points.getPointsS();

	vector<double> WeightsS = points.getWeightsS();

	// Inicijalizacija pomocnih promenljivih
	double* pomocvF = new double[3];
	double* pomocvS = new double[3];
	double* RWGf = new double[3];
	double* RWGs = new double[3];

	//#pragma omp parallel for

	for (int ele1 = 0; ele1 < Nt; ++ele1)
	{
		int* n1 = mesh.getNOvertex(ele1);
		//vector<double>  p1(3);
		double*p1 = mesh.getCoord(*n1);
		//vector<double>  p2(3);
		double* p2 = mesh.getCoord(*(n1 + 1));
		// vector<double>  p3(3);
		double* p3 = mesh.getCoord(*(n1 + 2));


		vector<int> rwg1 = mesh.getRWG(ele1);

		Trianinfo trianF = Triangles[ele1];
		double det1 = trianF.getDeter();
		double AR1 = det1 / 2.0;



		for (int ele2 = 0; ele2 != Nt; ++ele2) {

			//vector<double> RWGf;
			//vector<double> RWGs;

			fill(alok1E.begin(), alok1E.end(), COMPLEX(0));
			fill(alok2E.begin(), alok2E.end(), COMPLEX(0));

			int* n2 = mesh.getNOvertex(ele2);
			// vector<double>  q1(3);
			double* q1 = mesh.getCoord(*n2);
			// vector<double>  q2(3);
			double* q2 = mesh.getCoord(*(n2 + 1));
			// vector<double>  q3(3);
			double* q3 = mesh.getCoord(*(n2 + 2));

			vector<int> rwg2 = mesh.getRWG(ele2);
			Trianinfo trianS = Triangles[ele2];
			double det2 = trianS.getDeter();
			double AR2 = det2 / 2.0;
			vector<double> nvecS = trianS.getnorm();


			if (is_regular(trianF, trianS)) {

				/// computation of regular submatrix*******************************************
				for (int nf = 0; nf != 4; ++nf) {

					double N1 = PointsNS[nf][0];
					double N2 = PointsNS[nf][1];
					double N0 = 1.0 - N1 - N2;

					double wf = WeightsNS[nf];

					fielpoinE[0] = *p1*N1 + *p2*N2 + *p3*N0;
					fielpoinE[1] = *(p1 + 1)*N1 + *(p2 + 1)*N2 + *(p3 + 1)*N0;
					fielpoinE[2] = *(p1 + 2)*N1 + *(p2 + 2)*N2 + *(p3 + 2)*N0;

					for (int ns = 0; ns != 4; ++ns) {

						N1 = PointsNS[ns][0];
						N2 = PointsNS[ns][1];
						N0 = 1.0 - N1 - N2;
						double ws = WeightsNS[ns];

						sourcepoinE[0] = *q1*N1 + *q2*N2 + *q3*(N0);
						sourcepoinE[1] = *(q1 + 1)*N1 + *(q2 + 1)*N2 + *(q3 + 1)*(N0);
						sourcepoinE[2] = *(q1 + 2)*N1 + *(q2 + 2)*N2 + *(q3 + 2)*(N0);


						double R = norm(&fielpoinE[0], &sourcepoinE[0]);

						COMPLEX Green = exp(-I * k*R) / R;

						for (unsigned int i = 0; i != 3; ++i) {

							//vector<double> p(3);
							double* p = mesh.getCoord(n1[i]);

							subtract(pomocvFE, &fielpoinE[0], p);
							double konstF = (1 / (2 * AR1));

							multconst(RWGf, pomocvF, konstF);

							for (unsigned int j = 0; j != 3; ++j) {

								//vector<double> q(3);
								double* q = mesh.getCoord(n2[j]);

								subtract(pomocvS, &sourcepoinE[0], q);

								double konstS = (1 / (2 * AR2));

								multconst(RWGs, pomocvS, konstS);


								alok1E[i * 3 + j] = alok1E[i * 3 + j] + wf * ws*det1*det2*(I*k*eta*Green*dot(RWGf, RWGs) - ((I*eta) / k)*Green*(double(1 / AR1))*(double(1 / AR2)));


							}
						}
					}
				}
			}


			else {

				//computation of singular submatrix

				for (int nf = 0; nf != 12; ++nf) {
					double N1 = PointsS[nf][0];
					double N2 = PointsS[nf][1];
					double N0 = 1.0 - N1 - N2;

					double wf = WeightsS[nf];

					fielpoinE[0] = *p1 * N1 + *p2 * N2 + *p3 * N0;
					fielpoinE[1] = *(p1 + 1) * N1 + *(p2 + 1) * N2 + *(p3 + 1) * N0;
					fielpoinE[2] = *(p1 + 2) * N1 + *(p2 + 2) * N2 + *(p3 + 2) * N0;


					for (int ns = 0; ns != 12; ++ns) {

						N1 = PointsS[ns][0];
						N2 = PointsS[ns][1];
						N0 = 1.0 - N1 - N2;

						double ws = WeightsS[ns];

						sourcepoinE[0] = *q1*N1 + *q2*N2 + *q3*(N0);
						sourcepoinE[1] = *(q1 + 1)*N1 + *(q2 + 1)*N2 + *(q3 + 1)*(N0);
						sourcepoinE[2] = *(q1 + 2)*N1 + *(q2 + 2)*N2 + *(q3 + 2)*(N0);

						double R = norm(&fielpoinE[0], &sourcepoinE[0]);
						// const REAL invR = REAL(1) / R;

						COMPLEX GreenNS;

						if (R < eps)
							GreenNS = -I * k;
						else
							GreenNS = ((exp(-I * k*R) / R) - 1.0 / R);

						for (unsigned int i = 0; i != 3; ++i) {

							//vector<double> p(3);
							double* p = mesh.getCoord(n1[i]);

							subtract(pomocvF, &fielpoinE[0], p);
							double konstF = (1.0 / (2.0 * AR1));

							multconst(RWGf, pomocvF, konstF);


							for (unsigned int j = 0; j != 3; ++j) {

								//vector<double> q(3);
								double* q = mesh.getCoord(n2[j]);

								subtract(pomocvS, &sourcepoinE[0], q);

								double konstS = (1.0 / (2.0 * AR2));

								multconst(RWGs, pomocvS, konstS);

								alok1E[i * 3 + j] = alok1E[i * 3 + j] + wf * ws*det1*det2*((I*k*eta*GreenNS*dot(RWGf, RWGs) - ((I*eta) / k)*GreenNS*(double(1 / AR1))*(double(1 / AR2))));

							}
						}
					}
				}


				//CALL singularity i zbrojiti te dvije matrice
				//vector<COMPLEX> alok2(9);
				//fill(alok2.begin(), alok2.end(), COMPLEX(0));


				singularityEFIE(det1, AR1, AR2, p1, p2, p3, q1, q2, q3, nvecS, PointsS, WeightsS, n1, n2, mesh, k, eta, alok2E);

				for (unsigned int i = 0; i != 3; ++i) {
					for (unsigned int j = 0; j != 3; ++j) {

						alok1E[i * 3 + j] = alok1E[i * 3 + j] + alok2E[i * 3 + j];
					}
				}
			}

			for (unsigned int i1 = 0; i1 != 3; ++i1) {
				const int tmp1 = rwg1[i1];
				const int s1 = (tmp1 < 0) ? -1 : 1;
				const int j1 = int(tmp1 * s1 - 1);
				for (unsigned int i2 = 0; i2 != 3; ++i2) {
					const int tmp2 = rwg2[i2];
					const int s2 = (tmp2 < 0) ? -1 : 1;
					const int j2 = int(tmp2 * s2 - 1);

					//	#pragma omp critical
					A[j1*maxele + j2] += double(1 / (4 * PI))*double(s1 * s2) * alok1E[i1 * 3 + i2];
				}
			}
		}
	}


}

}//end of namespace
