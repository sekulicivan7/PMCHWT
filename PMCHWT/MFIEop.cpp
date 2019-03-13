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


#define COMPLEX complex<double>

#define PI           double(3.14159265358979323846)  /* pi */
#define I             COMPLEX (0,1)
#define eps         double(2.2204e-16)
#define lam         double(1) // lambda equals to 1m


using namespace std;

namespace MFIE{

double* rho = new double[3];
double* rho0N = new double[3];
double* rho1 = new double[3];
double* rho2 = new double[3];
double* rho3 = new double[3];
double* l1 = new double[3];
double* l2 = new double[3];
double* l3 = new double[3];
double* RWGf = new double[3];
COMPLEX* Ikon = new COMPLEX[3];// finish integral
double* POM = new double[3];
double CONST1;
double CONST2;
COMPLEX KONST1;
COMPLEX KONST2;
double LOGpl;
double LOGmn;
double ATANpl;
double ATANmn;
double Isca1R;
double Isca1R3;
double K;
double* KK = new double[3];
double* npom1 = new double[3];
double* npom2 = new double[3];
double* npom3 = new double[3];
double* POM1 = new double[3];
double* POM2 = new double[3];

vector<double> fielpoin(3);
vector<double> sourcepoin(3);
double* u1 = new double[3];
double* u2 = new double[3];
double* u3 = new double[3];
vector<double*> Pnulvec(3);
vector<double*> U(3);

vector<double> Lpl(3);
vector<double> Lmn(3);
vector<double> Pnul(3);
vector<double> Ppl(3);
vector<double> Pmn(3);
vector<double> D(3);
vector<double> Rnul(3);
vector<double> Rpl(3);
vector<double> Rmn(3);

double* rhotmp = new double[3];
vector<COMPLEX> alok1(9), alok2(9), alok3(9);

double* vecpom = new double[3];
double *p1nulvec = new double[3];
double *p2nulvec = new double[3];
double *p3nulvec = new double[3];
double* pomVEC = new double[3];
double* pomocvF = new double[3];
double* pomNvec = new double[3];
double* pomVec1 = new double[3];
COMPLEX* pomVec2 = new COMPLEX[3];
COMPLEX* pomVec3 = new COMPLEX[3];
COMPLEX* pomVec4 = new COMPLEX[3];
double* pomVec5 = new double[3];
double* pomVec6 = new double[3];

double *Ivec1R = new double[3];
double *Ivec1R3 = new double[3];
double* Rivec = new double[3];
COMPLEX* Ikon1R = new COMPLEX[3];
double* Ikon1R3 = new double[3];

COMPLEX* KROSS1 = new COMPLEX[3];
double* KROSS2 = new double[3];



bool is_regularM(const Trianinfo &trianF, const Trianinfo &trianS)
{
	vector<double> cp1 = trianF.getcp();
	vector<double> cp2 = trianS.getcp();
	double len = trianF.getLen() + trianS.getLen();
	if (norm2(&cp1[0], &cp2[0]) < 0.25 * len * len) return false;
	return true;
}

bool is_sameT(const Trianinfo &trianF, const Trianinfo &trianS)
{
	vector<double> cp1 = trianF.getcp();
	vector<double> cp2 = trianS.getcp();
	double distance = norm(&cp1[0], &cp2[0]);
	if (abs(distance) < 2 * eps) return true;
	return false;
}


void singularityMFIE(double &det1, double &AR1, double &AR2, double* p1, double* p2, double* p3,
	double* q1, double* q2, double* q3, vector<double> &nvec2, vector<double> &nvec1,
	vector<vector<double>> &PointsS, vector<double> &WeightsS, int* n1, int* n2,
	Mesh &mesh, COMPLEX k, vector<COMPLEX> &alok2) {


	K = dot(&nvec2[0], &q1[0]);

	multconst(npom1, &nvec2[0], K);

	subtract(rho1, &q1[0], npom1);

	K = dot(&nvec2[0], &q2[0]);

	multconst(npom2, &nvec2[0], K);

	subtract(rho2, &q2[0], npom2);

	K = dot(&nvec2[0], &q3[0]);

	multconst(npom3, &nvec2[0], K);

	subtract(rho3, &q3[0], npom3);

	subtract(KK, rho2, rho1);
	double norm1 = norm(KK);

	subtract(rhotmp, rho2, rho1);

	dividecon(l1, rhotmp, norm1);

	subtract(rhotmp, rho3, rho2);
	norm1 = norm(rhotmp);

	subtract(rhotmp, rho3, rho2);
	dividecon(l2, rhotmp, norm1);

	subtract(rhotmp, rho1, rho3);
	norm1 = norm(rhotmp);
	dividecon(l3, rhotmp, norm1);



	cross(u1, l1, &nvec2[0]);
	cross(u2, l2, &nvec2[0]);
	cross(u3, l3, &nvec2[0]);

	U[0] = u1;
	U[1] = u2;
	U[2] = u3;


	for (int nf = 0; nf != 12; ++nf) {

		double N1 = PointsS[nf][0];
		double N2 = PointsS[nf][1];
		double wf = WeightsS[nf];

		double N0 = 1.0 - N1 - N2;

		fielpoin[0] = *p1*N1 + *p2*N2 + *p3*N0;
		fielpoin[1] = *(p1 + 1)*N1 + *(p2 + 1)*N2 + *(p3 + 1)*N0;
		fielpoin[2] = *(p1 + 2)*N1 + *(p2 + 2)*N2 + *(p3 + 2)*N0;

		K = dot(&nvec2[0], &fielpoin[0]);

		multconst(vecpom, &nvec2[0], K);

		subtract(rho, &fielpoin[0], vecpom);

		subtract(rhotmp, rho2, rho);
		double l1pl = dot(rhotmp, l1);
		subtract(rhotmp, rho1, rho);
		double l1mn = dot(rhotmp, l1);

		subtract(rhotmp, rho3, rho);
		double l2pl = dot(rhotmp, l2);
		subtract(rhotmp, rho2, rho);
		double l2mn = dot(rhotmp, l2);

		subtract(rhotmp, rho1, rho);
		double l3pl = dot(rhotmp, l3);
		subtract(rhotmp, rho3, rho);
		double l3mn = dot(rhotmp, l3);

		Lpl[0] = l1pl;
		Lpl[1] = l2pl;
		Lpl[2] = l3pl;


		Lmn[0] = l1mn;
		Lmn[1] = l2mn;
		Lmn[2] = l3mn;

		subtract(rhotmp, rho2, rho);
		double p1nul = abs(dot(rhotmp, u1));
		subtract(rhotmp, rho3, rho);
		double p2nul = abs(dot(rhotmp, u2));
		subtract(rhotmp, rho1, rho);
		double p3nul = abs(dot(rhotmp, u3));

		Pnul[0] = p1nul;
		Pnul[1] = p2nul;
		Pnul[2] = p3nul;

		subtract(rhotmp, rho2, rho);
		double p1pl = norm(rhotmp);
		subtract(rhotmp, rho1, rho);
		double p1mn = norm(rhotmp);

		subtract(rhotmp, rho3, rho);
		double p2pl = norm(rhotmp);
		subtract(rhotmp, rho2, rho);
		double p2mn = norm(rhotmp);

		subtract(rhotmp, rho1, rho);
		double p3pl = norm(rhotmp);
		subtract(rhotmp, rho3, rho);
		double p3mn = norm(rhotmp);

		Ppl[0] = p1pl;
		Ppl[1] = p2pl;
		Ppl[2] = p3pl;
		Pmn[0] = p1mn;
		Pmn[1] = p2mn;
		Pmn[2] = p3mn;

		if (p1nul < 50 * eps) {

			p1nulvec[0] = 0.0;
			p1nulvec[1] = 0.0;
			p1nulvec[2] = 0.0;
		}

		else
		{
			subtract(POM1, rho2, rho);
			multconst(POM2, l1, l1pl);
			subtract(POM, POM1, POM2);

			dividecon(p1nulvec, POM, p1nul);

		}

		if (p2nul < 50 * eps) {

			p2nulvec[0] = 0.0;
			p2nulvec[1] = 0.0;
			p2nulvec[2] = 0.0;
		}

		else
		{
			subtract(POM1, rho3, rho);
			multconst(POM2, l2, l2pl);
			subtract(POM, POM1, POM2);

			dividecon(p2nulvec, POM, p2nul);

		}

		if (p3nul < 50 * eps)
		{
			p3nulvec[0] = 0.0;
			p3nulvec[1] = 0.0;
			p3nulvec[2] = 0.0;
		}
		else
		{
			subtract(POM1, rho1, rho);
			multconst(POM2, l3, l3pl);
			subtract(POM, POM1, POM2);
			dividecon(p3nulvec, POM, p3nul);

		}


		Pnulvec[0] = p1nulvec;
		Pnulvec[1] = p2nulvec;
		Pnulvec[2] = p3nulvec;


		subtract(rhotmp, &fielpoin[0], q1);
		double d1 = dot(&nvec2[0], rhotmp);
		subtract(rhotmp, &fielpoin[0], q2);
		double d2 = dot(&nvec2[0], rhotmp);
		subtract(rhotmp, &fielpoin[0], q3);
		double d3 = dot(&nvec2[0], rhotmp);

		D[0] = d1;
		D[1] = d2;
		D[2] = d3;

		double R1nul = sqrt(pow(p1nul, 2) + pow(d1, 2));
		double R2nul = sqrt(pow(p2nul, 2) + pow(d2, 2));
		double R3nul = sqrt(pow(p3nul, 2) + pow(d3, 2));

		Rnul[0] = R1nul;
		Rnul[1] = R2nul;
		Rnul[2] = R3nul;

		double R1pl = sqrt(pow(p1pl, 2) + pow(d1, 2));
		double R1mn = sqrt(pow(p1mn, 2) + pow(d1, 2));

		double R2pl = sqrt(pow(p2pl, 2) + pow(d2, 2));
		double R2mn = sqrt(pow(p2mn, 2) + pow(d2, 2));

		double R3pl = sqrt(pow(p3pl, 2) + pow(d3, 2));
		double R3mn = sqrt(pow(p3mn, 2) + pow(d3, 2));

		Rpl[0] = R1pl;
		Rpl[1] = R2pl;
		Rpl[2] = R3pl;

		Rmn[0] = R1mn;
		Rmn[1] = R2mn;
		Rmn[2] = R3mn;



		Ivec1R3[0] = Ivec1R3[1] = Ivec1R3[2] = Ivec1R[0] = Ivec1R[1] = Ivec1R[2] = 0.0;


		Isca1R = Isca1R3 = 0.0; // scalar part of singular integral

								///// CALCULATION OF INTEGRALSS
								// vec part
		for (unsigned int iv = 0; iv != 3; ++iv) {

			if ((Rpl[iv] + Lpl[iv]) < 50 * eps)
				LOGpl = 0.0;
			else
				LOGpl = log(Rpl[iv] + Lpl[iv]);


			if ((Rmn[iv] + Lmn[iv]) < 50 * eps)
				LOGmn = 0.0;
			else
				LOGmn = log(Rmn[iv] + Lmn[iv]);


			CONST1 = 0.5*(pow(Rnul[iv], 2)*(LOGpl - LOGmn) + Rpl[iv] * Lpl[iv] - Rmn[iv] * Lmn[iv]);


			multconst(pomVEC, U[iv], CONST1);

			add(Ivec1R, Ivec1R, pomVEC);

		}


		// scalar part
		for (unsigned int is = 0; is != 3; ++is) {

			if ((Rpl[is] + Lpl[is]) < 50 * eps)
				LOGpl = 0.0;
			else
				LOGpl = log(Rpl[is] + Lpl[is]);

			if ((Rmn[is] + Lmn[is]) < 50 * eps)
				LOGmn = 0.0;
			else
				LOGmn = log(Rmn[is] + Lmn[is]);

			if ((pow(Rnul[is], 2) + abs(D[is])*Rpl[is]) < 50 * eps)
				ATANpl = 0.0;
			else
				ATANpl = atan((Pnul[is] * Lpl[is]) / (pow(Rnul[is], 2) + abs(D[is])*Rpl[is]));

			if ((pow(Rnul[is], 2) + abs(D[is])*Rmn[is]) < 50 * eps)
				ATANmn = 0.0;
			else
				ATANmn = atan((Pnul[is] * Lmn[is]) / (pow(Rnul[is], 2) + abs(D[is])*Rmn[is]));


			Isca1R = Isca1R + dot(Pnulvec[is], U[is])*(Pnul[is] * (LOGpl - LOGmn) - abs(D[is])*(ATANpl - ATANmn));

		}


		// vec part 1/R^3 singularity
		for (unsigned int iv = 0; iv != 3; ++iv) {

			if ((Rpl[iv] + Lpl[iv]) < 50 * eps)
				LOGpl = 0.0;
			else
				LOGpl = log(Rpl[iv] + Lpl[iv]);


			if ((Rmn[iv] + Lmn[iv]) < 50 * eps)
				LOGmn = 0;
			else
				LOGmn = log(Rmn[iv] + Lmn[iv]);

			CONST1 = (LOGpl - LOGmn);

			multconst(pomVEC, U[iv], CONST1);

			subtract(Ivec1R3, Ivec1R3, pomVEC);

		}


		// scalar part  1/R^3 singularity
		for (unsigned int is = 0; is != 3; ++is) {

			if (abs(D[is]) >= 10 * eps) {

				ATANpl = atan((Pnul[is] * Lpl[is]) / (pow(Rnul[is], 2) + abs(D[is])*Rpl[is]));

				ATANmn = atan((Pnul[is] * Lmn[is]) / (pow(Rnul[is], 2) + abs(D[is])*Rmn[is]));

				Isca1R3 = Isca1R3 + dot(Pnulvec[is], U[is])*((1 / abs(D[is]))*(ATANpl - ATANmn));
			}
			else if ((abs(D[is]) < 10 * eps) && (Pnul[is] >= 10 * eps)) {

				Isca1R3 = Isca1R3 - dot(Pnulvec[is], U[is])*(Lpl[is] / (Pnul[is] * Rpl[is]) - Lmn[is] / (Pnul[is] * Rmn[is]));


			}
			else
				Isca1R3 = 0.0;

		}


		// testing Galerkin
		for (unsigned int i = 0; i != 3; ++i) {


			//vector<double> p(3);
			double* p = mesh.getCoord(n1[i]);

			subtract(pomocvF, &fielpoin[0], p);
			double konstF = (1 / (2 * AR1));

			multconst(RWGf, pomocvF, konstF);



			for (unsigned int j = 0; j != 3; ++j) {

				//vector<double> q(3);
				double* q = mesh.getCoord(n2[j]);

				subtract(Rivec, &fielpoin[0], q);

				double K1 = dot(&nvec2[0], q);

				multconst(pomNvec, &nvec2[0], K1);

				subtract(rho0N, q, pomNvec);

				KONST1 = ((pow(k, 2)) / COMPLEX(2));
				KONST2 = KONST1*Isca1R;

				subtract(pomVec1, rho, rho0N);

				multconst1(pomVec2, Ivec1R, KONST1);

				multconst1(pomVec3, pomVec1, KONST2);

				add(pomVec4, pomVec2, pomVec3);

				CONST1 = (1 / (2 * AR2));

				multconst2(Ikon1R, pomVec4, CONST1);

				multconst(pomVec5, pomVec1, Isca1R3);

				add(pomVec6, Ivec1R3, pomVec5);

				multconst(Ikon1R3, pomVec6, CONST1);



				add(Ikon, Ikon1R3, Ikon1R);

				cross(KROSS1, Rivec, Ikon);

				KROSS2= RWGf;

				alok2[i * 3 + j] = alok2[i * 3 + j] + wf*det1*(dot(KROSS2, KROSS1));

			}
		}
	}
}



void assemble_system_matrixMFIE(vector<COMPLEX> &A, Mesh &mesh, vector<Trianinfo> &Triangles, Points &points, unsigned int Nt, unsigned int maxele, COMPLEX k)


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
	double* CROSS1 = new double[3];
	double* CROSS2 = new double[3];
	double* Rvec = new double[3];
	double* Rivec = new double[3];

	//#pragma omp parallel for

	for (int ele1 = 0; ele1 < Nt; ++ele1)
	{
		int* n1 = mesh.getNOvertex(ele1);
		
		double*p1 = mesh.getCoord(*n1);
		
		double* p2 = mesh.getCoord(*(n1 + 1));
		
		double* p3 = mesh.getCoord(*(n1 + 2));

		vector<int> rwg1 = mesh.getRWG(ele1);

		Trianinfo trianF = Triangles[ele1];

		vector<double> nvecF = trianF.getnorm();

		double det1 = trianF.getDeter();

		double AR1 = det1 / 2.0;



		for (int ele2 = 0; ele2 != Nt; ++ele2) {

		

			fill(alok1.begin(), alok1.end(), COMPLEX(0));
			fill(alok2.begin(), alok2.end(), COMPLEX(0));

			int* n2 = mesh.getNOvertex(ele2);
			
			double* q1 = mesh.getCoord(*n2);
			
			double* q2 = mesh.getCoord(*(n2 + 1));
			
			double* q3 = mesh.getCoord(*(n2 + 2));

			vector<int> rwg2 = mesh.getRWG(ele2);
			Trianinfo trianS = Triangles[ele2];
			double det2 = trianS.getDeter();
			double AR2 = det2 / 2.0;
			vector<double> nvecS = trianS.getnorm();


			if (is_regularM(trianF, trianS)) {

				/// computation of regular submatrix*******************************************
				for (int nf = 0; nf != 4; ++nf) {

					double N1 = PointsNS[nf][0];
					double N2 = PointsNS[nf][1];
					double N0 = 1.0 - N1 - N2;

					double wf = WeightsNS[nf];

					fielpoin[0] = *p1*N1 + *p2*N2 + *p3*N0;
					fielpoin[1] = *(p1 + 1)*N1 + *(p2 + 1)*N2 + *(p3 + 1)*N0;
					fielpoin[2] = *(p1 + 2)*N1 + *(p2 + 2)*N2 + *(p3 + 2)*N0;

					for (int ns = 0; ns != 4; ++ns) {

						N1 = PointsNS[ns][0];
						N2 = PointsNS[ns][1];
						N0 = 1.0 - N1 - N2;

						double ws = WeightsNS[ns];

						sourcepoin[0] = *q1*N1 + *q2*N2 + *q3*(N0);
						sourcepoin[1] = *(q1 + 1)*N1 + *(q2 + 1)*N2 + *(q3 + 1)*(N0);
						sourcepoin[2] = *(q1 + 2)*N1 + *(q2 + 2)*N2 + *(q3 + 2)*(N0);


						double R = norm(&fielpoin[0], &sourcepoin[0]);

						subtract(Rvec, &fielpoin[0], &sourcepoin[0]);

						COMPLEX GradGreen = (exp(-I*k*R) / (pow(R, 3)))*(COMPLEX(1) + I*k*R);

						for (unsigned int i = 0; i != 3; ++i) {

							//vector<double> p(3);
							double* p = mesh.getCoord(n1[i]);

							subtract(pomocvF, &fielpoin[0], p);
							double konstF = (1 / (2 * AR1));

							multconst(RWGf, pomocvF, konstF);

							for (unsigned int j = 0; j != 3; ++j) {

								//vector<double> q(3);
								double* q = mesh.getCoord(n2[j]);

								subtract(pomocvS, &sourcepoin[0], q);

								double konstS = (1 / (2 * AR2));

								multconst(RWGs, pomocvS, konstS);


								cross(CROSS1, Rvec, RWGs);
								CROSS2 = CROSS1;

								alok1[i * 3 + j] = alok1[i * 3 + j] + wf*ws*det1*det2*(GradGreen*dot(RWGf, CROSS2));
							}
						}
					}
				}
			}

			

			else if (!is_sameT) {

				/// computation of near singular submatrix*******************************************
				for (int nf = 0; nf != 12; ++nf) {

					double N1 = PointsS[nf][0];
					double N2 = PointsS[nf][1];
					double N0 = 1.0 - N1 - N2;

					double wf = WeightsS[nf];

					fielpoin[0] = *p1*N1 + *p2*N2 + *p3*N0;
					fielpoin[1] = *(p1 + 1)*N1 + *(p2 + 1)*N2 + *(p3 + 1)*N0;
					fielpoin[2] = *(p1 + 2)*N1 + *(p2 + 2)*N2 + *(p3 + 2)*N0;

					for (int ns = 0; ns != 12; ++ns) {

						N1 = PointsS[ns][0];
						N2 = PointsS[ns][1];
						N0 = 1.0 - N1 - N2;
						double ws = WeightsS[ns];

						sourcepoin[0] = *q1*N1 + *q2*N2 + *q3*(N0);
						sourcepoin[1] = *(q1 + 1)*N1 + *(q2 + 1)*N2 + *(q3 + 1)*(N0);
						sourcepoin[2] = *(q1 + 2)*N1 + *(q2 + 2)*N2 + *(q3 + 2)*(N0);


						double R = norm(&fielpoin[0], &sourcepoin[0]);

						subtract(Rvec, &fielpoin[0], &sourcepoin[0]);

					COMPLEX GradGreenNS;
					
					if(R<eps)
				GradGreenNS=COMPLEX(0);
					else
				GradGreenNS = (exp(-I*k*R) / (pow(R, 3)))*(COMPLEX(1) + I*k*R) - COMPLEX(1) / (pow(R, 3)) - ((pow(k, 2)) / COMPLEX(2))*(1 / R);

						for (unsigned int i = 0; i != 3; ++i) {

							//vector<double> p(3);
							double* p = mesh.getCoord(n1[i]);

							subtract(pomocvF, &fielpoin[0], p);
							double konstF = (1 / (2 * AR1));

							multconst(RWGf, pomocvF, konstF);

							for (unsigned int j = 0; j != 3; ++j) {

								//vector<double> q(3);
								double* q = mesh.getCoord(n2[j]);

								subtract(pomocvS, &sourcepoin[0], q);

								double konstS = (1 / (2 * AR2));

								multconst(RWGs, pomocvS, konstS);

								subtract(Rivec, &fielpoin[0], q);

								cross(CROSS1, Rivec, RWGs);
								CROSS2= RWGf;

								alok1[i * 3 + j] = alok1[i * 3 + j] + wf*ws*det1*det2*(GradGreenNS*dot(CROSS1, CROSS2));

							}
						}
					}
				}


				//CALL singularity i zbrojiti te dvije matrice

				singularityMFIE(det1, AR1, AR2, p1, p2, p3, q1, q2, q3, nvecS, nvecF, PointsS, WeightsS, n1, n2, mesh, k, alok2);

				for (unsigned int i = 0; i != 3; ++i) {
					for (unsigned int j = 0; j != 3; ++j) {

						alok1[i * 3 + j] = alok1[i * 3 + j] + alok2[i * 3 + j];
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
					A[j1*maxele + j2] += double(1 / (4 * PI))*double(s1 * s2) * alok1[i1 * 3 + i2];
				}
			}
		}
	}
}

}//end of namespace



