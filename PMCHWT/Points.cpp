
#include <iostream>
#include "Points.h"
#include <vector>
using namespace std;



vector<double> Points::getWeightsNS() {

	vector<double> w(7);
	
	w[0] = 0.112500000000000;
	
	w[1] = w[2] = w[3] = 0.066197076394255;
	
	w[4] = w[5] = w[6] = 0.062969590272415;
	
	return w;

}

 vector<vector<double>> Points::getPointsNS() {

	vector<vector<double>>P;
	double intvalue = 0;
	P.resize(7, vector<double>(2, intvalue));

	P[0][0] = P[0][1] = (1.0 / 3.0);

	P[1][0] = P[1][1] = P[2][0] = P[3][1] = 0.47014206410511;
	
	P[2][1] = P[3][0] = 0.05971587178977;
	
	P[4][0]= P[4][1]= P[5][0]= P[6][1]= 0.10128650732346;
	
	P[5][1] = P[6][0] = 0.79742698535309;
	

	return P;

}



 vector<double> Points::getWeightsS() {

	 vector<double> w(12);

	 w[0] = w[1] = w[2] = 5.839313786318950e-02;
	 w[3] = w[4] = w[5] = 2.542245318510350e-02;
	 for (int j = 6; j < 12; ++j) w[j] = 4.142553780918700e-02;

	 return w;

 }

 vector<vector<double>> Points::getPointsS() {

	 vector<vector<double>>P;
	 double intvalue = 0;
	 P.resize(12, vector<double>(2, intvalue));

	 P[0][0] = P[0][1] = P[1][0] = P[2][1]
		 = 2.492867451709110e-01;
	 P[1][1] = P[2][0] = 5.014265096581790e-01;
	 P[3][0] = P[3][1] = P[4][0] = P[5][1]
		 = 6.308901449150210e-02;

	 P[4][1] = P[5][0] = 8.738219710169960e-01;
	 P[6][1] = P[7][0] = P[9][0] = P[11][1]
		 = 6.365024991213990e-01;
	 P[7][1] = P[8][0] = P[10][1] = P[11][0]
		 = 5.314504984481700e-02;
	 P[6][0] = P[8][1] = P[9][1] = P[10][0]
		 = 3.103524510337840e-01;

	 return P;

 }
