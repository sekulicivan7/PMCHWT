#include "stdafx.h"
#include <iostream>
#include "Points.h"
#include <vector>
using namespace std;

vector<double> Points::getWeightsNS() {

	vector<double> w(4);
	w[0] = -2.812500000000000e-01;
	w[1] = w[2] = w[3] =2.604166666666670e-01;
	
	return w;

}

 vector<vector<double>> Points::getPointsNS() {

	vector<vector<double>>P;
	double intvalue = 0;
	P.resize(4, vector<double>(2, intvalue));

	P[0][0] = P[0][1] = (1.0 / 3.0);

	//P[1][0] = P[2][1] = 0.73333333333333;
		
	//P[1][1] = P[2][0] = P[3][0] =P[3][1]= 0.1333333333333;
	

	P[1][0] = P[2][1] = P[3][0] = P[3][1] = 0.2;
	P[1][1] = P[2][0] = 0.6;
	
	return P;

}

 vector<double> Points::getWeightsS() {

	 vector<double> w(12);

	 w[0] = w[2] = w[4] = 5.839313786318950e-02;
	 w[1] = w[3] = w[5] = 2.542245318510350e-02;
	 for (int j = 6; j < 12; ++j) w[j] = 4.142553780918700e-02;

	 return w;

 }

 vector<vector<double>> Points::getPointsS() {

	 vector<vector<double>>P;
	 double intvalue = 0;
	 P.resize(12, vector<double>(2, intvalue));

	 P[0][0] = P[2][1] = P[4][0] = P[4][1]
		 = 2.492867451709110e-01;
	 P[0][1] = P[2][0] = 5.014265096581790e-01;
	 P[1][0] = P[3][1] = P[5][0] = P[5][1]
		 = 6.308901449150210e-02;

	 P[1][1] = P[3][0] = 8.738219710169960e-01;
	 P[6][0] = P[9][1] = P[10][1] = P[11][0]
		 = 6.365024991213990e-01;
	 P[6][1] = P[7][1] = P[8][0] = P[9][0]
		 = 5.314504984481700e-02;
	 P[7][0] = P[8][1] = P[10][0] = P[11][1]
		 = 3.103524510337840e-01;

	 return P;

 }
