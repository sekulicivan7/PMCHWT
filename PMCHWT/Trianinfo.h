#pragma once
#include <vector>

using namespace std;

class Trianinfo {
private:
	vector<double> coord1;
	vector<double> coord2;
	vector<double> coord3;
	vector<double> nvec;
	vector<double> cp;
	double length, deter;

	void calculate_parameters();

public:
	Trianinfo(double* c1, double* c2, double* c3); // creating a constructor
	const vector<double>& getnorm() const { return nvec; }
	const vector<double>& getcp() const { return cp; }
	double getLen() const { return length; }
	double getDeter() const { return deter; }

};
