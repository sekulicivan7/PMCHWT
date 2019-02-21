#pragma once
#include <vector>
using namespace std;

class Points {

public:
	Points() {}
	vector<vector<double>> getPointsNS();
	vector<double> getWeightsNS();
	vector<vector<double>> getPointsS();
	vector<double> getWeightsS();
private:
	vector<double> w;
	vector<vector<double>> P;
};
