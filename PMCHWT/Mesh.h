#pragma once
#include <vector>
using namespace std;

class Mesh {
private:
	vector<double> coord;
	vector<int> topol;
	vector<int> trian;

public:
	Mesh(vector<double>, vector<int>, vector<int>);
	double* getCoord(const int vertex_number) ;
	int* getNOvertex(const int trian_number);
	vector<int> getRWG(const int trian_number);
};