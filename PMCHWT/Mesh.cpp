#include "stdafx.h"
#include <iostream>
#include "Mesh.h"
#include <vector>
using namespace std;

Mesh::Mesh(vector<double> co, vector<int> top, vector<int> tri):
	coord(co), topol(top), trian(tri)
	
{

}

double* Mesh::getCoord(const int vertex_number) 
{	

	return &(coord[vertex_number * 3]);
}


int* Mesh::getNOvertex(const int trian_number)
{
	
	return &(topol[trian_number*3]);
}

vector<int> Mesh::getRWG(const int trian_number)
{
	vector<int> NOrwg;

	for (int i = 0; i < 3; ++i) {

		NOrwg.push_back(trian[trian_number*3+i]);

	}

	return NOrwg;
}
