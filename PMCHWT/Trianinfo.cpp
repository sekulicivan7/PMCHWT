
#include <iostream>
#include "Trianinfo.h"
#include "Products.h"
#include <vector>
using namespace std;

Trianinfo::Trianinfo(double* c1, double* c2, double* c3) :
	 coord1(3), coord2(3), coord3(3), cp(3), nvec(3), deter(0.0), length(0.0)
{
	 coord1={ *c1, *(c1 + 1), *(c1 + 2) };
	 coord2={ *c2, *(c2 + 1), *(c2 + 2) };
	 coord3={ *c3, *(c3 + 1), *(c3 + 2) };

	calculate_parameters();
}


void Trianinfo::calculate_parameters() {

	
	for (unsigned int j = 0; j != 3; ++j) {
		cp[j] = (1.0 / 3.0) * (coord1[j] + coord2[j] + coord3[j]); 
	}
	
	vector<double> p21(3), p32(3), p13(3);
	
	for (unsigned int j = 0; j != 3; ++j) {
		p21[j] = coord2[j] - coord1[j];
		p13[j] = coord1[j] - coord3[j];
		p32[j] = coord3[j] - coord2[j];
	}

	double* nvector = new double[3];
		
    cross(nvector, &p13[0], &p21[0]);
	nvec = { *nvector, *(nvector + 1), *(nvector + 2) };
	deter = norm(nvector);


	nvec[0] = nvec[0]/deter;
	nvec[1] = nvec[1]/deter;
	nvec[2] = nvec[2]/deter;
    
	length = norm(&p21[0]) + norm(&p32[0]) + norm(&p13[0]);

}







