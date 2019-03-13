#pragma once
#include <vector>
#include "Mesh.h"
#include "Points.h"
#include "Trianinfo.h"
#include <complex>

#define COMPLEX complex<double>

namespace EFIE{

void assemble_system_matrixEFIE(vector<COMPLEX> &A, Mesh &mesh, vector<Trianinfo> &Triangles, Points &points, unsigned int Nt, unsigned int maxele, COMPLEX k, COMPLEX eta);

}
