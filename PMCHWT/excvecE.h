#pragma once
#ifndef _EXCVECE_HEADER_
#define _EXCVECE_HEADER_
#include <vector>
#include "Mesh.h"
#include "Points.h"
#include "Trianinfo.h"
#include <complex>

#define COMPLEX complex<double>



namespace EFIE{

namespace excEFIE{


void assemble_exic_vector(vector<COMPLEX> &C, Mesh &mesh, vector<Trianinfo> &Triangles, Points &points, int Nt, COMPLEX k);


}

}

#endif
