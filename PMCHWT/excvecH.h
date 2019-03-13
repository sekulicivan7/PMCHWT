#pragma once
#ifndef _EXCVECH_HEADER_
#define _EXCVECH_HEADER_
#include <vector>
#include "Mesh.h"
#include "Points.h"
#include "Trianinfo.h"
#include <complex>

#define COMPLEX complex<double>



namespace MFIE{

namespace excMFIE{

void assemble_exic_vector(vector<COMPLEX> &C, Mesh &mesh, vector<Trianinfo> &Triangles, Points &points, int Nt, COMPLEX k, COMPLEX eta);

}

}

#endif

