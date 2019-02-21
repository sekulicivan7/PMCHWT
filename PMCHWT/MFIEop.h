#pragma once
#ifndef _MFIEOP_HEADER_
#define _MFIEOP_HEADER_
#include <vector>
#include "Mesh.h"
#include "Points.h"
#include "Trianinfo.h"
#include <complex>

#define COMPLEX complex<double>

void assemble_system_matrixMFIE(vector<COMPLEX> &A, Mesh &mesh, vector<Trianinfo> &Triangles, Points &points, unsigned int Nt, unsigned int maxele, COMPLEX k);
#endif