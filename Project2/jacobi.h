#ifndef JACOBI_H
#define	JACOBI_H

#include <iostream>
#include <armadillo>
#include <cmath>
#include <time.h>

using namespace std;
using namespace arma;

void jacobiMethod(mat &A, mat &R, int k, int l, int n);
void offdiag(mat A, int& p, int& q, int n);

#endif /* JACOBI_H */
