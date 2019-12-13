#ifndef FUNCTION_H
#define FUNCTION_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
#include <omp.h>

using namespace  std;
using namespace arma;



void forward_step(double, rowvec &, rowvec &, int, bool);
void forward_step_2dim(double, mat &, mat &, int);
void forward_euler(double, mat &, int, int);
void forward_euler_2dim(double, cube &, int, int);
void tridiagSolver(rowvec &, rowvec, double, int, bool);
void backward_euler(double, mat &, int, int);
void crank_nicolson(double, mat &, int, int);
void analytic(mat &, int, int, vec, double);
void analytic_2D(mat &, double, double);
int JacobiSolver(mat &, double, double, double);


#endif
