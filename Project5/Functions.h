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

using namespace  std;
using namespace arma;



void forward_step(double, vec &, vec &, int);
void forward_euler(double, vec &, int, int);
void tridiag(double, vec &, int);
void backward_euler(double, vec &, int, int);
void crank_nicolson(double, vec &, int, int);
void g(vec &, int);
