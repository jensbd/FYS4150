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


// Function to initialise energy and magnetization
void InitializeLattice(int, mat &, double&, double&);

void InitializeLattice2(int, mat &, double&, double&);
// The metropolis algorithm including the loop over Monte Carlo cycles
void MetropolisSampling(int, int, double, vec &);

void MetropolisSampling2(int, int, double, vec &);
// prints to file the results of the calculations
void WriteResultstoFile(ofstream&, int, int, double, vec);

#endif // ISING_H
