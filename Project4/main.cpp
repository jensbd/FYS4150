/*
   Program to solve the two-dimensional Ising model
   with zero external field and no parallelization
   The coupling constant J is set to J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis aolgorithm  is used as well as periodic boundary conditions.
   The code needs an output file on the command line and the variables mcs, nspins,
   initial temp, final temp and temp step.
   Run as
   ./executable Outputfile numberof spins number of MC cycles initial temp final temp tempstep
   ./test.x Lattice 100 10000000 2.1 2.4 0.01
   Compile and link as
   c++ -O3 -std=c++11 -Rpass=loop-vectorize -o Ising.x IsingModel.cpp -larmadillo
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
#include "omp.h"
#include "Functions.h"
using namespace  std;
using namespace arma;
// output file
ofstream ofile;



// Main program begins here

int main(int argc, char* argv[])
{

  string filename;
  int NSpins, MCcycles;
  //double InitialTemp, FinalTemp, TempStep;
  double Temp;
  if (argc <= 3) {
    cout << "Bad Usage: " << argv[0] <<
      " read output file, Number of spins, MC cycles and temperature" << endl;
    exit(1);
  }
  if (argc > 1) {
    filename=argv[1];
    NSpins = atoi(argv[2]);
    MCcycles = atoi(argv[3]);
    Temp = atof(argv[4]);
    //FinalTemp = atof(argv[5]);
    //TempStep = atof(argv[6]);
  }
  // Declare new file name and add lattice size to file name
  string fileout = filename;
  string argument = to_string(NSpins);
  fileout.append(argument);
  ofile.open(fileout);

  ofile << "|   # Monte Carlo cycles  | Energy-Mean | Magnetization-Mean | # Accepted configurations |  Specific heat  | Susceptibility | Temperature |\n";

  int Nconfigs;

  vec Energies = zeros<mat>(400);
  vec counter = zeros<mat>(400);


  // Start Monte Carlo sampling by looping over the selcted Temperatures
  //for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature+=TempStep){
  vec ExpectationValues = zeros<mat>(5);
  // Start Monte Carlo computation and get expectation values
  //MetropolisSampling(NSpins, MCcycles, Temp, ExpectationValues, Nconfigs, false, Energies, counter);

  //WriteResultstoFile(ofile, NSpins, MCcycles, Temp, ExpectationValues, Nconfigs);
 // }
  ofile.close();  // close output file


/*

  cout << "Project Task 4c for ordered and unordered spin: " << endl;
  //cout << "Declare new file name : " << endl;
  //string file;
  //cin >> file;

  string file = "Ordered";

  ofile.open(file);
  //ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "|   # Monte Carlo cycles  | Energy-Mean | Magnetization-Mean | # Accepted configurations |  Specific heat  | Susceptibility | Temperature |\n";
  // Start Monte Carlo sampling by looping over the selcted Temperatures
  int N;
  long int MC;
  double T;
  cout << "Read in the number of spins" << endl;
  cin >> N;
  cout << "Read in the number of Monte Carlo cycles in times of 10" << endl;
  cin >> MC;
  cout << "Read in the given value for the Temperature" << endl;
  cin >> T;

  int iterations;
  //#pragma omp parallel for
  for (int i=1; i <= MC; i++){
    vec ExpectationValue = zeros<mat>(5);
    iterations = 10*i;
    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(N, iterations, T, ExpectationValue, Nconfigs, false);
    //
    WriteResultstoFile(ofile, N, iterations, T, ExpectationValue, Nconfigs);
  }
  ofile.close();  // close output file


  //cout << "Project Task 4c for unordered spin: " << endl;
  //cout << "Declare new file name : " << endl;
  //string file2;
  //cin >> file;

  string file2 = "Unordered";

  ofile.open(file2);
  //ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "| # Monte Carlo cycles  |  Energy-Mean   |  Magnetization-Mean  |  # Accepted configurations  |  Specific heat  |  Susceptibility   |  Temperature |\n";


  int iterations2;
  //#pragma omp parallel for
  for (int i=1; i <= MC; i++){
    vec ExpectationValue2 = zeros<mat>(5);
    iterations2 = 10*i;
    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(N, iterations2, T, ExpectationValue2, Nconfigs, true);
    //
    WriteResultstoFile(ofile, N, iterations2, T, ExpectationValue2, Nconfigs);
  }
  ofile.close();  // close output file

*/


/*
  cout << "Project Task 4d for Probability: " << endl;
  //cout << "Declare new file name : " << endl;
  //string file;
  //cin >> file;

  string file = "Probability";

  ofile.open(file);
  //ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "|  Energies | Energy counts |\n";
  // Start Monte Carlo sampling by looping over the selcted Temperatures
  int N;
  long int MC;
  double T;
  cout << "Read in the number of spins" << endl;
  cin >> N;
  cout << "Read in the number of Monte Carlo cycles in times of 10" << endl;
  cin >> MC;
  cout << "Read in the given value for the Temperature" << endl;
  cin >> T;

  int iterations;
  //#pragma omp parallel for
  for (int i=1; i <= MC; i++){
    vec ExpectationValue = zeros<mat>(5);
    iterations = 10*i;
    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(N, iterations, T, ExpectationValue, Nconfigs, false, Energies, counter);

    if (i == MC){
    Writeprobabilities(ofile, Energies, counter, N, iterations, ExpectationValue);
  }

  }
  ofile.close();  // close output file
*/

/*
cout << "Project Task 4e: " << endl;

string file = "Temperature";

ofile.open(file);
ofile << "| Temperature | Energy-Mean | Magnetization-Mean |  Specific heat  | Susceptibility |\n";

// Start Monte Carlo sampling by looping over the selcted Temperatures
int N_start, N_step, N_final;
long int MC;
double T_start, T_step, T_final, T;

// cout << "Read in the number of spins" << endl;
// cin >> N;
cout << "Read in the number of Monte Carlo cycles" << endl;
cin >> MC;


cout << "Read in the initial value for the Temperature" << endl;
cin >> T_start;
cout << "Read in the step size for the Temperature" << endl;
cin >> T_step;
cout << "Read in the final value for the Temperature" << endl;
cin >> T_final;


T_start = 2.2;
T_step = 0.025;
T_final = 2.4;

// Declare a matrix which stores the expectation values
mat values = zeros<mat>(5, 9);

N_start = 40;
N_step = 20;
N_final = 100;


// for (int N = N_start; N <= N_final; N += N_step){
// ofile << "\n";
// ofile << "\n";
// ofile << "Nspin =  " << N;
// ofile << "\n";

#pragma omp parallel for
// Start Monte Carlo sampling by looping over the selcted Temperatures
for (int i = 0; i <= 8; i++){
  vec ExpectationValue = zeros<mat>(5);

  T = T_start + T_step*i;

  // Start Monte Carlo computation and get expectation values
  MetropolisSampling(40, MC, T, ExpectationValue, Nconfigs, false, Energies, counter);
  //
  //WriteResultstoFile2(ofile, N, MC, T, ExpectationValue, Nconfigs);
  values.col(i) = ExpectationValue;
}

ofile << values;

// }
ofile.close();  // close output file

*/
  return 0;

}
