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

  ofile << "| Temperature | Energy-Mean | Magnetization-Mean|    Cv    | Susceptibility |\n";

  // Start Monte Carlo sampling by looping over the selcted Temperatures
  //for (double Temperature = InitialTemp; Temperature <= FinalTemp; Temperature+=TempStep){
  vec ExpectationValues = zeros<mat>(5);
  // Start Monte Carlo computation and get expectation values
  MetropolisSampling(NSpins, MCcycles, Temp, ExpectationValues);

  WriteResultstoFile(ofile, NSpins, MCcycles, Temp, ExpectationValues);
 // }
  ofile.close();  // close output file




  cout << "Project Task 4c for ordered spin: " << endl;
  cout << "Declare new file name : " << endl;
  string file;
  cin >> file;

  ofile.open(file);
  //ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "| Temperature | Energy-Mean | Magnetization-Mean|    Cv    | Susceptibility |\n";
  // Start Monte Carlo sampling by looping over the selcted Temperatures
  int N;
  long int MC;
  double T;
  cout << "Read in the number of spins" << endl;
  cin >> N;
  cout << "Read in the number of Monte Carlo cycles in times of 100" << endl;
  cin >> MC;
  cout << "Read in the given value for the Temperature" << endl;
  cin >> T;

  int iterations;
  for (int i=1; i < MC; i++){
    vec ExpectationValue = zeros<mat>(5);
    iterations = 100*i;
    // Start Monte Carlo computation and get expectation values
    MetropolisSampling(N, iterations, T, ExpectationValue);
    //
    WriteResultstoFile(ofile, N, iterations, T, ExpectationValue);
  }
  ofile.close();  // close output file


  cout << "Project Task 4c for unordered spin: " << endl;
  cout << "Declare new file name : " << endl;
  string file2;
  cin >> file2;

  ofile.open(file2);
  //ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "| Temperature | Energy-Mean | Magnetization-Mean|    Cv    | Susceptibility |\n";


  int iterations2;
  for (int i=1; i < MC; i++){
    vec ExpectationValue2 = zeros<mat>(5);
    iterations2 = 100*i;
    // Start Monte Carlo computation and get expectation values
    MetropolisSampling2(N, iterations2, T, ExpectationValue2);
    //
    WriteResultstoFile(ofile, N, iterations2, T, ExpectationValue2);
  }
  ofile.close();  // close output file

  return 0;



}
