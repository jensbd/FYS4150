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
using namespace std;
using namespace arma;

ofstream ofile;

int main(int argc, char* argv[]){

  cout << "\n" << "Project 5: Partial Differential Equations " << endl;
  cout << "\n" << "Which dimension of the diffusion equation do you want to run?: " << endl;
  cout << "\n" << "Project 5 - 1-dim: " <<  "Write 1 " << endl;
  cout << "\n" << "Project 5 - 2-dim: " <<  "Write 2 " << endl;

  double PI = 4*atan(1);

  cout << "\n" << "Write here " << endl;
  string Dimension;
  cin >> Dimension;

  // 1-Dimensional Diffusion Equation - Analytical vs Numerical Results for a 2x2 Lattice

  if (Dimension == "1"){



  cout << "\n" << "Choose step size for Delta x: " << endl;
  cout << "\n" << "Delta x = 0.1: " <<  "Write 0.1 " << endl;
  cout << "\n" << "Delta x= 0.01 " <<  "Write 0.01 " << endl;

  // Step length in position x
  double dx;
  cin >> dx;

  // Step length in time
  double dt = 0.5*dx*dx;


  // Defining alpha
  double alpha = dt/dx/dx;

  // Number of integration points along x-axis (inner points only)
  int N = int(1.0/(dx));

  // Number of time steps till final time
  int T = int(1.0/dt);

  // Defining u
  mat u = zeros<mat>(T,N+2);
  mat u_analytic = zeros<mat>(T,N+2);



  vec x = linspace<vec>(0,1,N+2);


  // Implement boundaries rigidly
  // Boundary condition for first endpoint is already zero, so just need to
  // implement boundary for last end point
  for (int t = 0; t < T; t++){
    u(t,N+1) = 1.0;
  }

  // Initial codition
  //g(u, N);
  //g(u_analytic, N);

  analytic(u_analytic, N, T, x, dt);



/*  cout << "\n" << "Which method do you want to use to solve the diffusion equation?: " << endl;
  cout << "\n" << "Forward Euler- Explicit Scheme: " <<  "Write FE " << endl;
  cout << "\n" << "Backward Euler- Implicit Scheme: " <<  "Write BE " << endl;
  cout << "\n" << "Crank-Nicolson- Implicit Scheme: " <<  "Write CN " << endl;


  cout << "\n" << "Write here " << endl;
  string Method;
  cin >> Method;
*/

string Method = "FE";
string file = "Analytic:"+to_string(dx);
//Remove trailing zeros from string
file.erase ( file.find_last_not_of('0') + 1, std::string::npos );
ofile.open(file);
ofile << u_analytic;
ofile.close();

if (Method == "FE"){
  mat uf = u;
  forward_euler(alpha, uf, N, T);
  file = "FE:"+to_string(dx);
  file.erase ( file.find_last_not_of('0') + 1, std::string::npos );
  ofile.open(file);
  ofile << uf;
  ofile.close();
  Method = "BE";
}

if (Method == "BE"){
  mat ub = u;
  backward_euler(alpha,ub,N,T);
  file = "BE:"+to_string(dx);
  file.erase ( file.find_last_not_of('0') + 1, std::string::npos );
  ofile.open(file);
  ofile << ub;
  ofile.close();
  Method = "CN";

}


if (Method == "CN"){
  mat uc = u;
  crank_nicolson(alpha, uc, N, T);
  file = "CN:"+to_string(dx);
  file.erase ( file.find_last_not_of('0') + 1, std::string::npos );
  ofile.open(file);
  ofile << uc;
  ofile.close();

}
cout << "Text files generated" << endl;



  }


  if (Dimension == "2"){





cout << "\n" << "Choose step size for Delta x,y: " << endl;
cout << "\n" << "Delta x,y = 0.1: " <<  "Write 0.1 " << endl;
cout << "\n" << "Delta x,y = 0.01 " <<  "Write 0.01 " << endl;

// Step length for both x and y
double dx;
cin >> dx;

// Step length in time
double dt = 0.2*dx*dx;


// Defining alpha
double alpha = dt/dx/dx;

// Number of integration points along x and y (inner points only)
int N = int(1.0/(dx));

// Number of time steps till final time
int T = int(1.0/dt);


cube u = zeros<cube>(N+2, N+2, T);
cube u_analytic = zeros<cube>(N+2,N+2,T);
//Boundary conditions
for (int t = 0; t < T; t++){
  for (int n = 0; n < N+2; n++){
    u(n,N+1,t) = 0.0;
    u(N+1,n,t) = 0.0;
  }
}

// Initial conditions
  for (double i = 1; i < N+1; i++) {
    for (double j =1; j < N+1; j++) {
      u(i,j,0) = sin(i*PI/(double)N)*sin(j*PI/(double)N);
    }
  }

forward_euler_2dim(alpha,u,N,T);

string file = "2dim_explicit:"+to_string(dx);
file.erase ( file.find_last_not_of('0') + 1, std::string::npos );
ofile.open(file);
for (int t = 0; t < T; t++){
  ofile << u.slice(t);
}
ofile.close();



mat u_implicit = zeros<mat>(N+2, N+2);


// Implement boundaries rigidly
// Boundary condition for endpoints is already zero
//Boundary conditions
for (int n = 0; n < N+2; n++){
  u_implicit(n, N+1) = 0.0;
  u_implicit(N+1, n) = 0.0;
}
// Initial conditions
  for (double i = 1; i < N+1; i++) {
    for (double j =1; j < N+1; j++) {
      u_implicit(i,j) = sin(i*PI/(double)N)*sin(j*PI/(double)N);
    }
  }


double ExactSolution;
double tolerance = 1.0e-8;

double start = omp_get_wtime();
int itcount = JacobiSolver(u_implicit, dx, dt, tolerance);
double end = omp_get_wtime();
double comptime = end-start;
cout << "Time used for Jacobis method: " << comptime << " s" << endl;

/*
// Testing against exact solution
double sum = 0.0;
for(int i = 0; i < N; i++){
  for(int j=0;j < N; j++){
    ExactSolution = -sin(PI*dx*i)*sin(PI*dx*j);
    sum += fabs((u_implicit(i,j) - ExactSolution));
  }
}

cout << setprecision(5) << setiosflags(ios::scientific);
cout << "Jacobi method with error " << sum/N<< " in " << itcount << " iterations" << endl;
*/



}

  return 0;
}
