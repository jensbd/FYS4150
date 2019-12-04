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



int main(int argc, char* argv[]){

  cout << "\n" << "Project 5: Partial Differential Equations " << endl;
  cout << "\n" << "Which dimension of the diffusion equation do you want to run?: " << endl;
  cout << "\n" << "Project 5 - 1-dim: " <<  "Write 1 " << endl;
  cout << "\n" << "Project 5 - 2-dim: " <<  "Write 2 " << endl;



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
  int T = N;

  // Defining u
  mat u = zeros<mat>(T,N+2);
  mat u_analytic = zeros<mat>(T,N+2);



  vec x = linspace<vec>(0,1,N+2);


  // Implement boundaries rigidly
  // Boundary condition for first endpoint is already zero, so just need to
  // implement boundary for last end point
  for (int t = 0; t < T; t++){
    u(t,N+1) = 1.0;
    u_analytic(t,N+1) = 1.0;
  }

  // Initial codition
  g(u, N);
  g(u_analytic, N);

  analytic(u_analytic, N, T, x, dt);

//  u = u.t();
//  u_analytic = u_analytic.t();

  cout << "\n" << "Which method do you want to use to solve the diffusion equation?: " << endl;
  cout << "\n" << "Forward Euler- Explicit Scheme: " <<  "Write FE " << endl;
  cout << "\n" << "Backward Euler- Implicit Scheme: " <<  "Write BE " << endl;
  cout << "\n" << "Crank-Nicolson- Implicit Scheme: " <<  "Write CN " << endl;


  cout << "\n" << "Write here " << endl;
  string Method;
  cin >> Method;

if (Method == "FE"){
  forward_euler(alpha, u, N, T);
  cout << u << endl;

}
if (Method == "BE"){
  backward_euler(alpha,u,N,T);
  cout << u << endl;
}
if (Method == "CN"){
  crank_nicolson(alpha, u, N, T);
  cout << u << endl;
}





  }


  if (Dimension == "2"){


    cout << "\n" << "In progress: " << endl;


  /* Simple program for solving the two-dimensional diffusion
     equation or Poisson equation using Jacobi's iterative method
     Note that this program does not contain a loop over the time
     dependence. It uses OpenMP to parallelize

*/

cout << "\n" << "Choose step size for Delta x: " << endl;
cout << "\n" << "Delta x = 0.1: " <<  "Write 0.1 " << endl;
cout << "\n" << "Delta x= 0.01 " <<  "Write 0.01 " << endl;

double PI = 4*atan(1);


// Step length in position x
double dx;
cin >> dx;

// Step length in time
double dt = 0.25*dx*dx;


// Defining alpha
double alpha = dt/dx/dx;

// Number of integration points along x-axis (inner points only)
int N = (1.0/(dx));

// Number of time steps till final time
int T = N;

double ExactSolution;
double tolerance = 1.0e-10;
mat A = zeros<mat>(N,N);
mat q = zeros<mat>(N,N);

// setting up an additional source term
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
      q(i,j) = -2.0*PI*PI*sin(PI*dx*i)*sin(PI*dx*j);
    }
}
  int itcount = JacobiSolver(N,dx,dt,A,q,tolerance);

  // Testing against exact solution
  double sum = 0.0;
  for(int i = 0; i < N; i++){
    for(int j=0;j < N; j++){
      ExactSolution = -sin(PI*dx*i)*sin(PI*dx*j);
      sum += fabs((A(i,j) - ExactSolution));
    }
  }
  cout << setprecision(5) << setiosflags(ios::scientific);
  cout << "Jacobi method with error " << sum/N << " in " << itcount << " iterations" << endl;


}
  return 0;

}
