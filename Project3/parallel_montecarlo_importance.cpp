#include "ran0.cpp"
#include "time.h"
#include <cmath>
#include <iostream>
#include <omp.h>

using namespace std;


double exponential_MC(double *);

int main(){
  double const  pi = 3.14159265359;

  int nImportance = pow(10,5);   //Number of simulations for importance sampling MC.
  double xMCimportance[6], fxImportance;
  double int_mcImportance = 0.; double varianceImportance = 0.;
  double sum_sigmaImportance = 0;long idum2 = -1;
  double jacobidetImportance = (pow(pi,4)/4);

  // Start timing importance sampling Monte Carlo:
  double start4, finish4;


  int numthreads = 4;
  omp_set_dynamic(0);
  omp_set_num_threads(4);
  nImportance = nImportance/numthreads;
  start4 = omp_get_wtime();
  #pragma omp parallel
  {
    //Different seed for each thread
    idum2 = omp_get_thread_num();
    #pragma omp for private(xMCimportance) reduction(+:fxImportance) reduction(+:sum_sigmaImportance)
    for (int i = 1; i <= nImportance; i++){

        xMCimportance[0] = (-1./4)*log(1-ran0(&idum2));
        xMCimportance[1] = (-1./4)*log(1-ran0(&idum2));
        xMCimportance[2] = pi*ran0(&idum2);
        xMCimportance[3] = pi*ran0(&idum2);
        xMCimportance[4] = 2*pi*ran0(&idum2);
        xMCimportance[5] = 2*pi*ran0(&idum2);
        fxImportance = exponential_MC(xMCimportance);
        int_mcImportance += fxImportance;
        sum_sigmaImportance += fxImportance*fxImportance;
    }
  }
  int_mcImportance = int_mcImportance/((double) nImportance);
  sum_sigmaImportance = sum_sigmaImportance/((double) nImportance);
  varianceImportance = sum_sigmaImportance - int_mcImportance*int_mcImportance;
  // Finish timing importance sampling Monte Carlo:
  finish4 = omp_get_wtime();
  double ComputationTimeMCimportance = (finish4-start4);

  cout << endl
       << "Threads: " << numthreads << endl
       << "N: " << nImportance*numthreads << endl
       << "Montecarlo (Importance): " << jacobidetImportance*int_mcImportance << endl
       << "Sigma: " << jacobidetImportance*sqrt(varianceImportance/((double) nImportance)) << endl
       << "Computation time: " << ComputationTimeMCimportance << endl;
}
// This function defines the function to integrate for importance sampling Monte Carlo:
double exponential_MC(double *x)
{
    double r1 = x[0]; double r2 = x[1]; double theta1 = x[2];
    double theta2 = x[3]; double phi1 = x[4]; double phi2 = x[5];
    double cosbeta = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1 - phi2);
    // Evaluate the denominator:
    double denom = sqrt(r1*r1 + r2*r2 - 2*r1*r2*cosbeta);
    if(denom <pow(10.,-6)||isnan(denom)) { return 0;}
    else return r1*r1*r2*r2*sin(theta1)*sin(theta2)/denom;
}
