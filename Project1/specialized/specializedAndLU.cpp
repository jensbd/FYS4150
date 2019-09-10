#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "time.h"
#include <algorithm>
#include <armadillo>
using namespace std;
using namespace arma;
ofstream ofile;
void LU(int n);
double solution(double x){
    return 1.0-(1-exp(-10))*x-exp(-10*x);
}
double f(double x){
    return 100*exp(-10*x);
}

int main(int argc, char* argv[]){
    char *filename;
    int n;

    if(argc <= 2){
        cout << "Too few arguments, 1 filename and 1 integer needed" << endl;
        exit(1);
    }
    else{
        filename = argv[1];
        n = atoi(argv[2]);
    }
    double h = 1.0/(n);
    double *x = new double[n+1];
    double *rhs = new double[n+1]; //Right hand side

    //Diagonal
    double *b = new double[n+1];


    //Analytic and numerical solution
    double *analytic = new double[n+1];
    double *v = new double[n+1];
    analytic[0] = 0;
    v[0] = v[n] = 0;
    b[0] = b[n] = 0;
    for (int i = 1; i < n; i++){
        b[i] = (i+1.0)/((double) i);
    }
    for (int i = 0; i <= n; i++){
        x[i] = i*h;
        rhs[i] = h*h*f(x[i]);
    }
    // Assigning values to right hand side, numerical and analytical
    for (int i = 1; i <= n; i++){
        analytic[i] = solution(x[i]);
    }
    clock_t start, finish;
    start = clock();
    cout << "Start: " << start << endl;
    //Forward substitution
    for (int i = 2; i < n; i++){
        rhs[i] = rhs[i] + rhs[i-1]/b[i-1];
    }

    //Backward substitution
    v[n-1] = rhs[n-1]/b[n-1];
    for (int i = n-2; i>0; i--){
        v[i] = (rhs[i] + v[i+1])/b[i];
    }
    finish = clock();

    //Print amount of time used by algorithm
    cout << "Finish: " << finish << endl;
    cout << "Specialized algorithm time: "<<(( double(finish - start))/(double(CLOCKS_PER_SEC) )) << endl;

    //Relative error
    double *e = new double[n+1];
    double maxerror = 0.0;
    for (int i = 1; i < n; i++){
        double error = log10(abs((v[i]-analytic[i])/analytic[i]));
        e[i] = error;
        if ( abs(error) > abs(maxerror)){
            maxerror = error;
        }
    }
    cout << "Max error: " << maxerror<< endl;
    //Output to file
    ofile.open(filename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "       x:             u(x):          v(x):  " << endl;
    for (int i = 1; i<n; i++){
        ofile << setw(15) << setprecision(8) <<x[i];
        ofile << setw(15) << setprecision(8) <<analytic[i];
        ofile << setw(15) << setprecision(8) <<v[i]<< endl;
    }
    ofile.close();
    delete [] x;
    delete [] rhs;
    delete [] b;
    delete [] analytic;
    delete [] v;
    //Run LU decomp
    LU(n);

    return 0;
}

void LU(int n){
    mat A(n+2,n+2, fill::zeros); //n+2 to include the boundary conditions
    vec v(n+2, fill::zeros);
    vec x(n+2);
    vec rhs(n+2);
    v[0] = v[n+2] = 0;
    double h = 1.0/double(n-1);
    //Fill in for a, b and c
    A.diag(0).fill(2);
    A.diag(-1).fill(-1);
    A.diag(1).fill(1);
    for (int i = 0; i <n+2; i++){
        x[i] = i*h;
        rhs[i] = h*h*f(x[i]);
    }
    clock_t lustart, lufinish;
    lustart = clock();
    cout << "LU start: " << lustart << endl;
    v = solve(A,rhs);
    lufinish = clock();
    cout << "LU finish: " << lufinish << endl;
    cout << "LU time: "<<(( double(lufinish - lustart))/(double(CLOCKS_PER_SEC) )) << endl;
}
