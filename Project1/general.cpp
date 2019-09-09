#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
using namespace std;
ofstream ofile;
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
    double h = 1.0/(n+1.0);
    double *x = new double[n+2];
    double *rhs = new double[n+1]; //Right hand side
    rhs[0] = 0;

    //Diagonals
    int *a = new int[n+1];
    int *b = new int[n+1];
    int *c = new int[n+1];

    //Analytic and numerical solution
    double *analytic = new double[n+2];
    double *v = new double[n+2];
    analytic[0] = 0;
    v[0] = 0;

    for (int i = 0; i < n+1; i++){
        x[i] = i*h;
    }
    // Assigning values to right hand side, numerical and analytical, as well ass a, b, c diagonals.
    for (int i = 1; i <= n; i++){
        rhs[i] = h*h*f(x[i]);
        analytic[i] = solution(x[i]);
        b[i] = 2;
        a[i] = -1;
        c[i] = -1;
    }
    c[n] = 0;
    a[1] = 0;

    double *diag_temp = new double[n+1];
    double b_temp = b[1];
    v[1] = rhs[1]/(b_temp);

    //Forward substitution
    for (int i = 2; i <= n; i++){
        diag_temp[i] = c[i-1]/b_temp;
        b_temp = b[i] - a[i]*diag_temp[i];
        v[i] = (rhs[i] - v[i-1]*a[i])/b_temp;
    }

    //Backward substitution
    for (int i = n; i>=1; i--){
        v[i] -= diag_temp[i+1]*v[i+1];
    }

    //Output to file
    ofile.open(filename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "       x:             u(x):          v(x):  " << endl;
    for (int i = 1; i<=n; i++){
        ofile << setw(15) << setprecision(8) <<x[i];
        ofile << setw(15) << setprecision(8) <<analytic[i];
        ofile << setw(15) << setprecision(8) <<v[i]<< endl;
    }
    ofile.close();

    delete [] x;
    delete [] rhs;
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] analytic;
    delete [] v;
    return 0;
}
