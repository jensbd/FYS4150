#include <iostream>
#include <fstream>
#include <cmath>
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
    return 0;
}
