#include <iostream>
#include <armadillo>
#include <cmath>
#include "time.h"
using namespace std;
using namespace arma;

void jacobiMethod(mat &A, mat &R, int k, int l, int n){
    double s, c;
    if (A(k,l) != 0.0){
        double t, tau;
        tau = (A(l,l) - A(k,k))/(2.0*A(k,l));

        if (tau >= 0){
            t = 1.0/(tau+sqrt(1.0 +tau*tau));
        }
        else {
            t = -1.0/(-tau + sqrt(1.0 + tau*tau));
        }
        c = 1.0/sqrt(1+t*t);
        s = c*t;
    }
    else{
        c = 1.0;
        s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    //Assigning new values to the matrix
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0;
    A(l,k) = 0.0;
    for ( int i = 0; i < n; i++ ) {
        if ( i != k && i != l ) {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }
        // And finally the new eigenvectors
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
    return;
}

//Function for finding max nondiagonal element
void offdiag(mat A, int& p, int& q, int n){
    double max;
    for (int i = 0; i < n; ++i){
        for (int j = i+1; j<n; ++j){
            double aij = fabs(A(i,j));
            if (aij > max){
                max = aij;
                p = i;
                q = j;
            }
        }
    }
    return;
}



int main()
{
    int n = 400;
    float rho_min = 0.0;
    float rho_max = 10.0;
    double h = (rho_max-rho_min)/(n+1);
    vec rho = vec(n);
    for (int i = 0; i < n; i++){
        rho(i) = (i+1)*h +rho_min;
    }
    //Element-wise multiplication
    vec V =  rho%rho;
    vec d = 2.0/(h*h) + V;
    double a = -1.0/(h*h);
    vec e = vec(n-1);
    e.fill(a);
    Mat<double> A(n, n, fill::zeros);

    //Filling the matrix
    A.diag() = d;
    A.diag(1) = e;
    A.diag(-1) = e;

    Mat<double> eigenvectors(n,n, fill::zeros);

    vec values;
    mat vectors;
    clock_t start, finish;
    start = clock();
    cout << "Start: " << start << endl;
    eig_sym(values, vectors, A);
    finish = clock();
    cout << "Finish: "<< finish << endl;
    cout << "Eig_sym time: "  << ( double(finish - start)/double(CLOCKS_PER_SEC) ) << endl;

    double tol = 1.0E-8;
    int iterations = 0;
    int maxiter = 1000000;

    int p, q;
    offdiag(A,p,q,n);
    double maxnondiag = fabs(a);

    //Time Jacobi's method
    clock_t start2, finish2;
    start2 = clock();
    cout << "Start: " << start2 << endl;
    //Perform Jacobi's method until the off-diagonals are almost zero
    while(maxnondiag > tol && iterations <= maxiter){
        //mat eigenvectors = eigenvec(A);
        jacobiMethod(A, eigenvectors, p, q, n);
        offdiag(A,p,q,n);
        maxnondiag = fabs(A(p,q));
        iterations ++;
    }
    finish2 = clock();
    cout << "Finish: "<< finish2 << endl;
    cout << "Jacobi's method time: "  << ( double(finish2 - start2)/double(CLOCKS_PER_SEC) ) << endl;

    cout << "No. of iterations:" << iterations << endl;
    cout << "Analytical eigenvalues" << endl;
    cout << 3 << endl << 7 << endl << 11 << endl << 15 << endl;
    cout << "Matrix diagonal sorted" << endl;
    vec diag = diagvec(A);
    diag = sort(diag);
    for (int i = 0; i < 4; i++){
        cout << diag(i) << endl;
    }


    return 0;
}
