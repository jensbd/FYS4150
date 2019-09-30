#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"

#include <armadillo>
#include "jacobi.h"
using namespace arma;

TEST_CASE("Testing max a(i,j)"){
    int n = 4;
    Mat<double> A(n, n, fill::zeros);
    A(1,2) = 10;
    A(4,3) = 100;
    int p,q;
    offdiag(A, p, q, n);

    REQUIRE(p == 4 );
    REQUIRE(q == 3);
    REQUIRE(A(p,q) == 100);
}

TEST_CASE("Testing eigenvalues after one jacobi rotation"){
    int n = 4;
    mat A(n, n, fill::randu);
    mat eigvec;
    vec values;
    eig_sym(values, eigvec, A);
    int p,q;
    offdiag(A, p, q, n);
    jacobiMethod(A, eigvec, p, q, n);
    mat eigvec2;
    vec values2;
    eig_sym(values2, eigvec2, A);

    for (int i = 0; i < values.n_elem; i++){
        REQUIRE(values(i) == values2(i));
    }
}
