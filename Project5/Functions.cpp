#include "Functions.h"
#include "omp.h"


// Code for solving the 1+1 dimensional diffusion equation
// du/dt = ddu/ddx on a rectangular grid of size L x (T*dt),
// with with L = 1, u(x,0) = g(x), u(0,t) = u(L,t) = 0



void forward_step(double alpha, vec &u, vec &uPrev, int N){
    /*
    Steps forward-euler algo one step ahead.
    Implemented in a separate function for code-reuse from crank_nicolson()
    */

    for (int i = 1; i < N+1; i++){ //loop from i=1 to i=N
        u(i) = alpha*uPrev(i-1) + (1.0-2*alpha)*uPrev(i) + alpha*uPrev(i+1);
      }
}

void forward_euler(double alpha, mat &u, int N, int T){
    /*
    Implements the forward Euler sheme, results saved to
    array u
    */

    // Skip boundary elements
    for (int t = 1; t < T; t++){
        forward_step(alpha,u(t),u(t-1),N);
      }
}

void tridiag(double alpha, vec &u, int N){
    /*
    Tridiagonal gaus-eliminator, specialized to diagonal = 1+2*alpha,
    super- and sub- diagonal = - alpha
    */

    vec d = zeros<mat>(N).fill(1+2*alpha);
    vec b = zeros<mat>(N-1).fill(- alpha);

    //Forward eliminate
    for (int i = 1; i < N; i++){
        //Normalize row i (i in u convention):
        b(i-1) /= d(i-1);
        u(i) /= d(i-1); //Note: row i in u = row i-1 in the matrix
        d(i-1) = 1.0;
        //Eliminate
        u(i+1) += u(i)*alpha;
        d(i) += b(i-1)*alpha;
    }
    //Normalize bottom row
    u(N) /= d(N-1);
    d(N-1) = 1.0;

    // Backward substitute
    for (int i = N; i > 1; i--){ // loop from i=N to i=2
        u(i-1) -= u(i)*b(i-2);
        // b(i-2) = 0.0 // This is never read, why bother...
      }

}


void backward_euler(double alpha, mat &u, int N, int T){
    /*
    Implements backward euler scheme by gaus-elimination of tridiagonal matrix.
    Results are saved to u.
    */

    for (int t = 1; t < T; t++){
        u(t) = u(t-1); // .copy() ?
        tridiag(alpha,u(t),N); //Note: Passing a pointer to row t, which is modified in-place
      }
}

void crank_nicolson(double alpha, mat &u, int N, int T){
    /*
    Implents crank-nicolson scheme, reusing code from forward- and backward euler
    */
    for (int t = 1; t < T; t++){
        forward_step(alpha/2,u(t),u(t-1),N);
        tridiag(alpha/2,u(t),N);
      }
}
void g(vec &u, int N){
    // Initial condition u(x,0) = g(x), x \in (0,1)

    double pi = 4*atan(1); // This is the constant pi
    double dx = 1.0/N;

    for (int i = 1; i < N; i++){
      u(i) = sin(pi*i*dx)
  }
}
