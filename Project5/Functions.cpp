#include "Functions.h"
#include "omp.h"


// Code for solving the 1+1 dimensional diffusion equation
// du/dt = ddu/ddx on a rectangular grid of size L x (T*dt),
// with with L = 1, u(x,0) = g(x), u(0,t) = u(L,t) = 0

double pi = 4*atan(1); // This is the constant pi

void forward_step(double alpha, rowvec &u, rowvec &uPrev, int N){
    /*
    Steps forward-euler algo one step ahead.
    Implemented in a separate function for code-reuse from crank_nicolson()
    */

    for (int i = 1; i < N+1; i++){ //loop from i=1 to i=N
        u(i) = alpha*uPrev(i-1) + (1.0-2*alpha)*uPrev(i) + alpha*uPrev(i+1);
      }
    //return u;
}

void forward_step_2dim(double alpha, mat &u, mat &uPrev, int N){
    /*
    2-dimensional version of the forward-euler step. Uniform mesh.
    */

    for (int i = 1; i < N+1; i++){
      for (int j = 1; j < N+1; j++){
        u(i,j) = uPrev(i,j) + alpha*(uPrev(i+1,j) + uPrev(i-1,j)
        +uPrev(i,j+1)  + uPrev(i,j-1)- 4*uPrev(i,j));
      }
    }
}
void forward_euler(double alpha, mat &u, int N, int T){
    /*
    Implements the forward Euler sheme, results saved to
    matrix u
    */

    // Skip boundary elements
    rowvec u_temp_row;
    rowvec u_temp_row1;
    for (int t = 1; t < T; t++){

        u_temp_row = u.row(t);
        u_temp_row1 = u.row(t-1);

        forward_step(alpha,u_temp_row,u_temp_row1,N);
        u.row(t) = u_temp_row;

      }
}
void forward_euler_2dim(double alpha, cube &u, int N, int T){
    /*
    2-dimensional version of forward Euler scheme, uniform mesh, results saved to cube u
    */

    // Skip boundary elements
    mat u_temp_mat;
    mat u_temp_mat1;
    for (int t = 1; t < T; t++){

        u_temp_mat = u.slice(t);
        u_temp_mat1 = u.slice(t-1);
        //cout << u.n_slices <<u.n_rows<<u.n_cols<< endl;

        forward_step_2dim(alpha, u_temp_mat,u_temp_mat1,N);
        u.slice(t) = u_temp_mat;

      }
}

void tridiag(double alpha, rowvec &u, int N){
    /*
    Tridiagonal gaus-eliminator, specialized to diagonal = 1+2*alpha,
    super- and sub- diagonal = - alpha
    */

    vec d(N);
    d.fill(1+2*alpha);
    vec b(N-1);
    b.fill(- alpha);

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
        //b(i-2) = 0.0; // Never read
      }
}
/*
void tridiag_solver(rowvec &x, rowvec y, int N, double alpha)
{
  x.set_size(N);
  double b = -alpha;
  double d = 1+2*alpha;
  double bb = b*b;

  vec diag(N); diag.fill(d);

  y(1) += -b*y(0);
  for (int i = 2; i < N - 1; ++i) {
    diag(i) +=       -bb/diag(i-1);
       y(i) += -b*y(i-1)/diag(i-1);
  }

  x(0) = y(0); x(N - 1) = y(N - 1);
  for (int i = N - 2; i > 0; --i) {
    x(i) = (y(i) - b*x(i + 1))/diag(i);
  }

  return;
}
*/
void backward_euler(double alpha, mat &u, int N, int T){
    /*
    Implements backward euler scheme by gaus-elimination of tridiagonal matrix.
    Results are saved to u.
    */
    rowvec u_temp_row;
    for (int t = 1; t < T; t++){
        //u.row(t) = u.row(t-1); // .copy() ?
        u_temp_row = u.row(t-1);
        tridiag(alpha,u_temp_row, N); //Note: Passing a pointer to row t, which is modified in-place
        u.row(t) = u_temp_row;

      }
}

void crank_nicolson(double alpha, mat &u, int N, int T){
    /*
    Implents crank-nicolson scheme, reusing code from forward- and backward euler
    */
    rowvec u_temp_row;
    rowvec u_temp_row1;
    for (int t = 1; t < T; t++){
        u_temp_row = u.row(t);
        u_temp_row1 = u.row(t-1);
        forward_step(alpha/2,u_temp_row,u_temp_row1,N);
        tridiag(alpha/2,u_temp_row,N);
        u.row(t) = u_temp_row;

      }
}
/*
void g(mat &u, int N){
    // Initial condition u(x,0) = g(x), x \in (0,1)
    double dx = 1.0/N;
    for (int i = 1; i < N; i++){
      u(0,i) = 0;
  }
}
*/
void analytic(mat &u, int N, int T, vec x, double dt){
  double L = 1.0;
  for (int t = 0; t < u.n_rows; t++ ){
    for (int i = 0; i < u.n_cols; i++){
      u(t,i) = x(i)/L - (2.0/(pi))*sin(x(i)*pi/L)*exp(-pi*pi*t*dt/(L*L));
    }
  }
}



// Function for setting up the iterative Jacobi solver
int JacobiSolver(int N, double dx, double dt, mat &A, mat &q, double abstol)
{
  int MaxIterations = 100000;
  mat Aold = zeros<mat>(N,N);

  double D = dt/(dx*dx);

  for(int i=1;  i < N-1; i++)
    for(int j=1; j < N-1; j++)
      Aold(i,j) = 1.0;

  // Boundary Conditions -- all zeros
  for(int i=0; i < N; i++){
    A(0,i) = 0.0;
    A(N-1,i) = 0.0;
    A(i,0) = 0.0;
    A(i,N-1) = 0.0;
  }
  // Start the iterative solver
  int i, j;
  double sum = 0.0;
  for(int k = 0; k < MaxIterations; k++){
    #pragma omp parallel default(shared) private (i,j) reduction(+:sum)
    #pragma omp for
    for( i = 1; i < N-1; i++){
      for(  j=1; j < N-1; j++){
	A(i,j) = dt*q(i,j) + Aold(i,j) +
	  D*(Aold(i+1,j) + Aold(i,j+1) - 4.0*Aold(i,j) +
	     Aold(i-1,j) + Aold(i,j-1));
      }
    }
     sum = 0.0;
    for(int i = 0; i < N;i++){
      for(int j = 0; j < N;j++){
	sum += (Aold(i,j)-A(i,j))*(Aold(i,j)-A(i,j));
	Aold(i,j) = A(i,j);
      }
    }
    if(sqrt (sum) <abstol){
      return k;
    }
  }
  cout << "Jacobi: Maximum Number of Interations Reached Without Convergence\n";
  return MaxIterations;

}
