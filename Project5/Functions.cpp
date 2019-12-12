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









void tridiagSolver(rowvec &u, rowvec u_prev, double alpha, int N) {
  /*
  * Thomas algorithm:
  * Solves matrix vector equation Au = b,
  * for A being a tridiagonal matrix with constant
  * elements diag on main diagonal and offdiag on the off diagonals.
  */
 double diag, offdiag;
 diag = 1+2*alpha; offdiag = -alpha;
 vec beta = zeros<vec>(N+1); beta(0) = diag;
 vec u_old = zeros<vec>(N+1); u_old(0) = u_prev(0);
 double btemp;


 // Forward substitution
 for(int i = 1; i < N+1; i++){
   btemp = offdiag/beta[i-1];
   beta(i) = diag - offdiag*btemp;
   u_old(i) = u_prev(i) - u_old(i-1)*btemp;
 }


 // Special case, boundary conditions
 u(0) = 0;
 u(N+1) = 1;

 // backward substitution
 for(int i=N; i>0; i--){
   u(i) = (u_old(i) - offdiag*u(i+1))/beta(i);
  }

}





void backward_euler(double alpha, mat &u, int N, int T){
    /*
    Implements backward euler scheme by gaus-elimination of tridiagonal matrix.
    Results are saved to u.
    */
    rowvec u_temp;
    for (int t = 1; t < T; t++){
        u_temp = u.row(t);
        tridiagSolver(u_temp, u.row(t-1), alpha, N); //Note: Passing a pointer to row t, which is modified in-place
        u.row(t) = u_temp;
      }
}

void crank_nicolson(double alpha, mat &u, int N, int T){
    /*
    Implents crank-nicolson scheme, reusing code from forward- and backward euler
    */
    rowvec u_temp;
    rowvec u_temp1;
    rowvec u_temp2;
    for (int t = 1; t < T; t++){
      u_temp = u.row(t-1);
      u_temp1 = u.row(t);
      forward_step(alpha/2, u_temp1, u_temp, N);
      u_temp2 = u.row(t);
      tridiagSolver(u_temp2, u_temp1, alpha, N);
      u.row(t) = u_temp2;
      }
}

void analytic(mat &u, int N, int T, vec x, double dt){
  double L = 1.0;
  for (int t = 0; t < u.n_rows; t++ ){
    for (int i = 0; i < u.n_cols; i++){
      u(t,i) = x(i)/L - (2.0/(pi))*sin(x(i)*pi/L)*exp(-pi*pi*t*dt/(L*L));
    }
  }
}


void analytic_2D(mat &u, double dx, double dt){
  /*
   * Analytic solution for the two dimensional diffusion equation at a
   * given time. Boundary conditions are all equal to zero. A "sine-paraboloid"
   * which is centered in the middle of the medium is the initial condition.
   * Boundary conditions are zeros.
   */
  double L = 1.0;
  int N = int(1/dx);   // Number of integration points along x & y -axis (inner points only)
  int T = int(1/dt);   // Number of time steps till final time
  vec x = linspace<vec>(0,1,N+2);


  for (int t = 0; t < u.n_rows; t++){
    for (int i = 0; i < u.n_cols; i++){
      u(t,i) = x(i)/L - (2.0/(pi))*sin(x(i)*pi/L)*exp(-pi*pi*t*dt/(L*L));
    }
  }
}



// Function for setting up the iterative Jacobi solver
int JacobiSolver(mat &u, double dx, double dt, double abstol){

  mat u_prev = u;

// Constants, variables
  int MaxIterations = 100000;
  int N = int(1/dx);   // Number of integration points along x & y -axis (inner points only)
  int T = int(1/dt);   // Number of time steps till final time

  double alpha = dt/(dx*dx);
  double factor = 1.0/(1.0 + 4*alpha);
  double factor_a = alpha*factor;
  double scale = (N+2)*(N+2);
  double delta;
//  #pragma omp parallel default(shared) private (i,j) reduction(+:sum)

// Time loop
//#pragma omp for
for (int t = 1; t < T; t++){
  int iterations = 0;
  double diff=1;
  mat u_guess = ones<mat>(N+2,N+2);

  while (iterations < MaxIterations && diff > abstol){
    diff = 0;
    // Loop over all inner elements, will converge towards solution
    for (int j = 1; j < N; j++) {
      for (int i = 1; i < N; i++) {
        // u_guess is the previous u, which also work for a random guess
        delta = (u_guess(i,j+1)+u_guess(i,j-1)+u_guess(i+1,j)+u_guess(i-1,j));
        u(i,j) = factor_a*delta + factor*u_prev(i,j);
        diff += fabs(u(i,j) - u_guess(i,j));
      }
    } // end of double for loop
    u_guess = u;
    diff /= scale;
    iterations++;
  } // end iteration loop
  u_prev = u;
} // end time loop

}
