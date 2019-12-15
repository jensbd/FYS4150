#include "Functions.h"
#include "omp.h"


// Code for solving the 1+1 dimensional diffusion equation
// du/dt = ddu/ddx on a rectangular grid of size L x (T*dt),
// with with L = 1, u(x,0) = g(x), u(0,t) = u(L,t) = 0

double pi = 4*atan(1); // This is the constant pi

void forward_step(double alpha, rowvec &u, rowvec &uPrev, int N, bool CN){
    /*
    Steps forward-euler algo one step ahead.
    Implemented in a separate function for code-reuse from crank_nicolson()
    */
    if (CN == false){
      for (int i = 1; i < N+1; i++){ //loop from i=1 to i=N
        u(i) = alpha*uPrev(i-1) + (1.0-2*alpha)*uPrev(i) + alpha*uPrev(i+1);
      }
  }
    if (CN == true){
      for (int i = 1; i < N+1; i++){ //loop from i=1 to i=N
        u(i) = alpha*uPrev(i-1) + (2.0-2*alpha)*uPrev(i) + alpha*uPrev(i+1);
      }
  }
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

        forward_step(alpha,u_temp_row,u_temp_row1,N,false);
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





void tridiagSolver(rowvec &u, rowvec u_prev, double alpha, int N, bool CN) {
  /*
  * Thomas algorithm:
  * Solves matrix vector equation Au = b,
  * for A being a tridiagonal matrix with constant
  * elements diag on main diagonal and offdiag on the off diagonals.
  */
 double diag, offdiag;
 if (CN == false){
   diag = 1+2*alpha;
}

 if (CN == true){
  diag = 2+2*alpha;
}

 offdiag = -alpha;
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
 for(int i = N; i > 0; i--){
   u(i) = (u_old(i) - offdiag*u(i+1))/beta(i);
  }

}





void backward_euler(double alpha, mat &u, int N, int T){
    /*
    Implements backward euler scheme by gaus-elimination of tridiagonal matrix.
    Results are saved to u.
    */
    rowvec u_temp;
    rowvec u_temp1;
    for (int t = 1; t < T; t++){
        u_temp = u.row(t);
        u_temp1 = u.row(t-1);
        tridiagSolver(u_temp, u_temp1, alpha, N, false); //Note: Passing a pointer to row t, which is modified in-place
        u_temp(0) = 0;
        u_temp(N+1) = 1;
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
      forward_step(alpha, u_temp1, u_temp, N, true);

      u_temp1(0) = 0;
      u_temp1(N+1) = 1;

      u_temp2 = u.row(t);
      tridiagSolver(u_temp2, u_temp1, alpha, N, true);
      u_temp2(0) = 0;
      u_temp2(N+1) = 1;
      u.row(t) = u_temp2;
      }
}

void analytic(mat &u, int N, int T, int infty) {
  /*
   * Analytic solution for the one dimensional diffusion equation,
   * with boundary conditions:
   * u(L,t) = 1 and u(0,t) = 0, for all t.
   * u(x,0) = 0 for x < L.
   */
  double L = 1.0; // scale the rod such that x goes from 0 to L=1.
  double sum;
  vec x = linspace<vec>(0,1,N+2);
  vec t = linspace<vec>(0,1,T);

  // time loop
  for (int i = 0; i < T; i++){
    // position x loop
    for (int j = 1; j < N+1; j++){
      // calculate the transient solution.
      sum = 0;
      for (int n = 1; n < infty; n++){
        sum +=  (2*(-1)*n*n)/(n*pi)*sin(n*pi*x(j)/L)*exp(-n*n*pi*pi*t(i)/L*L);
      }
      u(i,j) = x(j)/L + sum;
    }
  }
}


void analytic_2D(mat &u, int N, int T){
  /*
   * Analytic solution for the two dimensional diffusion equation at a
   * given time. Boundary conditions are all equal to zero. A "sine-paraboloid"
   * which is centered in the middle of the medium is the initial condition.
   * Boundary conditions are zeros.
   */
  double L = 1.0;
  vec x = linspace<vec>(0,1,N+2);
  vec y = linspace<vec>(0,1,N+2);
  vec t = linspace<vec>(0,1,T);


  for (int i = 0; i < N+1; i++){
    for (int j = 0; j < N+1; j++){
      u(i,j) = sin(pi*x(j))*sin(pi*y(j))*exp(-2*pi*pi*t(i));
    }
  }
}



// Function for setting up the iterative Jacobi solver
int JacobiSolver(mat &u, double dx, double dt, double abstol){
  ofstream ofile;
  string file = "2dim_implicit:"+to_string(dx);
  file.erase ( file.find_last_not_of('0') + 1, std::string::npos );
  ofile.open(file);
  ofile << u;

  mat u_prev = u;

// Constants, variables
  int MaxIterations = 100000; // Max number of iterations
  int N = int(1/dx);   // Number of integration points along x & y -axis (inner points only)
  int T = int(1/dt);   // Number of time steps till final time

  double alpha = dt/(dx*dx);
  double factor = 1.0/(1.0 + 4*alpha);
  double factor_a = alpha*factor;
  double scale = (N+2)*(N+2);
  double delta;
  double diff;
  int i,j;


// Time loop
// Parallelization using openmp
//#pragma omp parallel for
for (int t = 1; t < T; t++){
  int iterations = 0;
  diff=1;
  mat u_guess = ones<mat>(N+2,N+2);

  while (iterations < MaxIterations && diff > abstol){
    diff = 0;
    // Define parallel region
    // Loop over all inner elements, will converge towards solution
    for (j = 1; j < N+1; j++) {
      for (i = 1; i < N+1; i++) {
        // u_guess is the previous u, which also work for a random guess
        delta = (u_guess(i,j+1)+u_guess(i,j-1)+u_guess(i+1,j)+u_guess(i-1,j));
        u(i,j) = factor_a*delta + factor*u_prev(i,j);
        diff += fabs(u(i,j) - u_guess(i,j));
          }
        } // end of double for loop
    u_guess = u;
    diff /= scale;
    iterations++;
    }  // end iteration loop
  u_prev = u;
  ofile << u;
  } // end time loop
ofile.close();
}


void Lithosphere(int Case, double dx, double dt){
  /*
   * Solver throughout 1 Ga of the lithosphere in three zones:
   * Zone 1: 0-20km,
   * Zone 2: 20-40km,
   * Zone 3: 40-120km.
   * Case=1: no heat production
   * Case=2: natural heat production
   * Case=3: enriched mantle (zone 3) without decay
   * Case=4: enriched mantle (zone 3) witho decay
   */

   ofstream ofile;
   string file;
   ofile.open(file);


  double mantleQ = 0; // represents the constant Q term in the mantle (zone 3)
  if (Case == 1) {
    cout << "Solving lithosphere with no heat production\n";
    file = "No_Heat";
  } else if (Case == 2) {
    cout << "Solving lithosphere with natural heat production\n";
    file = "Heat";
    mantleQ = 0.05;
  } else if (Case == 3) {
    cout << "Solving lithosphere with enrichment without decay\n";
    file = "No_Decay";
    mantleQ = 0.55;
  } else if (Case == 4) {
    cout << "Solving lithosphere with enrichment with decay\n";
    file = "Decay";
    mantleQ = 0.05;
  } else {
    cout << "Error, Case must be 1,2,3,4.\n";
    exit(1);
  }

  // initialization
  int nx = int(1.25/dx); // 0-150 km wide
  int ny = int(1/dx); // 0-120 km deep
  int nt = int(1/dt); // number of time steps
  //double dx = 0.01; // step length [1.2 km]: 0.01*[1.2 km]*100 points = 120 km
  //double dt = 0.01; // 1 = [10^9 yr]: 100*dt = 1 Ga
  mat u = zeros<mat>(nx+1,ny+1); // Scaled from 0-1 -> 8-1300 deg Celsius
  if (Case==1) {
    NoHeat(u,nx,ny);
  } else {
    Heat(u,nx,ny);
  }
  mat u_old = u;

  // Constants and variables
  int maxiter = 10000; // Max no. of iterations each time step in Jacobi Method
  double delta, diff, Q, time; // temporary constant that will be used
  double scale = (nx+1)*(ny+1); // No. of points in matrix
  double tol = 1e-10; // tolerance for convergence in Jacobi method
  double alpha = dt/(dx*dx);
  double rho = 3.51*1e3; // density [kg/m^3]
  double cp = 1000; // heat capacity [J/kg/deg Celsius]
  double k = 2.5; // thermal conductivity [W/m/deg Celsius]
  double uscale = 1292; // temperature scale
  double tscale = 3.1557e16; // time scale
  double xscale = 120000; // distance scale
  double beta = tscale*k/(rho*cp*xscale*xscale);
  double Qscale = 1e-6*tscale*dt/(rho*cp*uscale);
  double factor1 = beta*alpha;
  double factor2 = 1/(1 + 4*factor1);

  // Construction of Qvec: stores constant heat production in zones:
  vec Qvec = zeros<vec>(ny+1);
  if (Case != 1) {
    double Q1 = 1.4*Qscale;
    double Q2 = 0.35*Qscale;
    double Q3 = mantleQ*Qscale;
    for (int i=0; i<17; i++) Qvec(i) = Q1; // 0-20 km depth
    for (int i=17; i<34; i++) Qvec(i) = Q2; // 20-40 km depth
    for (int i=34; i<=ny; i++) Qvec(i) = Q3; // 40+ km depth
  }

  // EConstruction of Qtime: stores time dependet heat production (due to decay):
  vec Qtime = zeros<vec>(nt+1);
  if (Case == 4) { // Qtime is simply 0 for cases 1-3
    double U_enrich, Th_enrich, K_enrich;
    for (int t=0; t<=nt; t++) {
      time = (double) t/nt;
      U_enrich = exp (-0.155*time);
      Th_enrich = exp (-0.0495*time);
      K_enrich = exp (-0.555*time);
      Qtime(t) = Qscale*(0.2*U_enrich + 0.2*Th_enrich + 0.1*K_enrich);
    }
  }

  // time loop
  for (int t = 1; t <= nt; t++)
  {
    int iter = 0;
    double diff=1;
    mat u_guess = u; // Guess matrix, can be anything
    // Jacobi method loop
    while (iter < maxiter && diff > tol)
    {
      diff = 0;
      for (int j=1; j<ny; j++) {
          // Case separation
          Q=Qvec(j); // constant Q for depth at index j
          if (j >= 34) Q += Qtime(t); // time dependent Q, only at depth j>=34
        for (int i=1; i<nx; i++) {
          // Jacobi method algorithm, u should converge towards solution
          delta = (u_guess(i,j+1)+u_guess(i,j-1)+u_guess(i+1,j)+u_guess(i-1,j));
          u(i,j) = factor2*(factor1*delta + Q + u_old(i,j));
          diff += fabs(u(i,j) - u_guess(i,j));
        }
      }
      u_guess = u;
      diff /= scale;
      iter++;
    } // end Jacobi method loop
    u_old = u;
    //cout << "timestep: " <<  t << " iterations: " << iter << endl;
    ofile << u;
  } // end time loop
  ofile.close();
}

void NoHeat(mat& u, int nx, int ny) {
  // Creates boundary conditions for the case of Q = 0 everywhere
  for (int i=0; i<=nx; i++) {
    u(i,0) = 0;
    u(i,ny) = 1;
  }
  // linear temperature from 0-1 (scaled)
  for (int i=0; i<=nx; i+=nx) {
    for (double j=0; j<=ny; j++) {
      u(i,j) = j/ny;
      u(i,j) = j/ny;
    }
  }
}

void Heat(mat& u, int nx, int ny) {
  // Analytic solution to steady state Temp(depth) with natural heat production
  // temperature is scaled [8-1300] -> [0-1292] -> [0,1]
  double y, temp;
  for (double j=0; j<17; j++) {
    y = j*1.2;
    temp = (-0.28*y*y + 23.66*y)/1292.0;
    for (int i=0; i<=nx; i++) u(i,j) = temp;
  }
  for (double j=17; j<34; j++) {
    y = j*1.2;
    temp = (-0.07*y*y + 15.26*y + 86)/1292.0;
    for (int i=0; i<=nx; i++) u(i,j) = temp;
  }
  for (double j=34; j<=ny; j++) {
    y = j*1.2;
    temp = (-0.01*y*y + 10.46*y + 180)/1292.0;
    for (int i=0; i<=nx; i++) u(i,j) = temp;
  }
}
