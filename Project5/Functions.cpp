#include "Functions.h"
#include "omp.h"


void forward_euler(double dx, mat &u, vec xlist, vec tlist){
  double C = dt/(dx*dx)
  double dt = 0.5*dx*dx //Stability criterion
  vec u1 = u(1) //Initial condition, zeros
  for(int j = 1; j < tlist.n_elem;  j++){
    for (int i = 1; i < xlist.n_elem-1; i++){
      u(j,i) = C*(u1(i+1) - 2*u1(i) + u1(i-1)) + u1(i)
    }
    // Boundaries need not be updated, they are already zeros by def. of u-matrix
    u1 = u(j)
  }
}

void backward_euler(int n, int tsteps, double delta_x, double alpha){
   double a, b, c;
   vec u(n+1); // This is u  of Au = y
   vec y(n+1); // Right side of matrix equation Au=y, the solution at a previous step

   // Initial conditions
   for (int i = 1; i < n; i++) {
      y(i) = u(i) = func(delta_x*i);
   }
   // Boundary conditions (zero here)
   y(n) = u(n) = u(0) = y(0);
   // Matrix A, only constants
   a = c = - alpha;
   b = 1 + 2*alpha;
   // Time iteration
   for (int t = 1; t <= tsteps; t++) {
      //  here we solve the tridiagonal linear set of equations,
      tridag(a, b, c, y, u, n+1);
      // boundary conditions
      u(0) = 0;
      u(n) = 0;
      // replace previous time solution with new
      for (int i = 0; i <= n; i++) {
	       y(i) = u(i);
      }
   }
}

void crank_nicolson(int n, int tsteps, double delta_x, double alpha){
   double a, b, c;
   vec u(n+1); // This is u in Au = r
   vec r(n+1); // Right side of matrix equation Au=r
   ....
   // setting up the matrix
   a = c = - alpha;
   b = 2 + 2*alpha;

   // Time iteration
   for (int t = 1; t <= tsteps; t++) {
      // Calculate r for use in tridag, right hand side of the Crank Nicolson method
      for (int i = 1; i < n; i++) {
	 r(i) = alpha*u(i-1) + (2 - 2*alpha)*u(i) + alpha*u(i+1);
      }
      r(0) = 0;
      r(n) = 0;
      //  Then solve the tridiagonal matrix
      tridiag(a, b, c, r, u, xsteps+1);
      u(0) = 0;
      u(n) = 0;
}

void tridiag(alpha,u,N){
    """
    Tridiagonal gaus-eliminator, specialized to diagonal = 1+2*alpha,
    super- and sub- diagonal = - alpha
    """
    d = numpy.zeros(N) + (1+2*alpha)
    b = numpy.zeros(N-1) - alpha

    //Forward eliminate
    for i in xrange(1,N):
        //Normalize row i (i in u convention):
        b[i-1] /= d[i-1];
        u[i] /= d[i-1] #Note: row i in u = row i-1 in the matrix
        d[i-1] = 1.0
        //Eliminate
        u[i+1] += u[i]*alpha
        d[i] += b[i-1]*alpha
    //Normalize bottom row
    u[N] /= d[N-1]
    d[N-1] = 1.0

    //Backward substitute
    for i in xrange(N,1,-1): #loop from i=N to i=2
        u[i-1] -= u[i]*b[i-2]
        #b[i-2] = 0.0
}
