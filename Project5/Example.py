# Code for solving the 1+1 dimensional diffusion equation
# du/dt = ddu/ddx on a rectangular grid of size L x (T*dt),
# with with L = 1, u(x,0) = g(x), u(0,t) = u(L,t) = 0

import numpy, sys, math
from  matplotlib import pyplot as plt
import numpy as np

def forward_step(alpha,u,uPrev,N):
    """
    Steps forward-euler algo one step ahead.
    Implemented in a separate function for code-reuse from crank_nicolson()
    """

    for x in range(1,N+1): #loop from i=1 to i=N
        u[x] = alpha*uPrev[x-1] + (1.0-2*alpha)*uPrev[x] + alpha*uPrev[x+1]

def forward_euler(alpha,u,N,T):
    """
    Implements the forward Euler sheme, results saved to
    array u
    """

    #Skip boundary elements
    for t in range(1,T):
        forward_step(alpha,u[t],u[t-1],N)

def tridiag(alpha,u,N):
    """
    Tridiagonal gaus-eliminator, specialized to diagonal = 1+2*alpha,
    super- and sub- diagonal = - alpha
    """
    d = numpy.zeros(N) + (1+2*alpha)
    b = numpy.zeros(N-1) - alpha

    #Forward eliminate
    for i in range(1,N):
        #Normalize row i (i in u convention):
        b[i-1] /= d[i-1];
        u[i] /= d[i-1] #Note: row i in u = row i-1 in the matrix
        d[i-1] = 1.0
        #Eliminate
        u[i+1] += u[i]*alpha
        d[i] += b[i-1]*alpha
    #Normalize bottom row
    u[N] /= d[N-1]
    d[N-1] = 1.0

    #Backward substitute
    for i in range(N,1,-1): #loop from i=N to i=2
        u[i-1] -= u[i]*b[i-2]
        #b[i-2] = 0.0 #This is never read, why bother...

def backward_euler(alpha,u,N,T):
    """
    Implements backward euler scheme by gaus-elimination of tridiagonal matrix.
    Results are saved to u.
    """
    for t in range(1,T):
        u[t] = u[t-1].copy()
        tridiag(alpha,u[t],N) #Note: Passing a pointer to row t, which is modified in-place

def crank_nicolson(alpha,u,N,T):
    """
    Implents crank-nicolson scheme, reusing code from forward- and backward euler
    """
    for t in range(1,T):
        forward_step(alpha/2,u[t],u[t-1],N)
        tridiag(alpha/2,u[t],N)

def g(x):
    """Initial condition u(x,0) = 0, x \in [0,1]"""
    return np.zeros(len(x))

# Step length in postion x
print("Type in the value for dx")
dx      =   float(input())
# Number of integration points along x-axis
N       =   int(1.0/dx)
# Step length in time
dt      =   0.5*dx*dx
# Number of time steps till final time
T       =   int(1.0/dt)
# Define method to use 1 = explicit scheme, 2= implicit scheme, 3 = Crank-Nicolson
print("Define method to use 1 = explicit scheme, 2= implicit scheme, 3 = Crank-Nicolson")
method  =   float(input())

#dx = 1/float(N+1)
u = numpy.zeros((T,N+2))
x = numpy.linspace (0,1,N+2)
alpha = dt/(dx**2)

#Initial codition
u[0,:] = g(x)
u[0,0] = 0.0

#Implement boundaries rigidly
for t in range(0,T):
    u[t,N+1] = 1.0;

if   method == 1:
    forward_euler(alpha,u,N,T)
elif method == 2:
    backward_euler(alpha,u,N,T)
elif method == 3:
    crank_nicolson(alpha,u,N,T)
else:
    print ("Please select method 1,2, or 3!")
    sys.exit(0)
# To do: add movie
print(u)
