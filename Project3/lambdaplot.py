import matplotlib.pyplot as plt
import numpy as np

lambdas = np.linspace(0,10,1000)
func = np.exp(-4*(lambdas))

plt.plot(lambdas,func)
plt.xlabel("$\lambda$")
plt.ylabel("$e^{-4\lambda}$")
plt.grid()
func3 = np.exp(-4*3)
plt.text(np.max(lambdas)*0.5,np.max(func)*0.7,"When $\lambda = 3$, \n$e^{-2\lambda} = %f \\approx 0$" % (func3))
plt.title("Single-particle wave function")
plt.show()
