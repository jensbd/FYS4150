import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib.pyplot import cm
from tqdm import tqdm

print("Do you want to run 1 or 2-dimensional?")
print("Write 1 or 2")

dim = input("Write here: ")

if dim == "1":
    L = 1.0

    x_analytic = np.linspace(0,1,1000)
    t_analytic = np.linspace(0,1,1000)
    dt = 1.0/len(t_analytic)
    u_analytic = np.zeros((len(x_analytic),len(t_analytic)))
    for i in tqdm(range(len(t_analytic))):
        u_analytic[i,:] = x_analytic/L
        for n in range(1,1000):
            u_analytic[i] += (2*(-1)**n)/(n*np.pi)*np.sin(n*np.pi*x_analytic/L)*np.exp(-n**2*np.pi**2*t_analytic[i]/L**2)



    for method in ["FE:", "BE:", "CN:"]:
        for dx in [0.1, 0.01]:
            """
            with open ("Analytic:"+str(dx)) as file:
                lines = file.readlines()
                u_analytic = np.zeros((len(lines),len(lines[0].split())))
                for i in range(len(lines)):
                    u_analytic[i,:] = lines[i].split()
            """
            dt = 0.5*dx*dx
            #Generate t-mesh
            t = np.zeros(int(1.0/dt))
            for l in range(len(t)):
                t[l] = l*dt
            #Generate x-mesh
            x = np.zeros(int(1.0/dx)+2)
            for k in range(len(x)):
                x[k] = k*dx

            with open (method+str(dx)) as file:
                lines = file.readlines()
                u = np.zeros((len(lines),len(lines[0].split())))
                for i in range(len(lines)):
                    u[i,:] = lines[i].split()

                fig = plt.figure();
                x,t = np.meshgrid(x,t)
                ax = fig.gca(projection='3d');
                # Plot the surface.
                surf = ax.plot_surface(x, t, u, cmap=cm.coolwarm,
                                   linewidth=0, antialiased=False);
                                   # Customize the z axis.
                #ax.set_zlim(-0.10, 1.40);
                for angle in range(0,230):
                    ax.view_init(40,angle)
                ax.zaxis.set_major_locator(LinearLocator(10));
                ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'));
                plt.xlabel("x")
                plt.ylabel("t")
                name = method+" dx = "+str(dx)
                plt.title(name)
                fig.savefig("plots/"+name+".png")


                if method == "FE:" and dx == 0.1:
                    fig = plt.figure()
                    x,t = np.meshgrid(x_analytic,t_analytic)
                    ax = fig.gca(projection='3d');
                    # Plot the surface.
                    surf = ax.plot_surface(x, t, u_analytic, cmap=cm.coolwarm,
                                       linewidth=0, antialiased=False);
                                       # Customize the z axis.
                    #ax.set_zlim(-0.10, 1.40);
                    for angle in range(0,230):
                        ax.view_init(40,angle)
                    ax.zaxis.set_major_locator(LinearLocator(10));
                    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'));
                    plt.xlabel("x")
                    plt.ylabel("t")
                    plt.title("Analytic")
                    fig.savefig("plots/Analytic.png")

                plt.show()

if dim == "2":
    for dx in [0.1,0.01]:
        dt = 0.2*dx*dx
        T = int(0.1/dt)
        #Generate t-mesh
        t = np.linspace(0,0.1,T)
        #Generate x- and y-mesh
        N = int(1.0/dx)
        x = np.linspace(0,1,N+2)
        y = np.linspace(0,1,N+2)
        filename = "2dim_explicit:"+str(dx)
        with open(filename) as file:
            lines = file.readlines()
            for t in tqdm(range(T)):
                u = np.zeros((len(x),len(y)))
                for i in range(len(y)):
                    data = lines[t*len(x)+i].split()

                    u[i] = data

                fig = plt.figure();
                x_,y_ = np.meshgrid(x,y)

                ax = fig.gca(projection='3d');
                # Plot the surface.
                surf = ax.plot_surface(x_, y_, u, cmap=cm.coolwarm,
                                   linewidth=0, antialiased=False);
                                   # Customize the z axis.
                #ax.set_zlim(-0.10, 1.40);
                for angle in range(0,230):
                    ax.view_init(40,angle)
                ax.zaxis.set_major_locator(LinearLocator(10));
                ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'));
                plt.xlabel("x")                                         
                plt.ylabel("y")
                name = "2dim_explicit: dx = "+str(dx)
                plt.title(name)
                fig.savefig("plots/2dim/"+str(dx)+"/"+name+","+str(t)+".png")

else:
    print("Please write either 1 or 2")
