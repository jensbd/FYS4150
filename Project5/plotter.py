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

    x_analytic = np.linspace(0,1,1002)
    t_analytic = np.linspace(0,1,1000)
    dt = 1.0/len(t_analytic)
    u_analytic = np.zeros((len(t_analytic),len(x_analytic)))
    for i in tqdm(range(len(t_analytic))):
        u_analytic[i,:] = x_analytic/L
        for n in range(1,100):
            u_analytic[i] += (2*(-1)**n)/(n*np.pi)*np.sin(n*np.pi*x_analytic/L)*np.exp(-n**2*np.pi**2*t_analytic[i]/L**2)
    u_analytic[0,len(x_analytic)-1] = 1.0 # Hard coding the boundary conditions


    u_list = []
    x_list = []
    t_list = []
    for method in ["FE:", "BE:", "CN:"]:
        for dx in [0.1, 0.01]:
            dt = 0.5*dx*dx
            #Generate t-mesh
            T = int(1.0/dt) #Number of time steps till final time
            t = np.zeros(T)
            for l in range(len(t)):
                t[l] = l*dt
            #Generate x-mesh
            N = int(1.0/dx)   #Number of integration points along x-axis (inner points only)
            x = np.zeros(N+2)
            for k in range(len(x)):
                x[k] = k/(N+1)
            if method == "FE:":
                x_list.append(x)
                t_list.append(t)

            with open (method+str(dx)) as file:
                lines = file.readlines()
                u = np.zeros((len(lines),len(lines[0].split())))
                for i in range(len(lines)):
                    u[i,:] = lines[i].split()

                u_list.append(u)

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
    dx = [0.1, 0.01]
    for i in range(2):
        fig = plt.figure();
        plt.title("Numerical vs Analytical solution for t = 0 \n dx = %g" % (dx[i]))
        plt.plot(x_list[i], u_list[i][0], ".")
        plt.plot(x_list[i], u_list[2+i][0], ".")
        plt.plot(x_list[i], u_list[4+i][0], ".")
        plt.plot(x_analytic, u_analytic[0])
        plt.legend(["FE", "BE", "CN", "Analytic"])
    plt.show()

    for i in range(2):
        fig = plt.figure();
        plt.title("Numerical vs analytical solution for t = 0.5 \n dx = %g" % (dx[i]))
        plt.plot(x_list[i], u_list[i][int(len(t_list[i])/2)], ".")
        plt.plot(x_list[i], u_list[2+i][int(len(t_list[i])/2)], ".")
        plt.plot(x_list[i], u_list[4+i][int(len(t_list[i])/2)], ".")
        plt.plot(x_analytic, u_analytic[int(len(t_analytic)/2)])
        plt.legend(["FE", "BE", "CN", "Analytic"])
    plt.show()

    for i in range(2):
        fig = plt.figure();
        plt.title("Numerical vs analytical solution for t = 1.0 \n dx = %g" % (dx[i]))
        plt.plot(x_list[i], u_list[i][len(t_list[i])-1], ".")
        plt.plot(x_list[i], u_list[2+i][len(t_list[i])-1], ".")
        plt.plot(x_list[i], u_list[4+i][len(t_list[i])-1], ".")
        plt.plot(x_analytic, u_analytic[len(t_analytic)-1])
        plt.legend(["FE", "BE", "CN", "Analytic"])
    plt.show()


elif dim == "2":
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
