import matplotlib.pyplot as plt
import numpy as np

print("Which Project Task do you want to run")
print("Task C - equilibrium plots, write c")
print("Task D - probability histogram, write d")
print("Task E & F - phase transitions and critical temperature, write e")
Task = input("Write here: ")

"""
-------------
Equilibrium
-------------
"""

if Task == "c":
    filenames = ["Ordered1","Ordered"]
    filenames2 = ["Unordered1", "Unordered"]

    MCcycles = []
    energyO = []
    energyO2 = []
    energyU = []
    energyU2 = []
    magO = []
    magO2 = []
    magU = []
    magU2 = []
    NconfigsU = []
    NconfigsU2 = []
    for i in (filenames):
        with open(i) as file:
            lines = file.readlines()
            #Skip the first two lines
            for j in range(2,len(lines)):
                line = lines[j]
                pieces = line.split()
                if i == "Ordered1":
                    MCcycles.append(float(pieces[0]))
                    energyO.append(float(pieces[1]))
                    magO.append(float(pieces[2]))
                else:
                    energyO2.append(float(pieces[1]))
                    magO2.append(float(pieces[2]))

    for i in (filenames2):
        with open(i) as file:
            lines = file.readlines()
            #Skip the first two lines
            for j in range(2,len(lines)):
                line = lines[j]
                pieces = line.split()
                if i == "Unordered1":
                    energyU.append(float(pieces[1]))
                    magU.append(float(pieces[2]))
                    NconfigsU.append(float(pieces[3]))

                else:
                    energyU2.append(float(pieces[1]))
                    magU2.append(float(pieces[2]))
                    NconfigsU2.append(float(pieces[3]))
    plt.figure()
    plt.title("Ordered")
    plt.plot(MCcycles, energyO)
    plt.plot(MCcycles,energyO2)
    plt.legend(["T = 1.0","T = 2.4"])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Energy expectation value $\langle$E$\\rangle$")

    plt.figure()
    plt.title("Unordered")
    plt.plot(MCcycles, energyU)
    plt.plot(MCcycles,energyU2)
    plt.legend(["T = 1.0","T = 2.4"])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Energy expectation value $\langle$E$\\rangle$")

    plt.figure()
    plt.title("Ordered")
    plt.plot(MCcycles, magO, "")
    plt.plot(MCcycles, magO2, "")
    plt.legend(["T = 1.0","T = 2.4"])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Magnetization expectation value $\langle$|M|$\\rangle$")

    plt.figure()
    plt.title("Unordered")
    plt.plot(MCcycles, magU, "")
    plt.plot(MCcycles, magU2, "")
    plt.legend(["T = 1.0","T = 2.4"])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Magnetization expectation value $\langle$|M|$\\rangle$")

    plt.figure()
    plt.title("Unordered")
    plt.plot(MCcycles, NconfigsU, "")
    plt.plot(MCcycles, NconfigsU2, "")
    plt.legend(["T = 1.0","T = 2.4"])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Accepted configurations (normalized)")

    Temp = []
    configs = []
    with open("Nconfig_vs_Temp") as file:
        lines = file.readlines()
        for i in range(2,len(lines)):

            pieces = lines[i].split()
            Temp.append(float(pieces[0]))
            configs.append(float(pieces[1]))
    plt.figure()
    plt.plot(Temp,configs)
    plt.xlabel("Temperature [kT/J]")
    plt.ylabel("Accepted number of configurations (normalized)")
    plt.title("Accepted number of configurations (normalized) as a function of T")
    plt.show()

"""
-------------
Probabilities
-------------
"""

if Task == "d":
    filenames = ["Probability_1","Probability_24"]


    for i in filenames:
        with open(i) as file:
            lines = file.readlines()
        Energies = []
        counts = []
        max_count = 0
        most_probable_energy = 0
        for j in range(1,len(lines)):
            line = lines[j]
            pieces = line.split()
            energy = float(pieces[0])
            count = float(pieces[1])
            Energies.append((energy))
            counts.append((count))
            if count > max_count:
                max_count = count
                most_probable_energy = energy
        plt.bar(Energies,counts,width = 4 if i == "Probability_1" else 3)
        plt.xlim(-805,-770) if i == "Probability_1" else plt.xlim(-705,-305)
        plt.xlabel("Energy")
        plt.ylabel("Energy counts")
        if i == "Probability_1":
            plt.title("T = 1.0")
        else:
            plt.title("T = 2.4")
        props = dict(boxstyle='round', facecolor='wheat', alpha=1)
        plt.text(0.05*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0] ,plt.ylim()[1]*0.85, "Most probable energy:\n" + str(most_probable_energy), bbox = props)
        plt.show()

if Task == "e":
    with open("Temperature") as file:
        lines = file.readlines()
    temps = []
    energylist = []
    maglist = []
    Cvlist = []
    Suscplist = []
    indeks = 0
    for i in range(1, len(lines)):
        pieces = lines[i].split()
        temps.append(float(pieces[0]))
        energylist.append(float(pieces[1]))
        maglist.append(float(pieces[2]))
        Cvlist.append(float(pieces[3]))
        Suscplist.append(float(pieces[4]))

    firstTemp = temps[0]
    for i in range(1,len(temps)):
        if temps[i] == firstTemp:
            temps = temps[0:i]
            break

    TCCv = []
    TCX = []
    for i in range(int(len(energylist)/len(temps))):
        max_temp = 0
        sublistCv = Cvlist[i*len(temps):len(temps)*(i+1)]
        sublistSuscp = Suscplist[i*len(temps):len(temps)*(i+1)]
        maxCv = max(sublistCv)
        maxSuscp = max(sublistSuscp)
        TCCv.append(temps[sublistCv.index(maxCv)])
        TCX.append(temps[sublistSuscp.index(maxSuscp)])
        print("TC for Cv =",temps[sublistCv.index(maxCv)])
        print("TC for X =",temps[sublistSuscp.index(maxSuscp)])

    plt.figure()
    plt.title("Mean Energy")
    plt.xlabel("T [kT/J]")
    plt.ylabel("Energy expectation value $\langle$E$\\rangle$")
    for i in range(int(len(energylist)/len(temps))):
        plt.plot(temps,energylist[i*len(temps):len(temps)*(i+1)])
    plt.legend(["L = 40","L = 60","L = 80","L = 100"])

    plt.figure()
    plt.title("Absolute mean Magnetization")
    plt.xlabel("T [kT/J]")
    plt.ylabel("Magnetization expectation value $\langle$|M|$\\rangle$")
    for i in range(int(len(energylist)/len(temps))):
        plt.plot(temps,maglist[i*len(temps):len(temps)*(i+1)])
    plt.legend(["L = 40","L = 60","L = 80","L = 100"])

    plt.figure()
    plt.title("Specific heat")
    plt.xlabel("T [kT/J]")
    plt.ylabel("Specific heat $\langle$$C_v$$\\rangle$")
    for i in range(int(len(energylist)/len(temps))):
        plt.plot(temps,Cvlist[i*len(temps):len(temps)*(i+1)])
    plt.legend(["L = 40","L = 60","L = 80","L = 100"])


    plt.figure()
    plt.title("Susceptibility")
    plt.xlabel("T [kT/J]")
    plt.ylabel("Susceptibility $\langle$$\chi$$\\rangle$")
    for i in range(int(len(energylist)/len(temps))):
        plt.plot(temps,Suscplist[i*len(temps):len(temps)*(i+1)])
    plt.legend(["L = 40","L = 60","L = 80","L = 100"])
    plt.show()

    """
    Task f)
    """
    #Performing a linear regression to find critical temp in thermodyn. limit
    meanTC = 0.5*(np.array(TCCv) +np.array(TCX))
    Llist = np.array([40,60,80,100])
    Llist = 1.0/Llist

    linreg = np.polyfit(Llist,meanTC,1)
    values = np.polyval(linreg, np.linspace(2.275,2.325,100))
    plt.plot(Llist,meanTC,"o")
    plt.plot(Llist,np.polyval(linreg,Llist))
    plt.legend(["Mean $T_C$","Linear fit = %gx + %g" % (linreg[0],linreg[1])])

    plt.show()
