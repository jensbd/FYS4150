import matplotlib.pyplot as plt
import numpy as np

print("Which Project Task do you want to run")
print("Task C - equilibrium plots, write c")
print("Task D - probability histogram, write d")

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
print(len(energyU))
plt.figure()
plt.title("Ordered")
plt.plot(MCcycles, energyO)
plt.plot(MCcycles,energyO2)
#plt.axis([0, 5000,-2.1,-1.8])
plt.legend(["T = 1.0","T = 2.4"])
plt.xlabel("# of Monte Carlo cycles")
plt.ylabel("Energy expectation value $\langle$E$\\rangle$")

plt.figure()
plt.title("Unordered")
plt.plot(MCcycles, energyU)
plt.plot(MCcycles,energyU2)
#plt.axis([0, 5000,-2.1,-1.8])
plt.legend(["T = 1.0","T = 2.4"])
plt.xlabel("# of Monte Carlo cycles")
plt.ylabel("Energy expectation value $\langle$E$\\rangle$")

plt.figure()
plt.title("Ordered")
plt.plot(MCcycles, magO, "")
plt.plot(MCcycles, magO2, "")
#plt.axis([0, 5000,-2.1,-1.8])
plt.legend(["T = 1.0","T = 2.4"])
plt.xlabel("# of Monte Carlo cycles")
plt.ylabel("Magnetization expectation value $\langle$M$\\rangle$")

plt.figure()
plt.title("Unordered")
plt.plot(MCcycles, magU, "")
plt.plot(MCcycles, magU2, "")
#plt.axis([0, 5000,-2.1,-1.8])
plt.legend(["T = 1.0","T = 2.4"])
plt.xlabel("# of Monte Carlo cycles")
plt.ylabel("Magnetization expectation value $\langle$M$\\rangle$")

plt.figure()
plt.title("Unordered")
plt.plot(MCcycles, NconfigsU, "")
plt.plot(MCcycles, NconfigsU2, "")
#plt.axis([0, 5000,-2.1,-1.8])
plt.legend(["T = 1.0","T = 2.4"])
plt.xlabel("# of Monte Carlo cycles")
plt.ylabel("Accepted configurations (normalized)")

"""
plt.figure()
plt.plot(MCcycles, Cv, "")
#plt.axis([0, 5000,-2.1,-1.8])
plt.xlabel("# of Monte Carlo cycles")
plt.ylabel("Specific heat $\langle$$C_v$$\\rangle$")
"""
"""
plt.figure()
plt.plot(MCcycles, susceptibility, "")
#plt.axis([0, 5000,-2.1,-1.8])
plt.xlabel("# of Monte Carlo cycles")
plt.ylabel("Susceptibility $\langle$$\chi$$\\rangle$")
"""
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
        props = dict(boxstyle='round', facecolor='wheat', alpha=1)
        plt.text(0.05*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0] ,plt.ylim()[1]*0.85, "Most probable energy:\n" + str(most_probable_energy), bbox = props)
        plt.show()
