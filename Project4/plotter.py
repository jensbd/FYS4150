import matplotlib.pyplot as plt
import numpy as np
filenames = ["c_unordered","c_ordered"]

"""
for i in (filenames):
    with open(i) as file:
        lines = file.readlines()

    MCcycles = []
    energy_mean = []
    mag_mean = []
    Nconfigs = []
    Cv = []
    susceptibility = []
    temperature = []
    #Skip the first two lines
    for j in range(2,len(lines)):
        line = lines[j]
        pieces = line.split()
        MCcycles.append(float(pieces[0]))
        energy_mean.append(float(pieces[1]))
        mag_mean.append(float(pieces[2]))
        Nconfigs.append(float(pieces[3]))
        Cv.append(float(pieces[4]))
        susceptibility.append(float(pieces[5]))
        temperature.append(float(pieces[6]))

    plt.figure()
    plt.plot(MCcycles, energy_mean, ".")
    #plt.axis([0, 5000,-2.1,-1.8])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Energy expectation value")

    plt.figure(),
    plt.plot(MCcycles, mag_mean, ".")
    #plt.axis([0, 5000,-2.1,-1.8])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Magnetization expectation value")

    plt.figure()
    plt.plot(MCcycles, Nconfigs, ".")
    #plt.axis([0, 5000,-2.1,-1.8])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("# of accepted configurations")

    plt.figure()
    plt.plot(MCcycles, Cv, ".")
    #plt.axis([0, 5000,-2.1,-1.8])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Specific heat")

    plt.figure()
    plt.plot(MCcycles, susceptibility, ".")
    #plt.axis([0, 5000,-2.1,-1.8])
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Susceptibility")
plt.show()
"""

"""
-------------
Probabilities
-------------
"""

filenames = ["Probability1","Probability24"]

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
    plt.bar(Energies,counts,width = 8 if i == "Probability1" else 3)
    plt.xlim(-805,-770) if i == "Probability1" else plt.xlim(-705,-305)
    plt.xlabel("Energy")
    plt.ylabel("Energy counts")
    props = dict(boxstyle='round', facecolor='wheat', alpha=1)
    plt.text(0.05*(plt.xlim()[1]-plt.xlim()[0])+plt.xlim()[0] ,plt.ylim()[1]*0.85, "Most probable energy:\n" + str(most_probable_energy), bbox = props)
    plt.show()
