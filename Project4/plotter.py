import matplotlib.pyplot as plt

filenames = ["Lattice2","c_ordered","c_unordered"]


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
    plt.plot(MCcycles, energy_mean, "o")
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Energy expectation value")
    plt.figure()
    plt.plot(MCcycles, mag_mean, "o")
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Magnetization expectation value")
    plt.figure()
    plt.plot(MCcycles, Nconfigs, "o")
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("# of accepted configurations")
    plt.figure()
    plt.plot(MCcycles, Cv, "o")
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Specific heat")
    plt.figure()
    plt.plot(MCcycles, susceptibility, "o")
    plt.xlabel("# of Monte Carlo cycles")
    plt.ylabel("Susceptibility")
plt.show()
