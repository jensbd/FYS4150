import matplotlib.pyplot as plt

for i in range(3):
    file = open(str(10**(i+1))+".txt")
    x_list = []
    analytic = []
    numeric = []
    file.readline()
    for line in file:
        contents = line.split()
        x_list.append(float(contents[0]))
        analytic.append(float(contents[1]))
        numeric.append(float(contents[2]))

    plt.plot(x_list, numeric,'bo', x_list, analytic, 'k')
    plt.title("General algorithm")
    plt.legend(["Numerical", "Analytic"])
    plt.xlabel("x")
    plt.ylabel("u(x)")
    plt.grid()
    plt.savefig(str(10**(i+1))+".jpg")
    plt.clf()
