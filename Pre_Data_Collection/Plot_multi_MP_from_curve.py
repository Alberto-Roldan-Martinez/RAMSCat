



import sys
import numpy as np
import matplotlib.pyplot as plt


site = "dist_" + str(sys.argv[1])                                         # site description, e.g. Mg
if sys.argv[1] == "O":
    x_min = 1.5
    x_max = 5
    y_min = -3
    y_max = 0.1
elif sys.argv[1] == "Mg":
    x_min = 2
    x_max = 5
    y_min = -1.5
    y_max = 0.1

x = np.linspace(x_min, x_max, 101)
y_mean = np.zeros(len(x))
#figure = plt.figure(figsize=(12, 12), clear=True)

for i in range(len(sys.argv)-2):
    n = i + 2
    curve = open(sys.argv[n], "r").readlines()[0].split()
    d = float(curve[1])
    a = float(curve[3])
    r = float(curve[5])
    y = []

    for j in x:
       y.append(d * (np.exp(-2 * a * (j - r)) - 2 * np.exp(-1 * a * (j - r))))
    y_mean = [(y_mean[j] + y[j]) for j in range(len(y))]

    plt.plot(x, y, label=str(i+1))

y_mean = [i/(len(sys.argv)-2) for i in y_mean]
plt.plot(x, y_mean, "k:", label="mean", lw=2)
#plt.plot(x, y_mean, "k", lw=10, alpha=0.15)

plt.xlabel("$dist_{" + str(sys.argv[1]) + "}$" +" ($\\AA$)", fontsize=14)
plt.ylabel("$E_{Adh}$ (eV)", fontsize=14)
#plt.xticks(np.arange(0,Xmax, 0.25))
plt.xlim([x_min, x_max])#[min(x)*0.9, max(x)*1.2])
plt.ylim([y_min, y_max]) #[min(y)*1.2, 0.1])
#plt.tick_params(axis='both', labelrotation=0, labelsize=16)               # custimise tick labels
#plt.grid(True)
plt.legend(loc='best')
plt.tight_layout()
#plt.title("$TM_{" + str(int(data[0,0])) + "}$ $dist_{" + str(sys.argv[2]) + "}$")
#plt.savefig("MP" + str(int(data[0,0])) + "_" + site +".png", figsize=(11.69, 16.53), clear=True, dpi=300, orientation='landscape', transparent=True)
plt.show()


